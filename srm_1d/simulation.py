"""
simulation.py — Simulation Driver (Compiled Time Loop)
========================================================

The public API is run_simulation(), which is unchanged from previous
versions. Internally, the time loop is now a single @njit compiled
function (_run_time_loop) that calls all sub-functions as direct C
function calls with zero Python dispatch overhead.

Architecture:
    run_simulation()           ← Public API (Python)
        → extract scalars/arrays from Python objects
        → pre-allocate output arrays
        → _run_time_loop(...)  ← Single @njit function (C speed)
            → piso_step, post_piso_update, etc. (Numba→Numba)
        → wrap results into dict + summary
        → return
"""

import numpy as np
import time as clock

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper

from .solver import (
    FLOOR_MIN_RAW_P_NONE,
    LIMIT_TEMPERATURE_FLOOR,
    LIMIT_TEMPERATURE_CEILING,
    LIMIT_PRESSURE_FLOOR,
    N_LIMIT_DIAGNOSTICS,
    _piso_step_with_energy_diagnostics,
    compute_dt_cfl,
    compute_dt_source_cap,
    _nozzle_boundary_flow,
)
from .burn_rate import (
    compute_burn_rates, haaland_friction, gnielinski_nusselt,
    resolve_erosive_model,
)
from .igniter_plenum import (
    _cyl_area_fraction,
    _step_plenum_ode,
    chamber_params,
    initial_plenum_state,
    pyrogen_params,
)
from .solid_thermal import (_step_goodman_ode, _step_peretz_ode,
                            _step_conduction, _surface_has_ignited,
                            build_solid_grid, CONDUCTION_DY0,
                            CONDUCTION_RATIO, CONDUCTION_AGE_FRACTION)
from .propellant import (
    R_UNIVERSAL,
    create_gas_properties, speed_of_sound, critical_flow_function,
)
from .grain_geometry import (
    update_cell_geometry, advance_endface_regression,
    advance_bore_regression,
)


# ================================================================
# Snapshot channel indices (for the 3D snapshot array)
# ================================================================
_SNAP_P = 0
_SNAP_U = 1
_SNAP_MACH = 2
_SNAP_T = 3
_SNAP_R_TOTAL = 4
_SNAP_R_EROSIVE = 5
_SNAP_D_PORT = 6
_SNAP_C_BURN = 7
_SNAP_ENDFACE = 8
_SNAP_IS_BURNING = 9
_SNAP_IS_GRAIN = 10
_SNAP_T_SURF = 11
_SNAP_MASS_SOURCE = 12
_SNAP_THERMAL_SOURCE = 13
_SNAP_MOMENTUM_SOURCE = 14
_SNAP_PYROGEN_SURFACE_HEAT_FLUX = 15
_SNAP_RADIATION_HEAT_FLUX = 16
_SNAP_REGRESS = 17
_SNAP_RHO = 18          # per-cell gas density [kg/m^3]; enables mass flux G = rho*u
N_SNAP_CHANNELS = 19

CAL_CM2_S_TO_W_M2 = 41840.0
STEFAN_BOLTZMANN = 5.670374419e-8

# v0.7.4 Phase Z: burn-rate floor [m/s] in the Z-N relaxation time
# tau = kappa_zn * alpha_solid / max(r_dyn, ZN_R_FLOOR)^2, guarding the
# divide-by-zero at ignition. Internal constant, not a calibration knob.
ZN_R_FLOOR = 1.0e-4

# ================================================================
# v0.7.1 species indices (used inside @njit kernels — Numba resolves
# module-level constants at compile time). Mirrors the run_simulation
# species_list ordering.
# ================================================================
_SPECIES_IGNITER = 0
_SPECIES_PROPELLANT = 1
_SPECIES_AMBIENT = 2


# ================================================================
# Fused per-step helpers (called from inside _run_time_loop)
# ================================================================

@njit(cache=True, fastmath=True)
def _post_piso_update(
    rho, u, P, T, D_hyd, Re, Mach, u_cell, f_darcy,
    N, mu_gas, gamma_arr, R_arr, roughness,
):
    """Post-PISO: velocities, Re, Mach, friction, a_max — single pass.

    v0.7.1 (Phase 3): the sound speed and CFL wavespeed seed now use
    per-cell ``gamma_arr[i] * R_arr[i] * T[i]``. The returned ``a_max``
    is the maximum local sound speed across all cells (instead of the
    sqrt at the global hottest T scaled by a constant γR).
    """
    a_max = 0.0
    for i in range(N):
        u_cell[i] = 0.5 * (u[i] + u[i + 1])
        gR_i = gamma_arr[i] * R_arr[i]
        a_local = (gR_i * T[i]) ** 0.5
        Mach[i] = u_cell[i] / a_local
        Re[i] = rho[i] * abs(u_cell[i]) * D_hyd[i] / mu_gas
        f_darcy[i] = haaland_friction(Re[i], roughness, D_hyd[i])
        if a_local > a_max:
            a_max = a_local
    return a_max


@njit(cache=True, fastmath=True)
def _bare_heat_transfer_coeff(
    Re_local, D_hyd, x_from_head, f_local, Pr, k_thermal,
    T_gas, T_wall, kappa,
):
    """Gas-side heat-transfer coefficient for unignited Goodman heating."""
    if D_hyd <= 1e-10 or k_thermal <= 0.0:
        return 0.0
    Nu = gnielinski_nusselt(
        Re_local, Pr, D_hyd, x_from_head, f_local,
        T_gas, T_wall, kappa,
    )
    return Nu * k_thermal / D_hyd


# ---------------------------------------------------------------------------
# IG2: opt-in ignition heat-transfer closures (docs/pyrogen_realism/IG2_REPORT.md).
# They replace ONLY the h used to heat the UNIGNITED propellant surface in the
# Goodman ignition kernel (and its matching convective wall-loss sink). The Ma
# 2020 erosive path (burn_rate.py) is untouched.
#   0 'gnielinski' (pre-IG3 default) -- the Gnielinski/Ma chain.
#   1 'peretz' (DEFAULT since IG3) -- Peretz, Caveny, Kuo & Summerfield, Princeton AMS Rep. 1100
#                  (Apr 1973), Eq. III-22..III-26: bore heat transfer for the
#                  starting transient of high-velocity-bore SRMs.
#   2 'keller'  -- Keller, Baer & Ryan, AIAA J. 4(8):1358 (1966), Eq. 11,
#                  extrapolated to the local bore diameter (see function doc).
IGN_H_GNIELINSKI = 0
IGN_H_PERETZ = 1
IGN_H_KELLER = 2
_IGN_H_MODES = {'gnielinski': IGN_H_GNIELINSKI, 'peretz': IGN_H_PERETZ,
                'keller': IGN_H_KELLER}


def _ignition_h_code(mode):
    try:
        return _IGN_H_MODES[mode]
    except (KeyError, TypeError):
        raise ValueError(
            f"unknown ignition_heat_transfer: {mode!r} "
            f"(expected one of {sorted(_IGN_H_MODES)})") from None


# IG3: which port area the Peretz (A_p/A_t)^0.4 factor uses. The sources
# (Rocketdyne duct data, Peretz's own motor) only had uniform-port geometry, so
# 'local', 'head' and 'aft' coincide there and none is stated for a non-uniform
# bore (docs/pyrogen_realism/R9_SOURCES.md). Only read in 'peretz' mode.
#   0 'local' (default) -- the cell's own port area (IG2's choice).
#   1 'head'  -- port area of the fore-most grain cell (instantaneous).
#   2 'aft'   -- port area of the aft-most grain cell (instantaneous).
PERETZ_AP_LOCAL = 0
PERETZ_AP_HEAD = 1
PERETZ_AP_AFT = 2
_PERETZ_AP_MODES = {'local': PERETZ_AP_LOCAL, 'head': PERETZ_AP_HEAD,
                    'aft': PERETZ_AP_AFT}


def _peretz_ap_code(mode):
    try:
        return _PERETZ_AP_MODES[mode]
    except (KeyError, TypeError):
        raise ValueError(
            f"unknown peretz_ap_reference: {mode!r} "
            f"(expected one of {sorted(_PERETZ_AP_MODES)})") from None


# SM1 (docs/pyrogen_realism/SM1_REPORT.md): solid-heating model of the ignition
# gate. 'goodman' (default, bit-identical to pre-SM1) = algebraic-in-h cubic
# heat-balance integral; 'peretz' = Peretz's T_s ODE (III-50); 'conduction' =
# exact 1-D transient conduction on a graded grid.
SOLID_GOODMAN = 0
SOLID_PERETZ = 1
SOLID_CONDUCTION = 2
_SOLID_MODES = {'goodman': SOLID_GOODMAN, 'peretz': SOLID_PERETZ,
                'conduction': SOLID_CONDUCTION}


def _solid_model_code(mode):
    try:
        return _SOLID_MODES[mode]
    except (KeyError, TypeError):
        raise ValueError(
            f"unknown solid_heating_model: {mode!r} "
            f"(expected one of {sorted(_SOLID_MODES)})") from None


@njit(cache=True, fastmath=True)
def _peretz_heat_transfer_coeff(
    P, u_cell, T_gas, T_wall, R_gas, Cp_gas, gamma_gas,
    D_hyd, x_from_head, A_port, A_throat,
):
    """Peretz et al. (1973) bore heat-transfer coefficient [W/(m^2 K)].

    Eq. III-22:  Nu = a6 (x_r/d_h)^-0.1 Re^0.8 Pr^0.4 (A_p/A_t)^0.4,
    a6 = 2.56e-2 (Peretz's value; his h-form, Eq. III-31, carries the same
    constant: (0.87e-6)^0.2 * a6 = 1.56e-3). All gas properties except the
    velocity are evaluated at the film temperature T_af = (T_gas+T_wall)/2:
      rho_af = P/(R T_af)                                     (Eq. III-26)
      mu_af  = 0.87e-6 W^0.5 T_af^0.65 poise, W = molar mass  (Eq. III-24)
      Pr     = gamma/(1.77 gamma - 0.45)                      (Eq. III-25)
      k_af   = Cp mu_af / Pr
    x_r is the distance from the igniter-jet reattachment point (Peretz:
    ~one port diameter downstream of the igniter exit). Here
    x_r = max(x_from_head, d_h): a 0.1-power factor, so the choice moves h
    by <~ 10 %. Nu is floored at 3.66 (fully developed laminar limit, the
    same floor the default Gnielinski chain carries); Peretz's correlation
    is a turbulent one (his Table 4: Re_d 1.1e4-2.1e4 during induction, up
    to 4.5e5 during chamber filling; A_p/A_t 1.06-2.0).
    """
    if D_hyd <= 1e-10 or A_throat <= 0.0 or A_port <= 0.0 or R_gas <= 0.0:
        return 0.0
    T_af = 0.5 * (T_gas + T_wall)
    if T_af <= 0.0 or P <= 0.0:
        return 0.0
    W_gmol = 8314.462618 / R_gas
    mu = 0.87e-7 * W_gmol ** 0.5 * T_af ** 0.65
    rho_af = P / (R_gas * T_af)
    Re = rho_af * abs(u_cell) * D_hyd / mu
    Pr = gamma_gas / (1.77 * gamma_gas - 0.45)
    k_af = Cp_gas * mu / Pr
    x_r = x_from_head if x_from_head > D_hyd else D_hyd
    Nu = (2.56e-2 * (x_r / D_hyd) ** (-0.1) * Re ** 0.8 * Pr ** 0.4
          * (A_port / A_throat) ** 0.4)
    if Nu < 3.66:
        Nu = 3.66
    return Nu * k_af / D_hyd


# Keller, Baer & Ryan test-section hydraulic diameter: 0.500 x 0.125 in
# channel, D_h = 4A/P = 0.200 in.
_KELLER_D_H = 5.08e-3


@njit(cache=True, fastmath=True)
def _keller_heat_transfer_coeff(rho, u_cell, T_gas, D_hyd):
    """Keller-Baer-Ryan (1966) Eq. 11 [W/(m^2 K)], air/O2/N2 shock-tube channel.

    h = 1.435e-4 Tg^0.3 G^0.905 cal/(cm^2 s K), Tg [K], G [g/(cm^2 s)] (6 %
    average deviation; Tg 1000-2600 K, 14-25 atm, 50-800 m/s, in a 5.08 mm
    hydraulic-diameter channel; G itself is not tabulated, ~9-680 g/cm^2 s).
    Extrapolation to the bore: with Nu ~ Re^0.905 at fixed properties,
    h ~ G^0.905 D^-0.095, so h is scaled by (D_h/D_Keller)^-0.095. The
    property difference between Keller's N2/air and the bore gas is NOT
    corrected (the paper gives no property data); Tg is the local static
    gas temperature rather than Keller's stagnation temperature.
    """
    if D_hyd <= 1e-10 or T_gas <= 0.0:
        return 0.0
    G = rho * abs(u_cell) * 0.1          # kg/(m^2 s) -> g/(cm^2 s)
    if G <= 0.0:
        return 0.0
    h = 1.435e-4 * T_gas ** 0.3 * G ** 0.905 * 41840.0
    return h * (D_hyd / _KELLER_D_H) ** (-0.095)


@njit(cache=True, fastmath=True)
def _orifice_exit_velocity(P_ig, T_ig, P_main, gamma, M):
    """Ideal-gas pyrogen orifice exit velocity from plenum state."""
    if P_ig <= 0.0 or T_ig <= 0.0 or P_main >= P_ig:
        return 0.0
    if gamma <= 1.0 or M <= 0.0:
        return 0.0

    pressure_ratio = P_main / P_ig
    if pressure_ratio < 0.0:
        pressure_ratio = 0.0

    crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    if pressure_ratio < crit:
        T_exit = T_ig * 2.0 / (gamma + 1.0)
    else:
        T_exit = T_ig * pressure_ratio ** ((gamma - 1.0) / gamma)

    R_specific_ig = R_UNIVERSAL / M
    v2 = 2.0 * gamma / (gamma - 1.0) * R_specific_ig * (T_ig - T_exit)
    if v2 <= 0.0:
        return 0.0
    return v2 ** 0.5


@njit(cache=True, fastmath=True)
def _cal_cm2_s_to_w_m2(heat_flux_cal_cm2_s):
    """Convert cal/(cm^2*s) to W/m^2."""
    return heat_flux_cal_cm2_s * CAL_CM2_S_TO_W_M2


@njit(cache=True, fastmath=True)
def _pyrogen_surface_heat_power(
    mdot_igniter, T_ig, T_surf, C_burn, dx, Cp_pyrogen,
    measured_heat_flux_w_m2,
):
    """
    Delivered pyrogen surface heating, capped by available sensible power.

    Returns ``(power_w, heat_flux_w_m2)`` for the target cell.

    v0.7.1 (Phase 3.5): the gas carrying enthalpy from the plenum to the
    propellant surface IS pyrogen products. The sensible-power cap uses
    the pyrogen species's Cp, not the cell mixture's Cp and not the
    propellant Cp_gas (the v0.7.0 / Phase 3 placeholder). For BPNV that
    ratio is ~1385/2060 ≈ 0.67, so the sensible cap binds ~33% lower
    than under v0.7.0 — physically correct, and a recognized source of
    ignition-transient drift vs the v0.7.0 baseline (will be absorbed
    by Phase 5 re-LHS).
    """
    if (mdot_igniter <= 0.0 or T_ig <= T_surf or C_burn <= 0.0 or
            dx <= 0.0 or Cp_pyrogen <= 0.0 or measured_heat_flux_w_m2 <= 0.0):
        return 0.0, 0.0

    contact_area = C_burn * dx
    if contact_area <= 1.0e-16:
        return 0.0, 0.0

    measured_power = measured_heat_flux_w_m2 * contact_area
    sensible_power = mdot_igniter * Cp_pyrogen * (T_ig - T_surf)
    if sensible_power <= 0.0:
        return 0.0, 0.0

    delivered_power = measured_power
    if sensible_power < delivered_power:
        delivered_power = sensible_power
    if delivered_power <= 0.0:
        return 0.0, 0.0
    return delivered_power, delivered_power / contact_area


@njit(cache=True, fastmath=True)
def _pyrogen_surface_thermal_sink(
    surface_heat_power_w, dx, pyrogen_enthalpy_source_w_per_m,
):
    """Enthalpy-source sink matching solid heating power [W/m].

    v0.7.1 (Phase 3): ``thermal_source`` carries W/m (enthalpy injection
    per unit length), so the sink is ``surface_heat_power_w / dx``. The
    output is clamped to the available pyrogen enthalpy injection at this
    cell so we never extract more enthalpy than the pyrogen just added.
    """
    if surface_heat_power_w <= 0.0 or dx <= 0.0:
        return 0.0
    sink = surface_heat_power_w / dx
    if pyrogen_enthalpy_source_w_per_m <= 0.0:
        return 0.0
    if sink > pyrogen_enthalpy_source_w_per_m:
        return pyrogen_enthalpy_source_w_per_m
    return sink


@njit(cache=True, fastmath=True)
def _gas_sensible_energy(rho, T, A_port, dx, Cp_arr, N):
    """Discrete gas sensible energy used by diagnostics [J].

    v0.7.1 (Phase 3): ``Cp_arr`` is per-cell so cells with different
    mixtures contribute their own Cp to the total.
    """
    total = 0.0
    for i in range(N):
        total += rho[i] * A_port[i] * dx * Cp_arr[i] * T[i]
    return total


@njit(cache=True, fastmath=True)
def _gas_mass_energy(rho, T, A_port, dx, Cp_arr, N):
    """Gas mass [kg] and sensible energy [J] of a state (energy audit, T6).

    Read-only bookkeeping helper: ``(sum rho*A*dx, sum rho*A*dx*Cp*T)``.
    Used to evaluate the storage before/after the operations that act on the
    carried state between PISO steps (mixture refresh, geometry refresh).
    """
    mass = 0.0
    energy = 0.0
    for i in range(N):
        m_i = rho[i] * A_port[i] * dx
        mass += m_i
        energy += m_i * Cp_arr[i] * T[i]
    return mass, energy


@njit(cache=True, fastmath=True)
def _thermal_source_power(thermal_source, dx, N):
    """Sum solver enthalpy-source units to thermal power [W].

    v0.7.1 (Phase 3): ``thermal_source`` is already in W/m units (enthalpy
    injection per unit length), so the sum is ``thermal_source[i] * dx``.
    """
    total = 0.0
    for i in range(N):
        total += thermal_source[i] * dx
    return total


@njit(cache=True, fastmath=True)
def _compute_mixture_cell(Y_row, species_params):
    """
    Ideal-gas mass-fraction mixing for one cell.

    Standard textbook formulas (e.g. Kuo, *Principles of Combustion* §1.6;
    Bird/Stewart/Lightfoot §16.1):

        Cp_mix     = sum_s  Y[s] * Cp[s]                        (mass-weighted)
        1/M_mix    = sum_s  Y[s] / M[s]                         (harmonic)
        R_mix      = R_universal / M_mix
        gamma_mix  = Cp_mix / (Cp_mix - R_mix)                  (ideal gas)

    Parameters
    ----------
    Y_row : np.ndarray[S]
        Mass fractions in this cell. Sum must be ~1; the caller is
        responsible for renormalization. Each element is in [0, 1].
    species_params : np.ndarray[S, 4]
        Row layout: (gamma, Cp, molecular_weight, T_flame). gamma is not
        used here (it derives from Cp_mix and R_mix); T_flame is not
        used (it is a source-injection property, not a bulk-mixture
        property).

    Returns
    -------
    (gamma_mix, Cp_mix, R_mix, M_mix) : 4 floats
    """
    S = species_params.shape[0]
    Cp_mix = 0.0
    inv_M_mix = 0.0
    for s in range(S):
        y = Y_row[s]
        if y <= 0.0:
            continue
        Cp_s = species_params[s, 1]
        M_s = species_params[s, 2]
        Cp_mix += y * Cp_s
        inv_M_mix += y / M_s
    if inv_M_mix <= 0.0 or Cp_mix <= 0.0:
        # Degenerate: fall back to species 0 to keep the kernel total.
        gamma0 = species_params[0, 0]
        Cp0 = species_params[0, 1]
        M0 = species_params[0, 2]
        R0 = R_UNIVERSAL / M0
        return gamma0, Cp0, R0, M0
    M_mix = 1.0 / inv_M_mix
    R_mix = R_UNIVERSAL / M_mix
    denom = Cp_mix - R_mix
    if denom <= 0.0:
        # Physically pathological: Cp <= R implies gamma -> infinity.
        # Clamp to species-0 thermo. (Should not occur for real fuel
        # combustion products + air.)
        gamma0 = species_params[0, 0]
        Cp0 = species_params[0, 1]
        M0 = species_params[0, 2]
        R0 = R_UNIVERSAL / M0
        return gamma0, Cp0, R0, M0
    gamma_mix = Cp_mix / denom
    return gamma_mix, Cp_mix, R_mix, M_mix


@njit(cache=True, fastmath=True)
def _refresh_mixture_arrays(
    Y, species_params, gamma_arr, Cp_arr, R_arr, M_arr, N,
):
    """
    Refresh per-cell (gamma, Cp, R, M) arrays from Y[N, S].

    Called every time-step after _advect_species updates Y. The output
    arrays are then available for the solver to consume (Phase 3 will
    wire them into _piso_step).

    For Phase 2 this is a "ready but not yet consumed" hook: it runs
    every step so we can confirm zero-overhead correctness, but the
    arrays sit unused by the solver until Phase 3 lands.
    """
    for i in range(N):
        gamma_i, Cp_i, R_i, M_i = _compute_mixture_cell(Y[i, :], species_params)
        gamma_arr[i] = gamma_i
        Cp_arr[i] = Cp_i
        R_arr[i] = R_i
        M_arr[i] = M_i


@njit(cache=True, fastmath=True)
def _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N):
    """v0.7.2 Phase A — exponential-decay axial weights for pyrogen injection.

    Returns per-cell weights summing to exactly 1.0:

        w[i] = exp(-x_centers[i] / L_jet) * dx[i] / sum_j(exp(-x_centers[j] / L_jet) * dx[j])

    Pyrogen mass / enthalpy / momentum are split across cells by `w[i]`
    so the energy deposits along the bore over the plume's physical
    reach rather than concentrating in cell 0. ``L_jet`` is the
    characteristic decay length set by ``kappa_jet * d_throat`` at the
    Python boundary.

    Edge cases:
    - ``L_jet <= 0``: all weight goes to cell 0 (recovers v0.7.1
      behavior byte-for-byte). Used as the regression gate.
    - Pathological total < 0: defensive fallback to cell-0-only.

    Conservation guarantee: ``sum(w) == 1`` to within floating-point
    precision; verified in tests/test_pyrogen_axial_weights.py.

    See ``srm_1d/docs/v0_7_2/candidates/03_pyrogen_spatial_distribution.md``.
    """
    w = np.zeros(N)
    if L_jet <= 0.0:
        w[0] = 1.0
        return w
    total = 0.0
    for i in range(N):
        w[i] = np.exp(-x_centers[i] / L_jet) * dx[i]
        total += w[i]
    if total <= 0.0:
        w[:] = 0.0
        w[0] = 1.0
        return w
    inv_total = 1.0 / total
    for i in range(N):
        w[i] *= inv_total
    return w


@njit(cache=True, fastmath=True)
def _compute_uncontained_pyrogen_mdot_cells(
    P_bore, a, n, rho_p, A_cell0, m_cell, m_cell0, t, dt,
    i_start, i_end, N, mdot_arr, gas_mass_fraction, P_plateau, t_light,
    shape_LD=1.0,
):
    """OB1 (opt-in) — per-cell uncontained pyrogen mdot with the burning
    area following the pellet web consumption.

    Each cartridge cell i keeps its own unburnt mass ``m_cell[i]`` (initially
    ``m_cell0``) and burns it at its host-cell bore pressure:

        A_i     = A_cell0 * g(m_cell[i] / m_cell0; L/D)     [m^2]
        r_b[i]  = a * min(P_bore[i], P_plateau)^n           [m/s]  (cap if P_plateau > 0)
        mdot[i] = rho_p * r_b[i] * A_i                      [kg/s, SOLID]

    ``g`` is the exact right-cylinder form function
    (``igniter_plenum._cyl_area_fraction``: all faces regress at the same
    rate, D = D0-2s, H = H0-2s); for L/D = 1 (spheres, cubes, L = D cylinders)
    it is exactly ``f^(2/3)``. A cell that
    has not yet lit (``t < t_light[i]``) is inert. The last step of a cell is
    capped so it consumes exactly its remaining mass. Returns the total unburnt
    mass after the step (the caller stores it in ``plenum_state[0]``).
    ``mdot_arr`` receives the GAS mass rate (solid rate * gas_mass_fraction),
    exactly as in the pooled function.
    """
    for i in range(N):
        mdot_arr[i] = 0.0
    lo = i_start if i_start >= 0 else 0
    hi = i_end if i_end < N else N - 1
    total_remaining = 0.0
    if lo > hi or dt <= 0.0 or A_cell0 <= 0.0 or m_cell0 <= 0.0:
        for i in range(N):
            total_remaining += m_cell[i]
        return total_remaining
    for i in range(lo, hi + 1):
        mi = m_cell[i]
        if mi <= 0.0 or t < t_light[i]:
            total_remaining += mi if mi > 0.0 else 0.0
            continue
        P = P_bore[i]
        if P < 0.0:
            P = 0.0
        if P_plateau > 0.0 and P > P_plateau:
            P = P_plateau
        frac = mi / m_cell0
        A_i = A_cell0 * _cyl_area_fraction(frac, shape_LD)
        mdot_solid = rho_p * a * (P ** n) * A_i
        if mdot_solid * dt >= mi:
            mdot_solid = mi / dt
            m_cell[i] = 0.0
        else:
            m_cell[i] = mi - mdot_solid * dt
        total_remaining += m_cell[i]
        mdot_arr[i] = mdot_solid * gas_mass_fraction
    return total_remaining


@njit(cache=True, fastmath=True)
def _compute_uncontained_pyrogen_mdot(
    P_bore, a, n, rho_p, A_burn_per_cell, m_pyrogen_remaining, dt,
    i_start, i_end, N, mdot_arr, gas_mass_fraction=1.0, P_plateau=0.0,
):
    """v0.7.3 Phase A — per-cell pyrogen mdot for uncontained
    submerged-igniter topologies (head_basket, aft_basket).

    Each pyrogen pellet (cell in [i_start, i_end] inclusive) burns at
    its host cell's LOCAL bore pressure — no plenum chamber, no
    orifice, no separate P_ig. The physics is:

        r_b[i] = a * max(P_bore[i], 0)^n        [m/s]
        mdot[i] = rho_p * r_b[i] * A_burn_per_cell   [kg/s]

    where ``A_burn_per_cell = chamber.A_burn_initial / n_cartridge_cells``
    is the initial burning surface area distributed uniformly across
    the cartridge's axial extent.

    Mass conservation: if ``sum(mdot[i]) * dt`` would deplete more
    pyrogen than ``m_pyrogen_remaining``, all per-cell mdots are
    uniformly scaled down so the last burn step consumes exactly the
    remaining pyrogen. After depletion, subsequent steps return
    all-zero mdot.

    Returns the updated ``m_pyrogen_remaining`` (caller persists it
    in plenum_state[0] across steps; this mirrors the existing
    forward_plenum convention even though the uncontained model has
    no plenum-gas state per se).

    Edge cases:
    - ``m_pyrogen_remaining <= 0``: all-zero mdot, returns 0.0
    - ``A_burn_per_cell <= 0`` or ``dt <= 0``: all-zero mdot,
      returns ``m_pyrogen_remaining`` unchanged
    - ``i_start > i_end`` or fully out-of-range: same as above
    - ``P_bore[i] < 0`` (numerical artifact): clamped to 0

    See srm_1d/docs/v0_7_2/candidates_post_phaseA.md (v0.7.3 design)
    and the PyrogenChamber docstring in igniter_plenum.py for the
    full architecture (uncontained vs plenum split).
    """
    # Zero the full array each step (cheap O(N))
    for i in range(N):
        mdot_arr[i] = 0.0

    # Defensive early exits
    if m_pyrogen_remaining <= 0.0:
        return 0.0
    if A_burn_per_cell <= 0.0 or dt <= 0.0:
        return m_pyrogen_remaining
    if i_start > i_end:
        return m_pyrogen_remaining

    # Clamp range to [0, N-1]
    lo = i_start if i_start >= 0 else 0
    hi = i_end if i_end < N else N - 1
    if lo > hi:
        return m_pyrogen_remaining

    # Compute provisional per-cell mdot
    total_mdot = 0.0
    for i in range(lo, hi + 1):
        P = P_bore[i]
        if P < 0.0:
            P = 0.0
        if P_plateau > 0.0 and P > P_plateau:
            P = P_plateau
        r_b = a * (P ** n)
        m = rho_p * r_b * A_burn_per_cell
        mdot_arr[i] = m
        total_mdot += m

    # Mass conservation: cap SOLID consumption at m_pyrogen_remaining.
    # mdot_arr currently holds the SOLID burn rate (sets pyrogen depletion).
    if total_mdot * dt > m_pyrogen_remaining:
        if total_mdot > 0.0:
            scale = m_pyrogen_remaining / (total_mdot * dt)
            for i in range(lo, hi + 1):
                mdot_arr[i] *= scale
            new_remaining = 0.0
        else:
            new_remaining = m_pyrogen_remaining
    else:
        new_remaining = m_pyrogen_remaining - total_mdot * dt

    # Condensed-phase split: the solid depletes at the full burn rate (above),
    # but the bore receives only gas_mass_fraction of it as pressurizing gas
    # (the rest is condensed B/BN/Al2O3). Mirrors the forward_plenum split.
    if gas_mass_fraction != 1.0:
        for i in range(lo, hi + 1):
            mdot_arr[i] *= gas_mass_fraction
    return new_remaining


@njit(cache=True, fastmath=True)
def _compute_uniform_band_weights(dx, i_start, i_end, N, w):
    """v0.7.3 Phase A — uniform top-hat axial weights for submerged-igniter
    pyrogen injection.

    Fills ``w[i]`` with mass-conservative weights summing to 1.0:

        w[i] = dx[i] / sum(dx[j] for j in [i_start, i_end])  if i_start <= i <= i_end
        w[i] = 0                                              otherwise

    This is the load-bearing kernel for submerged-igniter topologies
    (head-end basket, aft-cavity cartridge, etc.). It distributes
    pyrogen mass / enthalpy / species uniformly across the cells the
    igniter physically occupies, in the same way Phase A's exponential
    decay distributes forward-plenum injection across the head-end
    impingement region.

    By design, this kernel unifies four topologies:
    - Head-end basket (4a): caller passes i_start=0, i_end = i_cartridge_end
    - Aft-cavity zero-axial cartridge (4b zero-mom): caller passes
      i_start = N - n_cartridge_cells, i_end = N - 1, and skips
      momentum injection entirely
    - Aft-cavity fore-firing cartridge (4b fore-firing): same i_start /
      i_end as zero-axial, but caller adds upstream-directed momentum
      at face i_start
    - Mid-bore submerged basket: caller passes arbitrary i_start, i_end

    For the existing forward-plenum default (Phase A pre-existing
    behavior), callers should NOT use this kernel — they use
    _compute_pyrogen_axial_weights with L_jet = kappa_jet * d_throat
    for the exponential-decay shape.

    Edge cases:
    - i_start > i_end (invalid): all weight goes to cell 0 (defensive
      fallback to match _compute_pyrogen_axial_weights's L_jet<=0
      convention).
    - i_start < 0 or i_end >= N: clamped to [0, N-1] internally.
    - All-zero dx in the range: defensive fallback to cell 0.

    Conservation: sum(w) == 1 to floating-point precision; verified
    in tests/test_uniform_band_weights.py.

    See srm_1d/docs/v0_7_2/candidates_post_phaseA.md (v0.7.3 design).
    """
    for i in range(N):
        w[i] = 0.0
    # Clamp range
    lo = i_start
    if lo < 0:
        lo = 0
    hi = i_end
    if hi >= N:
        hi = N - 1
    if lo > hi:
        w[0] = 1.0
        return
    total = 0.0
    for i in range(lo, hi + 1):
        total += dx[i]
    if total <= 0.0:
        w[0] = 1.0
        return
    inv_total = 1.0 / total
    for i in range(lo, hi + 1):
        w[i] = dx[i] * inv_total


# v0.7.3 Phase B.4: heat-delivery mode codes for uncontained pyrogen
# topologies. Mutually exclusive at the implementation level to avoid
# double-counting DeMar's already-included radiative component.
HEAT_DELIVERY_NONE      = 0  # No pyrogen surface heat flux applied
HEAT_DELIVERY_DEMAR     = 1  # Pyrogen.heat_flux_cal_cm2_s, distributed across cartridge cells
HEAT_DELIVERY_RADIATION = 2  # Stefan-Boltzmann pellet emission + view factor + absorption


def _heat_delivery_code(mode):
    if mode == 'none':
        return HEAT_DELIVERY_NONE
    if mode == 'demar':
        return HEAT_DELIVERY_DEMAR
    if mode == 'radiation':
        return HEAT_DELIVERY_RADIATION
    raise ValueError(f"unknown heat_delivery_mode: {mode!r}")


@njit(cache=True, fastmath=True)
def _compute_pyrogen_heat_flux_arr(
    topology_code, heat_delivery_mode_code,
    # Common
    is_grain, has_ignited, T_surf, C_burn, mdot_igniter, dx, N,
    pyrogen_heat_flux_arr,
    # Forward plenum
    head_grain_cell, fwd_plenum_flux_w_m2,
    # Uncontained DEMAR
    mdot_uncontained_arr, demar_flux_w_m2,
    Cp_pyrogen, T_ig,
    cart_i_start, cart_i_end,
    # Uncontained RADIATION
    T_flame_pyrogen, pellet_emissivity, radiation_absorption_length_m,
    x_centers, A_port,
):
    """v0.7.3 Phase B.4 — fill per-cell pyrogen surface heat-flux array.

    Three mutually-exclusive heat-delivery modes for uncontained
    topologies, plus the existing forward-plenum DeMar path (topology
    code 0). The output ``pyrogen_heat_flux_arr[i]`` is the per-cell
    surface heat flux [W/m²] applied to the propellant surface via the
    Goodman ignition kernel.

    Mode semantics:
    - ``HEAT_DELIVERY_NONE``: array stays at 0. Recovers v0.7.3-phaseA
      uncontained behavior (no pyrogen surface heat delivery).
    - ``HEAT_DELIVERY_DEMAR``: DeMar 2021 time-averaged flux applied
      uniformly across cartridge cells with grain. Per-cell sensible
      cap: flux_capped = min(demar_flux,
                              mdot_local · Cp_pyrogen · (T_ig - T_surf)
                              / (C_burn · dx))
      where mdot_local is the per-cell pyrogen mdot from
      _compute_uncontained_pyrogen_mdot. Empirically grounded for
      pellet pyrogens (BKNO3, MTV).
    - ``HEAT_DELIVERY_RADIATION``: Stefan-Boltzmann pellet emission
      from cartridge cells to all unignited grain cells via geometric
      view factor + exponential absorption-length attenuation:
          q[j] = sum over e in cartridge cells (
              σ · ε · (T_flame^4 - T_surf[j]^4) · F_ij · exp(-d/L)
          )
      where F_ij = A_port[j] / (4π·d² + A_port[j]) saturates to 1 for
      adjacent cells and falls as ~1/d² far field. Emission ceases per
      cell when its pyrogen mass is depleted (mdot_uncontained[e]==0).
      Physically modeled; extensible to powder and chunks via T_flame.

    For forward_plenum (topology_code == 0): cell-0 dispatch is
    mode-aware (v0.7.3.3, building on v0.7.3-phaseB's mode-dispatch
    architecture for uncontained topologies):
    - DEMAR (default): DeMar 2021 lumped flux with sensible cap based
      on total plenum vent mdot. Preserves v0.7.0-v0.7.3.2 behavior.
    - RADIATION: Stefan-Boltzmann emission at cell-0 from the
      pyrogen plume,
        q = σ · pellet_emissivity · (T_flame_pyrogen^4 - T_surf[0]^4)
      No sensible cap (radiative emission is energy-balanced by
      Stefan-Boltzmann, not by the mass-advection enthalpy budget).
      View factor and gas absorption assumed unity at the cell-0
      impingement zone (gas-path length << radiation_absorption_length_m
      for forward_plenum's compact plume).
    - NONE: no surface flux; useful for diagnostic isolation runs.

    Mutually exclusive at the implementation level so DeMar's
    already-radiative contribution doesn't double-count.

    See srm_1d/docs/v0_7_3/PHASE_B_SCOPE.md §B.4 for the
    double-counting analysis between 'demar' and 'radiation' modes.
    """
    for i in range(N):
        pyrogen_heat_flux_arr[i] = 0.0

    # Forward plenum: mode-aware cell-0 dispatch.
    if topology_code == 0:
        if not (head_grain_cell >= 0 and mdot_igniter > 0.0):
            return
        if heat_delivery_mode_code == HEAT_DELIVERY_NONE:
            return
        if heat_delivery_mode_code == HEAT_DELIVERY_DEMAR:
            if fwd_plenum_flux_w_m2 <= 0.0:
                return
            contact_area = C_burn[head_grain_cell] * dx
            if contact_area <= 1.0e-16:
                return
            dT = T_ig - T_surf[head_grain_cell]
            if dT <= 0.0:
                return
            measured_power = fwd_plenum_flux_w_m2 * contact_area
            sensible_power = mdot_igniter * Cp_pyrogen * dT
            delivered_power = (sensible_power if sensible_power < measured_power
                               else measured_power)
            if delivered_power > 0.0:
                pyrogen_heat_flux_arr[head_grain_cell] = (
                    delivered_power / contact_area
                )
            return
        if heat_delivery_mode_code == HEAT_DELIVERY_RADIATION:
            if pellet_emissivity <= 0.0:
                return
            T_em4 = T_flame_pyrogen ** 4
            T_rec4 = T_surf[head_grain_cell] ** 4
            if T_em4 <= T_rec4:
                return
            # σ · ε · (T_em^4 - T_rec^4); F_view ≈ 1 and gas
            # absorption negligible for the forward_plenum compact
            # plume → cell-0 impingement zone.
            q = 5.670374419e-8 * pellet_emissivity * (T_em4 - T_rec4)
            pyrogen_heat_flux_arr[head_grain_cell] = q
            return
        # Unknown mode → no flux (defensive).
        return

    # Uncontained: branch on heat delivery mode
    if heat_delivery_mode_code == HEAT_DELIVERY_NONE:
        return

    if heat_delivery_mode_code == HEAT_DELIVERY_DEMAR:
        if demar_flux_w_m2 <= 0.0:
            return
        lo = cart_i_start if cart_i_start >= 0 else 0
        hi = cart_i_end if cart_i_end < N else N - 1
        for i in range(lo, hi + 1):
            if not is_grain[i] or has_ignited[i]:
                continue
            if mdot_uncontained_arr[i] <= 0.0:
                continue
            contact_area = C_burn[i] * dx
            if contact_area <= 1.0e-16:
                continue
            dT = T_ig - T_surf[i]
            if dT <= 0.0:
                continue
            # Per-cell sensible cap using cell-local pyrogen mdot
            measured_power = demar_flux_w_m2 * contact_area
            sensible_power = mdot_uncontained_arr[i] * Cp_pyrogen * dT
            delivered_power = (sensible_power if sensible_power < measured_power
                               else measured_power)
            if delivered_power > 0.0:
                pyrogen_heat_flux_arr[i] = delivered_power / contact_area
        return

    if heat_delivery_mode_code == HEAT_DELIVERY_RADIATION:
        if pellet_emissivity <= 0.0:
            return
        lo = cart_i_start if cart_i_start >= 0 else 0
        hi = cart_i_end if cart_i_end < N else N - 1
        T_em4 = T_flame_pyrogen ** 4
        # Receivers: ALL unignited grain cells (including those inside
        # the cartridge range — pellets radiate to their own host cell's
        # propellant surface too with F_view = 1 since distance = 0).
        for j in range(N):
            if not is_grain[j] or has_ignited[j]:
                continue
            if A_port[j] <= 0.0:
                continue
            T_rec4 = T_surf[j] ** 4
            if T_em4 <= T_rec4:
                continue
            q_total = 0.0
            for e in range(lo, hi + 1):
                if mdot_uncontained_arr[e] <= 0.0:
                    continue
                dxc = x_centers[j] - x_centers[e]
                d_sq = dxc * dxc
                F_view = A_port[j] / (4.0 * np.pi * d_sq + A_port[j])
                if radiation_absorption_length_m > 0.0:
                    d_abs = np.sqrt(d_sq)
                    atten = np.exp(-d_abs / radiation_absorption_length_m)
                else:
                    atten = 1.0
                q_em = (5.670374419e-8 * pellet_emissivity *
                        (T_em4 - T_rec4))
                q_total += q_em * F_view * atten
            pyrogen_heat_flux_arr[j] = q_total


@njit(cache=True, fastmath=True)
def _compute_flame_front_augment(
    is_burning, has_ignited, ignition_time, t,
    tau_window, augment_value, N, augment_arr,
):
    """v0.7.2 Phase B-v2 — flame-front-marker h_c augmentation.

    Reformulation of the Phase B Kashiwagi/Han cumulative-G approach
    (commit 065d193) which was found to double-count with PISO's local-Re
    tracking and amplify the spike rather than smooth it. The new
    approach is strictly sequential: an unignited cell receives the
    augmentation ONLY if its immediate upstream neighbor has ignited
    within the last ``tau_window`` seconds.

    Physics: each ignited cell creates a transient flame-front that
    impinges on the next unignited cell, locally enhancing convective
    heat transfer for a brief window. Once that next cell ignites, the
    front passes to the following cell, and so on. Cells far from any
    recent ignition see no augmentation (default Bartz h_c).

    Algorithm:
        for i in 0..N-1:
            augment_arr[i] = 1.0   # default no augmentation
        for j in 0..N-2:
            if is_burning[j] AND (t - ignition_time[j]) < tau_window
                    AND NOT has_ignited[j + 1]:
                augment_arr[j + 1] = augment_value

    Strict-sequential property: cell 0 receives no augmentation (no
    upstream cell). Cell 1 receives augmentation only after cell 0
    ignites, within tau_window. Cell 2 receives augmentation only
    after cell 1 ignites, within tau_window. And so on.

    Knob defaults (see Propellant.flame_spread_tau / flame_spread_boost):
    - tau_window ~ 1 ms: time scale for one cell to ignite under
      enhanced h_c (calibrated to give physical flame-spread velocities
      of 10-100 m/s on typical motor grids).
    - augment_value ~ 3.0: peak h_c multiplier when boost is active.

    See srm_1d/docs/v0_7_2/candidates/02_spatial_ignition_front_coupling.md.
    """
    for i in range(N):
        augment_arr[i] = 1.0
    for j in range(N - 1):
        if not is_burning[j]:
            continue
        if (t - ignition_time[j]) >= tau_window:
            continue
        if has_ignited[j + 1]:
            continue
        augment_arr[j + 1] = augment_value


@njit(cache=True, fastmath=True)
def _advance_flame_front(
    is_burning, x_centers, flame_front_velocity,
    x_front, front_direction, front_seed_idx,
    cart_i_start, cart_i_end, dt, N, ignitable,
):
    """v0.7.4 Phase F — advance the flame-spread front at a bounded
    literature velocity, fill the ``ignitable`` exposure mask, and return
    the new ``x_front``.

    The front advances at ``flame_front_velocity`` [m/s] — a single
    physical constant (AP/HTPB lateral flame spread is ~1-10 m/s; Peretz-
    Kuo-Caveny-Summerfield 1973; Kumar & Kuo 1984), held across motors and
    decoupled from the acoustic/fill speed (~300-1000 m/s) that otherwise
    over-speeds ignition. The earlier per-step ``q''/(rho*Cps*dT)`` form
    was a burn/regression velocity (~mm/s) — the wrong quantity — and is
    abandoned; see docs/v0_7_4/.

    Once ANY cell is burning (a front exists), ``x_front`` advances by
    ``front_direction * flame_front_velocity * dt`` each step. A grain
    cell is ``ignitable`` if (a) it lies in the cartridge/igniter range
    ``[cart_i_start, cart_i_end]`` (induction-privileged — lit directly by
    the igniter), (b) it is the seed cell, or (c) the front has reached its
    center (``x_centers[i] <= x_front`` fore→aft / ``>=`` aft→fore). The
    caller withholds ALL surface heating (convective + pyrogen + radiation)
    from cells that are not ignitable, so ignition follows the front rather
    than the bulk-fill / broad pyrogen radiation.

    ``flame_front_velocity`` is a velocity, so the physical spread speed is
    grid-independent. Before any cell burns the front sits at the seed and
    only the cartridge region is exposed (induction).
    """
    any_burning = False
    for i in range(N):
        if is_burning[i]:
            any_burning = True
            break
    if any_burning and flame_front_velocity > 0.0:
        x_front = x_front + front_direction * flame_front_velocity * dt

    lo = cart_i_start if cart_i_start >= 0 else 0
    hi = cart_i_end if cart_i_end < N else N - 1
    for i in range(N):
        if i == front_seed_idx or (lo <= i <= hi):
            ignitable[i] = True
        elif front_direction > 0:
            ignitable[i] = x_centers[i] <= x_front
        else:
            ignitable[i] = x_centers[i] >= x_front
    return x_front


@njit(cache=True, fastmath=True)
def _advance_zn_burn_rate(r_dyn, r_qs, is_burning, alpha_solid,
                          kappa_zn, r_floor, dt, N):
    """v0.7.4 Phase Z — lumped Zeldovich-Novozhilov dynamic burn-rate
    relaxation. Each burning cell's rate relaxes toward the quasi-steady
    Ma value ``r_qs`` over the condensed-phase thermal-wave time:

        tau   = kappa_zn * alpha_solid / max(r_dyn, r_floor)^2
        r_dyn = r_qs + (r_dyn - r_qs) * exp(-dt / tau)

    The analytic exponential update is unconditionally stable and never
    overshoots. ``tau ~ 1/r^2`` self-attenuates as the burn rate climbs, so
    the high-r plateau is preserved while the low-r ignition transient
    (tau ~ ms) is smoothed. Greatrix 2008 lag form (Eq. 14); timescale per
    Lengelle 1975 / Greatrix Eq. 19. ``r_qs`` is the existing Ma-2020 total
    — Z-N lags it, it does NOT replace the erosive model.

    First-touch: a just-ignited cell (``r_dyn == 0``) seeds directly to
    ``r_qs`` (no ramp-from-zero artifact). Cells that are not burning reset
    to 0 so a re-ignition seeds cleanly.
    """
    for i in range(N):
        if not is_burning[i]:
            r_dyn[i] = 0.0
            continue
        rq = r_qs[i]
        if r_dyn[i] <= 0.0 or alpha_solid <= 0.0 or kappa_zn <= 0.0:
            r_dyn[i] = rq
            continue
        r_eff = r_dyn[i]
        if r_eff < r_floor:
            r_eff = r_floor
        tau = kappa_zn * alpha_solid / (r_eff * r_eff)
        if tau <= 0.0:
            r_dyn[i] = rq
            continue
        r_dyn[i] = rq + (r_dyn[i] - rq) * np.exp(-dt / tau)


@njit(cache=True, fastmath=True)
def _compute_T_ceiling_arr(
    Y, species_params, T_ceiling_arr, N, T_initial_gas,
    Y_min=0.05,
):
    """Refresh per-cell temperature ceiling from species mass fractions.

    DESIGN §5 specifies a per-cell ceiling tied to the active species
    in each cell:

        T_ceiling[i] = max(T_flame[s] for s with Y[i, s] > Y_min) * 1.01

    v0.7.1 (Phase 3 → strict-form follow-up, 2026-05-23): this is the
    strict DESIGN §5 form WITH an initial-condition guard. v0.7.1's
    documented IC (DESIGN §3) seeds ``T = T_flame_propellant`` while
    ``Y[:, ambient] = 1.0`` for v0.7.0 numerical-stability parity, so a
    naive strict §5 would clip the bore gas to ``T_ambient · 1.01`` on
    step 0 (only ambient passes the Y > Y_min filter; its T_flame is
    T_initial). The guard

        T_ceiling[i] = max(T_ceiling[i], T_initial_gas · 1.01)

    keeps the ceiling above the IC gas temperature during the pre-fill
    window. Once ambient purges (typically within ~1 ms in chamber
    fill), the per-species max becomes the binding bound.

    The strict ceiling tightens overshoot detection in three regimes:
    - Pyrogen-only cells get ceiling = T_flame_pyrogen · 1.01 once
      pyrogen displaces ambient (cell-0 during the early igniter
      pulse). Under the previous relaxed (max-of-all-species) form,
      these cells could climb ~9% above pyrogen T_flame.
    - Pure-ambient cells far from the igniter (last few cells late in
      ignition transient) cap at T_initial_gas · 1.01 — but only if
      that exceeds T_ambient · 1.01, which it does for the IC's hot
      seed temperature.
    - v0.8.0 multi-grain configurations don't broadcast the hottest
      species's T_flame to cells that don't contain it.
    """
    S = species_params.shape[0]
    ic_guard = T_initial_gas * 1.01
    for i in range(N):
        T_flame_max = 0.0
        for s in range(S):
            if Y[i, s] > Y_min:
                Tf = species_params[s, 3]
                if Tf > T_flame_max:
                    T_flame_max = Tf
        ceiling = T_flame_max * 1.01
        if ceiling < ic_guard:
            ceiling = ic_guard
        T_ceiling_arr[i] = ceiling


@njit(cache=True, fastmath=True)
def _advect_species(
    Y, rho_old, rho_new, u, A_port,
    nozzle_mdot, dx, dt,
    mass_source_by_species, N, S,
):
    """
    v0.7.1 — mass-fraction-conservative upwind passive-scalar advection.

    Updates ``Y[N, S]`` in-place after the PISO step has produced
    ``rho_new`` (cell densities) and ``u`` (face velocities). The
    interior-face density uses central averaging
    ``0.5 * (rho_old[j-1] + rho_old[j])`` to match the PISO mass balance;
    face area uses the same arithmetic mean of adjacent ``A_port``; Y at
    each face is taken from the upwind cell. The nozzle face (j=N)
    carries ``Y[N-1, :]`` outward at rate ``nozzle_mdot``.

    After flux + source accumulation, each cell's new species mass is
    divided by the new cell mass ``rho_new[i] * A_port[i] * dx`` to get
    Y_new; tiny FP drift is clamped via per-cell renormalization so
    ``sum_s Y[i, s] = 1`` and ``0 <= Y[i, s] <= 1``.

    Notes
    -----
    For Phase 1, mass closure with PISO is exact to O(round-off) for
    interior cells. The nozzle cell's closure depends on the
    ``nozzle_mdot`` argument matching the value PISO used to update
    ``rho_new[N-1]``; small inconsistencies (O(dt^2)) are absorbed by
    the renormalization step.
    """
    # Per-face species mass flux (mass-of-species crossing face j, signed).
    face_species_flux = np.zeros((N + 1, S))

    # Interior faces 1..N-1: rho_face and A_face from central averages
    # of adjacent cell quantities; upwind on Y per face velocity sign.
    for j in range(1, N):
        rho_face = 0.5 * (rho_old[j - 1] + rho_old[j])
        A_face_local = 0.5 * (A_port[j - 1] + A_port[j])
        mass_flux_face = rho_face * u[j] * A_face_local * dt
        if mass_flux_face >= 0.0:
            upwind = j - 1
        else:
            upwind = j
        for s in range(S):
            face_species_flux[j, s] = mass_flux_face * Y[upwind, s]

    # Nozzle face N: outflow carries cell N-1 composition.
    for s in range(S):
        face_species_flux[N, s] = nozzle_mdot * dt * Y[N - 1, s]

    # Face 0 (head-end wall) — flux remains zero (already initialized).

    # Per-cell update. NOTE on units: mass_source_by_species[i, s] is
    # a cell rate per unit axial length [kg/(m*s)] (mirrors the existing
    # mass_source convention used by PISO's continuity residual at
    # solver.py: b_rhs[i] = mass_source[i] * dx - ...). The mass added
    # to cell i per step is therefore source * dx * dt, NOT source * V * dt.
    for i in range(N):
        V_i = A_port[i] * dx
        m_old = rho_old[i] * V_i
        m_new = rho_new[i] * V_i
        if m_new < 1.0e-12:
            # Degenerate cell (no mass); leave Y unchanged.
            continue
        total = 0.0
        for s in range(S):
            new_mass_s = (m_old * Y[i, s]
                          + face_species_flux[i, s]      # in through west face
                          - face_species_flux[i + 1, s]  # out through east face
                          + mass_source_by_species[i, s] * dx * dt)
            if new_mass_s < 0.0:
                new_mass_s = 0.0
            Y[i, s] = new_mass_s / m_new
            total += Y[i, s]
        # Renormalize to sum=1; clamps fp drift over thousands of steps.
        if total > 1.0e-12:
            inv = 1.0 / total
            for s in range(S):
                Y[i, s] *= inv


@njit(cache=True, fastmath=True)
def _goodman_ignition_sources_and_mass(
    P, T, T_surf, delta, has_ignited, is_burning, is_grain, ignition_time,
    r_total, r_erosive, mass_source, thermal_source,
    C_burn, endface_msource, pyrogen_surface_heat_flux,
    radiation_heat_flux, radiation_sink_power, radiation_emitter,
    x_centers, Re, D_hyd, f_darcy,
    t, dt, rho_propellant, T_flame, T_initial,
    Pr, k_thermal, roughness, kappa, solid_alpha, k_solid,
    T_ignition, N, dx, mdot_igniter, T_ig,
    Cp_propellant, Cp_pyrogen,
    pyrogen_surface_heat_flux_w_m2, radiation_emissivity,
    diagnostic_disable_radiation_gas_sink,
    tau_establishment,
    mass_source_by_species,  # v0.7.1: [N, S] per-species rates [kg/s/m]
    # v0.7.2 Phase B-v2: flame-front h_c augmentation
    flame_spread_augment,    # [N] per-cell h_c multiplier (1.0 = no boost)
    flame_spread_enabled,    # bool — gate the augmentation off entirely
    # v0.7.3 Phase B.2: extend radiation_emitter to pyrogen-hot cells
    Y_species,               # [N, S] per-cell species mass fractions
    # v0.7.3 Phase B.4: pre-computed per-cell pyrogen surface heat flux
    pyrogen_heat_flux_arr_in, # [N] W/m², pre-capped per topology+mode
    # v0.7.4 Phase F: flame-spread front exposure gate
    flame_front_enabled,     # bool — gate ALL surface heating to the front
    ignitable,               # [N] bool — cell exposed (cartridge or behind front)
    topology_code,           # 0=forward_plenum (single-cell DeMar induction exempt)
    # Solver-verification T10: start-of-step local gas state, for the debit of
    # the gas that fills the volume the regressing grain opens
    rho_gas,                 # [N] cell gas density [kg/m^3]
    Cp_mix_gas,              # [N] cell mixture Cp [J/(kg K)]
    # IG2: opt-in ignition heat-transfer closure (0 = default Gnielinski)
    ignition_h_code,         # 0 gnielinski, 1 peretz, 2 keller
    u_cell, A_port, A_throat, R_mix_gas, gamma_mix_gas,
    # IG3: which port area Peretz's (A_p/A_t)^0.4 uses (0 local, 1 head, 2 aft)
    peretz_ap_code,
    # SM1: solid-heating model (0 goodman, 1 peretz, 2 conduction) + conduction state
    solid_model_code, solid_y, solid_T, solid_age, solid_nact, solid_work,
    solid_age_fraction,
    # SM2: ignition-event diagnostics (always on, bookkeeping only) and the opt-in
    # open-loop surface-forcing recorder
    ign_diag, rec_on, rec_h, rec_Td, rec_t, rec_dt, rec_k,
):
    """Goodman surface-temperature ignition and propellant source assembly.

    ``tau_establishment`` (seconds) is a post-ignition burn-establishment
    timescale. If positive, each cell's effective burn rate ramps linearly
    from 0 to its steady value over ``tau_establishment`` after the cell
    crosses ``T_ignition``. r_total and r_erosive are scaled in place so
    the next ``advance_bore_regression`` sees the same effective rate. Set
    to 0.0 to disable (no ramp; pure step at ignition, matching Peretz/
    Pardue/Cavallini's instantaneous-ignition convention).

    v0.7.1 (Phase 3.5): the scalar `Cp_gas` arg has split into two
    species-specific arguments. Propellant combustion sources (grain
    sidewall + endface) multiply T_flame contributions by
    ``Cp_propellant``; the pyrogen-to-surface heat-transfer cap inside
    ``_pyrogen_surface_heat_power`` uses ``Cp_pyrogen``. Each source
    species injects its OWN enthalpy, which is the physically correct
    multi-species accounting. Phase 3's previous use of a single
    ``Cp_gas`` (= propellant Cp) for both was a behavior-preserving
    placeholder during the unit shift.
    """
    n_burning = 0
    n_ignited = 0
    mass_sum = 0.0
    pyrogen_surface_heat_power = 0.0
    radiation_heat_power = 0.0
    radiation_sink_total_power = 0.0
    normal_sidewall_thermal_power = 0.0
    erosive_sidewall_thermal_power = 0.0
    endface_thermal_power = 0.0
    # Solver-verification T6: power the v0.7.4 convective wall heat-loss sink
    # removes from the gas [W] (sum of h_conv*dT*C_burn*dx). Previously folded
    # into thermal_source with no history of its own.
    convective_wall_sink_power = 0.0
    # Solver-verification T10: debit of the local gas that fills the volume
    # the regressing grain opens, as mass [kg/s] and sensible enthalpy [W].
    fill_debit_massflow = 0.0
    fill_debit_power = 0.0
    # v0.7.3 Phase B.2: a cell radiates meaningfully if EITHER its
    # propellant is burning (existing criterion) OR its bore gas is
    # majority-pyrogen species (new criterion). Pyrogen-hot cells with
    # T_gas ≈ T_flame_pyrogen emit Stefan-Boltzmann radiation just like
    # propellant-burning cells; the previous `is_burning[i]`-only gating
    # missed pyrogen-driven radiation entirely under uncontained
    # topologies. The Y > 0.5 threshold keeps the test cheap and
    # self-deactivates once propellant gas displaces pyrogen in cells
    # far from the cartridge.
    # v0.7.3 Phase B.4: pyrogen surface heat flux is now pre-computed
    # by _compute_pyrogen_heat_flux_arr (passed in as
    # pyrogen_heat_flux_arr_in), replacing the head_grain_cell /
    # pyrogen_heat_target single-cell special case.
    Y_emit_threshold = 0.5
    # IG3: fixed reference port area for the Peretz A_p/A_t factor (head / aft
    # grain cell, instantaneous); 'local' (code 0) uses the cell's own area.
    ap_ref_idx = -1
    if peretz_ap_code == 1:
        for i in range(N):
            if is_grain[i]:
                ap_ref_idx = i
                break
    elif peretz_ap_code == 2:
        for i in range(N - 1, -1, -1):
            if is_grain[i]:
                ap_ref_idx = i
                break
    for i in range(N):
        pyrogen_surface_heat_flux[i] = 0.0
        radiation_heat_flux[i] = 0.0
        radiation_sink_power[i] = 0.0
        radiation_emitter[i] = (is_burning[i] or
                                Y_species[i, 0] > Y_emit_threshold)

    # v0.7.1: reset per-species mass source rows. Each step recomputes
    # contributions from grain (s=1) inside this loop; pyrogen (s=0) is
    # set after this function returns. Other species remain zero.
    S_local = mass_source_by_species.shape[1]
    for i in range(N):
        for s in range(S_local):
            mass_source_by_species[i, s] = 0.0

    rec_idx = rec_k[0]
    rec_active = rec_on and rec_idx < rec_h.shape[0]
    if rec_active:
        rec_t[rec_idx] = t
        rec_dt[rec_idx] = dt
        for i in range(N):
            rec_h[rec_idx, i] = -1.0
            rec_Td[rec_idx, i] = 0.0
    for i in range(N):
        mass_source[i] = 0.0
        thermal_source[i] = 0.0

        if not is_grain[i]:
            is_burning[i] = False
            r_total[i] = 0.0
            r_erosive[i] = 0.0
        else:
            # v0.7.4 Phase F: withhold bulk-gas surface heating from grain
            # cells the flame front has not yet reached. A gated-out
            # unignited cell skips the whole heating + ignition sub-block →
            # T_surf stays frozen at T_initial (no pre-heat, no ignition
            # until the front arrives). flame_front_enabled=False →
            # heat_cell is always True → byte-for-byte the prior behaviour.
            #
            # The seed grain cell is always in `ignitable` (set by
            # _advance_flame_front), so it lights from the head-end pyrogen
            # plume + bore gas, and the front then propagates from it. DeMar
            # surface flux (applied below) is delivered only at ignitable
            # cells, so it cannot ignite cells ahead of the front.
            #
            # BUGFIX (2026-06-16): the previous gate also exempted ANY cell
            # carrying forward_plenum DeMar flux (pyrogen_heat_flux_arr_in[i]
            # > 0). But that flux targets `head_grain_cell` = the FIRST
            # UNIGNITED grain cell, which marches down the grain as cells
            # ignite — a self-propagating DeMar cascade that bypassed the
            # front entirely (flame_front_velocity had ~no effect on
            # forward_plenum; full ignition was fixed at ~4.6 ms for v=15..100
            # m/s on Chunc). Removing the exemption makes the front
            # authoritative: ignition spread now tracks grain_length / v.
            if not flame_front_enabled:
                heat_cell = True
            elif ignitable[i]:
                heat_cell = True
            else:
                heat_cell = False
            if rec_active and has_ignited[i]:
                # SM2: the forcing this cell WOULD see if it were still unlit (bare
                # convective h of the local gas state), so an open-loop replay can run
                # past the reference run's own ignition time.
                if ignition_h_code == 0:
                    h_v = _bare_heat_transfer_coeff(
                        Re[i], D_hyd[i], x_centers[i], f_darcy[i],
                        Pr, k_thermal, T[i], T_surf[i], kappa,
                    )
                elif ignition_h_code == 1:
                    if ap_ref_idx >= 0:
                        A_apv = A_port[ap_ref_idx]
                    else:
                        A_apv = A_port[i]
                    h_v = _peretz_heat_transfer_coeff(
                        P[i], u_cell[i], T[i], T_surf[i], R_mix_gas[i],
                        Cp_mix_gas[i], gamma_mix_gas[i], D_hyd[i],
                        x_centers[i], A_apv, A_throat,
                    )
                else:
                    h_v = _keller_heat_transfer_coeff(
                        rho_gas[i], u_cell[i], T[i], D_hyd[i],
                    )
                rec_h[rec_idx, i] = h_v
                rec_Td[rec_idx, i] = T[i]
            if (not has_ignited[i]) and heat_cell:
                if ignition_h_code == 0:
                    h_c = _bare_heat_transfer_coeff(
                        Re[i], D_hyd[i], x_centers[i], f_darcy[i],
                        Pr, k_thermal, T[i], T_surf[i], kappa,
                    )
                elif ignition_h_code == 1:
                    if ap_ref_idx >= 0:
                        A_ap = A_port[ap_ref_idx]
                    else:
                        A_ap = A_port[i]
                    h_c = _peretz_heat_transfer_coeff(
                        P[i], u_cell[i], T[i], T_surf[i], R_mix_gas[i],
                        Cp_mix_gas[i], gamma_mix_gas[i], D_hyd[i],
                        x_centers[i], A_ap, A_throat,
                    )
                else:
                    h_c = _keller_heat_transfer_coeff(
                        rho_gas[i], u_cell[i], T[i], D_hyd[i],
                    )
                # v0.7.2 Phase B-v2: flame-front augmentation. The
                # caller pre-fills flame_spread_augment[i] = boost only
                # when cell i is the immediate downstream neighbor of a
                # cell that ignited within tau_window seconds; otherwise
                # it stays at 1.0 (no-op). flame_spread_enabled=False
                # skips the multiply entirely → byte-for-byte Phase A.
                if flame_spread_enabled:
                    h_c = h_c * flame_spread_augment[i]
                h_total = h_c
                h_driver_num = h_c * T[i]
                # v0.7.4 energy-balance fix: capture the BARE convective
                # coefficient and the gas-wall ΔT (pre-Goodman-update) for
                # the convective wall heat-loss sink applied after the
                # update. Only the convective channel is debited here:
                # adjacent-cell propellant radiation has its own sink
                # (radiation_sink_power) and the pyrogen flux is the
                # pyrogen's energy budget (handled separately).
                h_conv = h_c
                dT_conv = T[i] - T_surf[i]
                # v0.7.3 Phase B.4: consume the pre-computed per-cell
                # pyrogen heat flux (already topology/mode-dispatched
                # and sensible-capped by _compute_pyrogen_heat_flux_arr).
                # Replaces the v0.7.0-v0.7.3-phaseA single-cell DeMar
                # special case at i == pyrogen_heat_target. The flux
                # may be zero (mode 'none', non-cartridge cells in
                # uncontained 'demar', etc.) — guarded below.
                flux_w_m2 = pyrogen_heat_flux_arr_in[i]
                if flux_w_m2 > 0.0:
                    contact_area = C_burn[i] * dx
                    if contact_area > 0.0:
                        delivered_power = flux_w_m2 * contact_area
                        pyrogen_surface_heat_power += delivered_power
                        pyrogen_surface_heat_flux[i] = flux_w_m2
                        h_ig = flux_w_m2 / max(T_ig - T_surf[i], 1.0e-9)
                        h_total += h_ig
                        h_driver_num += h_ig * T_ig

                if radiation_emissivity > 0.0 and C_burn[i] > 0.0:
                    # Adjacent-cell radiation uses the *local* gas
                    # temperature of each burning neighbor as the emitter
                    # temperature. Using the constant adiabatic T_flame
                    # here overstates the flux during a cold-start
                    # transient (a just-ignited neighbor whose gas is
                    # still ramping radiates as if fully developed),
                    # producing an unphysically fast radiative ignition
                    # chain. The gas-energy sink debits the same cell at
                    # the same temperature, so the exchange is
                    # self-consistent.
                    rad_flux = 0.0
                    rad_driver_num = 0.0
                    if (i > 0 and radiation_emitter[i - 1]
                            and T[i - 1] > T_surf[i]):
                        rad_left = radiation_emissivity * STEFAN_BOLTZMANN * (
                            T[i - 1] ** 4 - T_surf[i] ** 4
                        )
                        if rad_left > 0.0:
                            rad_flux += rad_left
                            rad_driver_num += rad_left * T[i - 1]
                            if not diagnostic_disable_radiation_gas_sink:
                                sink = rad_left * C_burn[i] * dx
                                radiation_sink_power[i - 1] += sink
                                radiation_sink_total_power += sink
                    if (i < N - 1 and radiation_emitter[i + 1]
                            and T[i + 1] > T_surf[i]):
                        rad_right = radiation_emissivity * STEFAN_BOLTZMANN * (
                            T[i + 1] ** 4 - T_surf[i] ** 4
                        )
                        if rad_right > 0.0:
                            rad_flux += rad_right
                            rad_driver_num += rad_right * T[i + 1]
                            if not diagnostic_disable_radiation_gas_sink:
                                sink = rad_right * C_burn[i] * dx
                                radiation_sink_power[i + 1] += sink
                                radiation_sink_total_power += sink

                    if rad_flux > 0.0:
                        radiation_heat_flux[i] = rad_flux
                        cell_power = rad_flux * C_burn[i] * dx
                        radiation_heat_power += cell_power
                        rad_driver_T = rad_driver_num / rad_flux
                        h_rad = rad_flux / max(rad_driver_T - T_surf[i], 1.0e-9)
                        h_total += h_rad
                        h_driver_num += h_rad * rad_driver_T

                h_driver = T[i]
                if h_total > 0.0:
                    h_driver = h_driver_num / h_total
                    h_c = h_total

                T_s_before = T_surf[i]
                if rec_active:
                    rec_h[rec_idx, i] = h_c
                    rec_Td[rec_idx, i] = h_driver
                # SM1: the solid-heating model. The surface boundary is the
                # same combined Robin flux for all three: h_c (h_driver - T_s),
                # with h_c / h_driver already folding in the convective,
                # pyrogen and adjacent-radiation channels above.
                if solid_model_code == 0:
                    new_delta, new_T_surf = _step_goodman_ode(
                        delta[i], T_surf[i], h_c, h_driver, T_initial,
                        solid_alpha, k_solid, dt,
                    )
                    delta[i] = new_delta
                elif solid_model_code == 1:
                    new_T_surf = _step_peretz_ode(
                        T_surf[i], h_c, h_driver, T_initial,
                        solid_alpha, k_solid, dt,
                    )
                else:
                    new_T_surf = T_surf[i]
                    if h_c > 0.0 or solid_age[i] > 0.0:
                        new_T_surf, solid_nact[i] = _step_conduction(
                            solid_T[i], solid_y, solid_nact[i], solid_age[i],
                            h_c, h_driver, T_initial, solid_alpha, k_solid,
                            dt, solid_age_fraction,
                            solid_work[0], solid_work[1],
                        )
                        solid_age[i] += dt
                T_surf[i] = new_T_surf

                # v0.7.4 energy-balance fix: convective wall heat-loss sink.
                # The bore gas convectively heated this UNIGNITED wall by
                # q = h_conv*(T_gas - T_surf)*C_burn [W/m]; debit it from the
                # gas so the system is NOT adiabatic-at-the-walls during the
                # ignition transient (the gas must cool as it pumps energy
                # into cold boundaries). Mirrors the existing radiation sink.
                # Clamped to heating (dT_conv > 0). Vanishes naturally once
                # the cell ignites (this block only runs while unignited) and
                # the wall becomes a transpiring source — no Boolean switch,
                # no temperature/ignition lever beyond the existing T_ignition
                # gate. thermal_source[i] may go negative for cells with hot
                # advected gas but no local source (a net energy sink) — PISO
                # handles this and the bore gas cools correctly.
                if h_conv > 0.0 and dT_conv > 0.0 and C_burn[i] > 0.0:
                    thermal_source[i] -= h_conv * dT_conv * C_burn[i]
                    convective_wall_sink_power += (
                        h_conv * dT_conv * C_burn[i] * dx
                    )

                if _surface_has_ignited(T_surf[i], T_ignition):
                    has_ignited[i] = True
                    is_burning[i] = True
                    ignition_time[i] = t
                    # SM2: what lit the cell (diagnostic only)
                    ign_diag[i, 0] = h_c
                    ign_diag[i, 1] = h_driver
                    ign_diag[i, 2] = h_c * (h_driver - T_s_before)
                    ign_diag[i, 3] = T_surf[i]
                    ign_diag[i, 4] = T[i]
                    ign_diag[i, 5] = u_cell[i]
                    ign_diag[i, 6] = abs(u_cell[i]) / np.sqrt(
                        gamma_mix_gas[i] * R_mix_gas[i] * T[i])
                    ign_diag[i, 7] = P[i]
                    ign_diag[i, 8] = t
                    ign_diag[i, 9] = h_c * (h_driver - T_surf[i])

            if has_ignited[i]:
                is_burning[i] = True
                n_ignited += 1
                n_burning += 1

        if is_burning[i] and C_burn[i] > 0.0:
            if tau_establishment > 0.0:
                dt_since_ign = t - ignition_time[i]
                if dt_since_ign < 0.0:
                    phi_est = 0.0
                elif dt_since_ign < tau_establishment:
                    phi_est = dt_since_ign / tau_establishment
                else:
                    phi_est = 1.0
                if phi_est < 1.0:
                    r_total[i] *= phi_est
                    r_erosive[i] *= phi_est
            prop_source = rho_propellant * r_total[i] * C_burn[i]
            r_normal = r_total[i] - r_erosive[i]
            if r_normal < 0.0:
                r_normal = 0.0
            normal_source = rho_propellant * r_normal * C_burn[i]
            erosive_source = prop_source - normal_source
            if erosive_source < 0.0:
                erosive_source = 0.0
            mass_source[i] += prop_source
            mass_source_by_species[i, _SPECIES_PROPELLANT] += prop_source
            thermal_source[i] += prop_source * T_flame * Cp_propellant
            normal_sidewall_thermal_power += normal_source * T_flame * Cp_propellant * dx
            erosive_sidewall_thermal_power += erosive_source * T_flame * Cp_propellant * dx

        if endface_msource[i] > 0.0:
            mass_source[i] += endface_msource[i]
            mass_source_by_species[i, _SPECIES_PROPELLANT] += endface_msource[i]
            thermal_source[i] += endface_msource[i] * T_flame * Cp_propellant
            endface_thermal_power += endface_msource[i] * T_flame * Cp_propellant * dx

        mass_sum += mass_source[i]

        # Solver-verification T10: the standard fixed-area formulation has a
        # net mass source (rho_p - rho_g) * dV/dt, because the burned
        # propellant first fills the volume it vacates. The solver advances
        # continuity at fixed A_port and update_cell_geometry then grows
        # A_port at fixed rho, which re-adds rho_g*dA*dx -- so the local gas
        # filling the opened volume is debited here (mass, sensible enthalpy,
        # per species), using the SAME effective regression rate x area as the
        # source above (sidewall r_total*C_burn, f_active already inside
        # C_burn; end face endface_msource/rho_p). Start-of-step cell state.
        # Momentum is deliberately left unchanged. mass_sum above stays the
        # GROSS generation (total_mass_produced / mass_balance_error).
        fill_vol = 0.0
        if is_burning[i] and C_burn[i] > 0.0:
            fill_vol += r_total[i] * C_burn[i]
        if endface_msource[i] > 0.0:
            fill_vol += endface_msource[i] / rho_propellant
        if fill_vol > 0.0:
            fill_mass = rho_gas[i] * fill_vol
            fill_enthalpy = fill_mass * Cp_mix_gas[i] * T[i]
            mass_source[i] -= fill_mass
            thermal_source[i] -= fill_enthalpy
            for sp in range(S_local):
                mass_source_by_species[i, sp] -= fill_mass * Y_species[i, sp]
            fill_debit_massflow += fill_mass * dx
            fill_debit_power += fill_enthalpy * dx

    # Radiation gas-energy sink debit (solver-verification T9). Done AFTER
    # the cell loop: the sink on the LEFT neighbour (radiation_sink_power
    # [i-1] += ...) is accumulated while processing cell i, i.e. after cell
    # i-1's iteration finished, so an in-loop debit never saw it (an energy
    # leak: recorded, never subtracted). No-op when emissivity == 0.
    if dx > 0.0:
        for i in range(N):
            if radiation_sink_power[i] > 0.0:
                thermal_source[i] -= radiation_sink_power[i] / dx

    if rec_active:
        rec_k[0] = rec_idx + 1
    return (n_burning, n_ignited, mass_sum,
            pyrogen_surface_heat_power, radiation_heat_power,
            radiation_sink_total_power,
            normal_sidewall_thermal_power, erosive_sidewall_thermal_power,
            endface_thermal_power, convective_wall_sink_power,
            fill_debit_massflow, fill_debit_power)


# ================================================================
# Compiled time loop
# ================================================================

@njit(cache=True, nogil=True, fastmath=True)
def _run_time_loop(
    # --- Cell arrays (N) ---
    rho, u, P, T,
    D_port, x_centers, A_port, C_burn, D_hyd,
    is_grain, endface_msource,
    is_burning, has_ignited, ignition_time,
    r_total, r_erosive,
    mass_source, thermal_source, momentum_source, pyrogen_surface_heat_flux,
    radiation_heat_flux, radiation_sink_power, radiation_emitter,
    f_darcy, Re, Mach, u_cell,
    T_surf, delta,
    regress,
    # --- Segment arrays (N_seg) ---
    seg_x_start, seg_length,
    seg_fwd_regression, seg_aft_regression,
    seg_inhibit_fwd, seg_inhibit_aft,
    seg_D_bore_fwd, seg_D_bore_aft,
    cell_D_bore_init, cell_segment_id,
    cell_wall_web, cell_segment_type, cell_fmm_idx,
    fmm_offset, fmm_reg_flat, fmm_perim_flat, fmm_port_flat,
    # --- Geometry scalars ---
    N, N_seg, dx, cell_D_outer,
    # --- Gas/propellant scalars ---
    gamma, R_specific, T_flame, Cp_gas, mu_gas, k_thermal, Pr,
    rho_propellant, Cps, T_surface, T_initial, T_initial_gas, k_solid,
    # --- Burn rate tabs (parallel arrays) ---
    tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        # --- Simulation parameters ---
    roughness, kappa,
    cfl_target, dt_max, burn_update_interval, geometry_update_interval,
    source_cfl_factor, port_mach_cap, pressure_floor_pa,
    piso_energy_coupling_iters, piso_corr2_storage,
    T_ignition, P_ambient, ambient_temperature,
    diagnostic_disable_erosive, diagnostic_disable_endfaces,
    diagnostic_disable_momentum, diagnostic_disable_pyrogen_surface_heating,
    diagnostic_disable_adjacent_radiation,
    diagnostic_disable_radiation_gas_sink,
    igniter_axial_momentum_fraction, pyrogen_surface_heat_flux_w_m2,
    radiation_emissivity,
    tau_establishment,
    t_max, P_cutoff,
    erosion_coeff, slag_coeff, throat_is_evolving,
    snapshot_interval,
    # --- Precomputed ---
    gamma_R, Gamma_crit, nozzle_denom,
    D_throat_init, A_throat_init,
    pyrogen_params_arr, chamber_params_arr, plenum_state,
    # --- Output: time history (pre-allocated) ---
    time_hist, P_head_hist, P_exit_hist, D_throat_hist,
    Kn_hist, massflow_hist, P_ig_hist, T_ig_hist, mdot_ig_hist,
    mdot_grain_hist,
    m_pyrogen_hist,
    gas_sensible_energy_before_hist, gas_sensible_energy_hist,
    gas_sensible_dE_dt_hist,
    normal_sidewall_thermal_power_hist, erosive_sidewall_thermal_power_hist,
    endface_thermal_power_hist, convective_scalar_flux_power_hist,
    clipping_correction_power_hist, pyrogen_enthalpy_power_hist,
    pyrogen_surface_heat_power_hist, gas_surface_heat_sink_power_hist,
    radiation_heat_power_hist, radiation_sink_power_hist,
    nozzle_enthalpy_power_hist, thermal_source_power_hist,
    energy_residual_hist,
    pyrogen_momentum_expected_hist, pyrogen_momentum_deposited_hist,
    pyrogen_momentum_residual_hist,
    dt_hist, n_burning_hist, n_ignited_hist,
    radiation_emitter_count_hist, radiation_receiver_count_hist,
    min_gas_temperature_hist, max_gas_temperature_hist,
    min_surface_temperature_hist, max_surface_temperature_hist,
    min_pressure_hist, max_pressure_hist, max_mach_hist,
    pressure_floor_cells_hist, pressure_floor_mass_hist,
    # --- T6 energy audit (accounting only) ---
    pressure_floor_energy_hist, convective_wall_sink_power_hist,
    fill_debit_massflow_hist, fill_debit_power_hist,
    nozzle_enthalpy_power_post_hist, nozzle_massflow_piso_hist,
    jump_eos_J_hist, jump_mixture_J_hist, jump_geometry_J_hist,
    jump_eos_kg_hist, jump_mixture_kg_hist, jump_geometry_kg_hist,
    jump_total_kg_hist,
    coupling_passes_hist, coupling_rel_dT_hist, eos_jump_p2_kg_hist,
    max_hist,
    # --- Output: snapshots (pre-allocated) ---
    snap_data, snap_times, max_snaps,
    snap_seg_fwd, snap_seg_aft,
    # --- v0.7.1: N-species state ---
    Y_species, species_params_arr, mass_source_by_species,
    gamma_mix_arr, Cp_mix_arr, R_mix_arr, M_mix_arr,
    T_ceiling_arr,
    # --- PB/P2: per-cell numerical-limit diagnostics (bookkeeping only) ---
    limit_activation_counts, limit_duration_s, limit_abs_energy_j,
    limit_first_time_s, limit_last_time_s,
    pressure_floor_max_deficit_pa,
    mach_limit_activation_counts, mach_limit_duration_s,
    mach_limit_first_time_s, mach_limit_last_time_s,
    # --- v0.7.2 Phase A: pyrogen axial distribution ---
    pyrogen_axial_weights,
    # --- v0.7.2 Phase B-v2: flame-front h_c augmentation ---
    flame_spread_augment, flame_spread_enabled,
    flame_spread_tau, flame_spread_boost,
    # --- v0.7.3 Phase A: igniter topology (uncontained vs plenum) ---
    topology_code,                # 0=forward_plenum, 1=head_basket, 2=aft_basket
    cart_i_start, cart_i_end,     # cartridge cell range (inclusive); both 0 for forward_plenum
    A_burn_per_cell,              # initial burning surface per cartridge cell (uncontained)
    mdot_uncontained_arr,         # [N] scratch for per-cell uncontained mdot
    # --- v0.7.3 Phase B.4: pyrogen-to-surface heat delivery ---
    heat_delivery_mode_code,      # 0=none, 1=demar, 2=radiation (uncontained only)
    demar_flux_w_m2_uncontained,  # heat_flux_cal_cm2_s converted to W/m² (uncontained DeMar)
    T_flame_pyrogen,              # adiabatic pyrogen flame T [K] (radiation emitter T)
    pellet_emissivity,            # Stefan-Boltzmann emissivity [-] (radiation only)
    radiation_absorption_length_m,# Beer-Lambert attenuation length [m] (radiation only)
    pyrogen_heat_flux_arr,        # [N] scratch buffer filled per step
    # --- v0.7.4 Phase F: flame-spread front gate ---
    flame_front_enabled,          # bool — opt-in front gate
    flame_front_velocity,         # [m/s] lateral flame-spread speed
    ignitable,                    # [N] bool scratch — front-exposure mask
    # --- v0.7.4 Phase Z: Z-N dynamic burn-rate relaxation ---
    zn_enabled,                   # bool — opt-in relaxation
    kappa_zn,                     # O(1) prefactor on tau = kappa_zn*alpha_s/r^2
    r_dyn,                        # [N] persistent relaxed burn rate
    r_qs_persist,                 # [N] held quasi-steady total target
    r_ero_qs_persist,             # [N] held quasi-steady erosive target
    # --- spike-reopen probe: Mukunda-Paul high-g erosive ceiling ---
    erosive_mp_cap,               # float — 0.0 off; else eta cap at value×MP Eq. 12
    erosive_model_code,           # int — registry code, 0=ma2020, 1=mukunda_paul, 2=saint_robert
    # --- IG2: opt-in ignition heat-transfer closure (0 = default Gnielinski) ---
    ignition_h_code,
    peretz_ap_code,               # IG3: 0 local, 1 head, 2 aft (peretz mode only)
    stop_if_unignitable,          # IG3: bool, end the run once ignition is impossible
    # --- SM1: solid-heating model + conduction state ---
    solid_model_code,             # 0 goodman, 1 peretz, 2 conduction
    solid_y, solid_T, solid_age, solid_nact, solid_work, solid_age_fraction,
    # --- SM2: ignition-event diagnostics + open-loop forcing recorder ---
    ign_diag, rec_on, rec_h, rec_Td, rec_t, rec_dt, rec_k,
    # --- v0.8.0 Phase 6: live progress / cooperative cancel ---
    progress_state,               # [2] float64: [0]=progress 0..1 (written),
                                  #              [1]=cancel flag (read; >0.5 aborts)
):
    """
    The complete simulation time loop, compiled to native code.

    All sub-function calls (piso_step, compute_burn_rates, etc.) are
    Numba→Numba direct C calls with zero Python dispatch.

    Returns
    -------
    n_steps : int
        Number of time steps completed.
    n_snaps : int
        Number of snapshots written.
    total_mass_produced : float
    total_mass_nozzle : float
    first_burnthrough_time : float
        -1.0 if no burnthrough occurred.
    D_throat_final : float
    termination_code : int
        0 = t_max reached, 1 = complete burnout, 2 = pressure cutoff,
        3 = history array full, 4 = numerical collapse aborted,
        5 = canceled, 6 = stopped early: ignition impossible (IG3).
    (... then pyrogen_mass_burned, pyrogen_duration, pyrogen_peak_P, and,
    appended by solver-verification T3 -- pressure-floor accounting:)
    floor_energy_total : float
        Run total of the sensible energy [J] the 1 kPa pressure-floor clamp
        added (sum over steps of the PISO ``floor_energy_added``).
    floor_min_raw_P_global : float
        Minimum raw pre-clamp pressure [Pa] over all clamped cells and
        steps; ``FLOOR_MIN_RAW_P_NONE`` (finite sentinel) if never clamped.
    (... and, appended by solver-verification T5 -- refresh-gap diagnostics:)
    max_burn_refresh_gap : float
        Largest simulated time [s] between two consecutive burn-rate
        refreshes (``compute_burn_rates``, cadence ``burn_update_interval``).
    max_geometry_refresh_gap : float
        Same for ``update_cell_geometry`` (cadence ``geometry_update_interval``).
    """
    PI = 3.141592653589793

    # Numerical-collapse trip: any one of (dt < 1e-9 s, max Mach > 100,
    # max pressure > 1 GPa) for N consecutive steps aborts the loop with
    # termination_code = 4. Chamber Mach should not approach 1, so >100 is
    # already pathology; 1 GPa is ~20x any physical motor; dt collapse to
    # ~1e-10 was the documented signature in 2026-05-11 radiation runs.
    COLLAPSE_DT_THRESHOLD = 1.0e-9
    COLLAPSE_MACH_THRESHOLD = 100.0
    COLLAPSE_PRESSURE_THRESHOLD = 1.0e9
    COLLAPSE_CONSECUTIVE_STEPS = 3
    collapse_consecutive = 0

    # Source-aware CFL cap, refreshed once per step from the prior
    # step's thermal_source. dt_max bootstrap (no constraint) until the
    # first source assembly produces a real cap.
    dt_source_cap = dt_max
    SOURCE_DT_FLOOR = 1.0e-8

    t = 0.0
    step = 0
    n_burning = 0
    n_ignited = 0
    total_mass_produced = 0.0
    total_mass_nozzle = 0.0
    # Pressure-floor accounting (solver-verification T3): run totals of the
    # sensible energy the 1 kPa clamp adds and the minimum raw pre-clamp
    # pressure. Finite sentinel (fastmath: no infs); gated on n_floor_cells.
    floor_energy_total = 0.0
    floor_min_raw_P_global = FLOOR_MIN_RAW_P_NONE
    # Refresh-gap diagnostics (solver-verification T5): simulated time between
    # consecutive refreshes of each kind. Both kinds refresh on step 0 (t=0),
    # so the first gap is measured from there. Bookkeeping only -- reads ``t``,
    # never feeds back into the state.
    t_last_burn_refresh = 0.0
    t_last_geometry_refresh = 0.0
    max_burn_refresh_gap = 0.0
    max_geometry_refresh_gap = 0.0
    # Energy audit (solver-verification T6): ledger storage at the END of the
    # previous step, to difference against this step's pre-PISO storage.
    prev_gas_energy_after = 0.0
    prev_gas_mass_after = 0.0
    first_bt_time = -1.0
    last_snapshot_t = -snapshot_interval
    termination_code = 0
    p_peak = 0.0  # v0.8.0: running max head pressure, for the taildown metric

    D_throat = D_throat_init
    A_throat = A_throat_init
    pyrogen_initial_mass = plenum_state[0]
    # OB1: per-cell basket bookkeeping (area follows web consumption via the exact
    # cylinder form function; optional sequential bed lighting). DEFAULT for the
    # baskets; basket_area_law='constant' selects the pooled constant-area path.
    use_basket_cells = (
        topology_code != 0 and
        (chamber_params_arr[5] > 0.5 or chamber_params_arr[6] > 0.0)
    )
    m_cell_basket = np.zeros(N)
    t_light_basket = np.zeros(N)
    m_cell0_basket = 0.0
    A_cell0_basket = A_burn_per_cell
    if use_basket_cells:
        lo_b = cart_i_start if cart_i_start >= 0 else 0
        hi_b = cart_i_end if cart_i_end < N else N - 1
        n_b = hi_b - lo_b + 1
        if n_b < 1:
            n_b = 1
        m_cell0_basket = pyrogen_initial_mass / n_b
        v_bed = chamber_params_arr[6]
        for i in range(lo_b, hi_b + 1):
            m_cell_basket[i] = m_cell0_basket
            if v_bed > 0.0:
                if topology_code == 2:
                    t_light_basket[i] = (hi_b - i) * dx / v_bed
                else:
                    t_light_basket[i] = (i - lo_b) * dx / v_bed
    pyrogen_done = False
    pyrogen_peak_P = 0.0
    pyrogen_duration = 0.0
    solid_alpha = 0.0
    if k_solid > 0.0 and rho_propellant > 0.0 and Cps > 0.0:
        solid_alpha = k_solid / (rho_propellant * Cps)

    # v0.7.4 Phase F: seed the flame-spread front from the igniter topology.
    #   forward_plenum(0)/head_basket(1): seed at the head of the cartridge
    #     range, spread fore→aft (front_direction +1).
    #   aft_basket(2): seed at the aft of the cartridge range, spread
    #     aft→fore (front_direction -1).
    # x_front starts at the seed cell centre; it only advances once a cell
    # is actually burning (induction lights the seed via pyrogen flux).
    if topology_code == 2:
        front_seed_idx = cart_i_end
        front_direction = -1
    else:
        front_seed_idx = cart_i_start
        front_direction = 1
    if front_seed_idx < 0:
        front_seed_idx = 0
    elif front_seed_idx >= N:
        front_seed_idx = N - 1
    # v0.7.4 Phase F: snap the seed to the first GRAIN cell in the
    # propagation direction. A head_basket cartridge often sits in a
    # non-grain head cavity (cart_i_start maps to the motor head, but the
    # grain starts a few cells aft); without this snap the seed + cartridge
    # cells are all non-grain, nothing can ignite, and the front never
    # starts (the whole grain stays gated). forward_plenum is unaffected —
    # its DeMar head-cell exemption already lights the real grain cell.
    if front_direction > 0:
        for i in range(front_seed_idx, N):
            if is_grain[i]:
                front_seed_idx = i
                break
    else:
        for i in range(front_seed_idx, -1, -1):
            if is_grain[i]:
                front_seed_idx = i
                break
    x_front = x_centers[front_seed_idx]
    for i in range(N):
        ignitable[i] = True  # default exposed; gate active only when enabled

    # v0.7.1: workspace for pre-PISO density snapshot. Used by
    # _advect_species after PISO updates rho in-place. Face areas are
    # computed inline inside the advection kernel from the per-step
    # A_port to handle regression-driven area evolution.
    S_species_local = Y_species.shape[1]
    rho_pre_step = np.zeros(N)

    # v0.7.1 Phase 3: seed the per-cell mixture arrays from the initial Y
    # so the first PISO call sees consistent thermo before the post-PISO
    # _refresh_mixture_arrays runs. Also seed T_ceiling_arr (a constant
    # array under the current relaxed ceiling formula; refreshed each
    # step in case species_params ever becomes time-varying).
    _refresh_mixture_arrays(
        Y_species, species_params_arr,
        gamma_mix_arr, Cp_mix_arr, R_mix_arr, M_mix_arr, N,
    )
    _compute_T_ceiling_arr(
        Y_species, species_params_arr, T_ceiling_arr, N, T_initial_gas,
    )

    # Initial a_max for CFL — use the hottest local sound speed across
    # cells (per-cell γ·R varies once species advect).
    a_max = 0.0
    for i in range(N):
        a_local_init = (gamma_mix_arr[i] * R_mix_arr[i] * T[i]) ** 0.5
        if a_local_init > a_max:
            a_max = a_local_init

    hist_idx = 0
    snap_idx = 0

    # Progress denominator: total burnable web over every real grain cell
    # (constant). Cells with negligible web — gaps, pinned buffers, near-closed
    # OD-taper tip stations — drop out, so the progress fraction is the bulk
    # web actually consumed rather than the single fastest-burning cell.
    web_total = 0.0
    for i in range(N):
        if cell_wall_web[i] > 1e-9:
            web_total += cell_wall_web[i]

    while t < t_max:
        # --- v0.8.0: publish progress + honor cooperative cancel ---
        # Composite, monotonic metric (plain array ops keep the loop nopython):
        #  • burn phase — web-consumed fraction = total web burned / total web
        #    available (mass-representative: weighted by each cell's web depth,
        #    so thin-web OD-taper tip cells can't snap the bar near 1 the
        #    instant they burn through), or simulated-time fraction, whichever
        #    leads;
        #  • tail phase (web > 0.9) — head-pressure decay toward P_cutoff fills
        #    the final 10%. Without this the bar stalls near 99% through the
        #    low-rate taildown, where regression (hence web_frac) barely moves.
        # Gated on web > 0.9 so mid-burn pressure dips can't make the bar jump.
        if P[0] > p_peak:
            p_peak = P[0]
        web_num = 0.0
        for i in range(N):
            if cell_wall_web[i] > 1e-9:
                r = regress[i]
                if r > cell_wall_web[i]:
                    r = cell_wall_web[i]
                web_num += r
        web_frac = web_num / web_total if web_total > 1e-12 else 0.0
        time_frac = t / t_max
        progress = web_frac if web_frac > time_frac else time_frac
        if web_frac > 0.9 and p_peak > 2.0 * P_cutoff:
            # Chamber blowdown is ~exponential (dP/dt ∝ -P once burning
            # stops), so a linear-in-pressure tail crawls as P→P_cutoff.
            # Using log-pressure linearizes it in time: for P≈P_peak·e^(-t/τ),
            # ln(P_peak/P)/ln(P_peak/P_cutoff) ≈ t/t_end, a steady fill.
            p_head = P[0] if P[0] > 1.0 else 1.0
            tail = np.log(p_peak / p_head) / np.log(p_peak / P_cutoff)
            if tail < 0.0:
                tail = 0.0
            elif tail > 1.0:
                tail = 1.0
            tail_prog = 0.9 + 0.1 * tail
            if tail_prog > progress:
                progress = tail_prog
        if progress < progress_state[0]:
            progress = progress_state[0]   # never let the bar regress
        progress_state[0] = progress
        if progress_state[1] > 0.5:
            termination_code = 5
            break

        # --- Termination: complete burnout ---
        if n_ignited > 0 and n_burning == 0 and pyrogen_done:
            termination_code = 1
            break

        # --- History array full ---
        if hist_idx >= max_hist:
            termination_code = 3
            break

        # --- Time step ---
        dt = compute_dt_cfl(u, a_max, dx, N + 1, cfl_target, dt_max)
        if source_cfl_factor > 0.0 and dt_source_cap < dt:
            dt = dt_source_cap

        # ============================================
        # STEP 1: GEOMETRY
        # ============================================
        P_for_endface = P[0] if P[0] > 1e4 else 101325.0
        if not diagnostic_disable_endfaces:
            advance_endface_regression(
                seg_fwd_regression, seg_aft_regression,
                seg_length, seg_x_start,
                seg_inhibit_fwd, seg_inhibit_aft,
                N_seg, P_for_endface,
                tab_min_p, tab_max_p, tab_a, tab_n, n_tabs, dt,
            )

        advance_bore_regression(
            regress, r_total, dt, N,
            cell_wall_web, cell_segment_id,
        )

        # T6 energy audit (bookkeeping only): the storage change the geometry
        # refresh (A_port) makes to the carried state, evaluated as sum
        # before/after. Zero on steps where no refresh runs.
        geo_dM = 0.0
        geo_dE = 0.0
        if step % geometry_update_interval == 0:
            gap_geometry = t - t_last_geometry_refresh
            if gap_geometry > max_geometry_refresh_gap:
                max_geometry_refresh_gap = gap_geometry
            t_last_geometry_refresh = t
            geo_mass_pre, geo_energy_pre = _gas_mass_energy(
                rho, T, A_port, dx, Cp_mix_arr, N,
            )
            update_cell_geometry(
                regress, D_port, x_centers, dx, N, N_seg, cell_D_outer,
                seg_x_start, seg_length,
                seg_fwd_regression, seg_aft_regression,
                seg_inhibit_fwd, seg_inhibit_aft,
                cell_segment_id,
                P, rho_propellant,
                tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
                A_port, C_burn, D_hyd, is_grain, endface_msource,
                cell_D_bore_init, cell_wall_web,
                cell_segment_type, cell_fmm_idx,
                fmm_offset, fmm_reg_flat, fmm_perim_flat, fmm_port_flat,
            )
            geo_mass_post, geo_energy_post = _gas_mass_energy(
                rho, T, A_port, dx, Cp_mix_arr, N,
            )
            geo_dM = geo_mass_post - geo_mass_pre
            geo_dE = geo_energy_post - geo_energy_pre
            if diagnostic_disable_endfaces:
                for i in range(N):
                    endface_msource[i] = 0.0

        # ============================================
        # STEP 2: BURN RATES
        # ============================================
        if step % burn_update_interval == 0:
            gap_burn = t - t_last_burn_refresh
            if gap_burn > max_burn_refresh_gap:
                max_burn_refresh_gap = gap_burn
            t_last_burn_refresh = t
            r_total_new, r_erosive_new = compute_burn_rates(
                P, Re, D_hyd, x_centers, is_burning, roughness,
                Pr, k_thermal, Cp_gas,
                T_flame, T_surface,
                rho_propellant, Cps, T_initial,
                tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
                kappa, mu_gas, erosive_mp_cap, erosive_model_code, N,
            )
            for i in range(N):
                if diagnostic_disable_erosive:
                    r_normal = r_total_new[i] - r_erosive_new[i]
                    if r_normal < 0.0:
                        r_normal = 0.0
                    r_total[i] = r_normal
                    r_erosive[i] = 0.0
                else:
                    r_total[i] = r_total_new[i]
                    r_erosive[i] = r_erosive_new[i]
            # v0.7.4 Phase Z: snapshot the quasi-steady targets for the
            # relaxation. Refreshed on each burn-rate update and held
            # constant between updates (same cadence the baseline already
            # holds r_total between updates).
            if zn_enabled:
                for i in range(N):
                    r_qs_persist[i] = r_total[i]
                    r_ero_qs_persist[i] = r_erosive[i]

        # v0.7.4 Phase Z: relax the dynamic burn rate toward the quasi-steady
        # target EVERY step (dt-accurate), then overwrite r_total / r_erosive
        # so both consumers — the mass-source assembly and the next step's
        # advance_bore_regression — see r_dyn with no further plumbing. The
        # held r_qs_persist (not r_total) is the relaxation target, so the
        # overwrite cannot poison it between burn updates. Disabled →
        # untouched (byte-for-byte the prior behaviour).
        if zn_enabled:
            _advance_zn_burn_rate(
                r_dyn, r_qs_persist, is_burning, solid_alpha,
                kappa_zn, ZN_R_FLOOR, dt, N,
            )
            for i in range(N):
                r_total[i] = r_dyn[i]
                if r_qs_persist[i] > 1.0e-12:
                    r_erosive[i] = r_dyn[i] * (r_ero_qs_persist[i]
                                               / r_qs_persist[i])
                else:
                    r_erosive[i] = 0.0

        # ============================================
        # STEP 3: IGNITION + SOURCE ASSEMBLY
        # ============================================

        # v0.7.3 Phase A: branch on igniter topology.
        #   - forward_plenum (code 0): plenum-with-orifice (v0.7.0+ path).
        #   - head_basket / aft_basket (codes 1/2): uncontained pyrogen
        #     at local bore P, no plenum dynamics, no momentum injection.
        if topology_code == 0:
            # Forward plenum (existing path, unchanged behavior).
            new_plenum_state, mdot_igniter, _mdot_generated, P_ig = _step_plenum_ode(
                plenum_state, pyrogen_params_arr, chamber_params_arr, dt, P[0]
            )
            plenum_state[0] = new_plenum_state[0]
            plenum_state[1] = new_plenum_state[1]
            plenum_state[2] = new_plenum_state[2]
            T_ig = plenum_state[2]
        else:
            # Uncontained pyrogen (head_basket / aft_basket): each pyrogen
            # pellet burns at its host cell's local bore P; mass enters
            # the bore directly. No plenum gas state, no choked vent,
            # no separate P_ig. plenum_state[0] is reused as
            # m_pyrogen_remaining; plenum_state[1, 2] are vestigial here
            # but left in place to keep the state-vector shape stable.
            a_pyro = pyrogen_params_arr[0]
            n_pyro = pyrogen_params_arr[1]
            rho_pyro = pyrogen_params_arr[2]
            T_flame_pyro = pyrogen_params_arr[3]
            if use_basket_cells:
                new_remaining = _compute_uncontained_pyrogen_mdot_cells(
                    P, a_pyro, n_pyro, rho_pyro, A_cell0_basket,
                    m_cell_basket, m_cell0_basket, t, dt,
                    cart_i_start, cart_i_end, N,
                    mdot_uncontained_arr, pyrogen_params_arr[6],
                    pyrogen_params_arr[7], t_light_basket,
                    chamber_params_arr[7],
                )
            else:
                new_remaining = _compute_uncontained_pyrogen_mdot(
                    P, a_pyro, n_pyro, rho_pyro, A_burn_per_cell,
                    plenum_state[0], dt, cart_i_start, cart_i_end, N,
                    mdot_uncontained_arr, pyrogen_params_arr[6],
                    pyrogen_params_arr[7],
                )
            plenum_state[0] = new_remaining
            # Diagnostics: aggregate per-cell mdot for the existing
            # mdot_igniter / pyrogen_enthalpy_power channels.
            mdot_igniter = 0.0
            for i in range(N):
                mdot_igniter += mdot_uncontained_arr[i]
            # Adiabatic flame T (no plenum gas dynamics in uncontained model).
            T_ig = T_flame_pyro
            # P_ig diagnostic: volume-avg bore P over cartridge cells.
            p_sum = 0.0
            p_vol = 0.0
            for i in range(cart_i_start, cart_i_end + 1):
                if 0 <= i < N:
                    p_sum += P[i] * A_port[i] * dx
                    p_vol += A_port[i] * dx
            if p_vol > 0.0:
                P_ig = p_sum / p_vol
            else:
                P_ig = P[0]
        if P_ig > pyrogen_peak_P:
            pyrogen_peak_P = P_ig
        if mdot_igniter > 1e-12:
            pyrogen_duration = t + dt
        pyrogen_done = plenum_state[0] <= 1e-12 and mdot_igniter <= 1e-9
        pyrogen_momentum_expected = 0.0
        pyrogen_momentum_deposited = 0.0
        for i in range(N + 1):
            momentum_source[i] = 0.0
        # v0.7.3 Phase A: momentum injection only for forward_plenum.
        # Uncontained topologies let PISO handle axial flow via the
        # high-P pyrogen-cell pressure gradient (no explicit momentum).
        if (topology_code == 0 and
                not diagnostic_disable_momentum and mdot_igniter > 0.0 and
                igniter_axial_momentum_fraction > 0.0 and N > 1):
            v_exit = _orifice_exit_velocity(
                P_ig, T_ig, P[0], pyrogen_params_arr[5], pyrogen_params_arr[4]
            )
            pyrogen_momentum_expected = mdot_igniter * v_exit * igniter_axial_momentum_fraction
            face_area = 0.5 * (A_port[0] + A_port[1])
            if face_area > 1e-12 and v_exit > 0.0:
                momentum_source[1] = (
                    pyrogen_momentum_expected / (face_area * dx)
                )
                pyrogen_momentum_deposited = momentum_source[1] * face_area * dx

        active_pyrogen_surface_heat_flux_w_m2 = pyrogen_surface_heat_flux_w_m2
        if diagnostic_disable_pyrogen_surface_heating:
            active_pyrogen_surface_heat_flux_w_m2 = 0.0

        # v0.7.1 Phase 3.5: each combustion source uses its own species Cp.
        Cp_propellant_species = species_params_arr[_SPECIES_PROPELLANT, 1]
        Cp_pyrogen_species = species_params_arr[_SPECIES_IGNITER, 1]

        # v0.7.3 Phase B.4: build the per-cell pyrogen surface heat
        # flux array. Forward_plenum uses single-cell DeMar at the head
        # grain cell; uncontained topologies dispatch on
        # Pyrogen.heat_delivery_mode (DeMar / Radiation / None).
        # Replaces the v0.7.0-v0.7.3-phaseA "disabled for uncontained"
        # logic above with a mode-aware path that finally lets
        # uncontained topologies deliver heat to the propellant
        # surface.
        # Find the head grain cell (smallest i with grain + C_burn > 0
        # AND not yet ignited) for forward_plenum DeMar targeting.
        head_grain_cell = -1
        for i in range(N):
            if (is_grain[i] and C_burn[i] > 0.0 and not has_ignited[i]):
                head_grain_cell = i
                break
        _compute_pyrogen_heat_flux_arr(
            topology_code, heat_delivery_mode_code,
            is_grain, has_ignited, T_surf, C_burn, mdot_igniter, dx, N,
            pyrogen_heat_flux_arr,
            head_grain_cell, active_pyrogen_surface_heat_flux_w_m2,
            mdot_uncontained_arr, demar_flux_w_m2_uncontained,
            Cp_pyrogen_species, T_ig,
            cart_i_start, cart_i_end,
            T_flame_pyrogen, pellet_emissivity,
            radiation_absorption_length_m,
            x_centers, A_port,
        )

        # v0.7.2 Phase B-v2: flame-front augmentation. Always computed
        # (cheap, O(N)); the multiply is gated INSIDE the Goodman kernel
        # by flame_spread_enabled so disabled runs reduce to Phase A
        # byte-for-byte. The kernel fills augment_arr[i]=boost for cells
        # immediately downstream of a recently-ignited cell, else 1.0.
        _compute_flame_front_augment(
            is_burning, has_ignited, ignition_time, t,
            flame_spread_tau, flame_spread_boost, N, flame_spread_augment,
        )

        # v0.7.4 Phase F: advance the flame-spread front and refresh the
        # `ignitable` exposure mask. Skipped when disabled — `ignitable`
        # then stays all-True and is never read (the Goodman heating gate
        # short-circuits on flame_front_enabled=False).
        if flame_front_enabled:
            x_front = _advance_flame_front(
                is_burning, x_centers, flame_front_velocity,
                x_front, front_direction, front_seed_idx,
                cart_i_start, cart_i_end, dt, N, ignitable,
            )

        (n_burning, n_ignited, mass_sum,
         pyrogen_surface_heat_power, radiation_heat_power,
         radiation_sink_total_power,
         normal_sidewall_thermal_power, erosive_sidewall_thermal_power,
         endface_thermal_power,
         convective_wall_sink_power,
         fill_debit_massflow, fill_debit_power) = _goodman_ignition_sources_and_mass(
            P, T, T_surf, delta, has_ignited, is_burning, is_grain,
            ignition_time, r_total, r_erosive,
            mass_source, thermal_source,
            C_burn, endface_msource, pyrogen_surface_heat_flux,
            radiation_heat_flux, radiation_sink_power, radiation_emitter,
            x_centers, Re, D_hyd, f_darcy,
            t, dt, rho_propellant, T_flame, T_initial,
            Pr, k_thermal, roughness, kappa, solid_alpha, k_solid,
            T_ignition, N, dx, mdot_igniter, T_ig,
            Cp_propellant_species, Cp_pyrogen_species,
            active_pyrogen_surface_heat_flux_w_m2, radiation_emissivity,
            diagnostic_disable_radiation_gas_sink,
            tau_establishment,
            mass_source_by_species,
            # v0.7.2 Phase B-v2: flame-front h_c augmentation
            flame_spread_augment, flame_spread_enabled,
            # v0.7.3 Phase B.2: pyrogen-hot cells emit radiation
            Y_species,
            # v0.7.3 Phase B.4: pre-computed per-cell pyrogen heat flux
            pyrogen_heat_flux_arr,
            # v0.7.4 Phase F: flame-spread front exposure gate
            flame_front_enabled, ignitable, topology_code,
            # Solver-verification T10: start-of-step gas state for the fill debit
            rho, Cp_mix_arr,
            # IG2: opt-in ignition heat-transfer closure
            ignition_h_code, u_cell, A_port, A_throat,
            R_mix_arr, gamma_mix_arr, peretz_ap_code,
            solid_model_code, solid_y, solid_T, solid_age, solid_nact,
            solid_work, solid_age_fraction,
            ign_diag, rec_on, rec_h, rec_Td, rec_t, rec_dt, rec_k,
        )

        # v0.7.1 Phase 3.5: pyrogen mass injection uses pyrogen Cp for the
        # enthalpy flow, not the scalar Cp_gas (= propellant Cp) that
        # Phase 3 used as a behavior-preserving placeholder. The enthalpy
        # injected per unit length is mdot/dx * Cp_pyrogen * T_ig — the
        # cell back-outs T_new = h_new / Cp_arr[i], where Cp_arr is the
        # mixture Cp (PISO).
        #
        # v0.7.2 Phase A: pyrogen mass / species mass / enthalpy are
        # distributed across cells via pyrogen_axial_weights (computed
        # once at sim init from L_jet = kappa_jet * d_throat_pyrogen).
        # Sum(weights) = 1 → total mass + enthalpy preserved. Momentum
        # stays at face 1 (head-end aperture) — distributed momentum
        # deferred to v0.7.3+ per design doc rationale. Pyrogen surface
        # heat sink also stays at cell 0 (Goodman surface heating acts on
        # leading-edge unignited cell; sink clamp uses cell 0's
        # distributed enthalpy share to avoid over-subtraction).
        pyrogen_enthalpy_power = 0.0
        pyrogen_surface_heat_sink_power = 0.0
        if mdot_igniter > 0.0:
            pyrogen_enthalpy_power = mdot_igniter * Cp_pyrogen_species * T_ig
            # v0.7.3 Phase A: branch on topology for mass / enthalpy
            # deposition. forward_plenum uses Phase A axial weights;
            # uncontained topologies use per-cell mdot directly.
            if topology_code == 0:
                # Forward plenum: Phase A exponential axial decay
                for i in range(N):
                    w_i = pyrogen_axial_weights[i]
                    if w_i <= 0.0:
                        continue
                    mass_source_i = w_i * mdot_igniter / dx
                    mass_source[i] += mass_source_i
                    mass_source_by_species[i, _SPECIES_IGNITER] += mass_source_i
                    thermal_source[i] += mass_source_i * T_ig * Cp_pyrogen_species
                    mass_sum += mass_source_i
                # Surface heat sink clamped against cell 0's distributed
                # enthalpy share (forward_plenum DeMar pyrogen-plume
                # convention; not applicable to uncontained).
                ign_enthalpy_cell0 = (
                    pyrogen_axial_weights[0] * mdot_igniter
                    * Cp_pyrogen_species * T_ig / dx
                )
                pyrogen_surface_heat_sink = _pyrogen_surface_thermal_sink(
                    pyrogen_surface_heat_power, dx, ign_enthalpy_cell0
                )
                thermal_source[0] -= pyrogen_surface_heat_sink
                pyrogen_surface_heat_sink_power = pyrogen_surface_heat_sink * dx
            else:
                # Uncontained (head_basket / aft_basket): per-cell
                # mdot directly from the cartridge cells.
                for i in range(N):
                    mdot_cell = mdot_uncontained_arr[i]
                    if mdot_cell <= 0.0:
                        continue
                    mass_source_i = mdot_cell / dx
                    mass_source[i] += mass_source_i
                    mass_source_by_species[i, _SPECIES_IGNITER] += mass_source_i
                    thermal_source[i] += mass_source_i * T_ig * Cp_pyrogen_species
                    mass_sum += mass_source_i
                # v0.7.4 energy-balance fix (item 2a): the pyrogen radiation
                # delivered to the propellant walls (pyrogen_surface_heat_power,
                # accumulated in the Goodman kernel) is energy LEAVING the
                # pyrogen products. The full pyrogen enthalpy was injected as
                # gas above AND the pellets radiate to the walls — that double-
                # counts the pyrogen energy unless the radiated portion is
                # debited from the gas. Distribute the debit across the
                # cartridge (emitter) cells proportional to their pyrogen mdot
                # and clamp the total to the injected enthalpy so we never
                # extract more than the pyrogen added.
                if pyrogen_surface_heat_power > 0.0 and mdot_igniter > 0.0:
                    sink_total = pyrogen_surface_heat_power
                    if sink_total > pyrogen_enthalpy_power:
                        sink_total = pyrogen_enthalpy_power
                    for i in range(N):
                        mdot_cell = mdot_uncontained_arr[i]
                        if mdot_cell <= 0.0:
                            continue
                        share = mdot_cell / mdot_igniter
                        thermal_source[i] -= (sink_total * share) / dx
                    pyrogen_surface_heat_sink_power = sink_total

        # Refresh the source-aware CFL cap from THIS step's complete
        # thermal_source. The next iteration's dt will use this as an
        # upper bound (with a one-step lag, which is fine because
        # thermal_source magnitudes change smoothly during ignition
        # cascades on a per-step basis).
        # v0.7.1 Phase 3: per-cell Cp_arr replaces the scalar Cp_gas.
        if source_cfl_factor > 0.0:
            dt_source_cap = compute_dt_source_cap(
                rho, A_port, thermal_source, mass_source, N,
                Cp_mix_arr, T_flame, ambient_temperature, source_cfl_factor,
                SOURCE_DT_FLOOR,
            )

        # ============================================
        # STEP 3b: THROAT EVOLUTION
        # ============================================
        if throat_is_evolving:
            P_throat_MPa = P[N - 1] / 1e6
            e_rate = erosion_coeff * 1e-6 * P_throat_MPa
            if P_throat_MPa > 0.01:
                s_rate = slag_coeff / P_throat_MPa
            else:
                s_rate = 0.0
            D_throat = D_throat + 2.0 * (e_rate - s_rate) * dt
            if D_throat < 1e-6:
                D_throat = 1e-6
            A_throat = PI / 4.0 * D_throat * D_throat

        # v0.7.1: snapshot rho before PISO mutates it (needed for
        # mass-conservative Y advection after PISO).
        for i in range(N):
            rho_pre_step[i] = rho[i]

        # ============================================
        # STEP 4: PISO  (v0.7.1 Phase 3: per-cell γ/R/Cp/T_ceiling)
        # ============================================
        (rho, u, P, T,
         gas_energy_before, gas_energy_after, gas_sensible_dE_dt,
         convective_scalar_flux_power, nozzle_enthalpy_power,
         thermal_power_before_piso, clipping_correction_power,
         energy_residual,
         n_floor_cells, floor_mass_added, floor_energy_added,
         floor_min_raw_P,
         gas_mass_before, gas_mass_after,
         eos_mass_jump, eos_energy_jump,
         nozzle_mdot_piso,
         n_coupling_passes, coupling_rel_dT,
         eos_mass_jump_p2) = _piso_step_with_energy_diagnostics(
            rho, u, P, T, A_port, D_hyd,
            mass_source, thermal_source, momentum_source, f_darcy,
            dx, dt, gamma_mix_arr, R_mix_arr, Cp_mix_arr, T_ceiling_arr,
            A_throat, P_ambient, ambient_temperature, N,
            port_mach_cap,
            limit_activation_counts, limit_duration_s, limit_abs_energy_j,
            limit_first_time_s, limit_last_time_s,
            pressure_floor_max_deficit_pa,
            mach_limit_activation_counts, mach_limit_duration_s,
            mach_limit_first_time_s, mach_limit_last_time_s,
            t,
            pressure_floor_pa,
            piso_energy_coupling_iters, 1.0e-6,
            piso_corr2_storage,
        )

        # T6 energy audit (bookkeeping only): the gap between the END of the
        # previous step (ledger storage) and the START of this step's PISO is
        # complete only now, so its geometry part and total mass jump are
        # written back into the PREVIOUS step's history slot. The last step of
        # a run therefore keeps 0.0 here (no following step; see the audit).
        if hist_idx > 0:
            jump_geometry_J_hist[hist_idx - 1] = geo_dE
            jump_geometry_kg_hist[hist_idx - 1] = geo_dM
            jump_total_kg_hist[hist_idx - 1] = gas_mass_before - prev_gas_mass_after
        prev_gas_energy_after = gas_energy_after
        prev_gas_mass_after = gas_mass_after

        # ============================================
        # STEP 5: POST-PISO
        # ============================================
        a_max = _post_piso_update(
            rho, u, P, T, D_hyd, Re, Mach, u_cell, f_darcy,
            N, mu_gas, gamma_mix_arr, R_mix_arr, roughness,
        )

        # ============================================
        # STEP 6: BOOKKEEPING
        # ============================================
        total_mass_produced += mass_sum * dx * dt
        # v0.7.1 Phase 3: nozzle uses cell-N-1 mixture thermo.
        nozzle_mdot, _dmdp_nozzle, nozzle_upstream_T, _nozzle_state = _nozzle_boundary_flow(
            P[N - 1], T[N - 1], A_throat,
            gamma_mix_arr[N - 1], R_mix_arr[N - 1],
            P_ambient, ambient_temperature,
        )
        total_mass_nozzle += nozzle_mdot * dt
        # T6: post-PISO nozzle enthalpy power (final P, final T) for the
        # PISO-vs-reconstruction nozzle ledger comparison.
        nozzle_enthalpy_power_post = (
            nozzle_mdot * Cp_mix_arr[N - 1] * nozzle_upstream_T
        )

        # v0.7.1: species advection using post-PISO (u, rho) and the
        # pre-PISO rho snapshot. mass_source_by_species was populated
        # earlier (pyrogen at cell 0, propellant grain at active cells).
        # T6 energy audit (bookkeeping only): storage before/after the
        # mixture refresh (Cp changes at fixed rho, T, A_port).
        mix_mass_pre, mix_energy_pre = _gas_mass_energy(
            rho, T, A_port, dx, Cp_mix_arr, N,
        )
        _advect_species(
            Y_species, rho_pre_step, rho, u, A_port,
            nozzle_mdot, dx, dt,
            mass_source_by_species, N, S_species_local,
        )
        # v0.7.1 Phase 3: refresh per-cell (γ, Cp, R, M) for the NEXT
        # step's PISO + post-PISO consumption. T_ceiling_arr is the
        # max-of-species cap and is therefore time-invariant under the
        # current relaxed formula, but we refresh it here for symmetry
        # in case species_params later grows time-varying entries.
        _refresh_mixture_arrays(
            Y_species, species_params_arr,
            gamma_mix_arr, Cp_mix_arr, R_mix_arr, M_mix_arr, N,
        )
        mix_mass_post, mix_energy_post = _gas_mass_energy(
            rho, T, A_port, dx, Cp_mix_arr, N,
        )
        _compute_T_ceiling_arr(
            Y_species, species_params_arr, T_ceiling_arr, N, T_initial_gas,
        )

        # Kn = total bore burning area / throat area
        Kn = 0.0
        for i in range(N):
            Kn += C_burn[i]
        Kn = Kn * dx / A_throat

        # Burnthrough detection
        if first_bt_time < 0.0 and n_ignited > 0:
            for i in range(N):
                if is_grain[i] and D_port[i] >= cell_D_outer[i]:
                    first_bt_time = t
                    break

        # Record time history
        radiation_emitter_count = 0
        radiation_receiver_count = 0
        min_gas_temperature = T[0]
        max_gas_temperature = T[0]
        min_surface_temperature = T_surf[0]
        max_surface_temperature = T_surf[0]
        min_pressure = P[0]
        max_pressure = P[0]
        max_mach_abs = abs(Mach[0])
        for i in range(N):
            if radiation_emitter[i]:
                radiation_emitter_count += 1
            if radiation_heat_flux[i] > 0.0:
                radiation_receiver_count += 1
            if T[i] < min_gas_temperature:
                min_gas_temperature = T[i]
            if T[i] > max_gas_temperature:
                max_gas_temperature = T[i]
            if T_surf[i] < min_surface_temperature:
                min_surface_temperature = T_surf[i]
            if T_surf[i] > max_surface_temperature:
                max_surface_temperature = T_surf[i]
            if P[i] < min_pressure:
                min_pressure = P[i]
            if P[i] > max_pressure:
                max_pressure = P[i]
            mach_abs = abs(Mach[i])
            if mach_abs > max_mach_abs:
                max_mach_abs = mach_abs

        step_collapsed = (
            dt < COLLAPSE_DT_THRESHOLD
            or max_mach_abs > COLLAPSE_MACH_THRESHOLD
            or max_pressure > COLLAPSE_PRESSURE_THRESHOLD
        )
        if step_collapsed:
            collapse_consecutive += 1
        else:
            collapse_consecutive = 0

        time_hist[hist_idx] = t
        P_head_hist[hist_idx] = P[0]
        P_exit_hist[hist_idx] = P[N - 1]
        D_throat_hist[hist_idx] = D_throat
        Kn_hist[hist_idx] = Kn
        massflow_hist[hist_idx] = nozzle_mdot
        P_ig_hist[hist_idx] = P_ig
        T_ig_hist[hist_idx] = T_ig
        mdot_ig_hist[hist_idx] = mdot_igniter
        # IG4: the GRAIN's own gas-generation rate [kg/s]: total mass source
        # (grain + igniter) minus the igniter's share (the igniter deposits
        # exactly mdot_igniter / dx per unit length summed over cells).
        mdot_grain_hist[hist_idx] = max(mass_sum * dx - mdot_igniter, 0.0)
        m_pyrogen_hist[hist_idx] = plenum_state[0]
        gas_sensible_energy_before_hist[hist_idx] = gas_energy_before
        gas_sensible_energy_hist[hist_idx] = gas_energy_after
        gas_sensible_dE_dt_hist[hist_idx] = gas_sensible_dE_dt
        normal_sidewall_thermal_power_hist[hist_idx] = normal_sidewall_thermal_power
        erosive_sidewall_thermal_power_hist[hist_idx] = erosive_sidewall_thermal_power
        endface_thermal_power_hist[hist_idx] = endface_thermal_power
        convective_scalar_flux_power_hist[hist_idx] = convective_scalar_flux_power
        clipping_correction_power_hist[hist_idx] = clipping_correction_power
        pyrogen_enthalpy_power_hist[hist_idx] = pyrogen_enthalpy_power
        pyrogen_surface_heat_power_hist[hist_idx] = pyrogen_surface_heat_power
        gas_surface_heat_sink_power_hist[hist_idx] = pyrogen_surface_heat_sink_power
        radiation_heat_power_hist[hist_idx] = radiation_heat_power
        radiation_sink_power_hist[hist_idx] = radiation_sink_total_power
        nozzle_enthalpy_power_hist[hist_idx] = nozzle_enthalpy_power
        thermal_source_power_hist[hist_idx] = thermal_power_before_piso
        energy_residual_hist[hist_idx] = energy_residual
        pyrogen_momentum_expected_hist[hist_idx] = pyrogen_momentum_expected
        pyrogen_momentum_deposited_hist[hist_idx] = pyrogen_momentum_deposited
        pyrogen_momentum_residual_hist[hist_idx] = pyrogen_momentum_expected - pyrogen_momentum_deposited
        dt_hist[hist_idx] = dt
        n_burning_hist[hist_idx] = n_burning
        n_ignited_hist[hist_idx] = n_ignited
        radiation_emitter_count_hist[hist_idx] = radiation_emitter_count
        radiation_receiver_count_hist[hist_idx] = radiation_receiver_count
        min_gas_temperature_hist[hist_idx] = min_gas_temperature
        max_gas_temperature_hist[hist_idx] = max_gas_temperature
        min_surface_temperature_hist[hist_idx] = min_surface_temperature
        max_surface_temperature_hist[hist_idx] = max_surface_temperature
        min_pressure_hist[hist_idx] = min_pressure
        max_pressure_hist[hist_idx] = max_pressure
        max_mach_hist[hist_idx] = max_mach_abs
        pressure_floor_cells_hist[hist_idx] = n_floor_cells
        pressure_floor_mass_hist[hist_idx] = floor_mass_added
        # T6 energy audit histories. Jump slot n = the gap between the end of
        # step n and the start of step n+1's PISO; the geometry part and the
        # total mass jump are filled in by step n+1 (0.0 until then).
        pressure_floor_energy_hist[hist_idx] = floor_energy_added
        convective_wall_sink_power_hist[hist_idx] = convective_wall_sink_power
        fill_debit_massflow_hist[hist_idx] = fill_debit_massflow
        fill_debit_power_hist[hist_idx] = fill_debit_power
        nozzle_enthalpy_power_post_hist[hist_idx] = nozzle_enthalpy_power_post
        nozzle_massflow_piso_hist[hist_idx] = nozzle_mdot_piso
        coupling_passes_hist[hist_idx] = n_coupling_passes
        coupling_rel_dT_hist[hist_idx] = coupling_rel_dT
        eos_jump_p2_kg_hist[hist_idx] = eos_mass_jump_p2
        jump_eos_J_hist[hist_idx] = eos_energy_jump
        jump_mixture_J_hist[hist_idx] = mix_energy_post - mix_energy_pre
        jump_geometry_J_hist[hist_idx] = 0.0
        jump_eos_kg_hist[hist_idx] = eos_mass_jump
        jump_mixture_kg_hist[hist_idx] = mix_mass_post - mix_mass_pre
        jump_geometry_kg_hist[hist_idx] = 0.0
        jump_total_kg_hist[hist_idx] = 0.0
        if n_floor_cells > 0:
            floor_energy_total += floor_energy_added
            if floor_min_raw_P < floor_min_raw_P_global:
                floor_min_raw_P_global = floor_min_raw_P
        hist_idx += 1

        # Snapshot
        if t - last_snapshot_t >= snapshot_interval and snap_idx < max_snaps:
            snap_times[snap_idx] = t
            for i in range(N):
                snap_data[snap_idx, _SNAP_P, i] = P[i]
                snap_data[snap_idx, _SNAP_U, i] = u_cell[i]
                snap_data[snap_idx, _SNAP_MACH, i] = Mach[i]
                snap_data[snap_idx, _SNAP_T, i] = T[i]
                snap_data[snap_idx, _SNAP_R_TOTAL, i] = r_total[i]
                snap_data[snap_idx, _SNAP_R_EROSIVE, i] = r_erosive[i]
                snap_data[snap_idx, _SNAP_D_PORT, i] = D_port[i]
                snap_data[snap_idx, _SNAP_C_BURN, i] = C_burn[i]
                snap_data[snap_idx, _SNAP_ENDFACE, i] = endface_msource[i]
                snap_data[snap_idx, _SNAP_IS_BURNING, i] = 1.0 if is_burning[i] else 0.0
                snap_data[snap_idx, _SNAP_IS_GRAIN, i] = 1.0 if is_grain[i] else 0.0
                snap_data[snap_idx, _SNAP_T_SURF, i] = T_surf[i]
                snap_data[snap_idx, _SNAP_MASS_SOURCE, i] = mass_source[i]
                snap_data[snap_idx, _SNAP_THERMAL_SOURCE, i] = thermal_source[i]
                snap_data[snap_idx, _SNAP_MOMENTUM_SOURCE, i] = momentum_source[i + 1]
                snap_data[snap_idx, _SNAP_PYROGEN_SURFACE_HEAT_FLUX, i] = pyrogen_surface_heat_flux[i]
                snap_data[snap_idx, _SNAP_RADIATION_HEAT_FLUX, i] = radiation_heat_flux[i]
                snap_data[snap_idx, _SNAP_REGRESS, i] = regress[i]
                snap_data[snap_idx, _SNAP_RHO, i] = rho[i]
            # Per-segment end-face regression (axial face burnback) — small
            # N_seg arrays kept parallel to the per-cell snapshot.
            for k in range(N_seg):
                snap_seg_fwd[snap_idx, k] = seg_fwd_regression[k]
                snap_seg_aft[snap_idx, k] = seg_aft_regression[k]
            snap_idx += 1
            last_snapshot_t = t

        # Pressure cutoff (only after pyrogen is consumed and vented)
        if n_ignited > 0 and pyrogen_done and P[0] < P_cutoff:
            termination_code = 2
            break

        # IG3: ignition-impossible early stop (a physical bound, not a tuned
        # timeout). No grain cell has lit, and the pyrogen charge is fully
        # consumed with no mass leaving the igniter (pyrogen_done), so no
        # energy source is left -- the propellant gives none until it lights --
        # and the gas in EVERY grain cell is colder than T_ignition. A wall
        # heated by convection (the only channel left: nothing burns, so no
        # radiation) can approach but never exceed the temperature of the gas
        # heating it, and that gas cannot get hotter without a source, so no
        # unignited surface can reach T_ignition. Residual assumption: gas in
        # the non-grain head-end cells (the cooled-by-nothing igniter pocket,
        # which stays at plenum temperature in this model) is not advected
        # into the grain; the grain gas keeps cooling beside it (measured).
        if (stop_if_unignitable and n_ignited == 0 and pyrogen_done
                and step % 16 == 0):
            T_gas_max = 0.0
            for i in range(N):
                if is_grain[i] and T[i] > T_gas_max:
                    T_gas_max = T[i]
            if T_gas_max < T_ignition:
                termination_code = 6
                break

        # Numerical-collapse trip: classified abort to avoid burning the
        # full history budget on a doomed run. Records the collapse step
        # in history (so the user sees the trip in diagnostics) and exits.
        if collapse_consecutive >= COLLAPSE_CONSECUTIVE_STEPS:
            termination_code = 4
            break

        t += dt
        step += 1

    # Publish terminal progress (full unless the run was canceled) so the GUI
    # bar lands at 100% on a normal finish rather than wherever the metric was.
    if termination_code != 5:
        progress_state[0] = 1.0

    return (hist_idx, snap_idx,
            total_mass_produced, total_mass_nozzle,
            first_bt_time, D_throat, termination_code,
            pyrogen_initial_mass - plenum_state[0],
            pyrogen_duration, pyrogen_peak_P,
            floor_energy_total, floor_min_raw_P_global,
            max_burn_refresh_gap, max_geometry_refresh_gap)


# ================================================================
# Public API
# ================================================================

# History pre-allocation. The per-step history arrays are sized up front from
# an estimated step count; a run that outgrows them stops early with
# termination "history array full". The time loop's acoustic CFL step is
# dt = cfl * dx / max(|u| + a), so the estimate bounds the wave speed by the
# hottest species' sound speed times a margin for the flow velocity (plateau
# port Mach ~0.3-0.4; the supersonic ignition fill lasts a few ms and adds few
# steps). The former fixed 1000 m/s under-estimated hot-gas sound speed
# (~1100-1200 m/s) and truncated fine-grid runs at t ~ 0.83 * t_max
# (solver-verification T8, 2026-10-06).
_HISTORY_WAVE_SPEED_MARGIN = 1.5


def _shrink_hist(a, n):
    """MEM1: trim a pre-allocated history array to exactly ``n`` rows IN PLACE
    (``ndarray.resize``), instead of ``a[:n].copy()`` which held the full buffer
    and the copy at once (a 2x peak across ~60 arrays). Values are untouched."""
    a.resize(int(n), refcheck=False)
    return a


def _estimate_history_steps(t_max, cfl_target, dx, dt_max, species_params_arr):
    """Estimated flow-step count for a run (history pre-allocation only)."""
    a_hot = 0.0
    for s in range(species_params_arr.shape[0]):
        gamma, _cp, molar_mass, T_flame = species_params_arr[s]
        a_hot = max(a_hot, (gamma * R_UNIVERSAL / molar_mass * T_flame) ** 0.5)
    dt_est = cfl_target * dx / (_HISTORY_WAVE_SPEED_MARGIN * a_hot)
    if dt_max > 0.0:
        dt_est = min(dt_est, dt_max)
    return int(t_max / max(dt_est, 1e-12)) + 1000


def _build_numerical_limits(
    x_centers, T_ambient, T_ceiling_arr, pressure_floor_pa, port_mach_cap,
    counts, duration, abs_energy, first_t, last_t, floor_deficit,
    mach_counts, mach_duration, mach_first_t, mach_last_t,
):
    """Assemble ``result['numerical_limits']`` from the loop's tally arrays
    (copies; -1.0 first/last sentinels become NaN)."""
    def _times(a):
        return np.where(a >= 0.0, a, np.nan)

    def _block(row):
        return {
            'activation_count_by_cell': counts[row].copy(),
            'duration_s_by_cell': duration[row].copy(),
            'first_activation_time_s_by_cell': _times(first_t[row]),
            'last_activation_time_s_by_cell': _times(last_t[row]),
        }

    tf = _block(LIMIT_TEMPERATURE_FLOOR)
    tf['threshold_k'] = max(1.0, float(T_ambient))
    tf['absolute_correction_energy_j_by_cell'] = (
        abs_energy[LIMIT_TEMPERATURE_FLOOR].copy())
    tc = _block(LIMIT_TEMPERATURE_CEILING)
    tc['final_threshold_k_by_cell'] = T_ceiling_arr.copy()
    tc['absolute_correction_energy_j_by_cell'] = (
        abs_energy[LIMIT_TEMPERATURE_CEILING].copy())
    pf = _block(LIMIT_PRESSURE_FLOOR)
    pf['threshold_pa'] = float(pressure_floor_pa)
    pf['maximum_raw_deficit_pa_by_cell'] = floor_deficit.copy()
    x_face = 0.5 * (x_centers[:-1] + x_centers[1:])
    return {
        'x_m': x_centers.copy(),
        'temperature_floor': tf,
        'temperature_ceiling': tc,
        'pressure_floor': pf,
        'port_mach_cap': {
            'threshold_mach': float(port_mach_cap),
            'x_m': x_face,
            'activation_count_by_face': mach_counts.copy(),
            'duration_s_by_face': mach_duration.copy(),
            'first_activation_time_s_by_face': _times(mach_first_t),
            'last_activation_time_s_by_face': _times(mach_last_t),
        },
    }


def _run_simulation_once(
    geo,
    propellant,
    nozzle,
    pyrogen_chamber,
    # --- Environment ---
    P_ambient=101325.0,
    ambient_temperature=None,
    # --- Surface / erosive burning ---
    # v0.7.5 cross-motor re-LHS optimum (fired-motor set; docs/v0_7_5/RESULT.md):
    # roughness 32um, kappa 0.44 (was 50um / 0.45). All physical.
    roughness=32e-6,
    kappa=0.44,
    # Spike-reopen probe (docs/v0_7_4/SPIKE_REOPEN_BRIEF.md): opt-in ceiling
    # on the Ma erosive enhancement at `erosive_mp_cap × eta_MP`, the
    # Mukunda-Paul (1997) universal law (published constants, no free
    # parameters). 0.0 = off (default). 1.0 = the law itself; 1.1 = MP's
    # published ±10% accuracy band, binding only where Ma exceeds what the
    # MP data cloud can support. Ma and MP agree within ~6% at plateau g;
    # Ma exceeds MP ~25-35% at ignition-transient high g — the measured
    # high-L/D over-spike factor. Default OFF: enabling it is a deliberate
    # model-fidelity tradeoff (second closure bounding Ma), a user-visible
    # decision, not a shipped default.
    erosive_mp_cap=0.0,
    # Erosive-closure registry (docs/erosive_model_registry/DESIGN.md):
    # selects the erosive-burning closure. 'ma2020' (default, no-tuning
    # shipping model) — the Ma et al. (2020) heat-transfer bisection
    # above, optionally bounded by erosive_mp_cap. 'mukunda_paul' — the
    # Mukunda-Paul (1997) Eq. 12 universal law, closed form; roughness/
    # kappa stop influencing the burn rate under this model (they keep
    # their unrelated ignition-gate h_c role in simulation.py's Goodman
    # heating, unchanged). 'saint_robert' — erosive burning off (r_total
    # = r0). Model choice is solver config, not propellant data. Never
    # "recommended" — an explicit user tradeoff, not a shipped default.
    erosive_model='ma2020',
    # IG2/IG3 (docs/pyrogen_realism/IG2_REPORT.md, IG3_TASK.md): closure for
    # the convective h that heats the UNIGNITED propellant surface in the
    # Goodman ignition kernel. 'peretz' (default since IG3) = Peretz et al.
    # (1973) bore correlation for starting transients. 'gnielinski' = the
    # Ma/Gnielinski chain used for erosive burning (the pre-IG3 default).
    # 'keller' = Keller-Baer-Ryan (1966) Eq. 11 extrapolated to the bore. The
    # Ma erosive path is not affected by this switch.
    ignition_heat_transfer='peretz',
    # IG3: port area used in Peretz's (A_p/A_t)^0.4 factor ('peretz' mode
    # only): 'local' (default, the cell's own port area), 'head' (fore-most
    # grain cell) or 'aft' (aft-most grain cell); instantaneous values.
    peretz_ap_reference='local',
    # SM1 (docs/pyrogen_realism/SM1_REPORT.md): solid-phase heating model of
    # the ignition gate. 'goodman' (default; algebraic-in-h cubic heat-balance
    # integral, bit-identical to pre-SM1), 'peretz' (Peretz's integral-method
    # T_s ODE, III-50: T_s is carried, not recomputed from h) or 'conduction'
    # (exact 1-D transient conduction in the solid, graded grid, implicit).
    solid_heating_model='goodman',
    # SM2: record the per-cell surface forcing (h_total, T_driver per step) of the first
    # n steps into result['forcing'] for an offline open-loop replay of the solid
    # models (0 = off). The per-cell ignition-event diagnostics (result['ignition_diag'])
    # are always recorded. Neither changes the trajectory.
    record_forcing_steps=0,
    # --- Solver ---
    # v0.7.3.2 (2026-05-27): cfl_target tightened 0.5 → 0.3 and
    # source_cfl_factor tightened 0.10 → 0.05 to absorb Phase B.0
    # cold-bore IC transient dynamics. Sutton-default Hasegawa A was
    # silently collapsing under the looser 0.5/0.10 settings because
    # the pyrogen-plenum-into-cold-cell-0 transient happens faster
    # than those dt caps allow. The tighter settings match the
    # calibrated test config (test_yns_phase4_validation
    # _short_hasegawa_a_run) which has been stable since v0.7.1.
    cfl_target=0.3,
    dt_max=0.002,
    source_cfl_factor=0.05,
    # Aerodynamic-choking velocity limiter: cap interior port face velocities
    # at this Mach (~local sound speed). Bounds the grid-divergent supersonic
    # fill artifact that the pressure-based PISO produces at the cold/hot
    # contact during ignition (no interior choking limit, not shock-capturing;
    # docs/v0_7_4/IGNITION_SPIKE_REOPENED.md §6). 0.0 disables (default; prior
    # behaviour). ~1.0 is the physical choking limit. No-op at the subsonic
    # plateau.
    port_mach_cap=0.0,
    # T11 (default ON since the T11 adoption, 2026-10-11): energy-coupled PISO
    # outer iteration. 0 = the pre-T11 step. K > 0 = up to K extra pressure-
    # correction + energy passes with T-updated compressibility
    # (solver._piso_step_with_energy_diagnostics docstring). Default K = 1.
    piso_energy_coupling_iters=1,
    # T11 (default ON): the storage term the pre-T11 second pressure
    # corrector omitted (it re-applied the first corrector's storage, so every
    # transient pressure increment was applied ~twice). False = pre-T11 step.
    piso_corr2_storage=True,
    burn_update_interval=None,
    geometry_update_interval=None,
    # --- Ignition ---
    T_ignition=756.0,  # v0.7.5 cross-motor re-LHS (was 850); docs/v0_7_5/RESULT.md
    tau_establishment=0.0,
    # --- Diagnostics ---
    initial_gas_temperature=None,
    diagnostic_disable_erosive=False,
    diagnostic_disable_endfaces=False,
    diagnostic_disable_momentum=False,
    diagnostic_disable_pyrogen_surface_heating=False,
    diagnostic_disable_adjacent_radiation=False,
    diagnostic_disable_radiation_gas_sink=False,
    history_capacity=None,
    igniter_axial_momentum_fraction=1.0,
    # --- Termination ---
    t_max=10.0,
    P_cutoff=0.5e6,
    # --- Output ---
    print_interval=0.2,
    snapshot_interval=0.2,
    verbose=True,
    # --- v0.8.0 Phase 6: live progress / cooperative cancel ---
    progress_state=None,
    # Limiter-sensitivity knob (solver-verification PB/P1). The PISO step
    # clamps the raw post-correction pressure to this floor [Pa]. The
    # default 1.0e3 is the historical constant and is bit-identical to it;
    # any other value CHANGES the physics of the run (diagnostic only).
    # Kept last for positional API compatibility.
    pressure_floor_pa=1.0e3,
    # IG3: end the run (termination_code 6) as soon as ignition is physically
    # impossible (see the check in _run_time_loop). Never affects a motor that
    # lights; only shortens a run that would otherwise idle to t_max.
    stop_if_unignitable=True,
):
    """
    Run a complete transient simulation (one attempt; see ``run_simulation``
    for the auto-sized-charge retry wrapper).

    Parameters
    ----------
    geo : MotorGeometry
        Motor geometry (grain side) from grain_geometry.py.
    propellant : Propellant
        Propellant properties from propellant.py.
    nozzle : Nozzle
        Nozzle configuration (D_throat, D_exit, efficiency, div_angle,
        conv_angle, throat_length, erosion_coeff, slag_coeff). See
        srm_1d.nozzle.Nozzle.
    pyrogen_chamber : PyrogenChamber
        Hot-gas pyrogen igniter chamber configuration.
    P_ambient : float
        Ambient pressure [Pa]. Default sea level (101325). Sim-environment
        config; openMotor calls this `motor.config.ambPressure`.
    ambient_temperature : float or None
        Ambient reservoir temperature [K] for reverse nozzle inflow.
        ``None`` uses ``propellant.T_initial``.
    roughness : float
        Propellant surface roughness [m]. Default: 50 μm.
    kappa : float
        Gnielinski temperature-ratio exponent [-]. Default: 0.45.
    erosive_mp_cap : float
        Opt-in (default 0.0 = off): ceiling the Ma erosive enhancement
        eta = r_total/r0 at `erosive_mp_cap × eta_MP`, the Mukunda-Paul
        (1997) Eq. 12 universal law (published constants; no free
        parameters). 1.0 = the law itself; 1.1 = MP's published ±10%
        accuracy band. Spike-reopen probe — cuts the high-L/D
        ignition-transient over-spike where Ma exceeds MP at high g.
    erosive_model : str
        Erosive-closure registry key. Default ``'ma2020'`` (the Ma et
        al. (2020) heat-transfer bisection above — no-tuning shipping
        default). ``'mukunda_paul'`` selects the closed-form
        Mukunda-Paul (1997) Eq. 12 universal law instead (roughness/
        kappa then stop influencing the burn rate — they keep their
        ignition-gate role elsewhere). ``'saint_robert'`` disables
        erosive burning entirely (r_total = r0). ``erosive_mp_cap`` is a
        ``ma2020``-only option; setting it nonzero with any other
        ``erosive_model`` raises ``ValueError``. See
        docs/erosive_model_registry/DESIGN.md.
    ignition_heat_transfer : str
        Closure for the convective heat-transfer coefficient that heats the
        unignited propellant surface toward ``T_ignition`` (and the matching
        gas wall-loss sink). Default ``'peretz'`` (IG3; was ``'gnielinski'``
        before): Peretz et al. 1973 Eq. III-22 (film-temperature properties,
        A_p/A_t and entrance factors). ``'gnielinski'``: the Ma/Gnielinski
        erosive-burning chain. ``'keller'``: Keller-Baer-Ryan 1966 Eq. 11 with
        a (D/D_K)^-0.095 bore scaling. The Ma 2020 erosive burn-rate path is
        never affected.
    peretz_ap_reference : str
        Port area in Peretz's (A_p/A_t)^0.4 factor, ``'peretz'`` mode only:
        ``'local'`` (default; the cell's own port area, bit-identical to IG2),
        ``'head'`` (port area of the fore-most grain cell / throat area) or
        ``'aft'`` (aft-most grain cell). Instantaneous values (they follow the
        burnback and the throat). The sources only had uniform ports, so none
        of the three is stated for a non-uniform bore (R9_SOURCES.md).
    solid_heating_model : str
        How the unignited propellant surface temperature is advanced toward
        ``T_ignition``. ``'goodman'`` (default): Goodman cubic heat-balance
        integral with an algebraic T_s(h, delta) -- a drop in h instantly
        lowers T_s. ``'peretz'``: Peretz's T_s ODE (thesis Eq. III-50), same
        profile, T_s carried as state (held when h falls); identical to
        goodman at constant h. ``'conduction'``: exact 1-D transient
        conduction in the solid (graded grid of ``solid_thermal.
        CONDUCTION_DY0/CONDUCTION_RATIO``, backward Euler sub-stepped against
        the heating age; uses only k_solid and alpha). All three see the same
        combined surface flux (convective + pyrogen + radiation). T_ignition /
        k_solid were calibrated with goodman. See SM1_REPORT.md.
    record_forcing_steps : int
        SM2. If > 0, record the per-cell surface forcing (``h_total``, ``T_driver``; ``-1`` =
        not a heated grain cell) of the first n steps into ``result['forcing']`` (keys
        ``h``, ``T_drv``, ``t``, ``dt``) so the solid models can be replayed offline on
        identical forcing. Already-lit cells record the bare convective h of their local gas
        (the forcing they would see unlit). ``result['ignition_diag']`` (always) holds, per cell
        at its ignition step, [h_total, T_driver, q_pre, T_s, T_gas, u, Mach, P, t, q_post].
        Neither changes the trajectory.
    cfl_target : float
        Target CFL number for time stepping. Default: 0.5.
    dt_max : float
        Maximum allowed time step [s]. Default: 0.002.
    burn_update_interval : int or None
        Recompute burn rates (``compute_burn_rates``) every N flow steps.
        If None, auto-set to ``max(10, N_cells // 5)``.
    geometry_update_interval : int or None
        Recompute the per-cell port geometry (``update_cell_geometry``:
        ``A_port``, ``C_burn``, ``D_hyd``, end-face mass source) every N flow
        steps. If None (default), uses the resolved ``burn_update_interval``,
        i.e. one shared cadence exactly as before this kwarg existed.
        Regression itself advances every step regardless. The summary reports
        both resolved intervals and the largest simulated-time gap between
        consecutive refreshes of each kind (``max_burn_refresh_gap_s`` /
        ``max_geometry_refresh_gap_s``).
    T_ignition : float
        Per-cell solid surface ignition threshold [K]. Default: 850 K.
    initial_gas_temperature : float or None
        Optional diagnostic override for the initial bore gas temperature
        [K]. ``None`` preserves the historical behavior: fill the bore
        gas at propellant flame temperature.
    diagnostic_disable_erosive : bool
        If True, remove the Ma erosive increment from burn rates while
        preserving the Saint-Robert normal rate. Diagnostic only.
    diagnostic_disable_endfaces : bool
        If True, suppress end-face regression and end-face mass source
        terms. Diagnostic only.
    diagnostic_disable_momentum : bool
        If True, suppress pyrogen axial momentum injection while
        preserving igniter mass and enthalpy. Diagnostic only.
    diagnostic_disable_pyrogen_surface_heating : bool
        If True, suppress direct pyrogen-to-propellant Goodman heating
        while preserving igniter mass, enthalpy, and momentum. Diagnostic
        only.
    diagnostic_disable_adjacent_radiation : bool
        If True, suppress adjacent-burning-cell radiation while preserving
        the material emissivity setting in the result summary.
    diagnostic_disable_radiation_gas_sink : bool
        If True, keep adjacent-radiation Goodman receiver heating but do
        not debit the emitting gas cell. Diagnostic isolation only.
    pressure_floor_pa : float
        Pressure floor [Pa] applied by the PISO step (default 1.0e3 = the
        historical constant, bit-identical). Overriding it is a
        limiter-sensitivity diagnostic and changes the run's physics. Must
        be finite and positive.
    history_capacity : int or None
        Number of preallocated per-step history rows. ``None`` (default)
        sizes them from the hottest species' sound speed (see
        ``_estimate_history_steps``), with a 5,000,000-row floor. An integer
        sets the capacity exactly — raise it for very fine grids / low CFL /
        long runs, or lower it to make a probe stop early. It never changes
        the equations or the time stepping; a run that fills its history
        stops with termination "history array full" and a warning.
    igniter_axial_momentum_fraction : float
        Fraction of pyrogen orifice momentum projected downstream into
        the bore. Default 1.0 represents a head-end axial jet.
    t_max : float
        Maximum simulation time [s]. Default: 10.
    P_cutoff : float
        Head-end pressure below which simulation terminates [Pa].
    print_interval : float
        Print status every this many seconds of simulated time.
    snapshot_interval : float
        Store flow field snapshots at this interval [s].
    verbose : bool
        If True, print setup and summary blocks. Set False for large sweeps.

    Returns
    -------
    dict with 'time', 'P_head', 'P_exit', 'D_throat', 'snapshots', 'summary'.
    """
    if pyrogen_chamber is None:
        raise ValueError("pyrogen_chamber is required for v0.7.0 ignition")
    if ambient_temperature is not None and ambient_temperature <= 0.0:
        raise ValueError("ambient_temperature must be positive")
    if initial_gas_temperature is not None and initial_gas_temperature <= 0.0:
        raise ValueError("initial_gas_temperature must be positive")
    if not 0.0 <= igniter_axial_momentum_fraction <= 1.0:
        raise ValueError("igniter_axial_momentum_fraction must be between 0 and 1")
    if not np.isfinite(pressure_floor_pa) or pressure_floor_pa <= 0.0:
        raise ValueError("pressure_floor_pa must be finite and positive")
    if burn_update_interval is not None and (
            int(burn_update_interval) != burn_update_interval
            or burn_update_interval < 1):
        raise ValueError(
            "burn_update_interval must be None or a positive integer; "
            f"got {burn_update_interval!r}")
    if geometry_update_interval is not None and (
            int(geometry_update_interval) != geometry_update_interval
            or geometry_update_interval < 1):
        raise ValueError(
            "geometry_update_interval must be None or a positive integer; "
            f"got {geometry_update_interval!r}")
    # Erosive-closure registry resolution + validation (before any
    # sim setup, per docs/erosive_model_registry/DESIGN.md §3): unknown
    # name -> ValueError listing registered names; erosive_mp_cap is a
    # ma2020-only option, so a nonzero cap under another model fails
    # loudly rather than being silently ignored.
    _ero_code = resolve_erosive_model(erosive_model)
    _ign_h_code = _ignition_h_code(ignition_heat_transfer)
    _ap_code = _peretz_ap_code(peretz_ap_reference)
    _solid_code = _solid_model_code(solid_heating_model)
    if erosive_mp_cap > 0.0 and erosive_model != 'ma2020':
        raise ValueError(
            f"erosive_mp_cap={erosive_mp_cap} is a ma2020-only option; "
            f"got erosive_model={erosive_model!r}. Set erosive_mp_cap=0.0 "
            "or erosive_model='ma2020'."
        )
    pyrogen_heat_flux = pyrogen_chamber.pyrogen.heat_flux_cal_cm2_s
    if diagnostic_disable_pyrogen_surface_heating:
        pyrogen_surface_heat_flux_w_m2 = 0.0
    else:
        if pyrogen_heat_flux is None or pyrogen_heat_flux <= 0.0:
            raise ValueError(
                "pyrogen.heat_flux_cal_cm2_s must be positive when pyrogen "
                "surface heating is enabled; add heat_flux_cal_cm2_s to the "
                "custom pyrogen YAML or set "
                "diagnostic_disable_pyrogen_surface_heating=True"
            )
        pyrogen_surface_heat_flux_w_m2 = _cal_cm2_s_to_w_m2(float(pyrogen_heat_flux))
    active_radiation_emissivity = float(propellant.radiation_emissivity)
    if diagnostic_disable_adjacent_radiation:
        active_radiation_emissivity = 0.0

    erosion_coeff = nozzle.erosion_coeff
    slag_coeff = nozzle.slag_coeff

    # Pick the representative tab for sim-start gas thermo. With multi-tab
    # propellants, gamma/T_flame/MW are frozen at sim start (v0.3.x scope);
    # only a/n vary in the hot loop. See propellant.py docstring TODO note.
    rep_tab = propellant.representative_tab()
    tab_min_p, tab_max_p, tab_a, tab_n = propellant.tab_arrays()
    n_tabs = len(propellant.tabs)
    N = geo.N_cells
    dx = geo.dx

    # v0.7.1: N-species bore-gas registry. Indices:
    #   s=0  igniter (pyrogen)
    #   s=1  main grain combustion products
    #   s=2  ambient (pre-fill, no continuing source)
    # Higher indices reserved for future sources (v0.8.0 head-end motor,
    # ablation, ...). See docs/v0_7_1/DESIGN.md.
    from .propellant import (
        species_array as _build_species_array,
        ambient_air_species as _ambient_air_species,
    )
    _ambient_T = float(ambient_temperature) if ambient_temperature is not None \
                 else propellant.T_initial
    species_list = [
        pyrogen_chamber.pyrogen.species,           # s=0
        propellant.species(),                      # s=1
        _ambient_air_species(T_ambient=_ambient_T),# s=2
    ]
    species_params_arr = _build_species_array(species_list)
    S_species = species_params_arr.shape[0]
    SPECIES_IGNITER = 0
    SPECIES_PROPELLANT = 1
    SPECIES_AMBIENT = 2

    if burn_update_interval is None:
        burn_update_interval = max(10, N // 5)
    # Solver-verification T5: geometry and burn rates have separate cadences.
    # None -> share the (resolved) burn cadence, i.e. the pre-T5 behaviour.
    if geometry_update_interval is None:
        geometry_update_interval = burn_update_interval
    else:
        geometry_update_interval = int(geometry_update_interval)

    # ============================================================
    # SETUP — extract everything into scalars and arrays
    # ============================================================
    gas = create_gas_properties(
        rep_tab.gamma, rep_tab.molecular_weight, rep_tab.T_flame,
        propellant.mu_gas, propellant.k_gas, propellant.Cp_gas,
    )
    Gamma_crit = critical_flow_function(gas.gamma)
    gamma_R = gas.gamma * gas.R_specific
    nozzle_denom = 1.0 / (gas.R_specific * rep_tab.T_flame) ** 0.5
    T_ambient = propellant.T_initial
    if ambient_temperature is not None:
        T_ambient = float(ambient_temperature)

    D_throat_init = nozzle.D_throat
    A_throat_init = np.pi / 4.0 * D_throat_init ** 2
    throat_is_evolving = (erosion_coeff != 0.0 or slag_coeff != 0.0)

    # Compile geometry to arrays
    ga = geo.compile_geometry_arrays()
    # Capture the INITIAL cell->grain map (t=0, as-designed grain layout)
    # before the time loop mutates ga via end-face regression. Station-viz
    # classification (fore/mid/aft, Head/Grain/Gap/Aft) keys off this so the
    # markers sit on the original grain, not the final receded extent.
    cell_segment_id_init = ga['cell_segment_id'].copy()
    D_port = ga['D_port']
    x_centers = ga['x_centers']
    regress = ga['regress']

    # Working arrays
    A_port = np.zeros(N)
    C_burn = np.zeros(N)
    D_hyd = np.zeros(N)
    is_grain = np.zeros(N, dtype=np.bool_)
    endface_msource = np.zeros(N)

    # Initialize geometry
    update_cell_geometry(
        regress, D_port, x_centers, dx, N, ga['N_seg'], ga['cell_D_outer'],
        ga['seg_x_start'], ga['seg_length'],
        ga['seg_fwd_regression'], ga['seg_aft_regression'],
        ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'],
        ga['cell_segment_id'],
        np.full(N, 1e6),
        propellant.rho_propellant,
        tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        A_port, C_burn, D_hyd, is_grain, endface_msource,
        ga['cell_D_bore_init'], ga['cell_wall_web'],
        ga['cell_segment_type'], ga['cell_fmm_idx'],
        ga['fmm_offset'], ga['fmm_reg_flat'],
        ga['fmm_perim_flat'], ga['fmm_port_flat'],
    )

    # Flow state
    P = np.full(N, P_ambient)
    # v0.7.3 Phase B.0 (2026-05-25): default bore IC switched from
    # rep_tab.T_flame (v0.7.0 numerical-stability shortcut) to _ambient_T
    # (= propellant.T_initial unless overridden via ambient_temperature
    # kwarg). The previous IC short-circuited temperature-gradient flow
    # under uncontained-pyrogen topologies (the bore was already at flame T
    # so pyrogen mass injection at T_flame_pyrogen created no T gradient,
    # no density gradient, no pressure gradient → no flow). Override with
    # initial_gas_temperature= for backward compat or special studies.
    T_initial_gas = _ambient_T
    if initial_gas_temperature is not None:
        T_initial_gas = float(initial_gas_temperature)
    rho = P / (gas.R_specific * T_initial_gas)
    u = np.zeros(N + 1)
    T = np.full(N, T_initial_gas)

    # v0.7.1: per-cell species mass fractions. Y_species[i, s] is the
    # mass fraction of species s in cell i. Invariant: sum_s Y[i, s] = 1.
    # Initial condition: all cells start as 100% ambient (species 2).
    # v0.7.3 Phase B.0: this is now self-consistent — Y=ambient at
    # T=ambient. Previously Y=ambient at T=T_flame_propellant (mildly
    # unphysical "cold air composition at hot temperature" per v0.7.1
    # DESIGN §3 — that workaround is now superseded).
    Y_species = np.zeros((N, S_species))
    Y_species[:, SPECIES_AMBIENT] = 1.0

    # v0.7.1 Phase 2/3: per-cell mixture thermophysical arrays derived
    # from Y each step. Allocated here, seeded + refreshed inside the
    # time loop after _advect_species runs. As of Phase 3 these arrays
    # are CONSUMED by PISO + post-PISO + source-CFL + nozzle BC.
    gamma_mix_arr = np.empty(N)
    Cp_mix_arr = np.empty(N)
    R_mix_arr = np.empty(N)
    M_mix_arr = np.empty(N)
    # v0.7.1 Phase 3: per-cell temperature ceiling. Refreshed each step
    # by ``_compute_T_ceiling_arr`` from species_params; consumed by the
    # PISO energy equation's T_raw clip.
    T_ceiling_arr = np.empty(N)
    # PB/P2: per-cell numerical-limit diagnostics (bookkeeping only; rows
    # are LIMIT_TEMPERATURE_FLOOR / _CEILING / LIMIT_PRESSURE_FLOOR).
    # First/last-time arrays use -1.0 as the "never active" sentinel.
    limit_activation_counts = np.zeros((N_LIMIT_DIAGNOSTICS, N), dtype=np.int64)
    limit_duration_s = np.zeros((N_LIMIT_DIAGNOSTICS, N))
    limit_abs_energy_j = np.zeros((N_LIMIT_DIAGNOSTICS, N))
    limit_first_time_s = np.full((N_LIMIT_DIAGNOSTICS, N), -1.0)
    limit_last_time_s = np.full((N_LIMIT_DIAGNOSTICS, N), -1.0)
    pressure_floor_max_deficit_pa = np.zeros(N)
    mach_limit_activation_counts = np.zeros(max(N - 1, 0), dtype=np.int64)
    mach_limit_duration_s = np.zeros(max(N - 1, 0))
    mach_limit_first_time_s = np.full(max(N - 1, 0), -1.0)
    mach_limit_last_time_s = np.full(max(N - 1, 0), -1.0)

    # Ignition state
    is_burning = np.zeros(N, dtype=np.bool_)
    has_ignited = np.zeros(N, dtype=np.bool_)
    ignition_time = np.full(N, 1e10)

    # Burn rates
    r_total = np.zeros(N)
    r_erosive = np.zeros(N)

    # Working arrays
    mass_source = np.zeros(N)
    # v0.7.1: per-species mass-source accounting. mass_source_by_species[i, s]
    # is the mass source rate [kg/s/m] (cell rate, same units as mass_source)
    # for species s into cell i. Sums across s equal mass_source[i].
    # Populated by the source-application sites (pyrogen injection,
    # Goodman grain-mass sources). Used in Phase 1d-e for Y advection.
    mass_source_by_species = np.zeros((N, S_species))
    thermal_source = np.zeros(N)
    momentum_source = np.zeros(N + 1)
    pyrogen_surface_heat_flux = np.zeros(N)
    radiation_heat_flux = np.zeros(N)
    radiation_sink_power = np.zeros(N)
    radiation_emitter = np.zeros(N, dtype=np.bool_)
    f_darcy = np.zeros(N)
    Re = np.zeros(N)
    Mach = np.zeros(N)
    u_cell = np.zeros(N)
    T_surf = np.full(N, propellant.T_initial)
    delta = np.full(N, 1.0e-6)
    # SM1: conduction-model state (1x2 dummies unless solid_heating_model=
    # 'conduction', so the default path allocates nothing).
    if _solid_code == SOLID_CONDUCTION and propellant.k_solid > 0.0:
        _alpha_s = propellant.k_solid / (propellant.rho_propellant * propellant.Cps)
        solid_y = build_solid_grid(_alpha_s, t_max, CONDUCTION_DY0,
                                   CONDUCTION_RATIO)
        solid_T = np.full((N, solid_y.shape[0]), propellant.T_initial)
        solid_work = np.zeros((2, solid_y.shape[0]))
    else:
        solid_y = np.array([0.0, 1.0e-6])
        solid_T = np.full((1, 2), propellant.T_initial)
        solid_work = np.zeros((2, 2))
    solid_age = np.zeros(N)
    solid_nact = np.full(N, 3, dtype=np.int64)
    # SM2: what lit each cell (always recorded; bookkeeping only) and the opt-in
    # open-loop forcing recorder (record_forcing_steps > 0 = first n steps).
    ign_diag = np.full((N, 10), np.nan)
    _rec_n = int(record_forcing_steps)
    rec_h = np.full((max(_rec_n, 1), N if _rec_n > 0 else 1), -1.0)
    rec_Td = np.zeros_like(rec_h)
    rec_t = np.zeros(max(_rec_n, 1))
    rec_dt = np.zeros(max(_rec_n, 1))
    rec_k = np.zeros(1, dtype=np.int64)

    plenum_state = initial_plenum_state(
        pyrogen_chamber, P_initial=P_ambient, T_initial=propellant.T_initial
    )
    pyrogen_params_arr = pyrogen_params(pyrogen_chamber.pyrogen)
    chamber_params_arr = chamber_params(pyrogen_chamber)

    # v0.7.2 Phase A: pyrogen-injection axial weights. Computed once at
    # sim init from L_jet = kappa_jet * d_throat_pyrogen (geometry doesn't
    # change during pyrogen burn, so weights are constant over the
    # simulation). See srm_1d/docs/v0_7_2/candidates/03_pyrogen_spatial_distribution.md.
    # The dx array is uniform (geo.dx is scalar); building np.full(N, dx)
    # keeps the kernel signature general for future non-uniform grids.
    d_throat_pyrogen = np.sqrt(4.0 * pyrogen_chamber.A_throat / np.pi)
    L_jet = float(pyrogen_chamber.pyrogen.kappa_jet) * d_throat_pyrogen
    dx_arr = np.full(N, dx)
    pyrogen_axial_weights = _compute_pyrogen_axial_weights(
        x_centers, dx_arr, L_jet, N
    )

    # v0.7.3 Phase A: igniter topology resolution. forward_plenum
    # preserves v0.7.0-v0.7.2 behavior byte-for-byte (i_start=i_end=0,
    # A_burn_per_cell unused, mdot_uncontained_arr unused — the
    # _run_time_loop pyrogen-injection block branches on topology_code).
    # head_basket / aft_basket use the resolved cartridge cell range
    # and per-cell A_burn distribution.
    from srm_1d.igniter_plenum import _topology_code as _resolve_topology_code
    topology_code = _resolve_topology_code(pyrogen_chamber.injection_topology)
    cart_i_start, cart_i_end = pyrogen_chamber.resolve_injection_cells(
        x_centers, dx, N, A_port,
    )
    n_cart_cells = max(1, cart_i_end - cart_i_start + 1)
    A_burn_per_cell = pyrogen_chamber.A_burn_initial / n_cart_cells
    mdot_uncontained_arr = np.zeros(N)

    # v0.7.3 Phase B.4: pyrogen-to-surface heat delivery resolution.
    # The mode is read from the Pyrogen YAML; the DeMar flux conversion
    # cal/cm²/s → W/m² happens here (1 cal = 4.184 J; 1 cal/cm²/s =
    # 4.184e4 W/m²). For radiation mode, pellet_emissivity and
    # radiation_absorption_length_m come straight from the Pyrogen.
    heat_delivery_mode_code = _heat_delivery_code(
        pyrogen_chamber.pyrogen.heat_delivery_mode
    )
    if pyrogen_chamber.pyrogen.heat_flux_cal_cm2_s is not None:
        demar_flux_w_m2_uncontained = float(
            pyrogen_chamber.pyrogen.heat_flux_cal_cm2_s
        ) * 4.184e4
    else:
        demar_flux_w_m2_uncontained = 0.0
    T_flame_pyrogen_val = float(pyrogen_chamber.pyrogen.T_flame)
    pellet_emissivity_val = float(pyrogen_chamber.pyrogen.pellet_emissivity)
    radiation_absorption_length_val = float(
        pyrogen_chamber.pyrogen.radiation_absorption_length_m
    )
    pyrogen_heat_flux_arr = np.zeros(N)

    # v0.7.2 Phase B-v2: flame-front h_c augmentation working state.
    # flame_spread_augment[i] is refreshed each step by
    # _compute_flame_front_augment inside _run_time_loop; allocated here
    # as a hot-loop scratch buffer initialized to 1.0 (no-op for the
    # first PISO step before the kernel fills it).
    # Knobs default to Propellant attributes (flame_spread_tau ~ 1 ms,
    # flame_spread_boost ~ 3.0). See
    # srm_1d/docs/v0_7_2/candidates/02_spatial_ignition_front_coupling.md
    # and the negative finding doc'd at commit 065d193 that motivated
    # the v2 reformulation.
    flame_spread_augment = np.ones(N)
    flame_spread_enabled = bool(getattr(propellant, 'flame_spread_enabled', True))
    flame_spread_tau = float(getattr(propellant, 'flame_spread_tau', 1.0e-3))
    flame_spread_boost = float(getattr(propellant, 'flame_spread_boost', 3.0))

    # v0.7.4 Phase F: flame-spread front exposure mask. Refreshed each step
    # by _advance_flame_front inside _run_time_loop when enabled; allocated
    # here as a hot-loop scratch buffer (bool, init True = exposed).
    flame_front_enabled = bool(getattr(propellant, 'flame_front_enabled', False))
    flame_front_velocity = float(getattr(propellant, 'flame_front_velocity', 3.0))
    ignitable = np.ones(N, dtype=np.bool_)

    # v0.7.4 Phase Z: Z-N dynamic burn-rate relaxation state. r_dyn is the
    # persistent relaxed rate; r_qs_persist / r_ero_qs_persist hold the
    # quasi-steady Ma targets between burn-rate updates. tau_establishment
    # is forced off when Z-N is enabled (they would double-damp the spike).
    zn_enabled = bool(getattr(propellant, 'zn_enabled', False))
    kappa_zn = float(getattr(propellant, 'kappa_zn', 1.0))
    r_dyn = np.zeros(N)
    r_qs_persist = np.zeros(N)
    r_ero_qs_persist = np.zeros(N)
    if zn_enabled:
        tau_establishment = 0.0

    theoretical_propellant_mass = (
        geo.total_propellant_volume() * propellant.rho_propellant
    )

    # Pre-allocate output arrays
    est_steps = _estimate_history_steps(
        t_max, cfl_target, dx, dt_max, species_params_arr)
    # MEM1: no fixed floor (it was 5,000,000 rows x ~60 arrays = 2.4 GB for every run,
    # even a 6 ms test). A run that outgrows the estimate is re-run once with double
    # the capacity by run_simulation (never silently truncated); an explicit
    # history_capacity is honoured as given.
    max_hist = est_steps
    if history_capacity is not None:
        if int(history_capacity) != history_capacity or history_capacity < 1:
            raise ValueError(
                "history_capacity must be None or a positive integer; "
                f"got {history_capacity!r}")
        max_hist = int(history_capacity)
    time_hist = np.empty(max_hist)
    P_head_hist = np.empty(max_hist)
    P_exit_hist = np.empty(max_hist)
    D_throat_hist = np.empty(max_hist)
    Kn_hist = np.empty(max_hist)
    massflow_hist = np.empty(max_hist)
    P_ig_hist = np.empty(max_hist)
    T_ig_hist = np.empty(max_hist)
    mdot_ig_hist = np.empty(max_hist)
    mdot_grain_hist = np.empty(max_hist)
    m_pyrogen_hist = np.empty(max_hist)
    gas_sensible_energy_before_hist = np.empty(max_hist)
    gas_sensible_energy_hist = np.empty(max_hist)
    gas_sensible_dE_dt_hist = np.empty(max_hist)
    normal_sidewall_thermal_power_hist = np.empty(max_hist)
    erosive_sidewall_thermal_power_hist = np.empty(max_hist)
    endface_thermal_power_hist = np.empty(max_hist)
    convective_scalar_flux_power_hist = np.empty(max_hist)
    clipping_correction_power_hist = np.empty(max_hist)
    pyrogen_enthalpy_power_hist = np.empty(max_hist)
    pyrogen_surface_heat_power_hist = np.empty(max_hist)
    gas_surface_heat_sink_power_hist = np.empty(max_hist)
    radiation_heat_power_hist = np.empty(max_hist)
    radiation_sink_power_hist = np.empty(max_hist)
    nozzle_enthalpy_power_hist = np.empty(max_hist)
    thermal_source_power_hist = np.empty(max_hist)
    energy_residual_hist = np.empty(max_hist)
    pyrogen_momentum_expected_hist = np.empty(max_hist)
    pyrogen_momentum_deposited_hist = np.empty(max_hist)
    pyrogen_momentum_residual_hist = np.empty(max_hist)
    dt_hist = np.empty(max_hist)
    n_burning_hist = np.empty(max_hist)
    n_ignited_hist = np.empty(max_hist)
    radiation_emitter_count_hist = np.empty(max_hist)
    radiation_receiver_count_hist = np.empty(max_hist)
    min_gas_temperature_hist = np.empty(max_hist)
    max_gas_temperature_hist = np.empty(max_hist)
    min_surface_temperature_hist = np.empty(max_hist)
    max_surface_temperature_hist = np.empty(max_hist)
    min_pressure_hist = np.empty(max_hist)
    max_pressure_hist = np.empty(max_hist)
    max_mach_hist = np.empty(max_hist)
    # Pressure-floor accounting (solver-verification T3): per-step count of
    # cells clamped to the 1 kPa floor and the gas mass [kg] the clamp added.
    pressure_floor_cells_hist = np.empty(max_hist)
    pressure_floor_mass_hist = np.empty(max_hist)
    # Energy audit (solver-verification T6): per-step ledger completeness
    # histories and the between-step storage-jump decomposition. See
    # tools/ignition_diagnostics.sensible_enthalpy_audit.
    pressure_floor_energy_hist = np.empty(max_hist)
    convective_wall_sink_power_hist = np.empty(max_hist)
    fill_debit_massflow_hist = np.empty(max_hist)
    fill_debit_power_hist = np.empty(max_hist)
    nozzle_enthalpy_power_post_hist = np.empty(max_hist)
    nozzle_massflow_piso_hist = np.empty(max_hist)
    coupling_passes_hist = np.empty(max_hist)
    coupling_rel_dT_hist = np.empty(max_hist)
    eos_jump_p2_kg_hist = np.empty(max_hist)
    jump_eos_J_hist = np.empty(max_hist)
    jump_mixture_J_hist = np.empty(max_hist)
    jump_geometry_J_hist = np.empty(max_hist)
    jump_eos_kg_hist = np.empty(max_hist)
    jump_mixture_kg_hist = np.empty(max_hist)
    jump_geometry_kg_hist = np.empty(max_hist)
    jump_total_kg_hist = np.empty(max_hist)

    # Pre-allocate snapshot storage
    max_snaps = int(t_max / snapshot_interval) + 10
    snap_data = np.empty((max_snaps, N_SNAP_CHANNELS, N))
    snap_times = np.empty(max_snaps)
    # Per-segment end-face regression history (axial face burnback) — parallel
    # to the snapshots, for the longitudinal-slice viewer.
    snap_seg_fwd = np.zeros((max_snaps, ga['N_seg']))
    snap_seg_aft = np.zeros((max_snaps, ga['N_seg']))

    # ============================================================
    # STATUS PRINT
    # ============================================================
    if verbose:
        numba_status = "Numba JIT enabled" if HAS_NUMBA else "Pure Python (no Numba)"
        print(f"PISO Solver ({numba_status}): {propellant.name}")
        print(f"  Motor: L={geo.L_motor*1e3:.0f}mm  D_outer={geo.D_outer*1e3:.0f}mm  "
              f"D_throat={nozzle.D_throat*1e3:.1f}mm  segments={ga['N_seg']}")
        print(f"  Params: Ts={propellant.T_surface:.0f}K  Cps={propellant.Cps:.0f}  "
              f"roughness={roughness*1e6:.0f}um  kappa={kappa}")
        print(f"  Pyrogen: {pyrogen_chamber.pyrogen.name}  "
              f"m={pyrogen_chamber.m_pyrogen_initial*1e3:.1f}g  "
              f"A_t={pyrogen_chamber.A_throat*1e6:.2f}mm^2  "
              f"T_ignition={T_ignition:.0f}K")
        print()

    wall_start = clock.time()

    # v0.8.0 Phase 6: shared progress/cancel cell. [0] = progress 0..1 written
    # by the @njit loop each step; [1] = cancel flag a caller (e.g. the GUI
    # progress dialog poller) sets to request a cooperative abort. Allocated
    # here when no external array is supplied so the loop signature is uniform.
    if progress_state is None:
        progress_state = np.zeros(2, dtype=np.float64)

    # ============================================================
    # RUN COMPILED TIME LOOP
    # ============================================================
    (n_steps, n_snaps,
     total_mass_produced, total_mass_nozzle,
     first_burnthrough_time, D_throat_final,
     termination_code,
     pyrogen_mass_burned, pyrogen_duration,
     pyrogen_peak_P,
     floor_energy_total, floor_min_raw_P_global,
     max_burn_refresh_gap, max_geometry_refresh_gap) = _run_time_loop(
        # Cell arrays
        rho, u, P, T,
        D_port, x_centers, A_port, C_burn, D_hyd,
        is_grain, endface_msource,
        is_burning, has_ignited, ignition_time,
        r_total, r_erosive,
        mass_source, thermal_source, momentum_source, pyrogen_surface_heat_flux,
        radiation_heat_flux, radiation_sink_power, radiation_emitter,
        f_darcy, Re, Mach, u_cell,
        T_surf, delta,
        regress,
        # Segment arrays
        ga['seg_x_start'], ga['seg_length'],
        ga['seg_fwd_regression'], ga['seg_aft_regression'],
        ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'],
        ga['seg_D_bore_fwd'], ga['seg_D_bore_aft'],
        ga['cell_D_bore_init'], ga['cell_segment_id'],
        ga['cell_wall_web'], ga['cell_segment_type'], ga['cell_fmm_idx'],
        ga['fmm_offset'], ga['fmm_reg_flat'],
        ga['fmm_perim_flat'], ga['fmm_port_flat'],
        # Geometry scalars
        N, ga['N_seg'], dx, ga['cell_D_outer'],
        # Gas/propellant scalars (gas thermo from representative tab)
        gas.gamma, gas.R_specific, rep_tab.T_flame,
        gas.Cp, gas.mu, gas.k_thermal, gas.Pr,
        propellant.rho_propellant, propellant.Cps,
        propellant.T_surface, propellant.T_initial, float(T_initial_gas),
        propellant.k_solid,
        # Burn rate tabs
        tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        # Simulation parameters
        roughness, kappa,
        cfl_target, dt_max, burn_update_interval, geometry_update_interval,
        source_cfl_factor, float(port_mach_cap), float(pressure_floor_pa),
        int(piso_energy_coupling_iters), bool(piso_corr2_storage),
        T_ignition, P_ambient, T_ambient,
        bool(diagnostic_disable_erosive), bool(diagnostic_disable_endfaces),
        bool(diagnostic_disable_momentum),
        bool(diagnostic_disable_pyrogen_surface_heating),
        bool(diagnostic_disable_adjacent_radiation),
        bool(diagnostic_disable_radiation_gas_sink),
        float(igniter_axial_momentum_fraction),
        float(pyrogen_surface_heat_flux_w_m2),
        float(active_radiation_emissivity),
        float(tau_establishment),
        t_max, P_cutoff,
        erosion_coeff, slag_coeff, throat_is_evolving,
        snapshot_interval,
        # Precomputed
        gamma_R, Gamma_crit, nozzle_denom,
        D_throat_init, A_throat_init,
        pyrogen_params_arr, chamber_params_arr, plenum_state,
        # Output: time history
        time_hist, P_head_hist, P_exit_hist, D_throat_hist,
        Kn_hist, massflow_hist, P_ig_hist, T_ig_hist, mdot_ig_hist,
        mdot_grain_hist,
        m_pyrogen_hist,
        gas_sensible_energy_before_hist,
        gas_sensible_energy_hist, gas_sensible_dE_dt_hist,
        normal_sidewall_thermal_power_hist,
        erosive_sidewall_thermal_power_hist,
        endface_thermal_power_hist,
        convective_scalar_flux_power_hist,
        clipping_correction_power_hist,
        pyrogen_enthalpy_power_hist,
        pyrogen_surface_heat_power_hist, gas_surface_heat_sink_power_hist,
        radiation_heat_power_hist, radiation_sink_power_hist,
        nozzle_enthalpy_power_hist, thermal_source_power_hist,
        energy_residual_hist,
        pyrogen_momentum_expected_hist, pyrogen_momentum_deposited_hist,
        pyrogen_momentum_residual_hist,
        dt_hist, n_burning_hist, n_ignited_hist,
        radiation_emitter_count_hist, radiation_receiver_count_hist,
        min_gas_temperature_hist, max_gas_temperature_hist,
        min_surface_temperature_hist, max_surface_temperature_hist,
        min_pressure_hist, max_pressure_hist, max_mach_hist,
        pressure_floor_cells_hist, pressure_floor_mass_hist,
        pressure_floor_energy_hist, convective_wall_sink_power_hist,
        fill_debit_massflow_hist, fill_debit_power_hist,
        nozzle_enthalpy_power_post_hist, nozzle_massflow_piso_hist,
        jump_eos_J_hist, jump_mixture_J_hist, jump_geometry_J_hist,
        jump_eos_kg_hist, jump_mixture_kg_hist, jump_geometry_kg_hist,
        jump_total_kg_hist,
        coupling_passes_hist, coupling_rel_dT_hist, eos_jump_p2_kg_hist,
        max_hist,
        # Output: snapshots
        snap_data, snap_times, max_snaps,
        snap_seg_fwd, snap_seg_aft,
        # v0.7.1: N-species state
        Y_species, species_params_arr, mass_source_by_species,
        gamma_mix_arr, Cp_mix_arr, R_mix_arr, M_mix_arr,
        T_ceiling_arr,
        # PB/P2: per-cell numerical-limit diagnostics
        limit_activation_counts, limit_duration_s, limit_abs_energy_j,
        limit_first_time_s, limit_last_time_s,
        pressure_floor_max_deficit_pa,
        mach_limit_activation_counts, mach_limit_duration_s,
        mach_limit_first_time_s, mach_limit_last_time_s,
        # v0.7.2 Phase A: pyrogen axial distribution
        pyrogen_axial_weights,
        # v0.7.2 Phase B-v2: flame-front h_c augmentation
        flame_spread_augment, flame_spread_enabled,
        flame_spread_tau, flame_spread_boost,
        # v0.7.3 Phase A: igniter topology (uncontained vs plenum)
        topology_code, cart_i_start, cart_i_end,
        A_burn_per_cell, mdot_uncontained_arr,
        heat_delivery_mode_code, demar_flux_w_m2_uncontained,
        T_flame_pyrogen_val, pellet_emissivity_val,
        radiation_absorption_length_val, pyrogen_heat_flux_arr,
        # v0.7.4 Phase F: flame-spread front gate
        flame_front_enabled, flame_front_velocity, ignitable,
        # v0.7.4 Phase Z: Z-N dynamic burn-rate relaxation
        zn_enabled, kappa_zn, r_dyn, r_qs_persist, r_ero_qs_persist,
        # spike-reopen probe: Mukunda-Paul high-g erosive ceiling
        float(erosive_mp_cap),
        int(_ero_code),
        # IG2: opt-in ignition heat-transfer closure
        int(_ign_h_code),
        int(_ap_code),
        bool(stop_if_unignitable),
        # SM1: solid-heating model + conduction state
        int(_solid_code), solid_y, solid_T, solid_age, solid_nact, solid_work,
        float(CONDUCTION_AGE_FRACTION),
        # SM2: ignition-event diagnostics + open-loop forcing recorder
        ign_diag, bool(_rec_n > 0), rec_h, rec_Td, rec_t, rec_dt, rec_k,
        # v0.8.0 Phase 6: live progress / cooperative cancel
        progress_state,
    )

    wall_elapsed = clock.time() - wall_start

    # ============================================================
    # TRIM AND WRAP RESULTS
    # ============================================================
    time_arr = _shrink_hist(time_hist, n_steps)
    P_head_arr = _shrink_hist(P_head_hist, n_steps)
    P_exit_arr = _shrink_hist(P_exit_hist, n_steps)
    D_throat_arr = _shrink_hist(D_throat_hist, n_steps)
    Kn_arr = _shrink_hist(Kn_hist, n_steps)
    massflow_arr = _shrink_hist(massflow_hist, n_steps)
    P_ig_arr = _shrink_hist(P_ig_hist, n_steps)
    T_ig_arr = _shrink_hist(T_ig_hist, n_steps)
    mdot_ig_arr = _shrink_hist(mdot_ig_hist, n_steps)
    mdot_grain_arr = _shrink_hist(mdot_grain_hist, n_steps)
    m_pyrogen_arr = _shrink_hist(m_pyrogen_hist, n_steps)
    gas_sensible_energy_before_arr = _shrink_hist(gas_sensible_energy_before_hist, n_steps)
    gas_sensible_energy_arr = _shrink_hist(gas_sensible_energy_hist, n_steps)
    gas_sensible_dE_dt_arr = _shrink_hist(gas_sensible_dE_dt_hist, n_steps)
    normal_sidewall_thermal_power_arr = _shrink_hist(normal_sidewall_thermal_power_hist, n_steps)
    erosive_sidewall_thermal_power_arr = _shrink_hist(erosive_sidewall_thermal_power_hist, n_steps)
    endface_thermal_power_arr = _shrink_hist(endface_thermal_power_hist, n_steps)
    convective_scalar_flux_power_arr = _shrink_hist(convective_scalar_flux_power_hist, n_steps)
    clipping_correction_power_arr = _shrink_hist(clipping_correction_power_hist, n_steps)
    pyrogen_enthalpy_power_arr = _shrink_hist(pyrogen_enthalpy_power_hist, n_steps)
    pyrogen_surface_heat_power_arr = _shrink_hist(pyrogen_surface_heat_power_hist, n_steps)
    gas_surface_heat_sink_power_arr = _shrink_hist(gas_surface_heat_sink_power_hist, n_steps)
    radiation_heat_power_arr = _shrink_hist(radiation_heat_power_hist, n_steps)
    radiation_sink_power_arr = _shrink_hist(radiation_sink_power_hist, n_steps)
    nozzle_enthalpy_power_arr = _shrink_hist(nozzle_enthalpy_power_hist, n_steps)
    thermal_source_power_arr = _shrink_hist(thermal_source_power_hist, n_steps)
    energy_residual_arr = _shrink_hist(energy_residual_hist, n_steps)
    pyrogen_momentum_expected_arr = _shrink_hist(pyrogen_momentum_expected_hist, n_steps)
    pyrogen_momentum_deposited_arr = _shrink_hist(pyrogen_momentum_deposited_hist, n_steps)
    pyrogen_momentum_residual_arr = _shrink_hist(pyrogen_momentum_residual_hist, n_steps)
    dt_arr = _shrink_hist(dt_hist, n_steps)
    n_burning_arr = _shrink_hist(n_burning_hist, n_steps)
    n_ignited_arr = _shrink_hist(n_ignited_hist, n_steps)
    radiation_emitter_count_arr = _shrink_hist(radiation_emitter_count_hist, n_steps)
    radiation_receiver_count_arr = _shrink_hist(radiation_receiver_count_hist, n_steps)
    min_gas_temperature_arr = _shrink_hist(min_gas_temperature_hist, n_steps)
    max_gas_temperature_arr = _shrink_hist(max_gas_temperature_hist, n_steps)
    min_surface_temperature_arr = _shrink_hist(min_surface_temperature_hist, n_steps)
    max_surface_temperature_arr = _shrink_hist(max_surface_temperature_hist, n_steps)
    min_pressure_arr = _shrink_hist(min_pressure_hist, n_steps)
    max_pressure_arr = _shrink_hist(max_pressure_hist, n_steps)
    max_mach_arr = _shrink_hist(max_mach_hist, n_steps)
    pressure_floor_cells_arr = _shrink_hist(pressure_floor_cells_hist, n_steps)
    pressure_floor_mass_arr = _shrink_hist(pressure_floor_mass_hist, n_steps)

    # ---- Energy audit (solver-verification T6): trim + derived ledgers ----
    pressure_floor_energy_arr = _shrink_hist(pressure_floor_energy_hist, n_steps)
    convective_wall_sink_power_arr = _shrink_hist(convective_wall_sink_power_hist, n_steps)
    fill_debit_massflow_arr = _shrink_hist(fill_debit_massflow_hist, n_steps)
    fill_debit_power_arr = _shrink_hist(fill_debit_power_hist, n_steps)
    nozzle_enthalpy_power_post_arr = _shrink_hist(nozzle_enthalpy_power_post_hist, n_steps)
    nozzle_massflow_piso_arr = _shrink_hist(nozzle_massflow_piso_hist, n_steps)
    coupling_passes_arr = _shrink_hist(coupling_passes_hist, n_steps)
    coupling_rel_dT_arr = _shrink_hist(coupling_rel_dT_hist, n_steps)
    eos_jump_p2_kg_arr = _shrink_hist(eos_jump_p2_kg_hist, n_steps)
    jump_eos_J_arr = _shrink_hist(jump_eos_J_hist, n_steps)
    jump_mixture_J_arr = _shrink_hist(jump_mixture_J_hist, n_steps)
    jump_geometry_J_arr = _shrink_hist(jump_geometry_J_hist, n_steps)
    jump_eos_kg_arr = _shrink_hist(jump_eos_kg_hist, n_steps)
    jump_mixture_kg_arr = _shrink_hist(jump_mixture_kg_hist, n_steps)
    jump_geometry_kg_arr = _shrink_hist(jump_geometry_kg_hist, n_steps)
    jump_total_kg_arr = _shrink_hist(jump_total_kg_hist, n_steps)

    # Source ledger: sum of the per-channel thermal-source powers vs the
    # PISO's own sum(thermal_source * dx) [W]. Positive channels inject,
    # negative channels (radiation sink, pyrogen surface sink, convective
    # wall sink) are debited from the gas.
    source_ledger_mismatch_arr = (
        normal_sidewall_thermal_power_arr
        + erosive_sidewall_thermal_power_arr
        + endface_thermal_power_arr
        + pyrogen_enthalpy_power_arr
        - radiation_sink_power_arr
        - gas_surface_heat_sink_power_arr
        - convective_wall_sink_power_arr
        - fill_debit_power_arr
        - thermal_source_power_arr
    )
    # Nozzle ledger: the PISO's explicit nozzle enthalpy outflow (start-of-step
    # T, pre-clamp P) minus the post-PISO reconstruction (final P, final T).
    nozzle_ledger_mismatch_arr = (
        nozzle_enthalpy_power_arr - nozzle_enthalpy_power_post_arr
    )
    nozzle_mass_mismatch_arr = nozzle_massflow_piso_arr - massflow_arr

    # Between-step storage jumps. Slot n = gap between the END of step n
    # (ledger storage) and the START of step n+1's PISO. The last slot has no
    # following step: its total/geometry/remainder are 0 and the audit
    # excludes it from every statistic. J jumps: the total comes straight
    # from the two stored ledger energies, so the remainder is an independent
    # check that the four named operations are exhaustive.
    jump_total_J_arr = np.zeros(n_steps)
    if n_steps > 1:
        jump_total_J_arr[:-1] = (
            gas_sensible_energy_before_arr[1:] - gas_sensible_energy_arr[:-1]
        )
    jump_remainder_J_arr = jump_total_J_arr - (
        jump_eos_J_arr + pressure_floor_energy_arr
        + jump_mixture_J_arr + jump_geometry_J_arr
    )
    jump_remainder_kg_arr = jump_total_kg_arr - (
        jump_eos_kg_arr + pressure_floor_mass_arr
        + jump_mixture_kg_arr + jump_geometry_kg_arr
    )
    if n_steps > 0:
        jump_remainder_J_arr[-1] = 0.0
        jump_remainder_kg_arr[-1] = 0.0

    # Pressure-floor accounting summary (solver-verification T3). Times are
    # the per-step ``time`` history convention (step START time).
    floor_clamped_steps = pressure_floor_cells_arr > 0.0
    floor_n_steps = int(np.count_nonzero(floor_clamped_steps))
    floor_activations = int(round(float(np.sum(pressure_floor_cells_arr))))
    floor_mass_total = float(np.sum(pressure_floor_mass_arr))
    if floor_n_steps > 0:
        floor_first_time = float(time_arr[int(np.argmax(floor_clamped_steps))])
        floor_min_raw_P = float(floor_min_raw_P_global)
    else:
        floor_first_time = -1.0
        floor_min_raw_P = float("nan")

    # Convert snapshot 3D array back to list of dicts for compatibility
    snapshots = []
    for s in range(n_snaps):
        snapshots.append({
            't': snap_times[s],
            'P': snap_data[s, _SNAP_P, :].copy(),
            'u': snap_data[s, _SNAP_U, :].copy(),
            'Mach': snap_data[s, _SNAP_MACH, :].copy(),
            'T': snap_data[s, _SNAP_T, :].copy(),
            'r_total': snap_data[s, _SNAP_R_TOTAL, :].copy(),
            'r_erosive': snap_data[s, _SNAP_R_EROSIVE, :].copy(),
            'D_port': snap_data[s, _SNAP_D_PORT, :].copy(),
            'x': x_centers.copy(),
            'C_burn': snap_data[s, _SNAP_C_BURN, :].copy(),
            'endface_msource': snap_data[s, _SNAP_ENDFACE, :].copy(),
            'is_burning': snap_data[s, _SNAP_IS_BURNING, :] > 0.5,
            'is_grain': snap_data[s, _SNAP_IS_GRAIN, :] > 0.5,
            'T_surf': snap_data[s, _SNAP_T_SURF, :].copy(),
            'mass_source': snap_data[s, _SNAP_MASS_SOURCE, :].copy(),
            'thermal_source': snap_data[s, _SNAP_THERMAL_SOURCE, :].copy(),
            'momentum_source': snap_data[s, _SNAP_MOMENTUM_SOURCE, :].copy(),
            'pyrogen_surface_heat_flux': snap_data[s, _SNAP_PYROGEN_SURFACE_HEAT_FLUX, :].copy(),
            'radiation_heat_flux': snap_data[s, _SNAP_RADIATION_HEAT_FLUX, :].copy(),
            'regress': snap_data[s, _SNAP_REGRESS, :].copy(),
            'rho': snap_data[s, _SNAP_RHO, :].copy(),
        })

    peak_idx = np.argmax(P_head_arr) if len(P_head_arr) > 0 else 0
    cstar = np.sqrt(gas.R_specific * rep_tab.T_flame) / Gamma_crit

    if first_burnthrough_time < 0:
        first_burnthrough_time = None

    finite_ignition_cells = np.flatnonzero(ignition_time < 1.0e9)
    if finite_ignition_cells.size:
        first_ignition_local = int(np.argmin(ignition_time[finite_ignition_cells]))
        first_ignition_cell = int(finite_ignition_cells[first_ignition_local])
        first_ignition_time = float(ignition_time[first_ignition_cell])
    else:
        first_ignition_cell = -1
        first_ignition_time = float("nan")

    # ============================================================
    # SUMMARY
    # ============================================================
    termination_names = {
        0: "t_max reached", 1: "complete burnout",
        2: "pressure cutoff", 3: "history array full",
        4: "numerical collapse aborted", 5: "canceled by user",
        6: "ignition impossible (stopped early)",
    }
    term_str = termination_names.get(termination_code, "unknown")
    if termination_code == 3:
        import warnings
        t_end = float(time_arr[-1]) if len(time_arr) > 0 else 0.0
        suggested = int(max_hist * t_max / max(t_end, 1e-12) * 1.1) + 1
        warnings.warn(
            f"run_simulation stopped at t = {t_end:.4g} s (< t_max = "
            f"{t_max:g} s): the history arrays filled after {max_hist} steps. "
            f"Results after t = {t_end:.4g} s are missing. Pass "
            f"history_capacity={suggested} (~{suggested * 450 / 1e9:.1f} GB "
            "of history) to complete the run.", RuntimeWarning, stacklevel=2)

    if verbose:
        print(f"\n{'='*65}")
        print(f"SUMMARY: {propellant.name}")
        print(f"  Motor: L={geo.L_motor*1e3:.0f}mm  D_outer={geo.D_outer*1e3:.0f}mm  "
              f"D_throat={nozzle.D_throat*1e3:.1f}mm  segments={ga['N_seg']}")
        print(f"  Params: Ts={propellant.T_surface:.0f}K  Cps={propellant.Cps:.0f}  "
              f"roughness={roughness*1e6:.0f}um  kappa={kappa}")
        print(f"  t_burn={time_arr[-1]:.3f}s  steps={n_steps}  cells={N}")
        print(f"  Termination: {term_str}")
        if first_burnthrough_time is not None:
            print(f"  t_first_burnout={first_burnthrough_time:.3f}s")
        print(f"  P_peak={P_head_arr[peak_idx]/1e6:.2f}MPa @ t={time_arr[peak_idx]:.3f}s")
        if len(P_head_arr) > 10:
            P_mid = np.mean(P_head_arr[max(0, len(P_head_arr)//4):len(P_head_arr)//2])
            print(f"  P_mid_burn={P_mid/1e6:.2f}MPa  P_final={P_head_arr[-1]/1e6:.2f}MPa")
        if throat_is_evolving:
            delta_mm = (D_throat_final - nozzle.D_throat) * 1000
            direction = "eroded" if delta_mm > 0 else "slagged"
            print(f"  Throat: {nozzle.D_throat*1e3:.2f} -> {D_throat_final*1e3:.2f} mm "
                  f"({direction} {abs(delta_mm):.3f} mm)")
        print(f"  Wall time: {wall_elapsed:.1f}s ({n_steps/max(wall_elapsed,0.01):.0f} steps/s)")
        print(f"  Mass: propellant={theoretical_propellant_mass:.3f}kg  "
              f"produced={total_mass_produced:.3f}kg  nozzle={total_mass_nozzle:.3f}kg  "
              f"balance_err={abs(total_mass_produced - total_mass_nozzle)/max(theoretical_propellant_mass,0.001)*100:.1f}%")
        print(f"  Pyrogen: burned={pyrogen_mass_burned*1e3:.2f}g  "
              f"peak_P_ig={pyrogen_peak_P/1e6:.2f}MPa  "
              f"duration={pyrogen_duration*1000:.1f}ms")
        print(f"  c*={cstar:.1f}m/s")
        print(f"{'='*65}")

    # Structured summary
    P_mid = float(np.mean(
        P_head_arr[max(0, len(P_head_arr)//4):len(P_head_arr)//2]
    )) if len(P_head_arr) > 10 else 0.0
    mass_balance_err = abs(
        total_mass_produced - total_mass_nozzle
    ) / max(theoretical_propellant_mass, 0.001)

    summary = {
        'propellant_mass': theoretical_propellant_mass,
        'mass_produced': total_mass_produced,
        'mass_nozzle': total_mass_nozzle,
        'mass_balance_error': mass_balance_err,
        # Guard the empty-history case (run cancelled before the first step).
        'P_peak': float(P_head_arr[peak_idx]) if len(P_head_arr) > 0 else 0.0,
        't_peak': float(time_arr[peak_idx]) if len(time_arr) > 0 else 0.0,
        'P_mid': P_mid,
        't_burn': float(time_arr[-1]) if len(time_arr) > 0 else 0.0,
        't_first_burnout': first_burnthrough_time,
        'c_star': cstar,
        'wall_time': wall_elapsed,
        'steps': n_steps,
        'cells': N,
        'history_capacity': int(max_hist),
        'history_capacity_explicit': history_capacity is not None,
        'history_estimated_steps': int(est_steps),
        'termination_code': int(termination_code),
        'history_cap_reached': bool(termination_code == 3),
        'dt_min': float(np.min(dt_arr)) if len(dt_arr) > 0 else float("nan"),
        'dt_median': float(np.median(dt_arr)) if len(dt_arr) > 0 else float("nan"),
        'dt_final': float(dt_arr[-1]) if len(dt_arr) > 0 else float("nan"),
        'first_ignition_time_s': first_ignition_time,
        'first_ignition_cell': first_ignition_cell,
        'D_throat_initial': nozzle.D_throat,
        'D_throat_final': float(D_throat_final),
        'termination': term_str,
        'pyrogen_name': pyrogen_chamber.pyrogen.name,
        'pyrogen_mass_initial': float(pyrogen_chamber.m_pyrogen_initial),
        'pyrogen_mass_burned': float(pyrogen_mass_burned),
        'pyrogen_mass_remaining': float(m_pyrogen_arr[-1]) if len(m_pyrogen_arr) > 0 else float(plenum_state[0]),
        'pyrogen_duration': float(pyrogen_duration),
        'pyrogen_peak_P': float(pyrogen_peak_P),
        'initial_gas_temperature': float(T_initial_gas),
        'ambient_temperature': float(T_ambient),
        # Erosive-closure registry provenance (docs/erosive_model_registry/
        # DESIGN.md §4) — traces, probes, and CSV exports are otherwise
        # unlabeled as to which erosive closure produced them.
        'erosive_model': erosive_model,
        'ignition_heat_transfer': ignition_heat_transfer,
        'peretz_ap_reference': peretz_ap_reference,
        'ignited': bool(first_ignition_time == first_ignition_time),  # not NaN
        # IG3: grain cells that never lit (convective-only ignition model:
        # no igniter-jet impingement heating, burning-neighbour radiation off
        # by default; NOT an igniter-size indicator). See run_simulation.
        'grain_cells': int(np.count_nonzero(is_grain)),
        'grain_cells_ignited': int(np.count_nonzero(
            is_grain & (ignition_time < 1.0e9))),
        'grain_never_lit_mask': (is_grain & (ignition_time >= 1.0e9)).copy(),
        'grain_cells_never_lit': int(np.count_nonzero(
            is_grain & (ignition_time >= 1.0e9))),
        'grain_fraction_lit': float(np.count_nonzero(
            is_grain & (ignition_time < 1.0e9))
            / max(int(np.count_nonzero(is_grain)), 1)),
        'grain_never_lit_cell_range': (
            [int(np.flatnonzero(is_grain & (ignition_time >= 1.0e9))[0]),
             int(np.flatnonzero(is_grain & (ignition_time >= 1.0e9))[-1])]
            if np.any(is_grain & (ignition_time >= 1.0e9)) else None),
        'grain_never_lit_x_range_mm': (
            [float(x_centers[np.flatnonzero(
                is_grain & (ignition_time >= 1.0e9))[0]]) * 1e3,
             float(x_centers[np.flatnonzero(
                 is_grain & (ignition_time >= 1.0e9))[-1]]) * 1e3]
            if np.any(is_grain & (ignition_time >= 1.0e9)) else None),
        'stop_if_unignitable': bool(stop_if_unignitable),
        'solid_heating_model': solid_heating_model,
        'diagnostic_disable_erosive': bool(diagnostic_disable_erosive),
        'diagnostic_disable_endfaces': bool(diagnostic_disable_endfaces),
        'diagnostic_disable_momentum': bool(diagnostic_disable_momentum),
        'diagnostic_disable_pyrogen_surface_heating': bool(diagnostic_disable_pyrogen_surface_heating),
        'diagnostic_disable_adjacent_radiation': bool(diagnostic_disable_adjacent_radiation),
        'diagnostic_disable_radiation_gas_sink': bool(diagnostic_disable_radiation_gas_sink),
        'igniter_axial_momentum_fraction': float(igniter_axial_momentum_fraction),
        'pyrogen_heat_flux_cal_cm2_s': (
            None if pyrogen_heat_flux is None else float(pyrogen_heat_flux)
        ),
        'pyrogen_surface_heat_flux_w_m2': float(pyrogen_surface_heat_flux_w_m2),
        'radiation_emissivity': float(propellant.radiation_emissivity),
        'active_radiation_emissivity': float(active_radiation_emissivity),
        'tau_establishment': float(tau_establishment),
        'energy_residual_convention': (
            'gas_sensible_dE_dt - convective_scalar_flux_power - '
            'thermal_source_power - clipping_correction_power; '
            'convective power is positive into the gas'
        ),
        # Pressure-floor accounting (solver-verification T3). The PISO step
        # clamps P to 1 kPa and rebuilds rho = P/(R*T); the mass/energy that
        # clamp injects is otherwise outside every ledger. Cell-steps =
        # activations; mass is vs the EOS at the unclamped pressure.
        'pressure_floor_activations': floor_activations,
        'pressure_floor_steps': floor_n_steps,
        'pressure_floor_mass_added_kg': floor_mass_total,
        'pressure_floor_energy_added_J': float(floor_energy_total),
        'pressure_floor_first_time_s': floor_first_time,
        'pressure_floor_min_raw_P': floor_min_raw_P,
        'pressure_floor_pa': float(pressure_floor_pa),
        # Refresh cadence + staleness (solver-verification T5). Intervals are
        # in flow steps; gaps are the largest simulated time [s] between two
        # consecutive refreshes of that kind (burn rate / port geometry).
        'burn_refresh_interval_steps': int(burn_update_interval),
        'geometry_refresh_interval_steps': int(geometry_update_interval),
        'max_burn_refresh_gap_s': float(max_burn_refresh_gap),
        'max_geometry_refresh_gap_s': float(max_geometry_refresh_gap),
        # Energy audit (solver-verification T6). The balance is a
        # sensible-enthalpy (Cp*T) one, NOT total energy: see
        # tools/ignition_diagnostics.sensible_enthalpy_audit for the ledger
        # closure, the between-step storage-jump decomposition and the size
        # of the omitted terms.
        'energy_convention': 'sensible_enthalpy_CpT',
        'energy_omitted_terms': [
            'kinetic_energy', 'pressure_work',
            'formation_energy_implicit_via_T_flame',
            'solid_energy_storage',
            'nozzle_kinetic_energy_static_outflow',
        ],
    }

    # Per-grain summary from snapshots
    #
    # v0.8.x: regression/web are derived from the per-cell ``regress``
    # distance carried in the snapshots (the openMotor-style
    # ``(avg_D - D_bore_init)/2`` form is BATES-only and produces -web for
    # FMM grains whose ``cell_D_bore_init == D_outer``). srm_1d regression
    # is axially varying, so we collapse to physically meaningful scalars:
    #   * ``regression`` = the FORE (foremost) cell's regress distance — the
    #     least-regressed cell, last to burn out, so the legacy grain
    #     cross-section stays non-blank while the head retains web.
    #   * ``web`` (remaining) = MIN over the grain's cells of
    #     (wall_web - regress) — burnout is governed by the first cell to
    #     break through. See STATION_VIZ_DESIGN.md §8a.
    grain_data = []
    cell_wall_web = ga['cell_wall_web']
    for seg_idx in range(ga['N_seg']):
        seg_cells = ga['cell_segment_id'] == seg_idx
        seg_idx_arr = np.where(seg_cells)[0]
        regression_hist = []
        web_hist = []
        for snap in snapshots:
            if len(seg_idx_arr) > 0:
                regress_seg = snap['regress'][seg_idx_arr]
                fore_cell = seg_idx_arr[0]  # foremost (head) cell of the span
                reg = snap['regress'][fore_cell]
                web = float(np.min(cell_wall_web[seg_idx_arr] - regress_seg))
                if web < 0.0:
                    web = 0.0
            else:
                reg = 0.0
                web = 0.0
            regression_hist.append(reg)
            web_hist.append(web)
        grain_data.append({
            'segment': seg_idx,
            'D_bore_fwd': ga['seg_D_bore_fwd'][seg_idx],
            'D_bore_aft': ga['seg_D_bore_aft'][seg_idx],
            'regression': np.array(regression_hist),
            'web': np.array(web_hist),
        })

    result_dict = {
        'time': time_arr, 'P_head': P_head_arr, 'P_exit': P_exit_arr,
        'D_throat': D_throat_arr, 'Kn': Kn_arr, 'massflow': massflow_arr,
        'P_ig': P_ig_arr, 'T_ig': T_ig_arr, 'mdot_ig': mdot_ig_arr, 'mdot_grain': mdot_grain_arr,
        'm_pyrogen': m_pyrogen_arr,
        'gas_sensible_energy_before': gas_sensible_energy_before_arr,
        'gas_sensible_energy': gas_sensible_energy_arr,
        'gas_sensible_dE_dt': gas_sensible_dE_dt_arr,
        'normal_sidewall_thermal_power': normal_sidewall_thermal_power_arr,
        'erosive_sidewall_thermal_power': erosive_sidewall_thermal_power_arr,
        'endface_thermal_power': endface_thermal_power_arr,
        'pyrogen_gas_thermal_power': pyrogen_enthalpy_power_arr,
        'convective_scalar_flux_power': convective_scalar_flux_power_arr,
        'nozzle_scalar_flux_power': nozzle_enthalpy_power_arr,
        'clipping_correction_power': clipping_correction_power_arr,
        'pyrogen_enthalpy_power': pyrogen_enthalpy_power_arr,
        'pyrogen_surface_heat_power': pyrogen_surface_heat_power_arr,
        'gas_surface_heat_sink_power': gas_surface_heat_sink_power_arr,
        'radiation_heat_power': radiation_heat_power_arr,
        'radiation_sink_power': radiation_sink_power_arr,
        'nozzle_enthalpy_power': nozzle_enthalpy_power_arr,
        'thermal_source_power': thermal_source_power_arr,
        'energy_residual': energy_residual_arr,
        'pyrogen_momentum_expected': pyrogen_momentum_expected_arr,
        'pyrogen_momentum_deposited': pyrogen_momentum_deposited_arr,
        'pyrogen_momentum_residual': pyrogen_momentum_residual_arr,
        'dt': dt_arr,
        'n_burning': n_burning_arr,
        'n_ignited': n_ignited_arr,
        'radiation_emitter_count': radiation_emitter_count_arr,
        'radiation_receiver_count': radiation_receiver_count_arr,
        'min_gas_temperature': min_gas_temperature_arr,
        'max_gas_temperature': max_gas_temperature_arr,
        'min_surface_temperature': min_surface_temperature_arr,
        'max_surface_temperature': max_surface_temperature_arr,
        'min_pressure': min_pressure_arr,
        'max_pressure': max_pressure_arr,
        'max_mach': max_mach_arr,
        'pressure_floor_cells': pressure_floor_cells_arr,
        'pressure_floor_mass': pressure_floor_mass_arr,
        # Energy audit (solver-verification T6). Slot n of the *_jump_*
        # arrays = gap between the end of step n and the start of step n+1's
        # PISO (J = sensible Cp*T energy, kg = gas mass); the last slot is
        # incomplete and excluded by sensible_enthalpy_audit.
        'pressure_floor_energy': pressure_floor_energy_arr,
        'convective_wall_sink_power': convective_wall_sink_power_arr,
        'regression_fill_debit_massflow': fill_debit_massflow_arr,
        'regression_fill_debit_power': fill_debit_power_arr,
        'nozzle_enthalpy_power_post_piso': nozzle_enthalpy_power_post_arr,
        'nozzle_massflow_piso': nozzle_massflow_piso_arr,
        'piso_coupling_passes': coupling_passes_arr,
        'piso_coupling_rel_dT': coupling_rel_dT_arr,
        'storage_jump_eos_p2_kg': eos_jump_p2_kg_arr,
        'source_ledger_mismatch_W': source_ledger_mismatch_arr,
        'nozzle_ledger_mismatch_W': nozzle_ledger_mismatch_arr,
        'nozzle_mass_mismatch_kg_s': nozzle_mass_mismatch_arr,
        'storage_jump_eos_J': jump_eos_J_arr,
        'storage_jump_floor_J': pressure_floor_energy_arr,
        'storage_jump_mixture_J': jump_mixture_J_arr,
        'storage_jump_geometry_J': jump_geometry_J_arr,
        'storage_jump_remainder_J': jump_remainder_J_arr,
        'storage_jump_total_J': jump_total_J_arr,
        'mass_jump_eos_kg': jump_eos_kg_arr,
        'mass_jump_floor_kg': pressure_floor_mass_arr,
        'mass_jump_mixture_kg': jump_mixture_kg_arr,
        'mass_jump_geometry_kg': jump_geometry_kg_arr,
        'mass_jump_remainder_kg': jump_remainder_kg_arr,
        'mass_jump_total_kg': jump_total_kg_arr,
        'ignition_time_by_cell': ignition_time.copy(),
        # SM2: per grain cell at its ignition step: [h_total, T_driver, q_pre, T_s, T_gas, u,
        # Mach, P, t, q_post] (NaN = never lit); forcing = open-loop recorder (or None)
        'ignition_diag': ign_diag.copy(),
        'forcing': ({'h': rec_h[:int(rec_k[0])].copy(), 'T_drv': rec_Td[:int(rec_k[0])].copy(),
                     't': rec_t[:int(rec_k[0])].copy(), 'dt': rec_dt[:int(rec_k[0])].copy()}
                    if _rec_n > 0 else None),
        # v0.8.x station-viz data contract: per-cell axial geometry so the
        # GUI station model + axial-payload builder work from `result`
        # alone (cell->grain map; gap sentinel -1). x_cell mirrors each
        # snapshot's constant 'x'. INITIAL map (t=0) so station placement /
        # classification reflects the as-designed grain, not the final
        # face-receded extent (which would push 'fore' inward over the burn).
        'cell_segment_id': cell_segment_id_init,
        'x_cell': x_centers.copy(),
        # v0.8.x roadmap #2 (longitudinal motor-slice viewer): constant
        # per-cell geometry for drawing the burnback slice. All live OUTSIDE
        # the @njit loop (scalars + a geometry array), so no snapshot channel
        # is needed: R_outer = D_outer/2, axial cell extent = dx, %web =
        # 1 - regress/cell_wall_web (R_bore = the snapshot D_port/2, already
        # the hydraulic-equivalent bore for FMM/non-circular grains).
        'dx': float(dx),
        'D_outer': float(ga['D_outer']),
        # Per-cell casing (outer) diameter for the slice viewer's OD-tapered
        # casing. Constant == D_outer for non-OD motors (slice falls back to
        # the scalar when this key is absent in pre-OD results).
        'cell_D_outer': ga['cell_D_outer'].copy(),
        'cell_wall_web': ga['cell_wall_web'].copy(),
        # Per-segment end-face burnback over the snapshot frames: the grain
        # extends axially from seg_x_start+fwd_reg to seg_x_start+seg_length-
        # aft_reg, so the slice viewer can recede the faces continuously
        # (not snap whole cells). fwd_reg/aft_reg are (n_snaps, N_seg).
        'seg_geom': {
            'seg_x_start': ga['seg_x_start'].copy(),
            'seg_length': ga['seg_length'].copy(),
            'seg_fwd_reg': snap_seg_fwd[:n_snaps].copy(),
            'seg_aft_reg': snap_seg_aft[:n_snaps].copy(),
        },
        'snapshots': snapshots, 'grains': grain_data,
        'summary': summary,
        'P_ambient': P_ambient,
        # Per-cell numerical-limit activation maps (solver-verification
        # PB/P2). Bookkeeping only. Times are step-START based: first =
        # start of the first active step, last = END of the last active
        # step; NaN = never active. Energies are sum(|clip energy|) [J].
        'numerical_limits': _build_numerical_limits(
            x_centers, T_ambient, T_ceiling_arr, pressure_floor_pa,
            port_mach_cap, limit_activation_counts, limit_duration_s,
            limit_abs_energy_j, limit_first_time_s, limit_last_time_s,
            pressure_floor_max_deficit_pa, mach_limit_activation_counts,
            mach_limit_duration_s, mach_limit_first_time_s,
            mach_limit_last_time_s),
        # v0.7.1: N-species final state and species registry
        'Y_species_final': Y_species.copy(),
        'species_params': species_params_arr.copy(),
        'species_names': [sp.name for sp in species_list],
        'rho_final': rho.copy(),
        'A_port_final': A_port.copy(),
        # v0.7.1 Phase 2: per-cell mixture arrays (final state).
        # Not yet consumed by the solver — solver still uses the rep-tab
        # scalars (gas.gamma, gas.R_specific, etc.). Phase 3 will wire
        # these arrays into the PISO step.
        'gamma_mix_final': gamma_mix_arr.copy(),
        'Cp_mix_final': Cp_mix_arr.copy(),
        'R_mix_final': R_mix_arr.copy(),
        'M_mix_final': M_mix_arr.copy(),
    }

    # v0.8.0 Phase 1 (capstone): return the channel object as the primary
    # result. It proxies item access to the raw dict (so legacy
    # result['P_head'] / result['summary'] code is unchanged) while exposing
    # the unit-aware .channels / .axial API for the openMotor frontend.
    from .channels import build_channels
    return build_channels(result_dict)


# ---------------------------------------------------------------------------
# IG3: auto-sized-charge retry wrapper
# ---------------------------------------------------------------------------
class PyrogenAutoSizeError(RuntimeError):
    """An auto-sized pyrogen charge failed to light the motor within the
    allowed number of attempts (see :func:`run_simulation`)."""


def _self_sustaining(res):
    """True iff the GRAIN, not the igniter, drives the chamber once the
    igniter is spent.

    Success = after the pyrogen is consumed, the grain's own gas-generation rate
    exceeds, at some step, the igniter's PEAK gas-generation rate of this attempt
    (two flows compared; no tuned number). If the run ended before the pyrogen
    was spent the comparison is made at the last step. There is no extinction
    model, so a grain that lights one cell and stays lit still fails this test
    when that cell cannot out-produce the igniter (a hang-fire at a chamber
    pressure set by the igniter's dribble). The two rates are written into the
    summary (``igniter_peak_mdot_kg_s``, ``grain_mdot_after_igniter_kg_s``).
    """
    summ = res['summary']
    t = np.asarray(res['time'])
    if t.size == 0:
        summ['igniter_peak_mdot_kg_s'] = 0.0
        summ['grain_mdot_after_igniter_kg_s'] = 0.0
        return False
    mig = np.asarray(res['mdot_ig'])
    mgr = np.asarray(res['mdot_grain'])
    peak_ig = float(mig.max()) if mig.size else 0.0
    after = t > float(summ['pyrogen_duration'])
    grain = float(mgr[after].max()) if np.any(after) else float(mgr[-1])
    summ['igniter_peak_mdot_kg_s'] = peak_ig
    summ['grain_mdot_after_igniter_kg_s'] = grain
    return bool(summ['ignited'] and int(summ['grain_cells_ignited']) >= 1
                and grain > peak_ig)


def _run_once_growing_history(bound):
    """One simulation attempt. MEM1: if the AUTO-sized history (no explicit
    ``history_capacity``) fills before the run ends (termination code 3), the
    attempt is re-run ONCE with double the capacity and a ``RuntimeWarning``
    instead of returning a truncated result. An explicit ``history_capacity`` is
    honoured as given (the 'history array full' warning of the single run stays).
    """
    import warnings
    explicit = bound.arguments.get('history_capacity') is not None
    if explicit:
        return _run_simulation_once(*bound.args, **bound.kwargs)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        res = _run_simulation_once(*bound.args, **bound.kwargs)
    full = int(res['summary']['termination_code']) == 3
    for w in caught:                         # re-emit, except the superseded one
        if full and 'history arrays filled' in str(w.message):
            continue
        warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
    if not full:
        return res
    cap = int(res['summary']['history_capacity'])
    res = None                               # free the truncated run before the re-run
    warnings.warn(
        f"the auto-sized history ({cap} rows, estimated from the CFL wave-speed "
        f"bound) filled before the run ended; re-running once with {2 * cap} rows. "
        "Pass history_capacity explicitly to avoid the re-run.",
        RuntimeWarning, stacklevel=3)
    bound.arguments['history_capacity'] = 2 * cap
    try:
        return _run_simulation_once(*bound.args, **bound.kwargs)
    finally:
        bound.arguments['history_capacity'] = None


def run_simulation(*args, ignition_retry_factor=1.25,
                   ignition_retry_max_attempts=10, **kwargs):
    """Run a complete transient simulation (all parameters as
    :func:`_run_simulation_once`, documented below), with the IG3 ignition
    fallback for AUTO-SIZED pyrogen charges.

    SUCCESS = the grain, not the igniter, drives the chamber once the igniter
    is spent (``_self_sustaining``: after the pyrogen is consumed the grain's
    own gas-generation rate exceeds, at some step, the igniter's peak
    gas-generation rate of that attempt; compared at the last step if the run
    ended first). One lit cell that cannot out-produce the igniter is a
    hang-fire, not a success. It does NOT require every grain cell to light: cells the convective-only
    ignition model never heats (e.g. BALLSstick's stagnant head-end 10 %) are
    reported by a warning and by the ``grain_*`` summary keys instead, because
    a bigger igniter charge is not the remedy (it only inflates the startup
    spike).

    If the pyrogen charge was auto-sized (``PyrogenChamber.mass_auto``, i.e.
    ``build_pyrogen_chamber(pyrogen_mass=None)``) and the run is not a
    success, a run in which NO cell can light stops early
    (``stop_if_unignitable``: no energy source left and every grain-cell gas
    temperature below ``T_ignition``) and is retried at ``ignition_retry_factor``
    times the charge (default 1.25, the user-approved step), at most
    ``ignition_retry_max_attempts`` attempts in all (default 10, i.e. up to
    1.25**9 = 7.5x the initial charge). If it still fails,
    :class:`PyrogenAutoSizeError` is raised. A successful upsize is reported
    with a ``RuntimeWarning`` naming the final charge.

    An EXPLICITLY set charge is never resized: a run that is not a success
    returns normally with ``summary['ignited'] == False`` and a warning.

    At the end of EVERY run (auto or explicit) grain cells that never lit are
    reported with a ``RuntimeWarning`` (count, cell-index range, x range) and
    ``summary['grain_cells_never_lit']`` / ``['grain_never_lit_mask']`` /
    ``['grain_fraction_lit']`` / ``['grain_never_lit_cell_range']`` /
    ``['grain_never_lit_x_range_mm']``.

    ``summary`` also gains ``ignition_attempts``, ``ignition_attempt_masses_kg``,
    ``ignition_attempt_wall_s``, ``pyrogen_mass_auto``, ``pyrogen_mass_final_kg``
    and ``pyrogen_resized``; ``summary['ignited']`` is the success flag.
    """
    import inspect
    import warnings
    bound = inspect.signature(_run_simulation_once).bind(*args, **kwargs)
    chamber = bound.arguments['pyrogen_chamber']
    if ignition_retry_factor <= 1.0:
        raise ValueError("ignition_retry_factor must be > 1")
    if int(ignition_retry_max_attempts) < 1:
        raise ValueError("ignition_retry_max_attempts must be >= 1")
    auto = bool(getattr(chamber, 'mass_auto', False))
    progress = bound.arguments.get('progress_state')
    masses = []
    t_wall = []
    attempt = 0
    while True:
        attempt += 1
        masses.append(float(getattr(chamber, 'm_pyrogen_initial', 0.0)))
        bound.arguments['pyrogen_chamber'] = chamber
        t0 = clock.time()
        res = _run_once_growing_history(bound)
        t_wall.append(clock.time() - t0)
        summ = res['summary']
        ignited = _self_sustaining(res)
        if ignited or int(summ['termination_code']) in (3, 4, 5):
            break
        if not auto or attempt >= int(ignition_retry_max_attempts):
            break
        res = summ = None          # free the failed attempt's history arrays (~GB) before the next
        chamber = chamber.resize(masses[-1] * float(ignition_retry_factor))
        if progress is not None:       # restart the live progress bar
            progress[0] = 0.0
    summ['ignited'] = bool(ignited)
    summ['ignition_attempts'] = attempt
    summ['ignition_attempt_masses_kg'] = list(masses)
    summ['ignition_attempt_wall_s'] = list(t_wall)
    summ['pyrogen_mass_auto'] = auto
    summ['pyrogen_mass_final_kg'] = masses[-1]
    summ['pyrogen_resized'] = bool(auto and attempt > 1)
    name = chamber.pyrogen.name
    if auto and not ignited and int(summ['termination_code']) not in (3, 4, 5):
        raise PyrogenAutoSizeError(
            f"auto-sized {name} charge did not light the motor in {attempt} "
            f"attempts (charges tried: "
            f"{', '.join(f'{m * 1e3:.3g}' for m in masses)} g, factor "
            f"{ignition_retry_factor:g}). Check the igniter topology / "
            f"propellant, or set the charge explicitly.")
    if ignited and auto and attempt > 1:
        warnings.warn(
            f"auto-sized {name} charge raised from {masses[0] * 1e3:.3g} g "
            f"to {masses[-1] * 1e3:.3g} g ({attempt} attempts, factor "
            f"{ignition_retry_factor:g}) to light the motor. Set the "
            f"charge explicitly to fix it.", RuntimeWarning, stacklevel=2)
    if not ignited and not auto:
        hang = (" The run ended by the pressure cutoff with the chamber still "
                "sitting at the igniter's dribble pressure (a hang-fire)."
                if int(summ['termination_code']) == 2 else "")
        warnings.warn(
            f"the motor did not reach self-sustaining combustion with the "
            f"explicitly set {name} charge of {masses[-1] * 1e3:.3g} g "
            f"({summ['grain_cells_ignited']} of {summ['grain_cells']} grain "
            f"cells lit; grain gas generation "
            f"{summ['grain_mdot_after_igniter_kg_s']:.3g} kg/s after the "
            f"igniter vs the igniter's peak {summ['igniter_peak_mdot_kg_s']:.3g} "
            f"kg/s). Explicit charges are never resized. Termination: "
            f"{summ['termination']}.{hang}", RuntimeWarning, stacklevel=2)
    n_never = int(summ['grain_cells_never_lit'])
    if n_never > 0 and summ['grain_cells_ignited'] > 0:
        c0, c1 = summ['grain_never_lit_cell_range']
        x0, x1 = summ['grain_never_lit_x_range_mm']
        warnings.warn(
            f"{n_never} of {summ['grain_cells']} grain cells never ignited "
            f"(cells {c0}-{c1}, x = {x0:.0f}-{x1:.0f} mm; "
            f"{100.0 * summ['grain_fraction_lit']:.0f} % of the grain lit). "
            "This is the convective-only ignition model (no igniter-jet "
            "impingement heating; burning-neighbour radiation is off by "
            "default), not an igniter-size problem: a larger charge would "
            "mostly inflate the startup spike. Their burning is missing from "
            "this trace.", RuntimeWarning, stacklevel=2)
    return res


run_simulation.__doc__ = (run_simulation.__doc__ + "\n\n"
                          + (_run_simulation_once.__doc__ or ""))
import inspect as _inspect  # noqa: E402

# Expose the single-run signature (plus the two retry kwargs) so
# inspect.signature(run_simulation).parameters[...] defaults keep working (the
# openMotor plugin's solver-config defaults are generated from it, I4).
_once_sig = _inspect.signature(_run_simulation_once)
run_simulation.__signature__ = _once_sig.replace(parameters=list(
    _once_sig.parameters.values()) + [
        _inspect.Parameter('ignition_retry_factor',
                           _inspect.Parameter.KEYWORD_ONLY, default=1.25),
        _inspect.Parameter('ignition_retry_max_attempts',
                           _inspect.Parameter.KEYWORD_ONLY, default=10)])
