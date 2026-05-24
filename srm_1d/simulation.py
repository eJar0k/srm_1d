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
    _piso_step_with_energy_diagnostics,
    compute_dt_cfl,
    compute_dt_source_cap,
    _nozzle_boundary_flow,
)
from .burn_rate import compute_burn_rates, haaland_friction, gnielinski_nusselt
from .igniter_plenum import (
    _step_plenum_ode,
    chamber_params,
    initial_plenum_state,
    pyrogen_params,
)
from .solid_thermal import _step_goodman_ode, _surface_has_ignited
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
N_SNAP_CHANNELS = 17

CAL_CM2_S_TO_W_M2 = 41840.0
STEFAN_BOLTZMANN = 5.670374419e-8

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

@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
def _cal_cm2_s_to_w_m2(heat_flux_cal_cm2_s):
    """Convert cal/(cm^2*s) to W/m^2."""
    return heat_flux_cal_cm2_s * CAL_CM2_S_TO_W_M2


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
def _gas_sensible_energy(rho, T, A_port, dx, Cp_arr, N):
    """Discrete gas sensible energy used by diagnostics [J].

    v0.7.1 (Phase 3): ``Cp_arr`` is per-cell so cells with different
    mixtures contribute their own Cp to the total.
    """
    total = 0.0
    for i in range(N):
        total += rho[i] * A_port[i] * dx * Cp_arr[i] * T[i]
    return total


@njit(cache=True)
def _thermal_source_power(thermal_source, dx, N):
    """Sum solver enthalpy-source units to thermal power [W].

    v0.7.1 (Phase 3): ``thermal_source`` is already in W/m units (enthalpy
    injection per unit length), so the sum is ``thermal_source[i] * dx``.
    """
    total = 0.0
    for i in range(N):
        total += thermal_source[i] * dx
    return total


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
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


@njit(cache=True)
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
    pyrogen_heat_target = -1
    head_grain_cell = -1

    for i in range(N):
        pyrogen_surface_heat_flux[i] = 0.0
        radiation_heat_flux[i] = 0.0
        radiation_sink_power[i] = 0.0
        radiation_emitter[i] = is_burning[i]
        if head_grain_cell < 0 and is_grain[i] and C_burn[i] > 0.0:
            head_grain_cell = i

    if (head_grain_cell >= 0 and mdot_igniter > 0.0 and
            pyrogen_surface_heat_flux_w_m2 > 0.0 and
            (not has_ignited[head_grain_cell])):
        pyrogen_heat_target = head_grain_cell

    # v0.7.1: reset per-species mass source rows. Each step recomputes
    # contributions from grain (s=1) inside this loop; pyrogen (s=0) is
    # set after this function returns. Other species remain zero.
    S_local = mass_source_by_species.shape[1]
    for i in range(N):
        for s in range(S_local):
            mass_source_by_species[i, s] = 0.0

    for i in range(N):
        mass_source[i] = 0.0
        thermal_source[i] = 0.0

        if not is_grain[i]:
            is_burning[i] = False
            r_total[i] = 0.0
            r_erosive[i] = 0.0
        else:
            if not has_ignited[i]:
                h_c = _bare_heat_transfer_coeff(
                    Re[i], D_hyd[i], x_centers[i], f_darcy[i],
                    Pr, k_thermal, T[i], T_surf[i], kappa,
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
                if i == pyrogen_heat_target:
                    power_w, flux_w_m2 = _pyrogen_surface_heat_power(
                        mdot_igniter, T_ig, T_surf[i], C_burn[i], dx,
                        Cp_pyrogen, pyrogen_surface_heat_flux_w_m2,
                    )
                    if power_w > 0.0 and flux_w_m2 > 0.0:
                        pyrogen_surface_heat_power = power_w
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

                new_delta, new_T_surf = _step_goodman_ode(
                    delta[i], T_surf[i], h_c, h_driver, T_initial,
                    solid_alpha, k_solid, dt,
                )
                delta[i] = new_delta
                T_surf[i] = new_T_surf

                if _surface_has_ignited(T_surf[i], T_ignition):
                    has_ignited[i] = True
                    is_burning[i] = True
                    ignition_time[i] = t

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

        if radiation_sink_power[i] > 0.0 and dx > 0.0:
            thermal_source[i] -= radiation_sink_power[i] / dx

        mass_sum += mass_source[i]

    return (n_burning, n_ignited, mass_sum,
            pyrogen_surface_heat_power, radiation_heat_power,
            radiation_sink_total_power,
            normal_sidewall_thermal_power, erosive_sidewall_thermal_power,
            endface_thermal_power)


# ================================================================
# Compiled time loop
# ================================================================

@njit(cache=True)
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
    N, N_seg, dx, D_outer,
    # --- Gas/propellant scalars ---
    gamma, R_specific, T_flame, Cp_gas, mu_gas, k_thermal, Pr,
    rho_propellant, Cps, T_surface, T_initial, T_initial_gas, k_solid,
    # --- Burn rate tabs (parallel arrays) ---
    tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        # --- Simulation parameters ---
    roughness, kappa,
    cfl_target, dt_max, burn_update_interval,
    source_cfl_factor,
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
    max_hist,
    # --- Output: snapshots (pre-allocated) ---
    snap_data, snap_times, max_snaps,
    # --- v0.7.1: N-species state ---
    Y_species, species_params_arr, mass_source_by_species,
    gamma_mix_arr, Cp_mix_arr, R_mix_arr, M_mix_arr,
    T_ceiling_arr,
    # --- v0.7.2 Phase A: pyrogen axial distribution ---
    pyrogen_axial_weights,
    # --- v0.7.2 Phase B-v2: flame-front h_c augmentation ---
    flame_spread_augment, flame_spread_enabled,
    flame_spread_tau, flame_spread_boost,
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
        3 = history array full, 4 = numerical collapse aborted.
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
    first_bt_time = -1.0
    last_snapshot_t = -snapshot_interval
    termination_code = 0

    D_throat = D_throat_init
    A_throat = A_throat_init
    pyrogen_initial_mass = plenum_state[0]
    pyrogen_done = False
    pyrogen_peak_P = 0.0
    pyrogen_duration = 0.0
    solid_alpha = 0.0
    if k_solid > 0.0 and rho_propellant > 0.0 and Cps > 0.0:
        solid_alpha = k_solid / (rho_propellant * Cps)

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

    while t < t_max:
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

        if step % burn_update_interval == 0:
            update_cell_geometry(
                regress, D_port, x_centers, dx, N, N_seg, D_outer,
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
            if diagnostic_disable_endfaces:
                for i in range(N):
                    endface_msource[i] = 0.0

        # ============================================
        # STEP 2: BURN RATES
        # ============================================
        if step % burn_update_interval == 0:
            r_total_new, r_erosive_new = compute_burn_rates(
                P, Re, D_hyd, x_centers, is_burning, roughness,
                Pr, k_thermal, Cp_gas,
                T_flame, T_surface,
                rho_propellant, Cps, T_initial,
                tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
                kappa, N,
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

        # ============================================
        # STEP 3: IGNITION + SOURCE ASSEMBLY
        # ============================================

        new_plenum_state, mdot_igniter, _mdot_generated, P_ig = _step_plenum_ode(
            plenum_state, pyrogen_params_arr, chamber_params_arr, dt, P[0]
        )
        plenum_state[0] = new_plenum_state[0]
        plenum_state[1] = new_plenum_state[1]
        plenum_state[2] = new_plenum_state[2]
        T_ig = plenum_state[2]
        if P_ig > pyrogen_peak_P:
            pyrogen_peak_P = P_ig
        if mdot_igniter > 1e-12:
            pyrogen_duration = t + dt
        pyrogen_done = plenum_state[0] <= 1e-12 and mdot_igniter <= 1e-9
        pyrogen_momentum_expected = 0.0
        pyrogen_momentum_deposited = 0.0
        for i in range(N + 1):
            momentum_source[i] = 0.0
        if (not diagnostic_disable_momentum and mdot_igniter > 0.0 and
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

        # v0.7.2 Phase B-v2: flame-front augmentation. Always computed
        # (cheap, O(N)); the multiply is gated INSIDE the Goodman kernel
        # by flame_spread_enabled so disabled runs reduce to Phase A
        # byte-for-byte. The kernel fills augment_arr[i]=boost for cells
        # immediately downstream of a recently-ignited cell, else 1.0.
        _compute_flame_front_augment(
            is_burning, has_ignited, ignition_time, t,
            flame_spread_tau, flame_spread_boost, N, flame_spread_augment,
        )

        (n_burning, n_ignited, mass_sum,
         pyrogen_surface_heat_power, radiation_heat_power,
         radiation_sink_total_power,
         normal_sidewall_thermal_power, erosive_sidewall_thermal_power,
         endface_thermal_power) = _goodman_ignition_sources_and_mass(
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
            # enthalpy share (not the full plenum enthalpy as in v0.7.1)
            # because cell 0 only receives w[0] * ign_enthalpy of the
            # bulk pyrogen enthalpy after distribution.
            ign_enthalpy_cell0 = (
                pyrogen_axial_weights[0] * mdot_igniter
                * Cp_pyrogen_species * T_ig / dx
            )
            pyrogen_surface_heat_sink = _pyrogen_surface_thermal_sink(
                pyrogen_surface_heat_power, dx, ign_enthalpy_cell0
            )
            thermal_source[0] -= pyrogen_surface_heat_sink
            pyrogen_surface_heat_sink_power = pyrogen_surface_heat_sink * dx

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
         energy_residual) = _piso_step_with_energy_diagnostics(
            rho, u, P, T, A_port, D_hyd,
            mass_source, thermal_source, momentum_source, f_darcy,
            dx, dt, gamma_mix_arr, R_mix_arr, Cp_mix_arr, T_ceiling_arr,
            A_throat, P_ambient, ambient_temperature, N,
        )

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

        # v0.7.1: species advection using post-PISO (u, rho) and the
        # pre-PISO rho snapshot. mass_source_by_species was populated
        # earlier (pyrogen at cell 0, propellant grain at active cells).
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
                if is_grain[i] and D_port[i] >= D_outer:
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
            snap_idx += 1
            last_snapshot_t = t

        # Pressure cutoff (only after pyrogen is consumed and vented)
        if n_ignited > 0 and pyrogen_done and P[0] < P_cutoff:
            termination_code = 2
            break

        # Numerical-collapse trip: classified abort to avoid burning the
        # full history budget on a doomed run. Records the collapse step
        # in history (so the user sees the trip in diagnostics) and exits.
        if collapse_consecutive >= COLLAPSE_CONSECUTIVE_STEPS:
            termination_code = 4
            break

        t += dt
        step += 1

    return (hist_idx, snap_idx,
            total_mass_produced, total_mass_nozzle,
            first_bt_time, D_throat, termination_code,
            pyrogen_initial_mass - plenum_state[0],
            pyrogen_duration, pyrogen_peak_P)


# ================================================================
# Public API
# ================================================================

def run_simulation(
    geo,
    propellant,
    nozzle,
    pyrogen_chamber,
    # --- Environment ---
    P_ambient=101325.0,
    ambient_temperature=None,
    # --- Surface / erosive burning ---
    roughness=50e-6,
    kappa=0.45,
    # --- Solver ---
    cfl_target=0.5,
    dt_max=0.002,
    source_cfl_factor=0.10,
    burn_update_interval=None,
    # --- Ignition ---
    T_ignition=850.0,
    tau_establishment=0.0,
    # --- Diagnostics ---
    initial_gas_temperature=None,
    diagnostic_disable_erosive=False,
    diagnostic_disable_endfaces=False,
    diagnostic_disable_momentum=False,
    diagnostic_disable_pyrogen_surface_heating=False,
    diagnostic_disable_adjacent_radiation=False,
    diagnostic_disable_radiation_gas_sink=False,
    diagnostic_history_capacity=None,
    igniter_axial_momentum_fraction=1.0,
    # --- Termination ---
    t_max=10.0,
    P_cutoff=0.5e6,
    # --- Output ---
    print_interval=0.2,
    snapshot_interval=0.2,
    verbose=True,
):
    """
    Run a complete transient simulation.

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
    cfl_target : float
        Target CFL number for time stepping. Default: 0.5.
    dt_max : float
        Maximum allowed time step [s]. Default: 0.002.
    burn_update_interval : int or None
        Recompute burn rates every N flow steps. If None, auto-set.
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
    diagnostic_history_capacity : int or None
        Optional diagnostic-only cap for preallocated history rows. This
        does not change equations or time stepping; it only allows probe
        runs to terminate earlier through the normal "history array full"
        path.
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
        regress, D_port, x_centers, dx, N, ga['N_seg'], ga['D_outer'],
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
    T_initial_gas = rep_tab.T_flame
    if initial_gas_temperature is not None:
        T_initial_gas = float(initial_gas_temperature)
    rho = P / (gas.R_specific * T_initial_gas)
    u = np.zeros(N + 1)
    T = np.full(N, T_initial_gas)

    # v0.7.1: per-cell species mass fractions. Y_species[i, s] is the
    # mass fraction of species s in cell i. Invariant: sum_s Y[i, s] = 1.
    # Initial condition: all cells start as 100% ambient (species 2).
    # The species composition tracks separately from temperature; the
    # numerical-stability shortcut of initializing T at the propellant
    # T_flame is preserved in T_initial_gas above.
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

    theoretical_propellant_mass = (
        geo.total_propellant_volume() * propellant.rho_propellant
    )

    # Pre-allocate output arrays
    # Conservative estimate: dt_min ~ cfl * dx / 1000 m/s
    est_steps = int(t_max / max(cfl_target * dx / 1000.0, 1e-8)) + 1000
    max_hist = max(est_steps, 5_000_000)
    if diagnostic_history_capacity is not None:
        diagnostic_history_capacity = int(diagnostic_history_capacity)
        if diagnostic_history_capacity < 1:
            raise ValueError("diagnostic_history_capacity must be positive")
        max_hist = min(max_hist, diagnostic_history_capacity)
    time_hist = np.empty(max_hist)
    P_head_hist = np.empty(max_hist)
    P_exit_hist = np.empty(max_hist)
    D_throat_hist = np.empty(max_hist)
    Kn_hist = np.empty(max_hist)
    massflow_hist = np.empty(max_hist)
    P_ig_hist = np.empty(max_hist)
    T_ig_hist = np.empty(max_hist)
    mdot_ig_hist = np.empty(max_hist)
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

    # Pre-allocate snapshot storage
    max_snaps = int(t_max / snapshot_interval) + 10
    snap_data = np.empty((max_snaps, N_SNAP_CHANNELS, N))
    snap_times = np.empty(max_snaps)

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

    # ============================================================
    # RUN COMPILED TIME LOOP
    # ============================================================
    (n_steps, n_snaps,
     total_mass_produced, total_mass_nozzle,
     first_burnthrough_time, D_throat_final,
     termination_code,
     pyrogen_mass_burned, pyrogen_duration,
     pyrogen_peak_P) = _run_time_loop(
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
        N, ga['N_seg'], dx, ga['D_outer'],
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
        cfl_target, dt_max, burn_update_interval,
        source_cfl_factor,
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
        max_hist,
        # Output: snapshots
        snap_data, snap_times, max_snaps,
        # v0.7.1: N-species state
        Y_species, species_params_arr, mass_source_by_species,
        gamma_mix_arr, Cp_mix_arr, R_mix_arr, M_mix_arr,
        T_ceiling_arr,
        # v0.7.2 Phase A: pyrogen axial distribution
        pyrogen_axial_weights,
        # v0.7.2 Phase B-v2: flame-front h_c augmentation
        flame_spread_augment, flame_spread_enabled,
        flame_spread_tau, flame_spread_boost,
    )

    wall_elapsed = clock.time() - wall_start

    # ============================================================
    # TRIM AND WRAP RESULTS
    # ============================================================
    time_arr = time_hist[:n_steps].copy()
    P_head_arr = P_head_hist[:n_steps].copy()
    P_exit_arr = P_exit_hist[:n_steps].copy()
    D_throat_arr = D_throat_hist[:n_steps].copy()
    Kn_arr = Kn_hist[:n_steps].copy()
    massflow_arr = massflow_hist[:n_steps].copy()
    P_ig_arr = P_ig_hist[:n_steps].copy()
    T_ig_arr = T_ig_hist[:n_steps].copy()
    mdot_ig_arr = mdot_ig_hist[:n_steps].copy()
    m_pyrogen_arr = m_pyrogen_hist[:n_steps].copy()
    gas_sensible_energy_before_arr = gas_sensible_energy_before_hist[:n_steps].copy()
    gas_sensible_energy_arr = gas_sensible_energy_hist[:n_steps].copy()
    gas_sensible_dE_dt_arr = gas_sensible_dE_dt_hist[:n_steps].copy()
    normal_sidewall_thermal_power_arr = normal_sidewall_thermal_power_hist[:n_steps].copy()
    erosive_sidewall_thermal_power_arr = erosive_sidewall_thermal_power_hist[:n_steps].copy()
    endface_thermal_power_arr = endface_thermal_power_hist[:n_steps].copy()
    convective_scalar_flux_power_arr = convective_scalar_flux_power_hist[:n_steps].copy()
    clipping_correction_power_arr = clipping_correction_power_hist[:n_steps].copy()
    pyrogen_enthalpy_power_arr = pyrogen_enthalpy_power_hist[:n_steps].copy()
    pyrogen_surface_heat_power_arr = pyrogen_surface_heat_power_hist[:n_steps].copy()
    gas_surface_heat_sink_power_arr = gas_surface_heat_sink_power_hist[:n_steps].copy()
    radiation_heat_power_arr = radiation_heat_power_hist[:n_steps].copy()
    radiation_sink_power_arr = radiation_sink_power_hist[:n_steps].copy()
    nozzle_enthalpy_power_arr = nozzle_enthalpy_power_hist[:n_steps].copy()
    thermal_source_power_arr = thermal_source_power_hist[:n_steps].copy()
    energy_residual_arr = energy_residual_hist[:n_steps].copy()
    pyrogen_momentum_expected_arr = pyrogen_momentum_expected_hist[:n_steps].copy()
    pyrogen_momentum_deposited_arr = pyrogen_momentum_deposited_hist[:n_steps].copy()
    pyrogen_momentum_residual_arr = pyrogen_momentum_residual_hist[:n_steps].copy()
    dt_arr = dt_hist[:n_steps].copy()
    n_burning_arr = n_burning_hist[:n_steps].copy()
    n_ignited_arr = n_ignited_hist[:n_steps].copy()
    radiation_emitter_count_arr = radiation_emitter_count_hist[:n_steps].copy()
    radiation_receiver_count_arr = radiation_receiver_count_hist[:n_steps].copy()
    min_gas_temperature_arr = min_gas_temperature_hist[:n_steps].copy()
    max_gas_temperature_arr = max_gas_temperature_hist[:n_steps].copy()
    min_surface_temperature_arr = min_surface_temperature_hist[:n_steps].copy()
    max_surface_temperature_arr = max_surface_temperature_hist[:n_steps].copy()
    min_pressure_arr = min_pressure_hist[:n_steps].copy()
    max_pressure_arr = max_pressure_hist[:n_steps].copy()
    max_mach_arr = max_mach_hist[:n_steps].copy()

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
        4: "numerical collapse aborted",
    }
    term_str = termination_names.get(termination_code, "unknown")

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
        'P_peak': float(P_head_arr[peak_idx]),
        't_peak': float(time_arr[peak_idx]),
        'P_mid': P_mid,
        't_burn': float(time_arr[-1]) if len(time_arr) > 0 else 0.0,
        't_first_burnout': first_burnthrough_time,
        'c_star': cstar,
        'wall_time': wall_elapsed,
        'steps': n_steps,
        'cells': N,
        'history_capacity': int(max_hist),
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
    }

    # Per-grain summary from snapshots
    grain_data = []
    for seg_idx in range(ga['N_seg']):
        seg_cells = ga['cell_segment_id'] == seg_idx
        # Average initial bore for this segment's cells
        D_bore_init_avg = np.mean(ga['cell_D_bore_init'][seg_cells]) if np.any(seg_cells) else geo.D_outer
        regression_hist = []
        web_hist = []
        for snap in snapshots:
            D_seg = snap['D_port'][seg_cells]
            if len(D_seg) > 0:
                avg_D = np.mean(D_seg)
                reg = (avg_D - D_bore_init_avg) / 2.0
                web = (geo.D_outer - avg_D) / 2.0
            else:
                reg = (geo.D_outer - D_bore_init_avg) / 2.0
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

    return {
        'time': time_arr, 'P_head': P_head_arr, 'P_exit': P_exit_arr,
        'D_throat': D_throat_arr, 'Kn': Kn_arr, 'massflow': massflow_arr,
        'P_ig': P_ig_arr, 'T_ig': T_ig_arr, 'mdot_ig': mdot_ig_arr,
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
        'ignition_time_by_cell': ignition_time.copy(),
        'snapshots': snapshots, 'grains': grain_data,
        'summary': summary,
        'P_ambient': P_ambient,
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
