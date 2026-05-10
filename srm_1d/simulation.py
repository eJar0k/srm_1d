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

from .solver import piso_step, compute_dt_cfl
from .burn_rate import compute_burn_rates, haaland_friction, gnielinski_nusselt
from .igniter_plenum import (
    _step_plenum_ode,
    chamber_params,
    initial_plenum_state,
    pyrogen_params,
)
from .solid_thermal import _step_goodman_ode, _surface_has_ignited
from .propellant import (
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
N_SNAP_CHANNELS = 14


# ================================================================
# Fused per-step helpers (called from inside _run_time_loop)
# ================================================================

@njit(cache=True)
def _post_piso_update(
    rho, u, P, T, D_hyd, Re, Mach, u_cell, f_darcy,
    N, mu_gas, gamma_R, roughness,
):
    """Post-PISO: velocities, Re, Mach, friction, a_max — single pass."""
    T_max = T[0]
    for i in range(N):
        u_cell[i] = 0.5 * (u[i] + u[i + 1])
        a_local = (gamma_R * T[i]) ** 0.5
        Mach[i] = u_cell[i] / a_local
        Re[i] = rho[i] * abs(u_cell[i]) * D_hyd[i] / mu_gas
        f_darcy[i] = haaland_friction(Re[i], roughness, D_hyd[i])
        if T[i] > T_max:
            T_max = T[i]
    return (gamma_R * T_max) ** 0.5


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
def _goodman_ignition_sources_and_mass(
    P, T, T_surf, delta, has_ignited, is_burning, is_grain, ignition_time,
    r_total, r_erosive, mass_source, thermal_source,
    C_burn, endface_msource,
    x_centers, Re, D_hyd, f_darcy,
    t, dt, rho_propellant, T_flame, T_initial,
    Pr, k_thermal, roughness, kappa, solid_alpha, k_solid,
    T_ignition, N,
):
    """Goodman surface-temperature ignition and propellant source assembly."""
    n_burning = 0
    n_ignited = 0
    mass_sum = 0.0

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
                new_delta, new_T_surf = _step_goodman_ode(
                    delta[i], T_surf[i], h_c, T[i], T_initial,
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
            prop_source = rho_propellant * r_total[i] * C_burn[i]
            mass_source[i] += prop_source
            thermal_source[i] += prop_source * T_flame

        if endface_msource[i] > 0.0:
            mass_source[i] += endface_msource[i]
            thermal_source[i] += endface_msource[i] * T_flame

        mass_sum += mass_source[i]

    return n_burning, n_ignited, mass_sum


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
    mass_source, thermal_source, f_darcy, Re, Mach, u_cell,
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
    rho_propellant, Cps, T_surface, T_initial, k_solid,
    # --- Burn rate tabs (parallel arrays) ---
    tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        # --- Simulation parameters ---
    roughness, kappa,
    cfl_target, dt_max, burn_update_interval,
    T_ignition,
    diagnostic_disable_erosive, diagnostic_disable_endfaces,
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
    max_hist,
    # --- Output: snapshots (pre-allocated) ---
    snap_data, snap_times, max_snaps,
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
        3 = history array full.
    """
    PI = 3.141592653589793

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

    # Initial a_max for CFL
    T_max = T[0]
    for i in range(N):
        if T[i] > T_max:
            T_max = T[i]
    a_max = (gamma_R * T_max) ** 0.5

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

        n_burning, n_ignited, mass_sum = _goodman_ignition_sources_and_mass(
            P, T, T_surf, delta, has_ignited, is_burning, is_grain,
            ignition_time, r_total, r_erosive,
            mass_source, thermal_source,
            C_burn, endface_msource,
            x_centers, Re, D_hyd, f_darcy,
            t, dt, rho_propellant, T_flame, T_initial,
            Pr, k_thermal, roughness, kappa, solid_alpha, k_solid,
            T_ignition, N,
        )

        if mdot_igniter > 0.0:
            ign_source = mdot_igniter / dx
            mass_source[0] += ign_source
            thermal_source[0] += ign_source * T_ig
            mass_sum += ign_source

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

        # ============================================
        # STEP 4: PISO
        # ============================================
        rho, u, P, T = piso_step(
            rho, u, P, T, A_port, D_hyd, mass_source, thermal_source, f_darcy,
            dx, dt, gamma, R_specific, T_flame,
            Cp_gas, A_throat, N,
        )

        # ============================================
        # STEP 5: POST-PISO
        # ============================================
        a_max = _post_piso_update(
            rho, u, P, T, D_hyd, Re, Mach, u_cell, f_darcy,
            N, mu_gas, gamma_R, roughness,
        )

        # ============================================
        # STEP 6: BOOKKEEPING
        # ============================================
        total_mass_produced += mass_sum * dx * dt
        nozzle_mdot = P[N - 1] * A_throat * Gamma_crit * nozzle_denom
        total_mass_nozzle += nozzle_mdot * dt

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
            snap_idx += 1
            last_snapshot_t = t

        # Pressure cutoff (only after pyrogen is consumed and vented)
        if n_ignited > 0 and pyrogen_done and P[0] < P_cutoff:
            termination_code = 2
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
    # --- Surface / erosive burning ---
    roughness=50e-6,
    kappa=0.45,
    # --- Solver ---
    cfl_target=0.5,
    dt_max=0.002,
    burn_update_interval=None,
    # --- Ignition ---
    T_ignition=850.0,
    # --- Diagnostics ---
    initial_gas_temperature=None,
    diagnostic_disable_erosive=False,
    diagnostic_disable_endfaces=False,
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
    if initial_gas_temperature is not None and initial_gas_temperature <= 0.0:
        raise ValueError("initial_gas_temperature must be positive")

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

    # Ignition state
    is_burning = np.zeros(N, dtype=np.bool_)
    has_ignited = np.zeros(N, dtype=np.bool_)
    ignition_time = np.full(N, 1e10)

    # Burn rates
    r_total = np.zeros(N)
    r_erosive = np.zeros(N)

    # Working arrays
    mass_source = np.zeros(N)
    thermal_source = np.zeros(N)
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

    theoretical_propellant_mass = (
        geo.total_propellant_volume() * propellant.rho_propellant
    )

    # Pre-allocate output arrays
    # Conservative estimate: dt_min ~ cfl * dx / 1000 m/s
    est_steps = int(t_max / max(cfl_target * dx / 1000.0, 1e-8)) + 1000
    max_hist = max(est_steps, 5_000_000)
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
        mass_source, thermal_source, f_darcy, Re, Mach, u_cell,
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
        propellant.T_surface, propellant.T_initial, propellant.k_solid,
        # Burn rate tabs
        tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        # Simulation parameters
        roughness, kappa,
        cfl_target, dt_max, burn_update_interval,
        T_ignition,
        bool(diagnostic_disable_erosive), bool(diagnostic_disable_endfaces),
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
        max_hist,
        # Output: snapshots
        snap_data, snap_times, max_snaps,
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
        })

    peak_idx = np.argmax(P_head_arr) if len(P_head_arr) > 0 else 0
    cstar = np.sqrt(gas.R_specific * rep_tab.T_flame) / Gamma_crit

    if first_burnthrough_time < 0:
        first_burnthrough_time = None

    # ============================================================
    # SUMMARY
    # ============================================================
    termination_names = {
        0: "t_max reached", 1: "complete burnout",
        2: "pressure cutoff", 3: "history array full",
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
        'diagnostic_disable_erosive': bool(diagnostic_disable_erosive),
        'diagnostic_disable_endfaces': bool(diagnostic_disable_endfaces),
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
        'snapshots': snapshots, 'grains': grain_data,
        'summary': summary,
        'P_ambient': P_ambient,
    }
