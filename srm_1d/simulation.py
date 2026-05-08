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

import math
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
from .burn_rate import compute_burn_rates, haaland_friction
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
N_SNAP_CHANNELS = 11


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
def _ignition_source_and_mass(
    P, has_ignited, is_burning, is_grain, ignition_time,
    r_total, r_erosive, mass_source,
    C_burn, endface_msource,
    ignition_ramp_tau, P_ignition,
    t, rho_propellant, N,
):
    """Hybrid Ignition: Gas-Dynamic Trigger + Thermal Soak Ramp."""
    n_burning = 0
    n_ignited = 0
    mass_sum = 0.0

    for i in range(N):
        if not is_grain[i]:
            is_burning[i] = False
            r_total[i] = 0.0
            r_erosive[i] = 0.0
            if has_ignited[i]:
                n_ignited += 1
        else:
            # 1. The Trigger: Pure local CFD pressure wave
            if not has_ignited[i]:
                if P[i] >= P_ignition:
                    has_ignited[i] = True
                    is_burning[i] = True
                    ignition_time[i] = t

            # 2. The Response: Thermal boundary layer development proxy
            if has_ignited[i]:
                n_ignited += 1
                elapsed = t - ignition_time[i]
                
                # Smooth asymptotic approach to 100% burn rate
                ramp = 1.0 - math.exp(-elapsed / ignition_ramp_tau)
                r_total[i] *= ramp
                r_erosive[i] *= ramp
                
                if is_burning[i]:
                    n_burning += 1

        if is_burning[i] and C_burn[i] > 0.0:
            mass_source[i] = rho_propellant * r_total[i] * C_burn[i]
        else:
            mass_source[i] = 0.0
            
        mass_source[i] += endface_msource[i]
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
    mass_source, f_darcy, Re, Mach, u_cell,
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
    rho_propellant, Cps, T_surface, T_initial,
    # --- Burn rate tabs (parallel arrays) ---
    tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        # --- Simulation parameters ---
    roughness, kappa,
    cfl_target, dt_max, burn_update_interval,
    igniter_mass_init, igniter_tau, n_ign_cells,  # <-- EDITED LINE
    P_ignition, ignition_ramp_tau,
    t_max, P_cutoff,
    erosion_coeff, slag_coeff, throat_is_evolving,
    snapshot_interval,
    # --- Precomputed ---
    gamma_R, Gamma_crit, nozzle_denom,
    D_throat_init, A_throat_init,
    # --- Output: time history (pre-allocated) ---
    time_hist, P_head_hist, P_exit_hist, D_throat_hist,
    Kn_hist, massflow_hist,
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
    igniter_mass_remaining = igniter_mass_init

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
        if n_ignited > 0 and n_burning == 0 and igniter_mass_remaining <= 0.0:
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
                r_total[i] = r_total_new[i]
                r_erosive[i] = r_erosive_new[i]

        # ============================================
        # STEP 3: IGNITION + SOURCE ASSEMBLY
        # ============================================

        n_burning, n_ignited, mass_sum = _ignition_source_and_mass(
            P, has_ignited, is_burning, is_grain, ignition_time,
            r_total, r_erosive, mass_source,
            C_burn, endface_msource,
            ignition_ramp_tau, P_ignition,
            t, rho_propellant, N,
        )

        # Add igniter: Decoupled Exponential Decay Model
        # mdot = (M_init / tau) * exp(-t / tau)
        if igniter_mass_remaining > 0.0:
            mdot_igniter = (igniter_mass_init / igniter_tau) * math.exp(-t / igniter_tau)

            # Don't consume more than remaining mass
            dm = mdot_igniter * dt
            if dm > igniter_mass_remaining:
                dm = igniter_mass_remaining
                mdot_igniter = dm / dt
            igniter_mass_remaining -= dm

            # Inject into head-end cells
            ign_per_cell = mdot_igniter / (n_ign_cells * dx)
            for j in range(n_ign_cells):
                mass_source[j] += ign_per_cell
            mass_sum += mdot_igniter / dx

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
            rho, u, P, T, A_port, D_hyd, mass_source, f_darcy,
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
            snap_idx += 1
            last_snapshot_t = t

        # Pressure cutoff (only after igniter consumed)
        if n_ignited > 0 and igniter_mass_remaining <= 0.0 and P[0] < P_cutoff:
            termination_code = 2
            break

        t += dt
        step += 1

    return (hist_idx, snap_idx,
            total_mass_produced, total_mass_nozzle,
            first_bt_time, D_throat, termination_code)


# ================================================================
# Public API
# ================================================================

def run_simulation(
    geo,
    propellant,
    nozzle,
    # --- Environment ---
    P_ambient=101325.0,
    # --- Surface / erosive burning ---
    roughness=50e-6,
    kappa=0.45,
    # --- Solver ---
    cfl_target=0.5,
    dt_max=0.002,
    burn_update_interval=None,
    # --- Igniter ---
    igniter_mass=0.010,
    igniter_tau=0.015,
    P_ignition=0.05e6,
    ignition_ramp_tau=0.010,
    # --- Termination ---
    t_max=10.0,
    P_cutoff=0.5e6,
    # --- Output ---
    print_interval=0.2,
    snapshot_interval=0.2,
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
    igniter_mass : float
        Total igniter propellant mass [kg]. Default: 1 g.
    igniter_tau : float
        Exponential decay time constant [s]. Default: 15 ms.
    P_ignition : float
        Pressure threshold for cell ignition [Pa]. Default: 0.05 MPa.
    ignition_ramp_tau : float
        Exponential ramp time constant [s]. Default: 10 ms.
    t_max : float
        Maximum simulation time [s]. Default: 10.
    P_cutoff : float
        Head-end pressure below which simulation terminates [Pa].
    print_interval : float
        Print status every this many seconds of simulated time.
    snapshot_interval : float
        Store flow field snapshots at this interval [s].

    Returns
    -------
    dict with 'time', 'P_head', 'P_exit', 'D_throat', 'snapshots', 'summary'.
    """
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
    RT = gas.R_specific * rep_tab.T_flame
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

    # Igniter: Decoupled exponential decay
    n_ign_cells = N # max(int(0.15 * N), 2)

    # Flow state
    P = np.full(N, 101325.0)
    rho = P / RT
    u = np.zeros(N + 1)
    T = np.full(N, rep_tab.T_flame)

    # Ignition state
    is_burning = np.zeros(N, dtype=np.bool_)
    has_ignited = np.zeros(N, dtype=np.bool_)
    ignition_time = np.full(N, 1e10)

    # Burn rates
    r_total = np.zeros(N)
    r_erosive = np.zeros(N)

    # Working arrays
    mass_source = np.zeros(N)
    f_darcy = np.zeros(N)
    Re = np.zeros(N)
    Mach = np.zeros(N)
    u_cell = np.zeros(N)

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

    # Pre-allocate snapshot storage
    max_snaps = int(t_max / snapshot_interval) + 10
    snap_data = np.empty((max_snaps, N_SNAP_CHANNELS, N))
    snap_times = np.empty(max_snaps)

    # ============================================================
    # STATUS PRINT
    # ============================================================
    numba_status = "Numba JIT enabled" if HAS_NUMBA else "Pure Python (no Numba)"
    print(f"PISO Solver ({numba_status}): {propellant.name}")
    print(f"  Motor: L={geo.L_motor*1e3:.0f}mm  D_outer={geo.D_outer*1e3:.0f}mm  "
          f"D_throat={nozzle.D_throat*1e3:.1f}mm  segments={ga['N_seg']}")
    print(f"  Params: Ts={propellant.T_surface:.0f}K  Cps={propellant.Cps:.0f}  "
          f"roughness={roughness*1e6:.0f}um  kappa={kappa}")
    print(f"  Igniter: m={igniter_mass*1e3:.1f}g  tau={igniter_tau*1000:.1f}ms")
    print()

    wall_start = clock.time()

    # ============================================================
    # RUN COMPILED TIME LOOP
    # ============================================================
    (n_steps, n_snaps,
     total_mass_produced, total_mass_nozzle,
     first_burnthrough_time, D_throat_final,
     termination_code) = _run_time_loop(
        # Cell arrays
        rho, u, P, T,
        D_port, x_centers, A_port, C_burn, D_hyd,
        is_grain, endface_msource,
        is_burning, has_ignited, ignition_time,
        r_total, r_erosive,
        mass_source, f_darcy, Re, Mach, u_cell,
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
        propellant.T_surface, propellant.T_initial,
        # Burn rate tabs
        tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
        # Simulation parameters
        roughness, kappa,
        cfl_target, dt_max, burn_update_interval,
        igniter_mass, igniter_tau, n_ign_cells,  # <-- EDITED LINE
        P_ignition, ignition_ramp_tau,
        t_max, P_cutoff,
        erosion_coeff, slag_coeff, throat_is_evolving,
        snapshot_interval,
        # Precomputed
        gamma_R, Gamma_crit, nozzle_denom,
        D_throat_init, A_throat_init,
        # Output: time history
        time_hist, P_head_hist, P_exit_hist, D_throat_hist,
        Kn_hist, massflow_hist,
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
        print(f"  Throat: {nozzle.D_throat*1e3:.2f} → {D_throat_final*1e3:.2f} mm "
              f"({direction} {abs(delta_mm):.3f} mm)")
    print(f"  Wall time: {wall_elapsed:.1f}s ({n_steps/max(wall_elapsed,0.01):.0f} steps/s)")
    print(f"  Mass: propellant={theoretical_propellant_mass:.3f}kg  "
          f"produced={total_mass_produced:.3f}kg  nozzle={total_mass_nozzle:.3f}kg  "
          f"balance_err={abs(total_mass_produced - total_mass_nozzle)/max(theoretical_propellant_mass,0.001)*100:.1f}%")
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
        'snapshots': snapshots, 'grains': grain_data,
        'summary': summary,
        'P_ambient': P_ambient,
    }
