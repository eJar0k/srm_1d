"""
solver.py — 1D Compressible PISO Solver on a Staggered Grid
=============================================================

PURPOSE:
    Solves the 1D unsteady Navier-Stokes equations (mass, momentum,
    energy) for compressible flow using the PISO algorithm. This module
    contains ONLY the numerical building blocks — it knows nothing about
    propellant, burn rates, or grain geometry. Those are handled by
    burn_rate.py and simulation.py.

THE PISO ALGORITHM (per time step):
    1. MOMENTUM PREDICTOR — Solve for u* using the old pressure field.
       Semi-implicit: new velocity, old pressure.

    2. PRESSURE CORRECTION 1 — u* does not satisfy continuity. Derive
       a correction P' from the continuity residual, yielding a
       tridiagonal system solved by the Thomas algorithm (TDMA).
       Update: P = P_old + P',  u = u* - (dt/ρ) · ∇P'

    3. PRESSURE CORRECTION 2 — Repeat correction on the updated field.
       This second pass captures coupling terms that the first correction
       misses, giving better temporal accuracy without outer iterations.
       This is PISO's key advantage over SIMPLE, which instead repeats
       the full predictor-corrector until convergence.

    4. ENERGY EQUATION — Advect temperature with the corrected velocity
       field. Mass source cells inject gas at the flame temperature.

    5. EQUATION OF STATE — Update density: ρ = P / (R·T).

GRID LAYOUT:
    Staggered (MAC) variable arrangement:

        Cell centers: |  0  |  1  |  2  | ... | N-1 |
        Cell faces:  0     1     2     3    ...  N-1    N
                     ^                                  ^
                 head end                            nozzle
                 (wall)                              (choked)

    Scalars (P, ρ, T) live at cell centers — arrays of length N.
    Velocities (u) live at cell faces — array of length N+1.
    Face j sits between cells j-1 and j.

    The staggered arrangement prevents pressure-velocity decoupling
    (checkerboarding). Each face's pressure gradient uses its two
    immediately adjacent cells, so odd and even cells cannot decouple.
    This eliminates the need for Rhie-Chow interpolation.

BOUNDARY CONDITIONS:
    Head end (face 0): Solid wall — u[0] = 0.
    Nozzle end (face N): Choked outflow — mass flux proportional to
    local pressure via the isentropic mass flow relation.

TIME STEPPING:
    Adaptive CFL-based: dt = CFL · dx / max(|u| + a_sound)

NUMBA JIT:
    All functions are compiled with @njit for C-level speed. This
    requires all data to be passed as arrays/scalars — no Python
    objects inside JIT boundaries.

REFERENCES:
    Issa, R.I. (1986). "Solution of the Implicitly Discretised Fluid
    Flow Equations by Operator-Splitting." J. Comput. Phys., 62.

    Moukalled, F., Mangani, L., Darwish, M. (2016). The Finite Volume
    Method in Computational Fluid Dynamics. Springer. Ch. 15.
"""

import numpy as np

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


# ================================================================
# TDMA (Thomas Algorithm)
# ================================================================

@njit(cache=True)
def thomas_solve(a, b, c, d, N):
    """
    Solve a tridiagonal system using the Thomas algorithm (TDMA).

    Solves: a[i]·x[i-1] + b[i]·x[i] + c[i]·x[i+1] = d[i]

    This is the workhorse of the pressure correction step. The pressure
    correction equation, when discretized on the staggered grid, produces
    a tridiagonal system because each cell's pressure correction couples
    only to its immediate east and west neighbors.

    Parameters
    ----------
    a : ndarray (N,)
        Sub-diagonal coefficients. a[0] is unused (no left neighbor
        for the first cell).
    b : ndarray (N,)
        Main diagonal coefficients.
    c : ndarray (N,)
        Super-diagonal coefficients. c[N-1] is unused.
    d : ndarray (N,)
        Right-hand side vector.
    N : int
        System size.

    Returns
    -------
    x : ndarray (N,)
        Solution vector.

    Notes
    -----
    TDMA is O(N) — it exploits the tridiagonal structure to avoid the
    O(N³) cost of general Gaussian elimination. For our 150-cell grid,
    this means ~300 operations instead of ~3.4 million.

    The algorithm has two phases:
    1. Forward sweep: eliminate the sub-diagonal by modifying the
       diagonal and RHS, working left to right.
    2. Back substitution: solve for x working right to left.

    Numerically stable as long as the system is diagonally dominant,
    which the pressure correction equation guarantees (the diagonal
    coefficient includes contributions from both faces plus the
    transient density term).
    """
    x = np.zeros(N)
    # Forward sweep
    c_prime = np.zeros(N)
    d_prime = np.zeros(N)
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    for i in range(1, N):
        m = a[i] / (b[i] - a[i] * c_prime[i - 1])
        c_prime[i] = c[i] / (b[i] - a[i] * c_prime[i - 1])
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / (b[i] - a[i] * c_prime[i - 1])
    # Back substitution
    x[N - 1] = d_prime[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return x


# ================================================================
# PISO Time Step
# ================================================================

@njit(cache=True)
def piso_step(
    rho, u, P, T, A_port, D_hyd, mass_source, f_darcy,
    dx, dt, gamma, R_specific, T_flame, Cp_gas, A_throat, N,
):
    """
    One complete PISO time step on a staggered grid.

    This function advances the flow field by one time step. It takes the
    current state (ρ, u, P, T) and source terms (mass_source, f_darcy)
    and returns the updated state. It has no knowledge of where the mass
    source or friction come from — that is the simulation driver's job.

    Parameters
    ----------
    rho : ndarray (N,)
        Cell-center density [kg/m³].
    u : ndarray (N+1,)
        Face velocities [m/s]. u[0] = head-end wall, u[N] = nozzle.
    P : ndarray (N,)
        Cell-center pressure [Pa].
    T : ndarray (N,)
        Cell-center temperature [K].
    A_port : ndarray (N,)
        Cell-center port cross-sectional area [m²].
    D_hyd : ndarray (N,)
        Cell-center hydraulic diameter [m].
    mass_source : ndarray (N,)
        Mass source per unit length [kg/(m·s)]. From propellant
        combustion and igniter.
    f_darcy : ndarray (N,)
        Darcy friction factor at cell centers [-].
    dx : float
        Cell width [m].
    dt : float
        Time step size [s].
    gamma : float
        Ratio of specific heats [-].
    R_specific : float
        Specific gas constant [J/(kg·K)].
    T_flame : float
        Flame temperature [K]. Temperature of injected mass.
    Cp_gas : float
        Specific heat at constant pressure [J/(kg·K)].
    A_throat : float
        Nozzle throat area [m²].
    N : int
        Number of cells.

    Returns
    -------
    rho_new : ndarray (N,)
        Updated density.
    u_new : ndarray (N+1,)
        Updated face velocities.
    P_new : ndarray (N,)
        Updated pressure.
    T_new : ndarray (N,)
        Updated temperature.
    """
    # Choked flow function: Γ = √γ · (2/(γ+1))^((γ+1)/(2(γ-1)))
    # Computed inline to keep solver.py dependency-free.
    Gamma_crit = np.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** (
        (gamma + 1.0) / (2.0 * (gamma - 1.0))
    )
    nozzle_coeff = A_throat * Gamma_crit / np.sqrt(R_specific * T[N - 1])

    # Face areas (interpolated from cell centers)
    A_face = np.zeros(N + 1)
    A_face[0] = A_port[0]
    for j in range(1, N):
        A_face[j] = 0.5 * (A_port[j - 1] + A_port[j])
    A_face[N] = A_port[N - 1]

    # -------------------------------------------------------
    # STEP 1: MOMENTUM PREDICTOR (at faces j = 1 .. N-1)
    # -------------------------------------------------------
    # For face j between cells j-1 and j:
    #   ρ_f · (u* - u^n)/dt = -(P[j] - P[j-1])/dx + convection + friction
    #
    # The pressure gradient (P[j] - P[j-1])/dx uses ADJACENT cells only.
    # This is the staggered grid's key property — no cell-skipping.

    u_star = np.zeros(N + 1)
    u_star[0] = 0.0  # Head-end wall BC

    for j in range(1, N):
        # Face j sits between cells j-1 and j
        rho_f = 0.5 * (rho[j - 1] + rho[j])
        A_f = A_face[j]

        if rho_f < 1e-10 or A_f < 1e-10:
            u_star[j] = 0.0
            continue

        # --- Pressure gradient (the clean staggered form) ---
        dPdx = (P[j] - P[j - 1]) / dx
        pres_force = -dPdx

        # --- Convective momentum flux ---
        # Momentum CV for face j spans from cell-center j-1 to cell-center j.
        # Fluxes cross at cell centers j-1 and j.
        # Velocity at cell center i ≈ 0.5*(u[i] + u[i+1])
        u_at_jm1 = 0.5 * (u[j - 1] + u[j])      # velocity at cell center j-1
        u_at_j = 0.5 * (u[j] + u[j + 1])          # velocity at cell center j

        # Mass flux at cell centers
        mdot_jm1 = rho[j - 1] * u_at_jm1 * A_port[j - 1]
        mdot_j = rho[j] * u_at_j * A_port[j]

        # Upwind velocity for momentum flux
        if mdot_jm1 >= 0:
            u_upwind_left = u[j - 1] if j > 1 else 0.0  # face j-1
        else:
            u_upwind_left = u[j]

        if mdot_j >= 0:
            u_upwind_right = u[j]
        else:
            u_upwind_right = u[j + 1] if j < N - 1 else u[j]

        conv = -(mdot_j * u_upwind_right - mdot_jm1 * u_upwind_left) / dx

        # --- Friction ---
        # Interpolate friction factor and hydraulic diameter to face
        f_f = 0.5 * (f_darcy[j - 1] + f_darcy[j])
        D_f = 0.5 * (D_hyd[j - 1] + D_hyd[j])
        friction = -f_f / (2.0 * D_f) * rho_f * abs(u[j]) * u[j]

        # --- Update ---
        RHS = pres_force + conv / A_f + friction
        u_star[j] = u[j] + dt / rho_f * RHS

    # Nozzle face: extrapolate for now (corrected by pressure step)
    u_star[N] = u_star[N - 1]

    # -------------------------------------------------------
    # STEP 2: PRESSURE CORRECTION 1
    # -------------------------------------------------------
    # Continuity at cell i:
    #   (ρ_new - ρ_old)·A·dx/dt + (ṁ*[i+1] - ṁ*[i]) = S·dx
    #
    # With ρ_new = (P + P')/(R·T) and velocity correction
    #   u'[j] = -dt/ρ_f[j] · (P'[j] - P'[j-1]) / dx
    #
    # This gives the tridiagonal system for P'.

    # Momentum equation inverse diagonal coefficient at each face
    d_face = np.zeros(N + 1)
    for j in range(1, N):
        rho_f = 0.5 * (rho[j - 1] + rho[j])
        d_face[j] = dt / max(rho_f, 1e-6)

    a_sub = np.zeros(N)
    a_diag = np.zeros(N)
    a_sup = np.zeros(N)
    b_rhs = np.zeros(N)

    for i in range(N):
        RT_local = R_specific * T[i]
        a_t = A_port[i] * dx / (RT_local * dt)  # Transient density term

        # West face coefficient (face i)
        if i > 0:
            coeff_w = A_face[i] * A_face[i] * d_face[i] / dx
        else:
            coeff_w = 0.0  # Wall: no flux correction

        # East face coefficient (face i+1)
        if i < N - 1:
            coeff_e = A_face[i + 1] * A_face[i + 1] * d_face[i + 1] / dx
        else:
            coeff_e = 0.0  # Nozzle: handled separately

        a_sub[i] = -coeff_w
        a_sup[i] = -coeff_e
        a_diag[i] = a_t + coeff_w + coeff_e

        # Nozzle BC: P' at last cell drives nozzle flow correction
        if i == N - 1:
            a_diag[i] += nozzle_coeff

        # Continuity residual using u*
        if i > 0:
            rho_w = 0.5 * (rho[i - 1] + rho[i])
            mdot_star_w = rho_w * u_star[i] * A_face[i]
        else:
            mdot_star_w = 0.0

        if i < N - 1:
            rho_e = 0.5 * (rho[i] + rho[i + 1])
            mdot_star_e = rho_e * u_star[i + 1] * A_face[i + 1]
        else:
            # Nozzle outflow
            mdot_star_e = P[i] * nozzle_coeff

        b_rhs[i] = mass_source[i] * dx - (mdot_star_e - mdot_star_w)

    P_prime = thomas_solve(a_sub, a_diag, a_sup, b_rhs, N)

    # Update pressure
    P_new = P + P_prime

    # Correct face velocities: u = u* - d_face · (P'[j] - P'[j-1])/dx
    u_new = np.copy(u_star)
    for j in range(1, N):
        u_new[j] = u_star[j] - d_face[j] * (P_prime[j] - P_prime[j - 1]) / dx
    u_new[0] = 0.0

    # -------------------------------------------------------
    # STEP 3: PRESSURE CORRECTION 2
    # -------------------------------------------------------
    rho_new_1 = np.zeros(N)
    for i in range(N):
        rho_new_1[i] = P_new[i] / (R_specific * T[i])

    # Recompute d_face with updated density
    for j in range(1, N):
        rho_f = 0.5 * (rho_new_1[j - 1] + rho_new_1[j])
        d_face[j] = dt / max(rho_f, 1e-6)

    for i in range(N):
        RT_local = R_specific * T[i]
        a_t = A_port[i] * dx / (RT_local * dt)

        if i > 0:
            coeff_w = A_face[i] * A_face[i] * d_face[i] / dx
        else:
            coeff_w = 0.0

        if i < N - 1:
            coeff_e = A_face[i + 1] * A_face[i + 1] * d_face[i + 1] / dx
        else:
            coeff_e = 0.0

        a_sub[i] = -coeff_w
        a_sup[i] = -coeff_e
        a_diag[i] = a_t + coeff_w + coeff_e
        if i == N - 1:
            a_diag[i] += nozzle_coeff

        if i > 0:
            rho_w = 0.5 * (rho_new_1[i - 1] + rho_new_1[i])
            mdot_w2 = rho_w * u_new[i] * A_face[i]
        else:
            mdot_w2 = 0.0

        if i < N - 1:
            rho_e = 0.5 * (rho_new_1[i] + rho_new_1[i + 1])
            mdot_e2 = rho_e * u_new[i + 1] * A_face[i + 1]
        else:
            mdot_e2 = P_new[i] * nozzle_coeff

        b_rhs[i] = mass_source[i] * dx - (mdot_e2 - mdot_w2)

    P_prime2 = thomas_solve(a_sub, a_diag, a_sup, b_rhs, N)

    P_new = P_new + P_prime2
    for j in range(1, N):
        u_new[j] = u_new[j] - d_face[j] * (P_prime2[j] - P_prime2[j - 1]) / dx
    u_new[0] = 0.0

    # -------------------------------------------------------
    # STEP 3b: ENERGY EQUATION
    # -------------------------------------------------------
    T_new = np.zeros(N)
    for i in range(N):
        rhoA = max(rho_new_1[i] * A_port[i], 1e-10)

        # West face flux (upwind using staggered face velocity)
        if i == 0:
            flux_w = 0.0
        else:
            rho_w = 0.5 * (rho_new_1[i - 1] + rho_new_1[i])
            mdot_w = rho_w * u_new[i] * A_face[i]
            if mdot_w >= 0:
                flux_w = mdot_w * T[i - 1]
            else:
                flux_w = mdot_w * T[i]

        # East face flux
        if i < N - 1:
            rho_e = 0.5 * (rho_new_1[i] + rho_new_1[i + 1])
            mdot_e = rho_e * u_new[i + 1] * A_face[i + 1]
            if mdot_e >= 0:
                flux_e = mdot_e * T[i]
            else:
                flux_e = mdot_e * T[i + 1]
        else:
            flux_e = P_new[i] * nozzle_coeff * T[i]

        conv_T = -(flux_e - flux_w) / dx
        source_T = mass_source[i] * T_flame

        T_new[i] = T[i] + dt / rhoA * (conv_T + source_T) * A_port[i]
        T_new[i] = max(T_new[i], 300.0)
        T_new[i] = min(T_new[i], T_flame * 1.01)

    # -------------------------------------------------------
    # STEP 4: UPDATE DENSITY
    # -------------------------------------------------------
    for i in range(N):
        P_new[i] = max(P_new[i], 1e3)
    rho_new = np.zeros(N)
    for i in range(N):
        rho_new[i] = P_new[i] / (R_specific * T_new[i])

    return rho_new, u_new, P_new, T_new


# ================================================================
# CFL Time Step
# ================================================================

@njit(cache=True)
def compute_dt_cfl(u, a_sound, dx, N, cfl_target, dt_max):
    """
    Compute the adaptive time step from the CFL condition.

    For compressible flow, information propagates at the speed of sound
    relative to the flow. The CFL condition requires that the numerical
    domain of dependence contains the physical domain of dependence:

        CFL = (|u| + a) · dt / dx ≤ CFL_target

    Rearranging: dt = CFL_target · dx / max(|u| + a)

    Parameters
    ----------
    u : ndarray (N,)
        Velocity field (face or cell-center, depending on caller).
    a_sound : float
        Speed of sound [m/s] (max across all cells).
    dx : float
        Cell width [m].
    N : int
        Length of the velocity array.
    cfl_target : float
        Target CFL number. 0.5 is stable for explicit schemes;
        0.3 is conservative.
    dt_max : float
        Maximum allowed time step [s]. Prevents excessively large
        steps during low-velocity phases (e.g., before ignition).

    Returns
    -------
    dt : float
        Time step [s].
    """
    max_wave_speed = a_sound  # At minimum, the sound speed
    for i in range(N):
        ws = abs(u[i]) + a_sound
        if ws > max_wave_speed:
            max_wave_speed = ws
    dt = cfl_target * dx / max_wave_speed
    return min(dt, dt_max)
