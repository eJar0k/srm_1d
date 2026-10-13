"""
solid_thermal.py -- solid-phase heating sub-solvers for the ignition gate
=========================================================================

Advances the surface temperature of an unignited propellant surface being
heated by the hot bore gas. Three models, selected by the
``solid_heating_model`` run kwarg (SM1, docs/pyrogen_realism/SM1_REPORT.md):

``goodman`` (default)
    Goodman cubic-polynomial heat-balance integral. State: penetration depth
    ``delta``; the surface temperature is ALGEBRAIC in the current h and delta
    (``_compute_T_surf``). ``_step_goodman_ode`` = RK4 on delta.
``peretz``
    Peretz's integral-method surface-temperature ODE (thesis Eq. III-50 /
    AIAA 72-1119 Eq. 19): same cubic profile and heat-balance integral, but T_s
    is the carried state and delta is diagnosed (``_step_peretz_ode``).
    Identical to Goodman for constant h and T_gas.
``conduction``
    Exact 1-D transient conduction in the solid on a graded grid, implicit
    (backward Euler, sub-stepped against the heating age), Robin surface
    boundary (``_step_conduction``). Uses only alpha and k_solid.

Wired into the time loop: ``simulation._goodman_ignition_sources_and_mass``
calls one of them for each unignited cell every step, then
``_surface_has_ignited`` (the ``T_surf > T_ignition`` gate). See
``docs/v0_7_0/DESIGN.md`` Eq. 8 for the Goodman surface-temperature closure.
"""

import math

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper


MIN_DELTA = 1.0e-12
GOODMAN_SUBSTEP_RATIO = 0.1


@njit(cache=True)
def _compute_T_surf(delta, h_c, T_gas, T_initial, k_solid):
    """
    Algebraic Goodman surface temperature relation.

    DESIGN.md Eq. 8:
        T_s = (3*k*T_i + h*delta*T_g) / (3*k + h*delta)
    """
    if h_c <= 0.0 or k_solid <= 0.0:
        return T_initial
    d = delta
    if d < MIN_DELTA:
        d = MIN_DELTA
    denom = 3.0 * k_solid + h_c * d
    if denom <= 0.0:
        return T_initial
    return (3.0 * k_solid * T_initial + h_c * d * T_gas) / denom


@njit(cache=True)
def _goodman_rhs(delta, h_c, alpha, k_solid):
    """Right-hand side of the Goodman penetration-depth ODE."""
    if h_c <= 0.0 or alpha <= 0.0 or k_solid <= 0.0:
        return 0.0
    d = delta
    if d < MIN_DELTA:
        d = MIN_DELTA
    denom = d * (6.0 * k_solid + h_c * d)
    if denom <= 0.0:
        return 0.0
    return 12.0 * alpha * (3.0 * k_solid + h_c * d) / denom


@njit(cache=True)
def _rk4_delta_step(delta, h_c, alpha, k_solid, dt):
    if dt <= 0.0:
        if delta < MIN_DELTA:
            return MIN_DELTA
        return delta

    d = delta
    if d < MIN_DELTA:
        d = MIN_DELTA

    k1 = _goodman_rhs(d, h_c, alpha, k_solid)
    k2 = _goodman_rhs(d + 0.5 * dt * k1, h_c, alpha, k_solid)
    k3 = _goodman_rhs(d + 0.5 * dt * k2, h_c, alpha, k_solid)
    k4 = _goodman_rhs(d + dt * k3, h_c, alpha, k_solid)

    d_new = d + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if d_new < MIN_DELTA:
        return MIN_DELTA
    return d_new


@njit(cache=True)
def _step_goodman_ode(delta, T_surf, h_c, T_gas, T_initial,
                      alpha, k_solid, dt):
    """
    RK4 step for the Goodman penetration-depth ODE.

    Returns ``(new_delta, new_T_surf)``. ``T_surf`` is accepted to match
    the Phase 2 task signature; the Goodman formulation recomputes
    surface temperature algebraically from ``new_delta``.
    """
    d = delta
    if d < MIN_DELTA:
        d = MIN_DELTA

    if h_c <= 0.0 or alpha <= 0.0 or k_solid <= 0.0 or dt <= 0.0:
        return d, _compute_T_surf(d, h_c, T_gas, T_initial, k_solid)

    rhs0 = _goodman_rhs(d, h_c, alpha, k_solid)
    ratio = 0.0
    if d > 0.0:
        ratio = abs(rhs0 * dt / d)

    n_sub = 1
    if ratio > GOODMAN_SUBSTEP_RATIO:
        n_sub = int(math.ceil(ratio / GOODMAN_SUBSTEP_RATIO))

    sub_dt = dt / n_sub
    for _ in range(n_sub):
        d = _rk4_delta_step(d, h_c, alpha, k_solid, sub_dt)

    return d, _compute_T_surf(d, h_c, T_gas, T_initial, k_solid)


@njit(cache=True)
def _surface_has_ignited(T_surf, T_ignition):
    """Strict ignition threshold helper for Phase 3 integration."""
    return T_surf > T_ignition


# ---------------------------------------------------------------------------
# SM1: Peretz integral-method surface-temperature ODE
# ---------------------------------------------------------------------------
# Peretz (1974 thesis Eq. III-50; AIAA 72-1119 Eq. 19), with u = T_s - T_i and
# a = T_g - T_i:
#     dT_s/dt = 4 alpha h^2 (T_g - T_s)^3 / [3 k^2 (T_s - T_i)(2 T_g - T_s - T_i)]
#             = 4 alpha h^2 (a - u)^3 / [3 k^2 u (2a - u)]
# derived with dh/dt and dT_g/dt neglected against dT_s/dt (his words: "completely
# justified for the entire induction interval"). The same cubic profile gives
# delta = 3 k (T_s - T_i) / (h (T_g - T_s)).  Initial condition T_s(0) = T_i + eps
# (eps = 0.1 K, Eq. III-51), to avoid the u = 0 singularity.
PERETZ_EPS = 0.1
PERETZ_SUBSTEP_FRACTION = 0.1


@njit(cache=True)
def _peretz_rhs(u, h_c, a, alpha, k_solid):
    """dT_s/dt of Peretz Eq. III-50 (zero when the gas is not hotter than the surface)."""
    if h_c <= 0.0 or u <= 0.0 or a <= u:
        return 0.0
    return (4.0 * alpha * h_c * h_c * (a - u) ** 3
            / (3.0 * k_solid * k_solid * u * (2.0 * a - u)))


@njit(cache=True)
def _step_peretz_ode(T_surf, h_c, T_gas, T_initial, alpha, k_solid, dt):
    """Advance T_s over dt with RK4 sub-stepping (each sub-step changes u or
    (a - u) by at most 10 %). T_s is the state: when h falls it is HELD (the ODE
    rate scales with h^2), it is not recomputed from h. If the gas is not hotter
    than the surface the ODE is outside its derivation (it was written for
    heating) and T_s is held. Returns the new T_s."""
    if h_c <= 0.0 or alpha <= 0.0 or k_solid <= 0.0 or dt <= 0.0:
        return T_surf
    a = T_gas - T_initial
    u = T_surf - T_initial
    if u < PERETZ_EPS:
        u = PERETZ_EPS
    if a <= u:
        return T_surf
    t_rem = dt
    for _ in range(100000):
        f0 = _peretz_rhs(u, h_c, a, alpha, k_solid)
        if f0 <= 0.0:
            break
        scale = u
        if (a - u) < scale:
            scale = a - u
        sub = PERETZ_SUBSTEP_FRACTION * scale / f0
        if sub > t_rem:
            sub = t_rem
        k1 = f0
        k2 = _peretz_rhs(u + 0.5 * sub * k1, h_c, a, alpha, k_solid)
        k3 = _peretz_rhs(u + 0.5 * sub * k2, h_c, a, alpha, k_solid)
        k4 = _peretz_rhs(u + sub * k3, h_c, a, alpha, k_solid)
        u_new = u + (sub / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if u_new > a:
            u_new = a
        u = u_new
        t_rem -= sub
        if t_rem <= 0.0:
            break
    return T_initial + u


# ---------------------------------------------------------------------------
# SM1: exact 1-D transient conduction (graded grid, implicit)
# ---------------------------------------------------------------------------
CONDUCTION_CUT_FACTOR = 7.0     # active depth = 7 sqrt(alpha (age + dt))  (erfc(3.5) ~ 7e-7)
# Resolution fixed by the SM1 convergence study (docs/pyrogen_realism/SM1_REPORT.md
# sec. 2.1), not tuned to sim results: T_s error <= 1 K (T_s <= 900 K) and
# time-to-756 K error <= 0.5 % against the exact Robin solution for
# h <= 20000 W/(m2 K) at sim-like and 10x coarser time steps.
CONDUCTION_AGE_FRACTION = 0.0125  # BE sub-step <= fraction * heating age
CONDUCTION_DY0 = 5.0e-8           # first node spacing [m]
CONDUCTION_RATIO = 1.1            # geometric growth of the node spacing


def build_solid_grid(alpha, t_window, dy0, ratio, depth_factor=8.0):
    """Node depths y[0..M-1] (m), y[0] = 0 = the surface, geometrically graded:
    y[j+1] - y[j] = dy0 * ratio**j, out to depth_factor * sqrt(alpha * t_window)
    (the insulated far boundary). Returns a float64 array."""
    import numpy as np
    L = depth_factor * math.sqrt(max(alpha, 1.0e-30) * max(t_window, 1.0e-3))
    ys = [0.0]
    d = dy0
    while ys[-1] < L:
        ys.append(ys[-1] + d)
        d *= ratio
    return np.asarray(ys, dtype=np.float64)


@njit(cache=True)
def _conduction_substep(T, y, m, h_c, T_drv, T_initial, alpha, k_solid, dt, cp, dp):
    """One backward-Euler step on nodes 0..m-1 (T absolute [K], in place).
    Node m (if m < M) is held at T_initial (Dirichlet; the active depth is
    several penetration depths so it is within ~1e-6 of that). m == M uses an
    insulated far end. Surface node: half control volume, Robin flux
    h (T_drv - T_s) implicit. Thomas algorithm, scratch cp/dp."""
    M = y.shape[0]
    bi = alpha * h_c / k_solid          # (alpha/k) h : flux / (rho c)
    for j in range(m):
        if j == 0:
            V = 0.5 * (y[1] - y[0])
            wl = 0.0
        elif j == M - 1:
            V = 0.5 * (y[j] - y[j - 1])
            wl = alpha / (y[j] - y[j - 1])
        else:
            V = 0.5 * (y[j + 1] - y[j - 1])
            wl = alpha / (y[j] - y[j - 1])
        if j == M - 1:
            wr = 0.0
        else:
            wr = alpha / (y[j + 1] - y[j])
        diag = V / dt + wl + wr
        rhs = V / dt * (T[j] - T_initial)           # work in theta = T - T_i
        if j == 0:
            diag += bi
            rhs += bi * (T_drv - T_initial)
        # forward elimination (lower = -wl, upper = -wr)
        if j == 0:
            cp[0] = -wr / diag
            dp[0] = rhs / diag
        else:
            den = diag + wl * cp[j - 1]
            cp[j] = -wr / den
            dp[j] = (rhs + wl * dp[j - 1]) / den
    th = dp[m - 1]
    T[m - 1] = T_initial + th
    for j in range(m - 2, -1, -1):
        th = dp[j] - cp[j] * th
        T[j] = T_initial + th


@njit(cache=True)
def _step_conduction(T, y, nact, age, h_c, T_drv, T_initial, alpha, k_solid,
                     dt, age_fraction, cp, dp):
    """Advance one cell's solid temperature profile T[0..M-1] over dt.

    ``age`` = time since heating began (the penetration-depth clock; the caller
    adds dt afterwards). The step is split into backward-Euler sub-steps no
    longer than age_fraction * max(age, dt), so the first-order
    time error stays a fixed small fraction of the heating age. The active
    node count grows with the penetration depth. Returns (T_surface, nact).
    """
    M = y.shape[0]
    if alpha <= 0.0 or k_solid <= 0.0 or dt <= 0.0:
        return T[0], nact
    ref = age
    if dt > ref:
        ref = dt
    n_sub = int(math.ceil(dt / (age_fraction * ref)))
    if n_sub < 1:
        n_sub = 1
    sub = dt / n_sub
    t_local = age
    m = nact
    for _ in range(n_sub):
        zcut = CONDUCTION_CUT_FACTOR * math.sqrt(alpha * (t_local + sub))
        while m < M and (m < 3 or y[m - 1] < zcut):
            m += 1
        _conduction_substep(T, y, m, h_c, T_drv, T_initial, alpha, k_solid,
                            sub, cp, dp)
        t_local += sub
    return T[0], m
