"""
solid_thermal.py -- Goodman solid-phase conduction sub-solver
=============================================================

Standalone Phase 2 implementation for v0.7.0. This module advances the
Goodman cubic-polynomial heat-balance integral state for an unignited
propellant surface. It is intentionally not wired into simulation.py yet.
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
