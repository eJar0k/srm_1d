"""
solid_verify.py -- exact solutions + prescribed-history drivers for the three
solid-heating models in ``srm_1d.solid_thermal`` (SM1).

The drivers feed each kernel a prescribed surface history ``h(t)``, ``T_g(t)``
at a fixed step ``dt`` and return the surface temperature history, so the
models can be checked against

* the exact semi-infinite Robin solution (constant h, T_g) -- Carslaw & Jaeger
  Sec. 2.7 / Eq. (2.7.?): ``T_s - T_i = (T_g - T_i) (1 - erfcx(beta))``,
  ``beta = h sqrt(alpha t) / k``;
* the exact constant-flux solution ``T_s - T_i = 2 q sqrt(alpha t / pi) / k``;
* an independent Volterra (Duhamel flux-superposition) solve for time-varying h,
  which does not use the conduction grid at all.
"""

import math

import numpy as np
from scipy.special import erfcx

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*a, **k):
        if len(a) == 1 and callable(a[0]):
            return a[0]
        return lambda f: f

from .. import solid_thermal as st


# ------------------------------------------------------------------ exact
def exact_robin(t, h, T_g, T_i, alpha, k):
    """Surface temperature, semi-infinite solid, Robin BC at constant h, T_g."""
    t = np.asarray(t, dtype=float)
    beta = h * np.sqrt(alpha * t) / k
    return T_i + (T_g - T_i) * (1.0 - erfcx(beta))


def exact_flux(t, q, T_i, alpha, k):
    """Surface temperature, semi-infinite solid, constant surface flux q [W/m^2]."""
    t = np.asarray(t, dtype=float)
    return T_i + 2.0 * q * np.sqrt(alpha * t / math.pi) / k


def time_to(T_hist, dt, T_target, T0):
    """First time T_hist (sampled at t = (n+1) dt, T(0) = T0) crosses T_target
    (linear interpolation within the step); NaN if never."""
    T = np.concatenate(([T0], np.asarray(T_hist, dtype=float)))
    above = np.nonzero(T >= T_target)[0]
    if above.size == 0:
        return float('nan')
    n = above[0]
    f = (T_target - T[n - 1]) / (T[n] - T[n - 1])
    return dt * (n - 1 + f)


# ---------------------------------------------------------------- drivers
@njit(cache=True)
def run_goodman(h_arr, Tg_arr, dt, T_i, alpha, k):
    n = h_arr.shape[0]
    out = np.empty(n)
    delta = 1.0e-6          # the sim's initial penetration depth
    Ts = T_i
    for i in range(n):
        delta, Ts = st._step_goodman_ode(delta, Ts, h_arr[i], Tg_arr[i], T_i,
                                         alpha, k, dt)
        out[i] = Ts
    return out


@njit(cache=True)
def run_peretz(h_arr, Tg_arr, dt, T_i, alpha, k):
    n = h_arr.shape[0]
    out = np.empty(n)
    Ts = T_i
    for i in range(n):
        Ts = st._step_peretz_ode(Ts, h_arr[i], Tg_arr[i], T_i, alpha, k, dt)
        out[i] = Ts
    return out


@njit(cache=True)
def run_conduction(h_arr, Tg_arr, dt, T_i, alpha, k, y, age_fraction):
    n = h_arr.shape[0]
    M = y.shape[0]
    out = np.empty(n)
    T = np.full(M, T_i)
    cp = np.zeros(M)
    dp = np.zeros(M)
    nact = 3
    age = 0.0
    started = False
    for i in range(n):
        if not started and h_arr[i] <= 0.0:
            out[i] = T_i
            continue
        started = True
        Ts, nact = st._step_conduction(T, y, nact, age, h_arr[i], Tg_arr[i], T_i,
                                       alpha, k, dt, age_fraction, cp, dp)
        age += dt
        out[i] = Ts
    return out


@njit(cache=True)
def volterra_flux_solve(h_arr, Tg_arr, dt, T_i, alpha, k):
    """Independent reference for time-varying h: semi-infinite solid, surface
    flux q(t) = h(t) (T_g(t) - T_s(t)), with T_s from Duhamel superposition of
    flux steps,
        T_s(t_n) - T_i = sqrt(alpha)/(k sqrt(pi)) * sum_j q_j * 2 (sqrt(t_n - t_{j-1})
                                                                 - sqrt(t_n - t_j)),
    q piecewise constant on each step, solved implicitly for the newest q_n.
    O(n^2); first order in dt (use Richardson)."""
    n = h_arr.shape[0]
    out = np.empty(n)
    q = np.empty(n)
    c = math.sqrt(alpha) / (k * math.sqrt(math.pi))
    for m in range(n):
        tm = (m + 1) * dt
        acc = 0.0
        for j in range(m):
            acc += q[j] * 2.0 * (math.sqrt(tm - j * dt) - math.sqrt(tm - (j + 1) * dt))
        w = 2.0 * math.sqrt(dt)
        # theta = c (acc + q_m w),  q_m = h (a - theta)  ->  solve linear
        a = Tg_arr[m] - T_i
        h = h_arr[m]
        theta = c * (acc + h * a * w) / (1.0 + c * h * w)
        q[m] = h * (a - theta)
        out[m] = T_i + theta
    return out


def volterra_richardson(h_fun, Tg_fun, t_end, dt, T_i, alpha, k):
    """Richardson-extrapolated Volterra reference sampled at multiples of dt.
    ``h_fun(t)``/``Tg_fun(t)`` are vectorised callables of time (step midpoint
    is NOT used: values are taken at the step END, like the kernels)."""
    n = int(round(t_end / dt))
    t1 = (np.arange(n) + 1) * dt
    r1 = volterra_flux_solve(h_fun(t1), Tg_fun(t1), dt, T_i, alpha, k)
    n2 = 2 * n
    d2 = dt / 2
    t2 = (np.arange(n2) + 1) * d2
    r2 = volterra_flux_solve(h_fun(t2), Tg_fun(t2), d2, T_i, alpha, k)
    return t1, 2.0 * r2[1::2] - r1


# ----------------------------------------------------------------- grid
def default_grid(alpha, t_window, dy0=None, ratio=None):
    """The grid the sim uses (constants fixed by the SM1 convergence study)."""
    return st.build_solid_grid(alpha, t_window,
                               st.CONDUCTION_DY0 if dy0 is None else dy0,
                               st.CONDUCTION_RATIO if ratio is None else ratio)
