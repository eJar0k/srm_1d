"""
solid_replay.py -- SM2 open-loop replay of the solid-heating models.

``run_simulation(record_forcing_steps=n)`` records, per step and per grain cell, the
surface forcing the ignition gate saw (``h_total`` and ``T_driver``; the flux is
``h_total (T_driver - T_s)``). This module pushes that SAME forcing through each solid
model offline, cell by cell, with no flow feedback: the per-cell ignition times then
differ only because of the solid model (and the ignition criterion, which is held fixed).

Conventions match the time loop: a cell lights at the step where the updated surface
temperature exceeds ``T_ignition`` (strictly), and its ignition time is the time at the
START of that step (``rec_t``). Forcing entries ``< 0`` (non-grain cells) are skipped.
"""

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*a, **k):
        if len(a) == 1 and callable(a[0]):
            return a[0]
        return lambda f: f

from .. import solid_thermal as st

NEVER = 1.0e10  # the sim's never-lit sentinel


@njit(cache=True)
def _replay_cell(model, h, Td, tt, dtt, T_i, T_ign, alpha, k, y, age_fraction):
    """Replay one cell. Returns (t_ignite, T_s_at_ign, q_pre, q_post, h_ign, Td_ign, step)."""
    n = h.shape[0]
    M = y.shape[0]
    T = np.full(M, T_i)
    cp = np.zeros(M)
    dp = np.zeros(M)
    nact = 3
    age = 0.0
    delta = 1.0e-6
    Ts = T_i
    t_c100 = np.nan
    Ts_max = T_i
    for j in range(n):
        hj = h[j]
        if hj < 0.0:
            continue
        Ts_before = Ts
        if model == 0:
            delta, Ts = st._step_goodman_ode(delta, Ts, hj, Td[j], T_i, alpha, k, dtt[j])
        elif model == 1:
            Ts = st._step_peretz_ode(Ts, hj, Td[j], T_i, alpha, k, dtt[j])
        else:
            if hj > 0.0 or age > 0.0:
                Ts, nact = st._step_conduction(T, y, nact, age, hj, Td[j], T_i, alpha, k,
                                               dtt[j], age_fraction, cp, dp)
                age += dtt[j]
        if Ts > Ts_max:
            Ts_max = Ts
        if np.isnan(t_c100) and Ts > T_ign - 100.0:
            t_c100 = tt[j]
        if Ts > T_ign:
            return (tt[j], Ts, hj * (Td[j] - Ts_before), hj * (Td[j] - Ts), hj, Td[j], j, t_c100, Ts_max)
    return (NEVER, Ts, np.nan, np.nan, np.nan, np.nan, -1, t_c100, Ts_max)


@njit(cache=True)
def _replay_all(model, h, Td, tt, dtt, T_i, T_ign, alpha, k, y, age_fraction):
    n, N = h.shape
    out = np.full((N, 9), np.nan)
    out[:, 0] = NEVER
    for i in range(N):
        col_h = np.ascontiguousarray(h[:, i])
        if col_h.max() < 0.0:
            continue
        col_T = np.ascontiguousarray(Td[:, i])
        r = _replay_cell(model, col_h, col_T, tt, dtt, T_i, T_ign, alpha, k, y, age_fraction)
        out[i, 0] = r[0]
        out[i, 1] = r[1]
        out[i, 2] = r[2]
        out[i, 3] = r[3]
        out[i, 4] = r[4]
        out[i, 5] = r[5]
        out[i, 6] = r[6]
        out[i, 7] = r[7]
        out[i, 8] = r[8]
    return out


MODEL_CODE = {'goodman': 0, 'peretz': 1, 'conduction': 2}


def replay(forcing, model, *, T_i=293.0, T_ign=756.0, k_solid=0.271, rho=None, cps=1500.0,
           alpha=None, grid=None):
    """Replay ``forcing`` (the dict from ``result['forcing']``) through ``model``.

    Returns an (N, 9) array: [t_ignite (NEVER if not lit within the recorded window), T_s,
    q_pre, q_post, h, T_driver at ignition, step index, time T_s first exceeded T_ign - 100 K
    (NaN if never), max T_s reached]. ``alpha`` defaults to
    ``k_solid / (rho cps)`` (give ``rho`` or ``alpha``).
    """
    if alpha is None:
        alpha = k_solid / (rho * cps)
    h = np.ascontiguousarray(forcing['h'])
    Td = np.ascontiguousarray(forcing['T_drv'])
    tt = np.ascontiguousarray(forcing['t'])
    dtt = np.ascontiguousarray(forcing['dt'])
    if grid is None:
        t_win = float(tt[-1] + dtt[-1])
        grid = st.build_solid_grid(alpha, t_win, st.CONDUCTION_DY0, st.CONDUCTION_RATIO)
    return _replay_all(MODEL_CODE[model], h, Td, tt, dtt, T_i, T_ign, alpha, k_solid, grid,
                       st.CONDUCTION_AGE_FRACTION)


def lighting_metrics(t_ign, grain_mask=None):
    """first / t10 / t50 / t90 / last light, number never lit, 90 %-lighting duration, for an
    array of per-cell ignition times (NEVER sentinel for unlit)."""
    t = np.asarray(t_ign, dtype=float)
    if grain_mask is not None:
        t = t[grain_mask]
    lit = t[t < 1.0e9]
    n = t.size
    out = dict(n_grain=int(n), n_lit=int(lit.size), n_never=int(n - lit.size))
    if lit.size == 0:
        return out
    out['first'] = float(lit.min())
    out['last'] = float(lit.max())
    for q in (10, 50, 90):
        # time by which q % of ALL grain cells are lit (NaN if fewer than q % ever light)
        need = int(np.ceil(q / 100.0 * n))
        out[f't{q}'] = float(np.sort(lit)[need - 1]) if lit.size >= need else float('nan')
    out['dur90'] = out['t90'] - out['first'] if np.isfinite(out['t90']) else float('nan')
    return out
