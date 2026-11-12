"""TP1: pathfin54 measured vs simulated thrust, before/after the stagnation-pressure fix.

Reads the output of ``tp1_motor.py pathfin OUT.json`` (+ ``OUT.json.npz``) and the measured
``static_fire_data/pathfin54.csv``.  The measured trace starts mid-ignition (80 N, 90 kPa at
t = 0) and is truncated in the tail (82 N, 1.08 MPa at the last sample), so the sim is time-
aligned to the measurement on the rising edge (time at which the thrust first crosses a
fraction of its own peak; reported for 10 %, 20 % and 30 %), and the integrals are compared
over the SAME span after alignment, as well as in full.

    python docs/solver_verification/probes/tp1_pathfin_compare.py OUT.json [FIG.png]
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def cross(t, y, level):
    i = np.where(y >= level)[0]
    if len(i) == 0:
        return None
    k = int(i[0])
    if k == 0:
        return float(t[0])
    return float(t[k - 1] + (level - y[k - 1]) * (t[k] - t[k - 1]) / (y[k] - y[k - 1]))


def stats(t, F, t0, t1):
    m = (t >= t0) & (t <= t1)
    tt, FF = t[m], F[m]
    imp = float(np.sum(0.5 * (FF[1:] + FF[:-1]) * np.diff(tt)))
    return imp, float(FF.max()), imp / (tt[-1] - tt[0])


def main():
    out = sys.argv[1]
    fig_path = sys.argv[2] if len(sys.argv) > 2 else None
    row = json.load(open(out))
    z = np.load(out + '.npz')
    t, Fb, Fa = z['t'], z['F_before'], z['F_after']
    Ph, P0 = z['P_head'], z['P0']
    d = np.loadtxt(ROOT / 'static_fire_data' / 'pathfin54.csv', delimiter=',', skiprows=1)
    te, Fe, Pe = d[:, 0], d[:, 1], d[:, 2]
    pke = Fe.max()
    res = {}
    for frac in (0.10, 0.20, 0.30):
        # measured rising edge: the first crossing (the trace starts at 80 N = 6.8 % of peak)
        tm = cross(te, Fe, frac * pke)
        ts_b = cross(t, Fb, frac * Fb.max())
        ts_a = cross(t, Fa, frac * Fa.max())
        # align on the measured crossing; sim shifted by (tm - ts)
        dt_b, dt_a = tm - ts_b, tm - ts_a
        span = (te[0], te[-1])
        out_f = {}
        for lbl, F, dt in (('before', Fb, dt_b), ('after', Fa, dt_a)):
            tsh = t + dt
            imp, pk, avg = stats(tsh, F, span[0], span[1])
            out_f[lbl] = dict(shift_s=dt, impulse_same_span=imp, peak=pk, avg_same_span=avg)
        res[f'{int(frac*100)}%'] = out_f
    # measured same span
    imp_e, pk_e, avg_e = stats(te, Fe, te[0], te[-1])
    summary = dict(measured=dict(impulse=imp_e, peak=float(pk_e), avg=avg_e, span=[te[0], te[-1]]),
                   aligned=res)
    # thrust-over-pressure at the plateau (80 % of peak thrust window), head-end pressure
    m80 = Fe > 0.8 * pke
    summary['F_over_P_N_per_MPa'] = dict(
        measured=float(np.median(Fe[m80] / (Pe[m80] / 1e6))),
        sim_before_over_Phead=float(np.median(Fb[Fb > 0.8 * Fb.max()] / (Ph[Fb > 0.8 * Fb.max()] / 1e6))),
        sim_after_over_Phead=float(np.median(Fa[Fa > 0.8 * Fa.max()] / (Ph[Fa > 0.8 * Fa.max()] / 1e6))),
        sim_after_over_P0=float(np.median(Fa[Fa > 0.8 * Fa.max()] / (P0[Fa > 0.8 * Fa.max()] / 1e6))))
    # peak pressure and pressure at the plateau
    summary['P_MPa'] = dict(measured_peak=float(Pe.max() / 1e6), sim_Phead_peak=float(Ph.max() / 1e6),
                            sim_P0_peak=float(P0.max() / 1e6),
                            measured_at_peak_thrust=float(Pe[np.argmax(Fe)] / 1e6),
                            sim_Phead_at_peak_thrust=float(Ph[np.argmax(Fa)] / 1e6))
    print(json.dumps(summary, indent=1))
    if fig_path:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        tm = cross(te, Fe, 0.2 * pke)
        for F, lbl, c in ((Fb, 'sim, P_head (before)', 'tab:orange'), (Fa, 'sim, P0_nozzle (after)', 'tab:blue')):
            ts = cross(t, F, 0.2 * F.max())
            ax[0].plot(t + tm - ts, F, c, lw=1.2, label=lbl)
        ax[0].plot(te, Fe, 'k', lw=1.5, label='measured')
        ax[0].set_ylabel('thrust [N]')
        ax[0].legend()
        ts = cross(t, Fa, 0.2 * Fa.max())
        ax[1].plot(t + tm - ts, Ph / 1e6, 'tab:green', lw=1.2, label='sim P_head')
        ax[1].plot(t + tm - ts, P0 / 1e6, 'tab:red', lw=1.0, label='sim P0_nozzle')
        ax[1].plot(te, Pe / 1e6, 'k', lw=1.5, label='measured pressure (location unknown)')
        ax[1].set_ylabel('pressure [MPa]')
        ax[1].set_xlabel('time, sim aligned on the 20 % thrust rising edge [s]')
        ax[1].set_xlim(0, 3.3)
        ax[1].legend()
        fig.tight_layout()
        fig.savefig(fig_path, dpi=130)


if __name__ == '__main__':
    main()
