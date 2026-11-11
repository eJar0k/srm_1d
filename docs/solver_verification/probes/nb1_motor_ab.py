"""NB1 motor A/B: nozzle_bc_stagnation 'static' vs 'isentropic' on the canonical motors.

Re-uses the three B0 rows (docs/v0_7_4/probes/mp_cap_crossmotor.py: committed .ric files,
canonical knobs, library BPNV igniter, N = 100 default) with the toggle added, and adds the
exit diagnostics at the plateau: grain-exit (cell N-1) Mach and the stagnation ratio P0/P.

  M_e is reconstructed from the run's own mass flow and cell state (NOT assumed from the
  area ratio):  mdot = A_e P M_e sqrt(gamma/(R T)),  P0/P = (1 + (gamma-1)/2 M_e^2)^(gamma/(gamma-1)).

Run from the checkout root through the guard, one (motor, mode) per process:
    python -m srm_1d.tools.guarded_run --mem-gb 8 --timeout-min 60 -- \
        python docs/solver_verification/probes/nb1_motor_ab.py chunc isentropic OUT.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import srm_1d  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "mp_cap_crossmotor", ROOT / "docs" / "v0_7_4" / "probes" / "mp_cap_crossmotor.py")
_probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_probe)

RIC = {'chunc': 'machbusterNew.ric', 'balls': 'BALLSstick.ric', 'haseg': 'hasegawa_a.ric'}
T_MAX = {'chunc': 1.0, 'balls': 2.0, 'haseg': 10.0}
P_CUT = {'chunc': 0.01e6, 'balls': 0.01e6, 'haseg': 0.05e6}
PLATEAU_T = {'chunc': (0.4, 0.9), 'balls': (1.0, 1.1), 'haseg': (2.0, 3.0)}


def exit_state(res, t_lo, t_hi):
    """Median grain-exit Mach and P0/P over [t_lo, t_hi] from the per-step history
    (massflow, P_exit) plus the exit-cell T / D_port / gamma / R of the snapshots."""
    t = np.asarray(res['time'])
    mdot = np.asarray(res['massflow'])
    Pe = np.asarray(res['P_exit'])
    sel = (t >= t_lo) & (t <= t_hi)
    snaps = [s for s in res['snapshots'] if t_lo <= s['t'] <= t_hi]
    if not snaps:
        snaps = [min(res['snapshots'], key=lambda s: abs(s['t'] - 0.5 * (t_lo + t_hi)))]
    # use the snapshot closest to the window centre for T, A_port, gamma, R
    snap = min(snaps, key=lambda s: abs(s['t'] - 0.5 * (t_lo + t_hi)))
    Te = float(snap['T'][-1])
    A_e = 0.25 * np.pi * float(snap['D_port'][-1]) ** 2
    gam = float(res['gamma_mix_final'][-1])
    Rg = float(res['R_mix_final'][-1])
    Me = mdot[sel] / (A_e * Pe[sel] * np.sqrt(gam / (Rg * Te)))
    Me = float(np.median(Me))
    g = (1.0 + 0.5 * (gam - 1.0) * Me ** 2) ** (gam / (gam - 1.0))
    A_t = 0.25 * np.pi * float(np.median(np.asarray(res['D_throat'])[sel])) ** 2
    return dict(M_exit=Me, P0_over_P=float(g), A_exit_over_A_throat=float(A_e / A_t),
                cell_Mach_snapshot=float(snap['Mach'][-1]), snap_t=float(snap['t']),
                T_exit=Te, P_exit_MPa=float(np.median(Pe[sel])) / 1e6,
                gamma=gam, D_port_exit_m=float(snap['D_port'][-1]),
                D_case_exit_m=float(np.asarray(res['cell_D_outer'])[-1]),
                A_throat_m2=float(A_t))


def main():
    motor, mode, out = sys.argv[1], sys.argv[2], sys.argv[3]
    ric = _probe._committed_ric_path(RIC[motor])
    print(f"[provenance] srm_1d from {Path(srm_1d.__file__).parent}")
    res = run_from_ric(str(ric), erosive_model='ma2020', erosive_mp_cap=0.0,
                       t_max=T_MAX[motor], P_cutoff=P_CUT[motor],
                       nozzle_bc_stagnation=mode, snapshot_interval=0.05, **_probe.CANON)[0]
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    summary = res['summary']
    row = dict(motor=motor, mode=mode, n_steps=int(summary['steps']),
               termination=summary.get('termination'),
               mass_balance_error=float(summary.get('mass_balance_error', float('nan'))),
               auto_charge_g=float(summary.get('pyrogen_mass_initial', float('nan'))) * 1e3,
               stamp=summary.get('nozzle_bc_stagnation'))
    if motor == 'chunc':
        k = int(P.argmax())
        row.update(peak=float(P[k]), t_peak=float(t[k]),
                   plateau=float(np.median(P[(t >= 0.4) & (t <= 0.9)])))
        row['ratio'] = row['peak'] / row['plateau']
    elif motor == 'balls':
        ms = t < 0.3
        kp = int(np.argmax(P[ms]))
        Ppk, tpk = float(P[ms][kp]), float(t[ms][kp])
        after = (t > tpk) & (t < 1.2)
        row.update(peak=Ppk, t_peak=tpk, trough=float(np.min(P[after])),
                   mid=float(P[int(np.argmin(np.abs(t - 1.03)))]))
        row['pk_tr'] = Ppk / row['trough']
        row['pk_mid'] = Ppk / row['mid']
    else:
        k = int(P.argmax())
        late = t > 1.5
        kl = int(np.argmax(np.where(late, P, -1.0)))
        early = t < 0.5
        ke = int(np.argmax(np.where(early, P, -1.0)))
        row.update(peak=float(P[k]), t_peak=float(t[k]),
                   late_peak=float(P[kl]), t_late_peak=float(t[kl]),
                   early_peak=float(P[ke]), t_early_peak=float(t[ke]))
    row['plateau_window'] = PLATEAU_T[motor]
    row['exit'] = exit_state(res, *PLATEAU_T[motor])
    # head-end vs exit static and stagnation (step 7 input): median over the window
    sel = (t >= PLATEAU_T[motor][0]) & (t <= PLATEAU_T[motor][1])
    Ph = np.asarray(res['P_head'])[sel]
    Pe = np.asarray(res['P_exit'])[sel]
    row['P_head_over_P_exit_static'] = float(np.median(Ph / Pe))
    row['P_head_over_P0_exit'] = float(np.median(Ph / (Pe * row['exit']['P0_over_P'])))
    # P_head history decimated for later inspection (memory-bounded: <= 2000 points)
    idx = np.linspace(0, len(t) - 1, min(2000, len(t))).astype(int)
    row['trace_t'] = t[idx].tolist()
    row['trace_P_MPa'] = P[idx].tolist()
    with open(out, 'w') as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: v for k, v in row.items() if not k.startswith('trace')}, indent=1))


if __name__ == '__main__':
    main()
