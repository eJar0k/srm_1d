"""IG3 Stage 1: port-area-to-throat ratios of the three Peretz Ap readings (local / head / aft).
Usage: ig3_ap.py <motor.ric> ...   Run from the worktree root with PYTHONPATH=$PWD.
Runs a short t_max with per-ms snapshots and reports, for the grain cells (is_grain):
head = fore-most grain cell, aft = aft-most grain cell, local = the cells that light first (the 8 first-lit),
at t=0 and at the first-ignition snapshot. Ap = pi/4 D_port^2 from the snapshots; At = nozzle throat area (initial)."""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric

for ric in sys.argv[1:]:
    res, perf, nozzle, geo, prop = run_from_ric(str(ROOT / ric), roughness=32e-6, kappa=0.44, T_ignition=756.0,
                                                print_interval=0.0, verbose=False, t_max=0.02,
                                                snapshot_interval=5e-4, P_cutoff=0.01e6)
    At = np.pi / 4 * nozzle.D_throat ** 2
    sn = res['snapshots']
    isg = sn[0]['is_grain']
    idx = np.where(isg)[0]
    ign = np.asarray(res['ignition_time_by_cell'])
    fin = np.where(isg & (ign < 1e8))[0]
    order = fin[np.argsort(ign[fin])][:8]
    tfirst = float(ign[order[0]]) if order.size else None
    def ratios(s):
        Ap = np.pi / 4 * np.asarray(s['D_port']) ** 2
        return dict(head=float(Ap[idx[0]] / At), aft=float(Ap[idx[-1]] / At),
                    local_firstlit=[round(float(Ap[i] / At), 3) for i in order],
                    local_firstlit_median=float(np.median(Ap[order] / At)) if order.size else None)
    t = np.array([s['t'] for s in sn])
    k1 = int(np.argmin(np.abs(t - tfirst))) if tfirst is not None else 0
    print(json.dumps(dict(ric=ric, At_mm2=At * 1e6, n_grain=int(idx.size), head_cell=int(idx[0]), aft_cell=int(idx[-1]),
                          first_lit_cells=[int(i) for i in order], t_first_ign_s=tfirst, t0=ratios(sn[0]),
                          at_first_ign=ratios(sn[k1]), snap_t=float(t[k1])), indent=1))
