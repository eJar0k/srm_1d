"""G2: run one motor at default settings and save the head-pressure trace + cost.

    python docs/solver_verification/probes/g2_trace.py <ric name in motors/> <out_prefix> [--cells N] [--t-max T]

Writes <out_prefix>.npz (time, P_head decimated to <= 4000 points, plus first-light / charge) and prints a JSON row
(N cells, dx, steps, wall time, peak, t_peak, burn time, charge, dark cells, termination, mass-balance error).
Run through srm_1d.tools.guarded_run from the checkout root; imports its own checkout.
"""
import argparse
import json
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

import srm_1d  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ric')
    ap.add_argument('out')
    ap.add_argument('--cells', type=int, default=100)
    ap.add_argument('--t-max', type=float, default=10.0)
    a = ap.parse_args()
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res, perf, nozzle, geo, prop = run_from_ric(
            str(ROOT / 'motors' / a.ric), t_max=a.t_max, verbose=False, print_interval=0.0,
            target_propellant_cells=a.cells, snapshot_interval=a.t_max / 20.0)
    wall = time.time() - t0
    s = res['summary']
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    idx = np.unique(np.linspace(0, len(t) - 1, min(len(t), 4000)).astype(int))
    ign = np.asarray(res['ignition_time_by_cell'], float)
    np.savez(a.out + '.npz', t=t[idx], P=P[idx])
    k = int(P.argmax())
    row = dict(ric=a.ric, cells=a.cells, N=int(geo.N_cells), dx_mm=float(geo.dx) * 1e3,
               steps=int(s['steps']), wall_s=round(wall, 1), peak=float(P[k]), t_peak=float(t[k]),
               t_end=float(t[-1]), termination=s['termination'],
               charge_g=float(s.get('pyrogen_mass_final_kg', float('nan'))) * 1e3,
               dark=int(s['grain_cells_never_lit']),
               first_light_ms=(float(ign[ign < 1e9].min() * 1e3) if (ign < 1e9).any() else None),
               mass_balance_error=float(s.get('mass_balance_error', float('nan'))),
               propellant_mass=float(s['propellant_mass']), srm_1d=str(Path(srm_1d.__file__).parent))
    print(json.dumps(row))
    with open(a.out + '.json', 'w') as f:
        f.write(json.dumps(row))


if __name__ == '__main__':
    main()
