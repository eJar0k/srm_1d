"""G1 stability probe: one short run of one committed motor at its default settings -> one JSON row.

    python docs/solver_verification/probes/g1_stability.py <ric file name in motors/> <out.jsonl>
        [--cells N] [--aft MODE] [--head-gap-mm X] [--t-max T] [--mass g]

Records what a numerical front instability would show: termination status (code), max Mach over the run
and where, number of steps and dt floor hits, pressure floor / limit activations, NaN/inf in the pressure trace, the
peak pressure, first light, cells never lit, the (auto) charge and attempts. Run from the checkout root; it imports
its own checkout (prints the path). Always run through srm_1d.tools.guarded_run.
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
    ap.add_argument('--aft', default=None)
    ap.add_argument('--head-gap-mm', type=float, default=None)
    ap.add_argument('--t-max', type=float, default=0.5)
    ap.add_argument('--mass', type=float, default=None, help='igniter charge [g]; default = the motor file / auto')
    ap.add_argument('--jet', default=None)
    ap.add_argument('--no-lump-gaps', action='store_true', help='G2: legacy gap cells + dx clamp')
    a = ap.parse_args()
    kw = dict(t_max=a.t_max, verbose=False, print_interval=0.0, target_propellant_cells=a.cells,
              snapshot_interval=a.t_max / 50.0, )
    if a.jet is not None:
        kw['head_end_jet_h'] = (False if a.jet == 'off' else a.jet)
    if a.no_lump_gaps:
        kw['lump_gaps'] = False
    if a.aft is not None:
        kw['aft_end'] = a.aft
    if a.head_gap_mm is not None:
        kw['head_gap'] = a.head_gap_mm * 1e-3
    if a.mass is not None:
        kw['pyrogen_mass'] = a.mass * 1e-3
        kw['pyrogen'] = 'bpnv'
    row = dict(ric=a.ric, cells=a.cells, aft=a.aft, head_gap_mm=a.head_gap_mm, t_max=a.t_max,
               srm_1d=str(Path(srm_1d.__file__).parent))
    t0 = time.time()
    try:
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter('always')
            res, perf, nozzle, geo, prop = run_from_ric(str(ROOT / 'motors' / a.ric), **kw)
        s = res['summary']
        t = np.asarray(res['time'])
        P = np.asarray(res['P_head'])
        mm = np.asarray(res['max_mach'])
        ign = np.asarray(res['ignition_time_by_cell'], float)
        lit = ign < 1e9
        row.update(
            status='ok', termination=s['termination'], steps=int(s['steps']), t_end=float(t[-1]),
            finite=bool(np.all(np.isfinite(P))), P_peak_MPa=float(P.max() / 1e6),
            mach_max=float(mm.max()), t_mach_max=float(t[int(mm.argmax())]),
            floor_activations=int(s.get('pressure_floor_activations', 0)),
            charge_g=float(s.get('pyrogen_mass_final_kg', float('nan'))) * 1e3,
            attempts=int(s.get('ignition_attempts', 1)),
            first_light_ms=(float(ign[lit].min() * 1e3) if lit.any() else None),
            dark=int(s['grain_cells_never_lit']), N=int(geo.N_cells), n_cone=int(geo.n_cone_cells) if hasattr(geo, 'n_cone_cells') else 0,
            head_gap_mm_used=float(getattr(geo, 'head_gap', 0.0)) * 1e3,
            mass_balance_error=float(s.get('mass_balance_error', float('nan'))),
            dt_min=float(np.min(res['dt'])) if 'dt' in res else None)
    except Exception as exc:
        row.update(status='exception', error=repr(exc)[:300])
    row['wall_s'] = time.time() - t0
    with open(a.out, 'a') as f:
        f.write(json.dumps(row, default=float) + '\n')
    print(json.dumps(row, default=float))


if __name__ == '__main__':
    main()
