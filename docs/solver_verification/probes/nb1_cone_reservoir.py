"""NB1 step 5: does G1's `bc_reservoir_area_ratio` still have a purpose?

Runs Chunc (committed machbusterNew.ric, canonical knobs, N = 100) with ``aft_end='cone'`` while
overriding the cone-cut area ratio A_BC/A_throat that `nozzle_convergent_cells` takes from
`grain_geometry.bc_reservoir_area_ratio` (monkeypatched in this process ONLY; no source change):

    python docs/solver_verification/probes/nb1_cone_reservoir.py MODE RATIO OUT.json [MOTOR]

    MODE  = static | isentropic          (nozzle_bc_stagnation)
    RATIO = default | <float>            (A_BC / A_throat used for the cut; default = G1's 1 % P0 value)
    MOTOR = chunc (default) | balls | haseg

Output: peak, plateau (median 0.4-0.9 s, chunc), n cells, cone cells, wall.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import srm_1d  # noqa: E402
from srm_1d import grain_geometry as gg  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "mp_cap_crossmotor", ROOT / "docs" / "v0_7_4" / "probes" / "mp_cap_crossmotor.py")
_probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_probe)

RIC = {'chunc': 'machbusterNew.ric', 'balls': 'BALLSstick.ric', 'haseg': 'hasegawa_a.ric'}
T_MAX = {'chunc': 1.0, 'balls': 2.0, 'haseg': 4.0}
WIN = {'chunc': (0.4, 0.9), 'balls': (1.0, 1.1), 'haseg': (2.0, 3.0)}


def main():
    mode, ratio, out = sys.argv[1], sys.argv[2], sys.argv[3]
    motor = sys.argv[4] if len(sys.argv) > 4 else 'chunc'
    orig = gg.bc_reservoir_area_ratio
    used = {}
    if ratio != 'default':
        r = float(ratio)
        gg.bc_reservoir_area_ratio = lambda p0_tol=None, gamma=None: r
    # record what the cut ratio is and how many cone cells result
    orig_ncc = gg.nozzle_convergent_cells

    def spy(*a, **k):
        res = orig_ncc(*a, **k)
        used['n_cone'] = int(res[0])
        used['lump_volume_m3'] = float(res[1])
        used['ratio'] = float(gg.bc_reservoir_area_ratio())
        return res
    gg.nozzle_convergent_cells = spy
    ric = _probe._committed_ric_path(RIC[motor])
    print(f"[provenance] srm_1d from {Path(srm_1d.__file__).parent}")
    res = run_from_ric(str(ric), erosive_model='ma2020', erosive_mp_cap=0.0, t_max=T_MAX[motor],
                       P_cutoff=0.01e6, nozzle_bc_stagnation=mode, aft_end='cone',
                       **_probe.CANON)[0]
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    lo, hi = WIN[motor]
    sel = (t >= lo) & (t <= hi)
    row = dict(motor=motor, mode=mode, ratio_arg=ratio, peak=float(P.max()),
               plateau=float(np.median(P[sel])), n_steps=int(res['summary']['steps']),
               cells=int(res['summary'].get('cells', -1)), **used)
    Path(out).write_text(json.dumps(row, indent=1))
    print(json.dumps(row, indent=1))


if __name__ == '__main__':
    main()
