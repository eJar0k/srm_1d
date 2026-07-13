"""Cross-motor CFL stability probe: does cfl 0.5 stay safe on the
violent-ignition motors (Chunc/Zerox/BALLSstick) the way it does on
Hasegawa A? Decides whether B1 is a simple default bump (0.3->0.5) or
must be phase-aware (conservative fill, aggressive plateau).

Metric: P_peak(0.5) / P_peak(0.3). ~1.0 + both stable => safe bump.
"""
import sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
REPO = Path(r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d")
sys.path.insert(0, str(REPO))
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric
MOTORS = REPO / 'motors'
BASE = dict(pyrogen='bpnv', P_cutoff=0.05e6, t_max=1.0,
            snapshot_interval=10.0, print_interval=10.0, verbose=False)
CASES = {'chunc': 'machbusterNew.ric', 'zerox': 'zerox.ric',
         'BALLSstick': 'BALLSstick.ric'}


def run(ric, cfl):
    t0 = time.perf_counter()
    try:
        result, *_ = run_from_ric(str(MOTORS / ric), cfl_target=cfl, **BASE)
    except Exception as e:
        return None, None, None, f"EXC:{type(e).__name__}:{str(e)[:50]}"
    wall = time.perf_counter() - t0
    P = np.asarray(result['max_pressure'])
    n = len(result['dt'])
    ppk = P.max() / 1e6
    sim_t = result['time'][-1]
    ok = np.isfinite(ppk) and sim_t > 0.3 and ppk < 200.0
    return ppk, n, wall, ("ok" if ok else f"SUSPECT(t={sim_t:.2f})")


# warmup
run('machbusterNew.ric', 0.3)
print(f"{'motor':>11} {'cfl':>4} {'P_peak':>9} {'n_steps':>9} {'wall':>6}  status")
for name, ric in CASES.items():
    res = {}
    for cfl in (0.3, 0.5):
        ppk, n, wall, status = run(ric, cfl)
        res[cfl] = ppk
        ps = f"{ppk:.4f}" if ppk is not None else "--"
        ns = f"{n}" if n is not None else "--"
        ws = f"{wall:.1f}" if wall is not None else "--"
        print(f"{name:>11} {cfl:>4.1f} {ps:>9} {ns:>9} {ws:>6}  {status}")
    if res.get(0.3) and res.get(0.5) and res[0.3] > 0:
        print(f"{name:>11}  -> P_peak ratio 0.5/0.3 = {res[0.5]/res[0.3]:.4f}\n")
