"""Probe the Courant headroom of the current PISO scheme.

Sweeps cfl_target on canonical Hasegawa A. The current default is 0.3
(tightened from 0.5 in v0.7.3.2). dt = cfl*dx/(|u|+a) includes the
acoustic speed a, so this maps how far the EXISTING scheme can be pushed
before it loses stability or distorts P_peak.

Decision:
  - If P_peak stays ~6.14 MPa and health passes up to cfl >> 0.3,
    the scheme has free Courant headroom -> cheap Lever B (raise/adapt CFL).
  - If P_peak distorts / run collapses near cfl ~0.3-0.5, the scheme is
    genuinely acoustic-CFL-limited -> need implicit-acoustic (IMEX/AUSM+-up).
"""
import sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
REPO = Path(r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d")
sys.path.insert(0, str(REPO))
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric
MOTOR = str(REPO / 'motors' / 'hasegawa_a.ric')
BASE = dict(roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None,
            T_ignition=756.0, P_cutoff=0.05e6, snapshot_interval=10.0,
            print_interval=10.0, verbose=False)

EXP_PEAK = 6.14  # MPa, canonical reference


def run(cfl):
    t0 = time.perf_counter()
    try:
        result, perf, nozzle, geo, prop = run_from_ric(
            MOTOR, target_propellant_cells=100, cfl_target=cfl, **BASE)
    except Exception as e:
        return None, None, None, f"EXC: {type(e).__name__}: {str(e)[:60]}"
    wall = time.perf_counter() - t0
    P = np.asarray(result['max_pressure'])
    n = len(result['dt'])
    ppk = P.max() / 1e6
    sim_t = result['time'][-1]
    # crude health: did it run to a sane sim time and finite P?
    healthy = np.isfinite(ppk) and sim_t > 1.0 and ppk < 100.0
    return ppk, n, wall, ("ok" if healthy else f"SUSPECT (t_end={sim_t:.2f}s)")


# warmup/compile
run(0.3)
print(f"{'cfl':>5} {'P_peak(MPa)':>12} {'dP%':>7} {'n_steps':>10} {'wall(s)':>8} {'steps/s':>8}  status")
base_ppk = None
for cfl in (0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 2.0):
    ppk, n, wall, status = run(cfl)
    if ppk is None:
        print(f"{cfl:>5.2f} {'--':>12} {'--':>7} {'--':>10} {'--':>8} {'--':>8}  {status}")
        continue
    if base_ppk is None:
        base_ppk = ppk
    dpct = 100.0 * (ppk - base_ppk) / base_ppk
    print(f"{cfl:>5.2f} {ppk:>12.4f} {dpct:>+6.2f}% {n:>10d} {wall:>8.2f} {n/wall:>8.0f}  {status}")
