"""Is the transient supersonic fill grid/CFL-converged (physical) or a
numerical artifact? Refine cells and tighten CFL; a physical blowdown Mach
should converge, a PISO over-acceleration artifact should drift."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"c:\Users\ejarocki\Documents\Rocketry\Code Stuff\Erosive Burning Solver\srm_1d")
sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric

BASE = dict(roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
            T_ignition=756.0, P_cutoff=0.01e6, t_max=0.05,
            snapshot_interval=0.5, print_interval=0.0, verbose=False)

def run(cells, cfl, dtmax=0.002):
    kw = dict(BASE); kw["cfl_target"]=cfl; kw["dt_max"]=dtmax
    if cells is not None: kw["target_propellant_cells"]=cells
    r = run_from_ric(str(ROOT/"motors"/"machbusterNew.ric"), **kw)[0]
    t=np.asarray(r["time"]); P=np.asarray(r["P_head"])/1e6
    M=np.asarray(r["max_mach"])
    fill = t < 0.003
    return M[fill].max(), float(P.max()), int(np.argmax(P)>0)

print(f"{'config':>22} {'maxMach(fill)':>13} {'P_peak_MPa':>11}")
print("--- grid refinement (cfl=0.3) ---")
for c in (50, 100, 200):
    mm, pk, _ = run(c, 0.3)
    print(f"{'cells='+str(c):>22} {mm:>13.2f} {pk:>11.2f}")
print("--- CFL tightening (cells=100) ---")
for cfl in (0.3, 0.15, 0.05):
    mm, pk, _ = run(100, cfl)
    print(f"{'cfl='+str(cfl):>22} {mm:>13.2f} {pk:>11.2f}")
print("\nIf maxMach drifts up with cells / down with cfl -> numerical (PISO over-accel).")
print("If converged -> physical blowdown within the model (no aero-choking limit).")
