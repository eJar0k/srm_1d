"""RT1: flow regime / recovery factor / driver-temperature rise at each cell's ignition.
   python rt1_regime.py  (Chunc, BALLSstick short runs, fixed charge, recovery driver)"""
import sys, warnings
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
warnings.simplefilter('ignore')
from srm_1d.openmotor_adapter import run_from_ric
for ric, mass, tmax in (('machbusterNew.ric', 2.6225688088052457e-3, 0.03), ('BALLSstick.ric', 4.829229907013092e-3, 0.03),
                        ('hasegawa_a.ric', 14.976521608656354e-3, 0.03)):
    for ef, mf in (('sensible', 'legacy'), ('total', 'conservative')):
        r = run_from_ric(str(ROOT / 'motors' / ric), t_max=tmax, pyrogen='bpnv', pyrogen_mass=mass,
                         roughness=32e-6, kappa=0.44, T_ignition=756.0, print_interval=0.0, verbose=False,
                         ignition_driver_temperature='recovery', energy_formulation=ef,
                         momentum_formulation=mf, history_capacity=150000)[0]
        d = np.asarray(r['ignition_diag']); lit = np.isfinite(d[:, 0])
        dT = d[lit, 10] - d[lit, 4]
        print(f"{ric:18s} {ef:8s} cells lit {lit.sum():3d} turbulent {int(d[lit,13].sum()):3d} laminar {int((1-d[lit,13]).sum()):3d} "
              f"Re {d[lit,12].min():.0f}..{d[lit,12].max():.0f} r {d[lit,11].min():.3f}..{d[lit,11].max():.3f} "
              f"T_aw-T {dT.min():.0f}..{dT.max():.0f} K (median {np.median(dT):.0f})  Mach at ign {np.nanmax(d[lit,6]):.2f}")
