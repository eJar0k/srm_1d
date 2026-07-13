"""Test: does correcting the bpnv igniter gas-generation reduce the spike?
(A) impetus match: M 30 -> 41 g/mol (W 6869 -> 5000 psi*in3/g, DeMar)
(B) condensed-phase: only a fraction of charge becomes gas (Ma2019 black powder ~56%)
    proxied here by reducing the injected pyrogen mass (gas mass) for the same charge.
Also sweep the as-fired 6 g vs Sutton 0.9 g.
"""
import numpy as np, copy
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
M=r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/motors/machbusterNew.ric"

def run(pyro, mass, label):
    r,_,_,_,_=run_from_ric(M,roughness=32e-6,kappa=0.44,pyrogen=pyro,pyrogen_mass=mass,
        T_ignition=756.0,P_cutoff=0.01e6,t_max=0.5,snapshot_interval=0.2,print_interval=0.0,verbose=False)
    s=r['summary']
    print(f"  {label:42} P_peak={s['P_peak']/1e6:6.2f} t_peak={s['t_peak']*1e3:6.1f}ms "
          f"P_mid={s['P_mid']/1e6:5.2f} ratio={s['P_peak']/s['P_mid']:.3f} "
          f"P_ig={s.get('pyrogen_peak_P',0)/1e6:5.1f}MPa")

base=load_pyrogen('bpnv')
hiM=copy.deepcopy(base); hiM.M=0.0412                       # impetus -> DeMar 5000
print(f"baseline bpnv M={base.M*1000:.0f} (W=6869); corrected M={hiM.M*1000:.1f} (W=5000 DeMar)\n")

print("--- Sutton-default mass (0.9 g) ---")
run(base, None, "baseline (M=30, 100% gas)")
run(hiM,  None, "impetus-corrected (M=41)")

print("\n--- as-fired charge 6 g (real Chunc) ---")
run(base, 6e-3, "baseline 6g (M=30)")
run(hiM,  6e-3, "impetus-corrected 6g (M=41)")

print("\n--- combined: impetus + ~56% gas proxy (inject 0.56*charge as gas) ---")
# proxy condensed phase by passing reduced *gas* mass for the same nominal charge
run(hiM, 0.56*0.9e-3, "M=41 + 0.56x gas, nominal 0.9g charge")
run(hiM, 0.56*6e-3,   "M=41 + 0.56x gas, nominal 6g charge")
print("\nExperimental: ratio ~1.0 (NO spike), plateau ~8.8 MPa; real charge ~6 g")
