"""
G-reconstruction discriminator (thread D) — Chunc / machbusterNew.

Compares the erosion-driving aft-port mass flux G three ways:
  (1) EXPERIMENTAL quasi-steady  : G_exp_qs = P_exp * A_t / (c* * A_port_aft)
  (2) SIM quasi-steady           : G_sim_qs = P_sim * A_t / (c* * A_port_aft)
  (3) SIM actual local rho*|u|   : max over grain cells (what Ma erosive sees)

(1 vs 2) isolates the higher-P effect; (2 vs 3) isolates the transient
gas-dynamic excess (supersonic fill); (1 vs 3) is the total over-production
of the erosion-driving G vs experiment.

High-res experimental trace: C:/Users/ejarocki/Downloads/ThomasMach5_edited.xlsx
"""
import numpy as np
import pandas as pd
from srm_1d.openmotor_adapter import run_from_ric

MOTOR = r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/motors/machbusterNew.ric"
XLSX  = r"C:/Users/ejarocki/Downloads/ThomasMach5_edited.xlsx"

# ---- run canonical Chunc, fine snapshots over the transient + early plateau
res, perf, nozzle, geo, prop = run_from_ric(
    MOTOR, roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None,
    T_ignition=756.0, P_cutoff=0.01e6, t_max=0.6,
    snapshot_interval=0.0005, print_interval=0.0, verbose=False,
)
summ = res['summary']
print("=== SIM baseline ===")
print("summary keys:", sorted(summ.keys()))
print(f"P_peak = {summ['P_peak']/1e6:.3f} MPa @ t={summ.get('t_peak', float('nan'))*1e3:.2f} ms")

# ---- geometry / nozzle / cstar
print("\n=== nozzle attrs ===")
for a in dir(nozzle):
    if not a.startswith('_'):
        v = getattr(nozzle, a)
        if isinstance(v, (int, float)):
            print(f"  nozzle.{a} = {v}")
# throat area
A_t = getattr(nozzle, 'throat_area', None)
if A_t is None:
    Dt = getattr(nozzle, 'throat_diameter', None) or getattr(nozzle, 'throat', None)
    A_t = np.pi/4*Dt**2 if Dt else None
print(f"\nA_t = {A_t} m^2")

# cstar
cstar = summ.get('cstar', None) or summ.get('c_star', None)
print(f"cstar (summary) = {cstar}")
# gas props for fallback
print("prop gas attrs:", [a for a in dir(prop) if not a.startswith('_')][:30])

snaps = res['snapshots']
print(f"\nn_snaps={len(snaps)}, n_cells={len(snaps[0]['x'])}")
sn0 = snaps[0]
print("grain cells:", int(sn0['is_grain'].sum()), " x-range",
      f"{sn0['x'][0]*1e3:.1f}..{sn0['x'][-1]*1e3:.1f} mm")
np.save(r"c:/tmp/_chunc_snaps.npy", np.array(snaps, dtype=object), allow_pickle=True)
np.save(r"c:/tmp/_chunc_meta.npy", np.array(
    {'A_t': A_t, 'cstar': cstar,
     'R': getattr(prop, 'R_specific', None),
     'summary': {k: summ[k] for k in summ if isinstance(summ[k], (int, float))}},
    dtype=object), allow_pickle=True)
print("\nsaved snapshots + meta to c:/tmp/")
