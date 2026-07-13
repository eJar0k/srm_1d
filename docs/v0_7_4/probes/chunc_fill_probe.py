"""Resolve the Chunc ignition fill: how does the surface reach T_ign across
the whole 1m grain in <1.5 ms? Dump T_gas / T_surf / P / is_burning profiles
at 0.25 ms cadence over the first 6 ms."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"c:\Users\ejarocki\Documents\Rocketry\Code Stuff\Erosive Burning Solver\srm_1d")
sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric

MOTOR = str(ROOT / "motors" / "machbusterNew.ric")
res = run_from_ric(MOTOR, roughness=32e-6, kappa=0.44, pyrogen="bpnv",
                   pyrogen_mass=None, T_ignition=756.0, P_cutoff=0.01e6,
                   t_max=0.03, snapshot_interval=0.00025, print_interval=0.0,
                   verbose=False)[0]
snaps = res["snapshots"]
g = snaps[0]["is_grain"]
x = snaps[0]["x"][g]
N = g.sum()
# sample 5 axial stations: head, 25%, 50%, 75%, aft
idx = [0, N//4, N//2, 3*N//4, N-1]
print(f"grain cells={N}, x stations (mm): {[round(x[i]*1e3) for i in idx]}")
print(f"\n{'t_ms':>5} {'P0_MPa':>7} | T_gas[head,25,50,75,aft] (K) | T_surf[...] (K) | nburn")
for s in snaps:
    if s["t"] > 0.006:
        break
    Tg = s["T"][g]; Ts = s["T_surf"][g]; P0 = s["P"][g][0]/1e6
    nb = int(s["is_burning"][g].sum())
    tg = " ".join(f"{Tg[i]:5.0f}" for i in idx)
    ts = " ".join(f"{Ts[i]:4.0f}" for i in idx)
    print(f"{s['t']*1e3:5.2f} {P0:7.2f} | {tg} | {ts} | {nb}")
