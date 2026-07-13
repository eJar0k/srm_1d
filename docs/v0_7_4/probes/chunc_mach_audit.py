"""Confirm the fill is supersonic (solver's own Mach field) and that the
fill velocity is continuity-consistent (not a PISO blowup). If the ignition
gate is driven by Mach>1 blowdown velocities fed into the steady-pipe-flow
Gnielinski correlation, that is a structural misuse outside the correlation's
validity -> over-predicted ignition heat transfer -> collapsed ignition seq."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"c:\Users\ejarocki\Documents\Rocketry\Code Stuff\Erosive Burning Solver\srm_1d")
sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric

res = run_from_ric(str(ROOT / "motors" / "machbusterNew.ric"),
                   roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
                   T_ignition=756.0, P_cutoff=0.01e6, t_max=0.02,
                   snapshot_interval=0.0001, print_interval=0.0, verbose=False)[0]
snaps = res["snapshots"]
g = snaps[0]["is_grain"]
dx = float(res["dx"])

print(f"{'t_ms':>5} {'nburn':>5} {'maxMach':>7} {'medMach_fill':>12} "
      f"{'max|u|':>7} {'maxRe':>9} {'P0_MPa':>7}")
MU = 9.125e-5
for s in snaps:
    if s["t"] > 0.004:
        break
    gg = g
    M = s["Mach"][gg]; u = s["u"][gg]; rho = s["rho"][gg]; Dp = s["D_port"][gg]
    burn = s["is_burning"][gg]
    Re = rho*np.abs(u)*Dp/MU
    medM = np.median(np.abs(M[burn])) if burn.any() else 0.0
    print(f"{s['t']*1e3:5.2f} {int(burn.sum()):5d} {np.abs(M).max():7.2f} "
          f"{medM:12.2f} {np.abs(u).max():7.0f} {Re.max():9.0f} {s['P'][gg][0]/1e6:7.2f}")

# Continuity check during the fill: rho*u*A at face i ~ cumulative mass_source
# upstream. Pick the mid-fill frame.
frame = None
for s in snaps:
    if 30 <= int(s["is_burning"][g].sum()) <= 70:
        frame = s; break
if frame is not None:
    gg = g
    u = frame["u"][gg]; rho = frame["rho"][gg]; Dp = frame["D_port"][gg]
    msrc = frame["mass_source"][gg]
    A = np.pi/4*Dp**2
    flux = rho*u*A                       # local rho*u*A
    cum = np.cumsum(msrc*dx)             # cumulative generated mdot (QS expectation)
    print(f"\nContinuity during fill (t={frame['t']*1e3:.2f} ms):")
    print(f"{'cell':>4} {'rho*u*A':>9} {'cum_src':>9} {'ratio':>6}")
    for i in range(0, 40, 4):
        r = flux[i]/cum[i] if abs(cum[i])>1e-9 else float('nan')
        print(f"{i:>4} {flux[i]:9.4f} {cum[i]:9.4f} {r:6.2f}")
