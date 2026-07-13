"""Acceleration-parameter (relaminarization) diagnostic for Chunc.

Tests whether the favorable-pressure-gradient relaminarization criterion
  K = (nu / u^2) * (du/dx)   [Kays-Crawford / Narasimha-Sreenivasan; crit ~3e-6]
is satisfied at the erosive-spike cells DURING the ignition spike but NOT at
the steady plateau. If yes, accel-relaminarization is a live, transient,
plateau-preserving lever to gate Ma's (turbulent) erosive term. If K exceeds
the threshold at the plateau too, the lever is dead (it would wrongly suppress
steady-state erosion).

nu = mu/rho (per cell), u = u_cell, du/dx = central diff of u along x.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"c:\Users\ejarocki\Documents\Rocketry\Code Stuff\Erosive Burning Solver\srm_1d")
sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric

MU = 9.125e-5; KCRIT = 3.0e-6
res = run_from_ric(str(ROOT/"motors"/"machbusterNew.ric"),
                   roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
                   T_ignition=756.0, P_cutoff=0.01e6, t_max=1.3,
                   snapshot_interval=0.002, print_interval=0.0, verbose=False)[0]
snaps = res["snapshots"]
xfull = snaps[0]["x"]

def K_profile(s):
    """Per-cell K over the whole domain (NaN where u too small)."""
    u = s["u"]; rho = s["rho"]; x = s["x"]
    nu = MU/np.maximum(rho, 1e-9)
    dudx = np.gradient(u, x)              # central difference
    with np.errstate(divide='ignore', invalid='ignore'):
        K = nu * dudx / np.maximum(u*u, 1e-12)
    return K, u, nu, dudx

def near(t):
    return min(snaps, key=lambda s: abs(s["t"]-t))

for label, tt in (("SPIKE", 0.012), ("PLATEAU", 0.6)):
    s = near(tt)
    K, u, nu, dudx = K_profile(s)
    g = s["is_grain"] & s["is_burning"]
    re = s["r_erosive"]
    # erosive-active cells: burning grain with meaningful erosion
    act = g & (re > 0.2*re[g].max() if re[g].max() > 0 else g)
    ipk = int(np.argmax(np.where(g, re, -1)))   # peak-erosive cell
    print(f"\n=== {label} (t={s['t']*1e3:.0f} ms) ===")
    print(f"  peak-erosive cell x={xfull[ipk]*1e3:.0f}mm: "
          f"u={u[ipk]:.0f} m/s, nu={nu[ipk]:.2e}, du/dx={dudx[ipk]:.0f} 1/s")
    print(f"     K = {K[ipk]:.2e}   (crit {KCRIT:.0e}; relaminarizing if K>crit) "
          f"-> {'RELAMINARIZING' if K[ipk]>KCRIT else 'turbulent'}")
    Kact = K[act]
    Kact = Kact[np.isfinite(Kact)]
    if Kact.size:
        frac = (Kact > KCRIT).mean()
        print(f"  erosive-active cells (n={act.sum()}): "
              f"median K={np.median(Kact):.2e}, "
              f"K range [{Kact.min():.1e}, {Kact.max():.1e}]")
        print(f"     fraction with K > {KCRIT:.0e} (relaminarizing): {frac*100:.0f}%")

# Time history of K at a mid-grain and aft cell through the spike
print("\n--- K(t) at fixed stations through the ignition transient ---")
midx = int(np.argmin(np.abs(xfull - 0.435)))   # ~mid grain
aftx = int(np.argmin(np.abs(xfull - 0.825)))   # ~peak-erosive station
print(f"{'t_ms':>6} {'K_mid':>10} {'rE_mid':>8} {'K_aft':>10} {'rE_aft':>8}")
for s in snaps:
    if s["t"] > 0.05: break
    K, u, nu, dudx = K_profile(s)
    rE = s["r_erosive"]
    print(f"{s['t']*1e3:6.0f} {K[midx]:10.2e} {rE[midx]*1e3:8.2f} "
          f"{K[aftx]:10.2e} {rE[aftx]*1e3:8.2f}")
