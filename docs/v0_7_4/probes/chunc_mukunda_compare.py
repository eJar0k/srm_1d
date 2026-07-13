"""Compare the sim's Ma-2020 erosive enhancement against the Mukunda-Paul
(1997) universal law at the Chunc spike vs plateau, per-cell.

Mukunda-Paul (Eq. 12):  eta = 1 + 0.023*(g^0.8 - 35^0.8) * H(g-35)
  g0  = G / (rho_p * r0)         r0 = NON-erosive (Saint-Robert) rate
  Re0 = rho_p * r0 * d0 / mu
  g   = g0 * (Re0/1000)^(-0.125)
Sim's Ma enhancement: eta_Ma = r_total / r0.

Tests: (a) does MP's threshold (g_th=35) suppress the transient? (g>>35 => no)
       (b) does Ma over-predict vs the universal law, and by how much?
       (c) would MP itself produce a spike (eta>1) at the transient condition?
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"c:\Users\ejarocki\Documents\Rocketry\Code Stuff\Erosive Burning Solver\srm_1d")
sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric

RHO_P = 1700.0; MU = 9.125e-5; GTH = 35.0
res = run_from_ric(str(ROOT/"motors"/"machbusterNew.ric"),
                   roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
                   T_ignition=756.0, P_cutoff=0.01e6, t_max=1.3,
                   snapshot_interval=0.002, print_interval=0.0, verbose=False)[0]
snaps = res["snapshots"]

def mp_eta(G, r0, d0):
    if r0 <= 0 or d0 <= 0: return 1.0, 0.0
    g0 = G/(RHO_P*r0)
    Re0 = RHO_P*r0*d0/MU
    g = g0*(Re0/1000.0)**(-0.125)
    if g <= GTH: return 1.0, g
    return 1.0 + 0.023*(g**0.8 - GTH**0.8), g

def near(t):
    return min(snaps, key=lambda s: abs(s["t"]-t))

for label, tt in (("SPIKE", 0.012), ("PLATEAU", 0.6)):
    s = near(tt); g = s["is_grain"] & s["is_burning"]
    rt = s["r_total"][g]; re = s["r_erosive"][g]; r0 = rt-re
    G = s["rho"][g]*np.abs(s["u"][g]); Dp = s["D_port"][g]
    eta_ma = np.where(r0>0, rt/np.maximum(r0,1e-9), 1.0)
    eta_mp = np.array([mp_eta(G[i], r0[i], Dp[i])[0] for i in range(len(r0))])
    gval   = np.array([mp_eta(G[i], r0[i], Dp[i])[1] for i in range(len(r0))])
    # focus on the peak-erosive cell (drives the spike)
    ipk = int(np.argmax(re))
    print(f"\n=== {label} (t={s['t']*1e3:.0f} ms) ===")
    print(f"  cells burning: {g.sum()}, g range [{gval.min():.0f}, {gval.max():.0f}], "
          f"cells with g>35: {(gval>GTH).sum()}/{len(gval)}")
    print(f"  peak-erosive cell: G={G[ipk]:.0f}, r0={r0[ipk]*1e3:.2f} mm/s, "
          f"D={Dp[ipk]*1e3:.1f} mm, g={gval[ipk]:.0f}")
    print(f"     eta_Ma (sim)   = {eta_ma[ipk]:.2f}")
    print(f"     eta_MukundaPaul= {eta_mp[ipk]:.2f}")
    print(f"     Ma / MP        = {eta_ma[ipk]/max(eta_mp[ipk],1e-9):.2f}")
    print(f"  domain max eta_Ma={eta_ma.max():.2f}, max eta_MP={eta_mp.max():.2f}")
