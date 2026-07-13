"""
Structural ignition-spike probe for Chunc / machbusterNew.

Re-opening the investigation per user direction: treat the spike as a
candidate STRUCTURAL artifact, not assumed-faithful Ma physics. This script
gathers the raw evidence needed to discriminate among structural causes:

  1. Ignition simultaneity   -> ignition_time_by_cell spread (grain cells).
  2. Spike = erosive or not  -> baseline vs diagnostic_disable_erosive.
  3. Transient mass-flux G    -> rho*u profile at the spike vs at the plateau,
     and whether G is consistent with locally-generated mass (continuity
     check) or transiently inflated (a fill/velocity overshoot).

Canonical Chunc knobs (examples/machbusterNew.py): bpnv, Sutton sizing,
roughness 32um, kappa 0.44, T_ign 756.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"c:\Users\ejarocki\Documents\Rocketry\Code Stuff\Erosive Burning Solver\srm_1d")
sys.path.insert(0, str(ROOT))

from srm_1d.openmotor_adapter import run_from_ric

MOTOR = str(ROOT / "motors" / "machbusterNew.ric")
KN = dict(roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
          T_ignition=756.0, P_cutoff=0.01e6, t_max=1.3,
          snapshot_interval=0.002, print_interval=0.0, verbose=False)


def nearest_snap(snaps, t):
    return min(snaps, key=lambda s: abs(s["t"] - t))


def summarize(tag, res):
    t = np.asarray(res["time"]); P = np.asarray(res["P_head"]) / 1e6
    ipk = int(np.argmax(P)); tpk = float(t[ipk]); Ppk = float(P[ipk])
    plat_mask = (t >= 0.4) & (t <= 0.9)
    Pplat = float(np.median(P[plat_mask])) if plat_mask.any() else float("nan")
    print(f"\n=== {tag} ===")
    print(f"  P_peak = {Ppk:.2f} MPa @ t = {tpk*1e3:.1f} ms")
    print(f"  P_plateau(0.4-0.9s median) = {Pplat:.2f} MPa")
    print(f"  spike ratio P_peak/P_plateau = {Ppk/Pplat:.2f}")
    return tpk, Ppk, Pplat


def main():
    print("Running Chunc baseline (this compiles numba on first call)...")
    res = run_from_ric(MOTOR, **KN)[0]
    tpk, Ppk, Pplat = summarize("BASELINE", res)

    # --- 1. Ignition simultaneity ---
    ign = np.asarray(res["ignition_time_by_cell"])
    snaps = res["snapshots"]
    is_grain0 = snaps[0]["is_grain"]
    g = is_grain0 & (ign < 1e8)
    igt = ign[g]
    print("\n--- Ignition simultaneity (grain cells) ---")
    print(f"  grain cells = {int(is_grain0.sum())}, ignited = {g.sum()}")
    if igt.size:
        print(f"  t_first = {igt.min()*1e3:.2f} ms, t_last = {igt.max()*1e3:.2f} ms")
        print(f"  spread (max-min) = {(igt.max()-igt.min())*1e3:.2f} ms")
        print(f"  p10-p90 spread   = {(np.percentile(igt,90)-np.percentile(igt,10))*1e3:.2f} ms")

    # --- 3. Mass-flux profile at spike vs plateau ---
    for label, tt in (("SPIKE", tpk), ("PLATEAU", 0.6)):
        s = nearest_snap(snaps, tt)
        gg = s["is_grain"]
        x = s["x"][gg]
        rho = s["rho"][gg]; u = s["u"][gg]
        G = rho * np.abs(u)                      # local mass flux [kg/m^2/s]
        rt = s["r_total"][gg]; re = s["r_erosive"][gg]
        Dp = s["D_port"][gg]
        msrc = s["mass_source"][gg]              # [kg/s/m]
        A = np.pi/4.0 * Dp**2                     # conical bore -> circular
        ef = np.where(rt > 0, re/np.maximum(rt, 1e-12), 0.0)
        # continuity-consistent mass flux: cumulative generated mass / A
        dx = float(res["dx"])
        cum = np.cumsum(msrc * dx)                # ~ rho*u*A expected (QS)
        G_qs = cum / np.maximum(A, 1e-12)
        ipk = int(np.argmax(re)) if re.size else 0
        print(f"\n--- {label} snapshot t={s['t']*1e3:.1f} ms ---")
        print(f"  burning grain cells = {int(s['is_burning'][gg].sum())}/{gg.sum()}")
        print(f"  max G (rho|u|)      = {G.max():.1f} kg/m^2/s")
        print(f"  max r_total         = {rt.max()*1e3:.3f} mm/s")
        print(f"  max r_erosive       = {re.max()*1e3:.3f} mm/s  (at x={x[ipk]*1e3:.0f} mm)")
        print(f"  erosive frac @ that cell = {ef[ipk]:.2f}")
        print(f"  mean erosive frac (burning) = {ef[s['is_burning'][gg]].mean() if s['is_burning'][gg].any() else float('nan'):.2f}")
        print(f"  G(rho*u) vs G_qs(cum src/A) at max-erosive cell: "
              f"{G[ipk]:.1f} vs {G_qs[ipk]:.1f}  ratio={G[ipk]/max(G_qs[ipk],1e-9):.2f}")
        print(f"  D_port at max-erosive cell = {Dp[ipk]*1e3:.1f} mm")

    # --- 2. Erosive on/off differential ---
    print("\nRunning Chunc with erosive DISABLED...")
    res2 = run_from_ric(MOTOR, diagnostic_disable_erosive=True, **KN)[0]
    summarize("NO-EROSIVE", res2)


if __name__ == "__main__":
    main()
