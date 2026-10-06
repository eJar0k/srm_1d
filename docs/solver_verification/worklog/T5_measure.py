"""T5 measurements: Chunc refresh gaps + default vs split vs every-step.

Run from the worktree root with PYTHONPATH="$PWD".
usage: python measure.py <t_max> [out.json]
"""
import json
import sys
import time as clock

import numpy as np

import srm_1d
from srm_1d.openmotor_adapter import run_from_ric

print("[provenance]", srm_1d.__file__)

t_max = float(sys.argv[1])
out = sys.argv[2] if len(sys.argv) > 2 else None

CANON = dict(roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
             T_ignition=756.0, print_interval=0.0, verbose=False, t_max=t_max)

CONFIGS = {
    "default/default": {},
    "1/default(geom=21)": dict(burn_update_interval=1, geometry_update_interval=21),
    "default(21)/1": dict(burn_update_interval=21, geometry_update_interval=1),
    "1/1": dict(burn_update_interval=1, geometry_update_interval=1),
}

rows = []
for name, kw in CONFIGS.items():
    a = dict(CANON)
    a.update(kw)
    w0 = clock.time()
    result, *_ = run_from_ric("motors/machbusterNew.ric", **a)
    wall = clock.time() - w0
    s = result["summary"]
    P = np.asarray(result["P_head"])
    t = np.asarray(result["time"])
    minP = np.asarray(result["min_pressure"])
    ip = int(np.argmax(P))
    row = dict(
        config=name, n_steps=int(s["steps"]), cells=int(s["cells"]),
        peak_MPa=float(P[ip] / 1e6), t_peak_s=float(t[ip]),
        global_min_P_Pa=float(minP.min()), t_min_P_s=float(t[int(np.argmin(minP))]),
        floor_activations=int(s["pressure_floor_activations"]),
        floor_min_raw_P=float(s["pressure_floor_min_raw_P"]),
        burn_k=s["burn_refresh_interval_steps"],
        geom_k=s["geometry_refresh_interval_steps"],
        max_burn_gap_us=float(s["max_burn_refresh_gap_s"] * 1e6),
        max_geom_gap_us=float(s["max_geometry_refresh_gap_s"] * 1e6),
        dt_max_us=float(np.max(result["dt"]) * 1e6),
        dt_median_us=float(np.median(result["dt"]) * 1e6),
        wall_s=wall, steps_per_s=s["steps"] / wall,
    )
    rows.append(row)
    print(json.dumps(row), flush=True)

base = rows[0]
print(f"\nt_max = {t_max}")
print(f"{'config':22s} {'steps':>8s} {'peak MPa':>10s} {'dPeak%':>8s} "
      f"{'minP Pa':>12s} {'dMinP%':>8s} {'t_minP s':>9s} "
      f"{'burn gap us':>11s} {'geom gap us':>11s} {'floor':>6s} {'rawP min':>10s}")
for r in rows:
    dpk = 100 * (r["peak_MPa"] - base["peak_MPa"]) / base["peak_MPa"]
    dmn = 100 * (r["global_min_P_Pa"] - base["global_min_P_Pa"]) / base["global_min_P_Pa"]
    print(f"{r['config']:22s} {r['n_steps']:8d} {r['peak_MPa']:10.5f} {dpk:+8.4f} "
          f"{r['global_min_P_Pa']:12.2f} {dmn:+8.3f} {r['t_min_P_s']:9.5f} "
          f"{r['max_burn_gap_us']:11.3f} {r['max_geom_gap_us']:11.3f} "
          f"{r['floor_activations']:6d} {r['floor_min_raw_P']:10.2f}")
if out:
    with open(out, "w") as f:
        json.dump(rows, f, indent=1)
