"""T3 step-4 measurement: pressure-floor activations vs grid resolution.

Chunc t_max=0.05 at target_propellant_cells in {100,200,400,800} and
BALLSstick t_max=0.05 at {100,400}; canonical knobs. Prints one JSON row per
run and writes <out>/measure_floor.json.

Run from the worktree root (SRM1D_OPENMOTOR_PATH exported):
    PYTHONPATH="$PWD" python docs/solver_verification/worklog/T3_scripts/measure_floor.py [out_dir]
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
import numpy as np  # noqa: E402

import srm_1d  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
print("[provenance] srm_1d from", srm_1d.__file__, flush=True)

CANON = dict(roughness=32e-6, kappa=0.44, T_ignition=756.0,
             pyrogen='bpnv', pyrogen_mass=None,
             print_interval=0.0, verbose=False)
CASES = [
    ('Chunc', 'motors/machbusterNew.ric', [100, 200, 400, 800]),
    ('BALLSstick', 'motors/BALLSstick.ric', [100, 400]),
]

rows = []
for name, ric, cell_list in CASES:
    for tpc in cell_list:
        t0 = time.time()
        res = run_from_ric(str(ROOT / ric), t_max=0.05,
                           target_propellant_cells=tpc, **CANON)[0]
        wall = time.time() - t0
        s = res['summary']
        mp = np.asarray(res['min_pressure'])
        cells = np.asarray(res['pressure_floor_cells'])
        mass = np.asarray(res['pressure_floor_mass'])
        mp_prod = float(s['mass_produced'])
        row = dict(
            motor=name, target_cells=tpc, cells=int(s['cells']),
            n_steps=int(s['steps']),
            activations=int(s['pressure_floor_activations']),
            floor_steps=int(s['pressure_floor_steps']),
            mass_kg=float(s['pressure_floor_mass_added_kg']),
            mass_produced=mp_prod,
            pyrogen_burned=float(s['pyrogen_mass_burned']),
            mass_frac=(float(s['pressure_floor_mass_added_kg']) / mp_prod
                       if mp_prod > 0 else float('nan')),
            energy_J=float(s['pressure_floor_energy_added_J']),
            first_time_s=float(s['pressure_floor_first_time_s']),
            min_raw_P=float(s['pressure_floor_min_raw_P']),
            global_min_pressure=float(np.min(mp)),
            t_of_global_min=float(res['time'][int(np.argmin(mp))]),
            peak_MPa=float(s['P_peak'] / 1e6), t_peak=float(s['t_peak']),
            termination=s['termination'], wall_s=wall,
            hist_sum_cells=int(round(cells.sum())),
            hist_sum_mass=float(mass.sum()),
        )
        rows.append(row)
        print(json.dumps(row), flush=True)

with open(OUT / 'measure_floor.json', 'w') as f:
    json.dump(rows, f, indent=1)
print("done")
