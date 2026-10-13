# SM1 Part 3.1: head-end dead zone (explicit BPNV charge, forward plenum, DeMar off) per solid-heating model.
# Usage (repo root, PYTHONPATH=$PWD): python sm1_deadzone_model.py <outjson_dir> <motor balls|haseg> <model> <charge_g> [t_max]
# Writes <outdir>/dz_<motor>_<model>_<charge>.json: ignition time per grain cell [s] (null = never lit), P history summary.
import sys, json, time, dataclasses
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
out, motor, model, charge = Path(sys.argv[1]), sys.argv[2], sys.argv[3], float(sys.argv[4])
RIC = {'balls': ('BALLSstick.ric', 2.0, 0.01e6), 'haseg': ('hasegawa_a.ric', 10.0, 0.05e6)}
ric, t_max, pc = RIC[motor]
if len(sys.argv) > 5: t_max = float(sys.argv[5])
pyro = dataclasses.replace(load_pyrogen('bpnv'), heat_delivery_mode='none')
t0 = time.time()
res = run_from_ric(str(ROOT / 'motors' / ric), pyrogen=pyro, pyrogen_mass=charge * 1e-3, t_max=t_max,
                   P_cutoff=pc, roughness=32e-6, kappa=0.44, T_ignition=756.0, snapshot_interval=0.01,
                   print_interval=0.0, verbose=False, solid_heating_model=model)[0]
wall = time.time() - t0
ign = np.asarray(res['ignition_time_by_cell']); x = np.asarray(res['snapshots'][0]['x'])
gi = np.where(res['snapshots'][0]['is_grain'])[0]
lit = {int(i): (float(ign[i]) if ign[i] < 1e9 else None) for i in gi}
t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6
late = t > 0.5
json.dump({'motor': motor, 'model': model, 'charge_g': charge, 'wall_s': wall, 'steps': int(res['summary']['steps']),
           'ign': lit, 'x_mm': {int(i): float(x[i] * 1e3) for i in gi},
           'P_max': float(P.max()), 't_Pmax': float(t[P.argmax()]),
           'P_peak_after_0p5s': float(P[late].max()) if late.any() else None,
           't_end': float(t[-1]), 'termination': res['summary'].get('termination')},
          open(out / f'dz_{motor}_{model}_{charge}.json', 'w'))
never = [i for i, v in lit.items() if v is None]
print(motor, model, charge, 'wall %.1f s' % wall, 'never lit:', never)
