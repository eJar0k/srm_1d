# SM1 Part 3.1: BALLSstick dead zone (4.8 g BPNV, forward plenum, DeMar off) per solid-heating model.
# Usage: python sm1_deadzone_model.py <outdir> <model> [charge_g] [t_max]
# <outdir>/BALLSstick.ric must exist (run sm1_deadzone_probe.py once; it writes the committed copy).
# Writes <outdir>/dz_<model>.json: ignition time per cell [s] (null = never lit) + wall time.
import sys, json, time, dataclasses, numpy as np
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
out, model = sys.argv[1], sys.argv[2]
charge = float(sys.argv[3]) if len(sys.argv) > 3 else 4.8
t_max = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
pyro = dataclasses.replace(load_pyrogen('bpnv'), heat_delivery_mode='none')
t0 = time.time()
res = run_from_ric(f'{out}/BALLSstick.ric', pyrogen=pyro, pyrogen_mass=charge * 1e-3, t_max=t_max,
                   P_cutoff=0.01e6, roughness=32e-6, kappa=0.44, T_ignition=756.0,
                   snapshot_interval=0.01, print_interval=0.0, verbose=False,
                   solid_heating_model=model)[0]
wall = time.time() - t0
ign = np.asarray(res['ignition_time_by_cell'])
gi = np.where(res['snapshots'][0]['is_grain'])[0]
lit = {int(i): (float(ign[i]) if ign[i] < 1e9 else None) for i in gi}
P = np.asarray(res['P_head']) if 'P_head' in res else None
json.dump({'model': model, 'charge_g': charge, 'wall_s': wall, 'ign': lit,
           'steps': int(len(res['time'])) if 'time' in res else None,
           'P_peak': float(np.max(P)) if P is not None else None}, open(f'{out}/dz_{model}_{charge}.json', 'w'))
never = [i for i, v in lit.items() if v is None]
print(model, 'wall %.1f s' % wall, 'never lit:', never[:3], '..', never[-3:] if never else '', len(never))
