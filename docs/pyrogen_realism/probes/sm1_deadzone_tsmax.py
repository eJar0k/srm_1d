# SM1: peak surface temperature reached by each head-end cell per solid model (first t_max s; snapshots every 0.5 ms).
# Usage: python sm1_deadzone_tsmax.py <balls|haseg> <charge_g> [t_max]
import sys, dataclasses
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
motor, ch = sys.argv[1], float(sys.argv[2]); tm = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
ric = {'balls': 'BALLSstick.ric', 'haseg': 'hasegawa_a.ric'}[motor]
pyro = dataclasses.replace(load_pyrogen('bpnv'), heat_delivery_mode='none')
rows = {}
for m in ('goodman', 'peretz', 'conduction'):
    r = run_from_ric(str(ROOT / 'motors' / ric), pyrogen=pyro, pyrogen_mass=ch * 1e-3, t_max=tm, P_cutoff=0.01e6,
                     roughness=32e-6, kappa=0.44, T_ignition=756.0, snapshot_interval=5e-4, print_interval=0.0,
                     verbose=False, solid_heating_model=m)[0]
    S = r['snapshots']; ts = np.array([s['t'] for s in S]); Tm = np.max([s['T_surf'] for s in S], axis=0)
    Tend = S[-1]['T_surf']; rows[m] = (Tm, Tend, S[-1]['t'])
print('cell | peak T_s [K] (goodman peretz conduction) | T_s at t=%.3f s' % rows['goodman'][2])
for i in range(3, 14):
    print(i, ' '.join('%6.0f' % rows[m][0][i] for m in rows), '|', ' '.join('%6.0f' % rows[m][1][i] for m in rows))
