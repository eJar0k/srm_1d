"""SM1 Part 3.2b: the three canonical B0 rows (Chunc, BALLSstick, Hasegawa A; library BPNV igniter, Sutton mass,
forward plenum, DeMar as shipped) per solid-heating model. Usage: sm1_canon.py <motor chunc|balls|haseg> <model> <out.jsonl>
Run from the worktree root with PYTHONPATH=$PWD. Same metrics as baseline_gate.py plus first light / t10 / t90."""
import sys, json, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric

CANON = dict(roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None, T_ignition=756.0,
             print_interval=0.0, verbose=False)
MOT = {'chunc': ('machbusterNew.ric', 1.0, 0.01e6), 'balls': ('BALLSstick.ric', 2.0, 0.01e6),
       'haseg': ('hasegawa_a.ric', 10.0, 0.05e6)}

motor, model, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
ric, t_max, pc = MOT[motor]
t0 = time.time()
res = run_from_ric(str(ROOT / 'motors' / ric), t_max=t_max, P_cutoff=pc, solid_heating_model=model, **CANON)[0]
wall = time.time() - t0
t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6
s = res['summary']
isg = res['snapshots'][0]['is_grain']
ig = np.asarray(res['ignition_time_by_cell'])[isg]
fin = ig[ig < 1e8]
k = int(P.argmax())
m = dict(motor=motor, solid=model, peak=float(P[k]), t_peak=float(t[k]), n_steps=int(s['steps']), wall_s=wall,
         frac_ignited=float(fin.size / ig.size), t_first=float(fin.min()), t10=float(np.percentile(ig, 10)),
         t50=float(np.percentile(ig, 50)), t90=float(np.percentile(ig, 90)) if fin.size >= .9 * ig.size else None,
         t_last=float(fin.max()))
if motor == 'chunc':
    m['plateau'] = float(np.median(P[(t >= 0.4) & (t <= 0.9)])); m['ratio'] = m['peak'] / m['plateau']
if motor == 'balls':
    ms = t < 0.3
    kp = int(np.argmax(P[ms])); tpk = float(t[ms][kp]); m['peak'] = float(P[ms][kp]); m['t_peak'] = tpk
    m['trough'] = float(np.min(P[(t > tpk) & (t < 1.2)])); m['mid'] = float(P[int(np.argmin(np.abs(t - 1.03)))])
with open(out, 'a') as f:
    f.write(json.dumps(m) + '\n')
print(json.dumps(m))
