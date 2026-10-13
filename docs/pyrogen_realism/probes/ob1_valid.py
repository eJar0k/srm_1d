"""OB1 Part 5.3 (post-rebase): canonical motors with their usual igniter (library BPNV, Sutton-mass, forward_plenum).
Usage: ob1_valid.py '<json spec>' <out.jsonl>; spec = {motor: chunc|balls|haseg, geom: lib|old32, demar: on|off,
law: lib|capZ_lib|zhu_mean}. Defaults = main's ignition h (Peretz), DeMar on, library I-A pellet. 'old32' = the pre-OB1
library pellet (3.2 mm, L/D 1; forward-plenum law identical to the pre-OB1 (m/m0)^(2/3)). Mirrors ig2_valid.py (same knobs and
windows). Run from the worktree root with PYTHONPATH=$PWD."""
import sys, json, time, dataclasses
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
from ob1_case import apply_law, GEOM

spec = json.loads(sys.argv[1])
out = Path(sys.argv[2])
tag = f"{spec['motor']}_{spec.get('geom', 'lib')}_{spec.get('demar', 'on')}_{spec.get('law', 'lib')}"
if out.exists() and any(json.loads(l).get('tag') == tag for l in out.read_text().splitlines() if l.strip()):
    print('skip', tag)
    sys.exit(0)
M = dict(chunc=('machbusterNew.ric', 1.0, 0.01e6), balls=('BALLSstick.ric', 2.0, 0.01e6),
         haseg=('hasegawa_a.ric', 10.0, 0.05e6))[spec['motor']]
over = dict(GEOM.get(spec.get('geom', 'lib'), {}))
if spec.get('demar', 'on') == 'off':
    over['heat_delivery_mode'] = 'none'
pyro = apply_law(dataclasses.replace(load_pyrogen('bpnv'), **over), spec.get('law', 'lib'))
kw = dict(roughness=32e-6, kappa=0.44, pyrogen=pyro, pyrogen_mass=None, T_ignition=756.0, print_interval=0.0,
          verbose=False, erosive_model='ma2020', erosive_mp_cap=0.0, t_max=M[1], P_cutoff=M[2])
t0 = time.time()
res = run_from_ric(str(ROOT / 'motors' / M[0]), **kw)[0]
wall = time.time() - t0
t = np.asarray(res['time'])
P = np.asarray(res['P_head']) / 1e6
s = res['summary']
k = int(P.argmax())
lit = float(P.max()) > 1.0
row = dict(tag=tag, **spec, wall_s=wall, steps=int(s['steps']), termination=s.get('termination'),
           first_ign_s=s.get('first_ignition_time_s'), pyro_mass_g=float(s['pyrogen_mass_initial']) * 1e3,
           pyro_duration_s=float(s['pyrogen_duration']), lit=lit)
if spec['motor'] == 'chunc':
    pl = float(np.median(P[(t >= 0.4) & (t <= 0.9)]))
    row.update(peak=float(P[k]), t_peak=float(t[k]), plateau=pl, ratio=float(P[k]) / pl if pl > 0 else None)
elif spec['motor'] == 'balls':
    ms = t < 0.3
    kp = int(np.argmax(P[ms]))
    Ppk, tpk = float(P[ms][kp]), float(t[ms][kp])
    after = (t > tpk) & (t < 1.2)
    Ptr = float(np.min(P[after]))
    Pmid = float(P[int(np.argmin(np.abs(t - 1.03)))])
    row.update(peak=Ppk, t_peak=tpk, trough=Ptr, mid=Pmid, pk_tr=Ppk / Ptr, pk_mid=Ppk / Pmid)
else:
    row.update(peak=float(P[k]), t_peak=float(t[k]))
with open(out, 'a') as f:
    f.write(json.dumps(row) + chr(10))
print(json.dumps(row))
