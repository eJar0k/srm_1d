"""IG3 stage-2 validation rows (library pyrogen as shipped: DeMar flux off, DeMar-sized auto charge + upsize fallback): canonical motors with their usual (library BPNV, Sutton-mass, forward_plenum) igniter.
Usage: ig3_valid.py '<json spec>' <out.jsonl>;  spec = {motor: chunc|balls|haseg, ign_h (default peretz), ap_ref (default local), demar}.
Mirrors docs/v0_7_4/probes/mp_cap_crossmotor.py rows (same knobs/windows). Run from the worktree root with
PYTHONPATH=$PWD."""
import sys, json, time, dataclasses
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen

spec = json.loads(sys.argv[1]); out = Path(sys.argv[2])
tag = f"s2_{spec['motor']}_{spec.get('pyrogen','bpnv')}"
if out.exists() and any(json.loads(l).get('tag') == tag for l in out.read_text().splitlines() if l.strip()):
    print('skip', tag); sys.exit(0)
M = dict(chunc=('machbusterNew.ric', 1.0, 0.01e6), balls=('BALLSstick.ric', 2.0, 0.01e6),
         haseg=('hasegawa_a.ric', 10.0, 0.05e6))[spec['motor']]
pyro = load_pyrogen(spec.get('pyrogen', 'bpnv'))
kw = dict(roughness=32e-6, kappa=0.44, pyrogen=pyro, pyrogen_mass=None, T_ignition=756.0, print_interval=0.0,
          verbose=False, erosive_model='ma2020', erosive_mp_cap=0.0, t_max=M[1], P_cutoff=M[2])
kw['ignition_heat_transfer'] = spec.get('ign_h', 'peretz')
if spec.get('ap_ref', 'local') != 'local':
    kw['peretz_ap_reference'] = spec['ap_ref']
t0 = time.time()
res = run_from_ric(str(ROOT / 'motors' / M[0]), **kw)[0]
wall = time.time() - t0
t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6; s = res['summary']
k = int(P.argmax())
row = dict(tag=tag, **spec, wall_s=wall, steps=int(s['steps']), termination=s.get('termination'),
           first_ign_s=s.get('first_ignition_time_s'), pyro_mass_g=float(s['pyrogen_mass_initial']) * 1e3,
           attempts=int(s['ignition_attempts']), attempt_masses_g=[m * 1e3 for m in s['ignition_attempt_masses_kg']],
           attempt_wall_s=s['ignition_attempt_wall_s'], ignited=bool(s['ignited']),
           never_lit=int(s['grain_cells_never_lit']), never_lit_cells=s['grain_never_lit_cell_range'],
           never_lit_x_mm=s['grain_never_lit_x_range_mm'], frac_lit=float(s['grain_fraction_lit']))
if spec['motor'] == 'chunc':
    pl = float(np.median(P[(t >= 0.4) & (t <= 0.9)]))
    row.update(peak=float(P[k]), t_peak=float(t[k]), plateau=pl, ratio=float(P[k]) / pl,
               P_at_ms={str(ms): float(np.interp(ms * 1e-3, t, P)) for ms in (10, 20, 30, 40, 60, 100, 200)})
elif spec['motor'] == 'balls':
    ms = t < 0.3; kp = int(np.argmax(P[ms])); Ppk, tpk = float(P[ms][kp]), float(t[ms][kp])
    after = (t > tpk) & (t < 1.2); Ptr = float(np.min(P[after])); Pmid = float(P[int(np.argmin(np.abs(t - 1.03)))])
    row.update(peak=Ppk, t_peak=tpk, trough=Ptr, mid=Pmid, pk_tr=Ppk / Ptr, pk_mid=Ppk / Pmid)
else:
    row.update(peak=float(P[k]), t_peak=float(t[k]), early_peak_lt_0p5s=float(P[t < 0.5].max()),
               late_peak_gt_1s=float(P[t > 1.0].max()), t_late_peak=float(t[t > 1.0][int(np.argmax(P[t > 1.0]))]))
with open(out, 'a') as f:
    f.write(json.dumps(row) + chr(10))
print(json.dumps(row))
