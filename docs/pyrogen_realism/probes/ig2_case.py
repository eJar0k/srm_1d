"""IG2 single-case runner (opt-in modes only). Usage: ig2_case.py '<json spec>' <out.jsonl>

spec keys: tag (unique), topo (head_basket|forward_plenum), g (pyrogen mass in GRAMS) or 'lb',
pyrogen (library name, default bpnv), cubes (bool: 3 mm cubes L/D=1), ign_h (gnielinski|peretz|keller),
demar (step|off), t_max, ric (default motors/machbusterNew.ric), snap (snapshot interval s).
Resumable: skips if tag already in the jsonl. Appends one JSON line. Run from the worktree root
with PYTHONPATH=$PWD (docs/solver_verification/TASKS.md)."""
import sys, json, time, dataclasses
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen



def _done(out, tag):
    if out.exists():
        for ln in out.read_text().splitlines():
            if ln.strip() and json.loads(ln).get('tag') == tag:
                return True
    return False


def run_case(spec):
    """Run one spec in-process and return the metrics dict (no file I/O)."""
    tag = spec['tag']
    LB = 0.45359237
    mass = spec['g'] * 1e-3 if 'g' in spec else spec['lb'] * LB
    kw = dict(roughness=32e-6, kappa=0.44, T_ignition=spec.get('T_ign', 756.0), print_interval=0.0,
              verbose=False, t_max=spec.get('t_max', 0.2), snapshot_interval=spec.get('snap', 1e-3),
              injection_topology=spec['topo'])
    if spec.get('ign_h', 'gnielinski') != 'gnielinski':
        kw['ignition_heat_transfer'] = spec['ign_h']
    if spec.get('erosive_model'):
        kw['erosive_model'] = spec['erosive_model']
    if spec.get('P_cutoff') is not None:
        kw['P_cutoff'] = spec['P_cutoff']
    over = {}
    if spec.get('demar', 'step') == 'off':
        over['heat_delivery_mode'] = 'none'
    if spec.get('cubes'):
        over['particle_diameter_m'] = 3e-3; over['particle_LD_ratio'] = 1.0
    pyro = dataclasses.replace(load_pyrogen(spec.get('pyrogen', 'bpnv')), **over)
    t0 = time.time()
    res, perf, nozzle, geo, prop = run_from_ric(str(ROOT / spec.get('ric', 'motors/machbusterNew.ric')),
                                                pyrogen=pyro, pyrogen_mass=mass, **kw)
    wall = time.time() - t0
    t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6
    s = res['summary']
    isg = res['snapshots'][0]['is_grain']
    ig = np.asarray(res['ignition_time_by_cell'])[isg]
    fin = ig[ig < 1e8]; frac = fin.size / max(ig.size, 1)
    k = int(P.argmax())
    tm = spec.get('plateau_win', (0.1, 0.2))
    w = (t >= tm[0]) & (t <= tm[1])
    m = dict(tag=tag, **{kk: spec[kk] for kk in spec if kk != 'tag'},
             peak_MPa=float(P[k]), t_peak_s=float(t[k]), P_end_MPa=float(P[-1]),
             plateau_MPa=float(np.median(P[w])) if w.any() else None,
             frac_ignited=float(frac), t_first_ign_s=float(fin.min()) if fin.size else None,
             t90_s=float(np.percentile(ig, 90)) if frac >= .9 else None,
             pyro_duration_s=float(s['pyrogen_duration']), pyro_mass_burned_g=float(s['pyrogen_mass_burned']) * 1e3,
             steps=int(s['steps']), wall_s=wall, n_grain=int(ig.size))
    return m


if __name__ == '__main__':
    spec = json.loads(sys.argv[1]); out = Path(sys.argv[2])
    if _done(out, spec['tag']):
        print('skip', spec['tag']); sys.exit(0)
    m = run_case(spec)
    with open(out, 'a') as f:
        f.write(json.dumps(m) + chr(10))
    print(json.dumps(m))
