"""OB1 single-case runner. Usage: ob1_case.py '<json spec>' <out.jsonl>

spec keys: tag, topo (head_basket|forward_plenum), g (grams) or lb, pyrogen (bpnv|mtv, default bpnv),
cubes (3 mm cubes L/D 1; default True for head_basket), combo ('current'|'forward'; forward = peretz ignition h +
DeMar off), law ('lib' | 'zhu_mean' | 'zhu_hi' | 'zhu_lo' | 'cap7_lib' | 'capZ_lib'), area ('constant'|'shrinking'),
bed_v (diagnostic bed flame speed m/s, needs area=shrinking), t_max, plateau_win, ric, snap, erosive_model.

Laws (all opt-in; 'lib' = shipped default):
  zhu_*   : Zhu 2008 piecewise  r = r0*P (P<=7 MPa), 7*r0 (P>7 MPa);  r0 mean/hi/lo of Zhu Table 2 (m/s/MPa)
  cap7_lib: library a,n capped at P = 7 MPa (rate flat above 7 MPa)
  capZ_lib: library a,n capped at the Zhu plateau RATE (mean 7*r0 = 42.5 mm/s) -> P_cap = (r/a)^(1/n)
Resumable (tag in jsonl -> skip). Run from the worktree root with PYTHONPATH=$PWD."""
import sys, json, time, dataclasses
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen

R0 = {'zhu_mean': 6.07375e-3, 'zhu_hi': 6.848e-3, 'zhu_lo': 5.271e-3}   # m/s/MPa, Zhu 2008 Table 2 (mean of 4)


def apply_law(pyro, law):
    if law == 'lib':
        return pyro
    if law in R0:
        return dataclasses.replace(pyro, a=R0[law] / 1e6, n=1.0, P_plateau=7e6)
    if law == 'cap7_lib':
        return dataclasses.replace(pyro, P_plateau=7e6)
    if law == 'capZ_lib':
        r_plat = 7.0 * R0['zhu_mean']                    # m/s
        return dataclasses.replace(pyro, P_plateau=(r_plat / pyro.a) ** (1.0 / pyro.n))
    raise ValueError(law)


def _done(out, tag):
    if out.exists():
        for ln in out.read_text().splitlines():
            if ln.strip() and json.loads(ln).get('tag') == tag:
                return True
    return False


def run_case(spec, return_res=False):
    tag = spec['tag']
    LB = 0.45359237
    mass = spec['g'] * 1e-3 if 'g' in spec else spec['lb'] * LB
    topo = spec.get('topo', 'head_basket')
    kw = dict(roughness=32e-6, kappa=0.44, T_ignition=spec.get('T_ign', 756.0), print_interval=0.0,
              verbose=False, t_max=spec.get('t_max', 0.2), snapshot_interval=spec.get('snap', 1e-3),
              injection_topology=topo)
    if spec.get('combo', 'current') == 'forward':
        kw['ignition_heat_transfer'] = 'peretz'
    if spec.get('erosive_model'):
        kw['erosive_model'] = spec['erosive_model']
    if spec.get('area', 'constant') != 'constant':
        kw['basket_area_law'] = spec['area']
    if spec.get('bed_v', 0.0) > 0:
        kw['bed_flame_speed_m_s'] = spec['bed_v']
    over = {}
    if spec.get('combo', 'current') == 'forward':
        over['heat_delivery_mode'] = 'none'
    if spec.get('cubes', topo == 'head_basket'):
        over['particle_diameter_m'] = 3e-3
        over['particle_LD_ratio'] = 1.0
    if spec.get('d_mm'):
        over['particle_diameter_m'] = spec['d_mm'] * 1e-3
    pyro = dataclasses.replace(load_pyrogen(spec.get('pyrogen', 'bpnv')), **over)
    pyro = apply_law(pyro, spec.get('law', 'lib'))
    t0 = time.time()
    res, perf, nozzle, geo, prop = run_from_ric(str(ROOT / spec.get('ric', 'motors/machbusterNew.ric')),
                                                pyrogen=pyro, pyrogen_mass=mass, **kw)
    wall = time.time() - t0
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    s = res['summary']
    isg = res['snapshots'][0]['is_grain']
    ig = np.asarray(res['ignition_time_by_cell'])[isg]
    fin = ig[ig < 1e8]
    frac = fin.size / max(ig.size, 1)
    k = int(P.argmax())
    tm = spec.get('plateau_win', (0.1, 0.2))
    w = (t >= tm[0]) & (t <= tm[1])
    mp = np.asarray(res['m_pyrogen'])
    m0 = mp[0] if mp.size else 0.0

    def t_frac(f):
        if m0 <= 0:
            return None
        idx = np.nonzero(mp <= (1.0 - f) * m0)[0]
        return float(t[idx[0]]) if idx.size else None
    plat = float(np.median(P[w])) if w.any() else None
    m = dict(tag=tag, **{kk: spec[kk] for kk in spec if kk != 'tag'},
             peak_MPa=float(P[k]), t_peak_s=float(t[k]), P_end_MPa=float(P[-1]), plateau_MPa=plat,
             peak_over_plateau=(float(P[k]) / plat if plat else None),
             frac_ignited=float(frac), t_first_ign_s=float(fin.min()) if fin.size else None,
             t90_s=float(np.percentile(ig, 90)) if frac >= .9 else None,
             pyro_duration_s=float(s['pyrogen_duration']), pyro_t50_s=t_frac(0.5), pyro_t90_s=t_frac(0.9),
             pyro_mass_burned_g=float(s['pyrogen_mass_burned']) * 1e3,
             steps=int(s['steps']), wall_s=wall, n_grain=int(ig.size))
    if return_res:
        return m, res
    return m


if __name__ == '__main__':
    spec = json.loads(sys.argv[1])
    out = Path(sys.argv[2])
    if _done(out, spec['tag']):
        print('skip', spec['tag'])
        sys.exit(0)
    m = run_case(spec)
    with open(out, 'a') as f:
        f.write(json.dumps(m) + chr(10))
    print(json.dumps(m))
