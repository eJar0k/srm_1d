"""SM1 single-case runner: the OB1 case runner plus spec key solid = goodman (default) | peretz | conduction
(run_simulation(solid_heating_model=...)). Usage: sm1_case.py '<json spec>' <out.jsonl>

spec keys: tag, topo (head_basket|forward_plenum), g (grams) or lb, pyrogen (bpnv|mtv, default bpnv),
geom ('lib' = library I-A pellet 3.175 mm L/D 1.504 | 'cube3' = 3 mm cubes L/D 1 | 'old32' = pre-OB1 library 3.2 mm L/D 1),
area ('shrinking' default | 'constant' = pre-OB1 pooled constant-area basket law), demar ('on' default = library
heat_delivery_mode 'demar' | 'off' = 'none'), law ('lib' | 'zhu_mean' | 'zhu_hi' | 'zhu_lo' | 'cap7_lib' | 'capZ_lib'),
bed_v (diagnostic), t_max, plateau_win, ric, snap, erosive_model, T_ign, ign_h (override).

Laws (all opt-in): zhu_*: Zhu 2008 piecewise r = r0 P (P <= 7 MPa), 7 r0 above; cap7_lib / capZ_lib: library a, n capped at
7 MPa / at the Zhu plateau RATE. Resumable (tag in jsonl -> skip). Run from the worktree root with PYTHONPATH=$PWD."""
import sys, json, time, dataclasses
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen

R0 = {'zhu_mean': 6.07375e-3, 'zhu_hi': 6.848e-3, 'zhu_lo': 5.271e-3}   # m/s/MPa, Zhu 2008 Table 2
GEOM = {'cube3': dict(particle_diameter_m=3.0e-3, particle_LD_ratio=1.0),
        'old32': dict(particle_diameter_m=3.2e-3, particle_LD_ratio=1.0)}


def apply_law(pyro, law):
    if law == 'lib':
        return pyro
    if law in R0:
        return dataclasses.replace(pyro, a=R0[law] / 1e6, n=1.0, P_plateau=7e6)
    if law == 'cap7_lib':
        return dataclasses.replace(pyro, P_plateau=7e6)
    if law == 'capZ_lib':
        r_plat = 7.0 * R0['zhu_mean']
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
    mass = spec['g'] * 1e-3 if 'g' in spec else (spec['lb'] * LB if 'lb' in spec else None)
    topo = spec.get('topo', 'head_basket')
    kw = dict(roughness=32e-6, kappa=0.44, T_ignition=spec.get('T_ign', 756.0), print_interval=0.0,
              verbose=False, t_max=spec.get('t_max', 0.2), snapshot_interval=spec.get('snap', 1e-3),
              injection_topology=topo)
    kw['solid_heating_model'] = spec.get('solid', 'goodman')
    if spec.get('ign_h'):
        kw['ignition_heat_transfer'] = spec['ign_h']
    if spec.get('erosive_model'):
        kw['erosive_model'] = spec['erosive_model']
    if spec.get('P_cutoff') is not None:
        kw['P_cutoff'] = spec['P_cutoff']
    if topo != 'forward_plenum' and spec.get('area', 'shrinking') != 'shrinking':
        kw['basket_area_law'] = spec['area']
    if spec.get('bed_v', 0.0) > 0:
        kw['bed_flame_speed_m_s'] = spec['bed_v']
    over = dict(GEOM.get(spec.get('geom', 'lib'), {}))
    if spec.get('demar', 'on') == 'off':
        over['heat_delivery_mode'] = 'none'
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
             pyro_mass_g=float(s['pyrogen_mass_initial']) * 1e3,
             pyro_mass_burned_g=float(s['pyrogen_mass_burned']) * 1e3,
             steps=int(s['steps']), wall_s=wall, n_grain=int(ig.size),
             t10_s=float(np.percentile(ig, 10)) if frac >= .9 else None,
             t50_s=float(np.percentile(ig, 50)) if frac >= .5 else None)
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
