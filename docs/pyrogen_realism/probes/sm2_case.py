"""SM2 single-case runner. Usage: sm2_case.py '<json spec>' <out.jsonl>
spec keys:
  tag, motor (chunc|chunc_hb|balls|haseg), solid (goodman|peretz|conduction),
  charge (grams; 'auto' = IG3 DeMar sizing + x1.25 fallback), demar ('off' default | 'on'),
  ign_h ('peretz' default | 'gnielinski'), k_solid (default 0.271), T_ign (default 756),
  cells (target_propellant_cells, default 100), t_max (override), rec (record_forcing_steps),
  forcing_dir (where to write <tag>_forcing.npz; default C:/tmp/sv/sm2_data).
Run from the repo root with PYTHONPATH=$PWD. Appends one JSON line; resumable (tag in jsonl -> skip)."""
import sys, json, time, dataclasses, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen

# motor -> (ric, topology, default fixed charge [g], t_max, P_cutoff, plateau window, label)
MOTORS = {
    'chunc':    ('machbusterNew.ric', 'forward_plenum', 2.62, 1.0, 0.01e6, (0.4, 0.9)),
    'chunc_hb': ('machbusterNew.ric', 'head_basket', 4.5, 0.2, 0.01e6, (0.1, 0.2)),
    'balls':    ('BALLSstick.ric', 'forward_plenum', 4.83, 2.0, 0.01e6, None),
    'haseg':    ('hasegawa_a.ric', 'forward_plenum', 14.98, 10.0, 0.05e6, None),
}
P_EXP_HASEG = 6.44


def _done(out, tag):
    if out.exists():
        for ln in out.read_text().splitlines():
            if ln.strip() and json.loads(ln).get('tag') == tag:
                return True
    return False


def run_case(spec, return_res=False):
    ric, topo, g0, tmax0, pc, plat = MOTORS[spec['motor']]
    t_max = spec.get('t_max', tmax0)
    charge = spec.get('charge', g0)
    over = {}
    if spec.get('demar', 'off') == 'on':
        over['heat_delivery_mode'] = 'demar'
    pyro = dataclasses.replace(load_pyrogen('bpnv'), **over)
    kw = dict(roughness=32e-6, kappa=0.44, T_ignition=spec.get('T_ign', 756.0), print_interval=0.0,
              verbose=False, t_max=t_max, P_cutoff=pc, snapshot_interval=spec.get('snap', 0.05),
              injection_topology=topo, solid_heating_model=spec['solid'],
              ignition_heat_transfer=spec.get('ign_h', 'peretz'),
              target_propellant_cells=int(spec.get('cells', 100)),
              record_forcing_steps=int(spec.get('rec', 0)))
    if spec.get('k_solid') is not None:
        kw['k_solid'] = spec['k_solid']
    mass = None if charge == 'auto' else charge * 1e-3
    t0 = time.time()
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        res, perf, nozzle, geo, prop = run_from_ric(str(ROOT / 'motors' / ric), pyrogen=pyro,
                                                    pyrogen_mass=mass, **kw)
    wall = time.time() - t0
    s = res['summary']
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    isg = np.asarray(res['snapshots'][0]['is_grain'])
    ig = np.asarray(res['ignition_time_by_cell'])[isg]
    diag = np.asarray(res['ignition_diag'])[isg]
    n = ig.size
    lit = ig[ig < 1e9]
    srt = np.sort(lit)

    def tq(q):
        need = int(np.ceil(q / 100.0 * n))
        return float(srt[need - 1]) if srt.size >= need else None
    k = int(P.argmax())
    m = dict(tag=spec['tag'], **{kk: v for kk, v in spec.items() if kk != 'tag'},
             charge_final_g=float(s['pyrogen_mass_final_kg']) * 1e3, resized=bool(s['pyrogen_resized']),
             n_grain=int(n), n_lit=int(lit.size), n_never=int(n - lit.size),
             never_cells=[int(i) for i in np.flatnonzero(ig >= 1e9)],
             t_first=float(srt[0]) if srt.size else None, t10=tq(10), t50=tq(50), t90=tq(90),
             t_last=float(srt[-1]) if srt.size else None,
             peak=float(P[k]), t_peak=float(t[k]), steps=int(s['steps']), wall_s=wall,
             termination=int(s['termination_code']))
    m['dur90'] = (m['t90'] - m['t_first']) if (m['t90'] is not None and m['t_first'] is not None) else None
    if plat:
        w = (t >= plat[0]) & (t <= plat[1])
        m['plateau'] = float(np.median(P[w])) if w.any() else None
        m['peak_over_plateau'] = (m['peak'] / m['plateau']) if m['plateau'] else None
    if spec['motor'] == 'balls':
        ms = t < 0.3
        kp = int(np.argmax(P[ms])); tpk = float(t[ms][kp])
        after = (t > tpk) & (t < 1.2)
        m['peak'] = float(P[ms][kp]); m['t_peak'] = tpk
        m['trough'] = float(P[after].min()); m['mid'] = float(P[int(np.argmin(np.abs(t - 1.03)))])
        m['pk_trough'] = m['peak'] / m['trough']; m['pk_mid'] = m['peak'] / m['mid']
    if spec['motor'] == 'haseg' and (t > 0.5).any():
        late = t > 0.5
        m['late_peak'] = float(P[late].max()); m['t_late_peak'] = float(t[late][P[late].argmax()])
        m['late_over_exp'] = m['late_peak'] / P_EXP_HASEG
        m['early_peak'] = float(P[t <= 0.5].max())
    elif spec['motor'] == 'haseg':
        m['early_peak'] = float(P.max())
    # fill Mach (thread-F confounder): max |Mach| over the bore at fixed times, + at each first-lit cell
    mm = np.asarray(res['max_mach'])
    m['max_mach_at_ms'] = {str(ms): float(mm[min(int(np.searchsorted(t, ms * 1e-3)), mm.size - 1)])
                           for ms in (1, 2, 5, 10, 20)}
    m['max_mach_first_50ms'] = float(mm[t <= 0.05].max()) if (t <= 0.05).any() else None
    # ignition-event diagnostics per grain cell (NaN -> None): h, T_drv, q_pre, T_s, T_gas, u, Mach, P, t, q_post
    m['diag'] = [[None if not np.isfinite(v) else float(v) for v in row] for row in diag]
    m['ign'] = [float(v) for v in ig]
    # first-lit cell summary
    if lit.size:
        j = int(np.argmin(ig))
        m['first_cell_diag'] = dict(zip(['h', 'T_drv', 'q_pre', 'T_s', 'T_gas', 'u', 'Mach', 'P', 't', 'q_post'],
                                        [None if not np.isfinite(v) else float(v) for v in diag[j]]))
    m['warnings'] = sorted({str(w.message)[:160] for w in wl if 'never' in str(w.message).lower()
                            or 'charge' in str(w.message).lower()})[:3]
    fo = res.get('forcing')
    if fo is not None and spec.get('forcing_dir', 'C:/tmp/sv/sm2_data'):
        d = Path(spec.get('forcing_dir', 'C:/tmp/sv/sm2_data')); d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / f"{spec['tag']}_forcing.npz", h=fo["h"][:, isg],
                            Td=fo["T_drv"][:, isg], t=fo['t'], dt=fo['dt'],
                            rho=float(prop.rho_propellant), ign=ig)
        m['forcing_steps'] = int(fo['t'].size)
        m['forcing_t_end'] = float(fo['t'][-1] + fo['dt'][-1])
    m['rho_propellant'] = float(prop.rho_propellant)
    if return_res:
        return m, res
    return m


if __name__ == '__main__':
    spec = json.loads(sys.argv[1]); out = Path(sys.argv[2])
    if _done(out, spec['tag']):
        print('skip', spec['tag']); sys.exit(0)
    m = run_case(spec)
    with open(out, 'a') as f:
        f.write(json.dumps(m) + '\n')
    print(spec['tag'], 'wall %.1f' % m['wall_s'], 'first', m['t_first'], 'peak', m['peak'])
