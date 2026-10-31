"""GS1: one (motor, charge, topology, grid) row. Appends a JSON line to the out file.

Usage (from the worktree root, PYTHONPATH=$PWD, under the guard):
  python -m srm_1d.tools.guarded_run --mem-gb 4 --timeout-min 60 -- \
      python docs/pyrogen_realism/probes/gs1_row.py '<json spec>' <out.jsonl>
spec keys: tag, motor (chunc|balls|haseg), g (charge, grams), topo (forward_plenum|head_basket), N (target
propellant cells), t_max, plateau_win=[t0,t1], extra (dict of run kwargs, optional).
Metrics: ignition (first light, t50/t90 of the grain, front speed, dark-zone end POSITION), startup peak /
t_peak / plateau, and the contact-front width vs time (10-90 % of the gas-T rise along the bore, from the
snapshots taken before first light).  Explicit charge: the ignition fallback never resizes it.
"""
import json
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np

from srm_1d.openmotor_adapter import run_from_ric

RIC = dict(chunc='motors/machbusterNew.ric', balls='motors/BALLSstick.ric', haseg='motors/hasegawa_a.ric')
SNAP = 2.0e-4


def _front_width(x, T, T_amb=293.0):
    """10-90 % width [m] of the monotone-from-the-head gas-T rise, plus x(90%) and x(10%)."""
    Tmax = float(T.max())
    if Tmax < T_amb + 50.0:
        return None
    th = (T - T_amb) / (Tmax - T_amb)
    i90 = np.argmax(th < 0.9)
    if th[i90] >= 0.9:
        return None
    # first crossing below 0.1 AFTER the 0.9 crossing
    rest = th[i90:]
    j = np.argmax(rest < 0.1)
    if rest[j] >= 0.1:
        return None
    i10 = i90 + j

    def cross(i, level):
        if i == 0:
            return float(x[0])
        a, b = th[i - 1], th[i]
        f = (a - level) / (a - b) if a != b else 0.0
        return float(x[i - 1] + f * (x[i] - x[i - 1]))
    x90, x10 = cross(i90, 0.9), cross(i10, 0.1)
    return x10 - x90, x90, x10


def run(spec):
    kw = dict(roughness=32e-6, kappa=0.44, T_ignition=756.0, print_interval=0.0, verbose=False,
              t_max=spec['t_max'], snapshot_interval=SNAP, injection_topology=spec['topo'],
              target_propellant_cells=int(spec['N']), pyrogen='bpnv', pyrogen_mass=spec['g'] * 1e-3,
              P_cutoff=spec.get('P_cutoff', 0.05e6))
    kw.update(spec.get('extra', {}))
    t0 = time.time()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        res = run_from_ric(str(ROOT / RIC[spec['motor']]), **kw)[0]
    wall = time.time() - t0
    s = res['summary']
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    snaps = res['snapshots']
    isg = np.asarray(snaps[0]['is_grain'])
    x_all = np.asarray(snaps[0]['x'])
    ign_all = np.asarray(res['ignition_time_by_cell'])
    ign = ign_all[isg]
    xg = x_all[isg]
    lit = ign < 1e8
    n_g = int(isg.sum())
    out = dict(spec=spec, wall_s=wall, steps=int(s['steps']), N_cells=int(len(x_all)), n_grain=n_g,
               dx_mm=float(res['dx']) * 1e3, termination=s['termination'], ignited=bool(s['ignited']),
               n_lit=int(lit.sum()), frac_lit=float(lit.mean()))
    if lit.any():
        tl = np.sort(ign[lit])
        out['t_first_ms'] = float(tl[0] * 1e3)
        for p in (50, 90):
            k = int(np.ceil(p / 100.0 * n_g)) - 1
            out[f't{p}_ms'] = float(tl[k] * 1e3) if k < tl.size else None
        order = np.argsort(ign[lit])
        xs, ts = xg[lit][order], ign[lit][order]
        # ignition-front speed: slope of x vs t over each quartile of the lit cells (by order of lighting)
        q = np.array_split(np.arange(xs.size), 4)
        sp = []
        for idx in q:
            if idx.size >= 3 and np.ptp(ts[idx]) > 0:
                sp.append(float(np.polyfit(ts[idx], xs[idx], 1)[0]))
            else:
                sp.append(None)
        out['front_speed_quartiles_m_s'] = sp
        k10, k90 = int(0.1 * (xs.size - 1)), int(0.9 * (xs.size - 1))
        out['front_speed_10_90_m_s'] = float((xs[k90] - xs[k10]) / max(ts[k90] - ts[k10], 1e-12))
        out['lit_x_mm'] = [float(xg[lit].min() * 1e3), float(xg[lit].max() * 1e3)]
    dark = ~lit
    out['n_dark'] = int(dark.sum())
    out['dark_x_mm'] = [float(xg[dark].min() * 1e3), float(xg[dark].max() * 1e3)] if dark.any() else None
    out['dark_end_mm'] = float(xg[dark].max() * 1e3) if dark.any() else None
    k = int(P.argmax())
    pw = spec.get('plateau_win', [0.6 * spec['t_max'], spec['t_max']])
    m = (t >= pw[0]) & (t <= pw[1])
    out.update(peak_MPa=float(P[k]), t_peak_ms=float(t[k] * 1e3),
               plateau_MPa=float(np.median(P[m])) if m.any() else None)
    if out['plateau_MPa']:
        out['peak_over_plateau'] = out['peak_MPa'] / out['plateau_MPa']
    # contact front: width of the gas-T rise along the bore before the first light
    cf = []
    t_first = out.get('t_first_ms', 1e9) * 1e-3
    for sn in snaps:
        if sn['t'] <= 0.0 or sn['t'] > t_first + 2e-3:
            continue
        r = _front_width(x_all, np.asarray(sn['T']))
        if r:
            cf.append([float(sn['t'] * 1e3), r[0] * 1e3, r[1] * 1e3, r[2] * 1e3, float(np.asarray(sn['T']).max())])
    out['contact_front'] = cf          # [t_ms, width_mm, x90_mm, x10_mm, Tmax]
    out['warn'] = [str(x.message)[:100] for x in w if issubclass(x.category, RuntimeWarning)][:3]
    if spec.get('save_ign'):
        out['ign_ms'] = [float(v * 1e3) if v < 1e8 else None for v in ign]
        out['x_grain_mm'] = [float(v * 1e3) for v in xg]
    return out


if __name__ == '__main__':
    spec = json.loads(sys.argv[1])
    out = Path(sys.argv[2])
    if out.exists() and any(json.loads(l)['spec']['tag'] == spec['tag'] for l in out.read_text().splitlines() if l.strip()):
        print('skip', spec['tag'])
        sys.exit(0)
    row = run(spec)
    with open(out, 'a') as f:
        f.write(json.dumps(row) + '\n')
    print(json.dumps({k: row[k] for k in ('wall_s', 'steps', 'n_lit', 'n_dark', 't_first_ms', 'peak_MPa') if k in row}))
