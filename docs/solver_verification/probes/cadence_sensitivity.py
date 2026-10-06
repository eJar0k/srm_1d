"""T7 -- timestep / refresh-cadence sensitivity sweep (measurement only).

3 motors x 6 configs = 18 rows. Config = (burn_update_interval,
geometry_update_interval, cfl_target); None = solver default (burn
max(10, N//5), geometry follows burn, cfl 0.3). Canonical knobs, ma2020,
cap 0, library BPNV, committed .ric files of THIS checkout.

Metric definitions follow docs/v0_7_4/probes/mp_cap_crossmotor.py
(Chunc peak/t_peak/plateau[0.4,0.9] median/ratio; BALLSstick peak in
t<0.3, trough = min over (t_pk,1.2), mid = P at t nearest 1.03,
pk/trough, pk/mid; Hasegawa peak/t_peak). Extra: global min of the
min_pressure history + time, summary pressure_floor_min_raw_P, fill Mach
= max(max_mach[t < t_peak]), T3 floor counters, T5 refresh gaps,
n_steps, wall, steps/s, plateau-window stationarity (linear-fit slope as
%/s of window mean, std/mean).

Resumable: each finished row is appended to a JSONL; finished
(motor, config) pairs are skipped.
Usage (repo root, PYTHONPATH=$PWD, SRM1D_OPENMOTOR_PATH set):
  python docs/solver_verification/probes/cadence_sensitivity.py
      [--motor chunc|balls|haseg] [--config NAME] [--out FILE.jsonl]
  python ... --report [--out FILE.jsonl]     # print the table/deltas
"""
import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

CANON = dict(roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None,
             T_ignition=756.0, print_interval=0.0, verbose=False)

CONFIGS = {  # name: (burn, geometry, cfl)
    'base':       (None, None, 0.3),
    'b1':         (1, 'default', 0.3),  # geometry pinned to the resolved default (None would FOLLOW burn=1)
    'g1':         (None, 1, 0.3),
    'b1g1':       (1, 1, 0.3),
    'cfl015':     (None, None, 0.15),
    'b1g1cfl015': (1, 1, 0.15),
}
MOTORS = {  # name: (ric, t_max, P_cutoff)
    'chunc': ('machbusterNew.ric', 1.0, 0.01e6),
    'balls': ('BALLSstick.ric', 2.0, 0.01e6),
    'haseg': ('hasegawa_a.ric', 10.0, 0.05e6),
}
DEFAULT_OUT = Path('C:/tmp/sv_t7/results.jsonl')


def stationarity(t, P, lo, hi):
    m = (t >= lo) & (t <= hi)
    if m.sum() < 3:
        return dict(n=int(m.sum()), slope_pct_per_s=None, std_over_mean=None)
    tt, pp = t[m], P[m]
    mean = float(pp.mean())
    slope = float(np.polyfit(tt, pp, 1)[0])
    return dict(n=int(m.sum()), mean=mean, slope_pct_per_s=100 * slope / mean,
                std_over_mean=float(pp.std() / mean))


def run_row(motor, cname, out=DEFAULT_OUT):
    from srm_1d.openmotor_adapter import run_from_ric
    ric, tmax, pcut = MOTORS[motor]
    b, g, cfl = CONFIGS[cname]
    kw = dict(CANON, t_max=tmax, P_cutoff=pcut, cfl_target=cfl,
              erosive_model='ma2020', erosive_mp_cap=0.0)
    if b is not None:
        kw['burn_update_interval'] = b
    if g == 'default':
        # geometry_update_interval=None FOLLOWS an explicit burn interval, so
        # 'burn every step, geometry default' needs the resolved default
        # (max(10, N//5)) passed explicitly; take it from the base row.
        base = [r for r in load(out) if r['motor'] == motor and r['config'] == 'base']
        if not base:
            raise RuntimeError(f"run {motor}/base first (need its resolved geometry interval)")
        g = base[0]['geom_int']
    if g is not None:
        kw['geometry_update_interval'] = g
    t0 = time.time()
    res = run_from_ric(str(ROOT / 'motors' / ric), **kw)[0]
    wall = time.time() - t0
    s = res['summary']
    t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6
    k = int(P.argmax())
    row = dict(motor=motor, config=cname, burn=b, geom=g, cfl=cfl,
               peak=float(P[k]), t_peak=float(t[k]))
    if motor == 'chunc':
        w = (t >= 0.4) & (t <= 0.9)
        pl = float(np.median(P[w]))
        row.update(plateau=pl, ratio=float(P[k]) / pl,
                   stat_plateau_0p4_0p9=stationarity(t, P, 0.4, 0.9))
    elif motor == 'balls':
        ms = t < 0.3
        kp = int(np.argmax(P[ms])); Ppk, tpk = float(P[ms][kp]), float(t[ms][kp])
        after = (t > tpk) & (t < 1.2)
        Ptr = float(np.min(P[after]))
        Pmid = float(P[int(np.argmin(np.abs(t - 1.03)))])
        row.update(peak=Ppk, t_peak=tpk, trough=Ptr, mid=Pmid,
                   pk_tr=Ppk / Ptr, pk_mid=Ppk / Pmid,
                   stat_mid_0p9_1p2=stationarity(t, P, 0.9, 1.2))
    tp = row['t_peak']
    mp = np.asarray(res['min_pressure'])
    j = int(mp.argmin())
    mm = np.asarray(res['max_mach'])
    pre = t < tp
    row.update(min_P_Pa=float(mp[j]), t_min_P=float(t[j]),
               floor_min_raw_P=float(s['pressure_floor_min_raw_P']),
               fill_mach=float(mm[pre].max()) if pre.any() else None,
               floor_activations=int(s['pressure_floor_activations']),
               floor_steps=int(s['pressure_floor_steps']),
               floor_mass_kg=float(s['pressure_floor_mass_added_kg']),
               floor_energy_J=float(s['pressure_floor_energy_added_J']),
               floor_first_t=float(s['pressure_floor_first_time_s']),
               max_burn_gap_us=1e6 * float(s['max_burn_refresh_gap_s']),
               max_geom_gap_us=1e6 * float(s['max_geometry_refresh_gap_s']),
               burn_int=int(s['burn_refresh_interval_steps']),
               geom_int=int(s['geometry_refresh_interval_steps']),
               n_steps=int(s['steps']), wall_s=wall,
               steps_per_s=int(s['steps']) / wall,
               termination=s.get('termination'))
    return row


def load(out):
    rows = []
    if Path(out).exists():
        for line in open(out):
            if line.strip():
                rows.append(json.loads(line))
    return rows


METRICS = {
    'chunc': ['peak', 't_peak', 'plateau', 'ratio', 'min_P_Pa', 'fill_mach', 'n_steps'],
    'balls': ['peak', 't_peak', 'trough', 'mid', 'pk_tr', 'pk_mid', 'min_P_Pa', 'fill_mach', 'n_steps'],
    'haseg': ['peak', 't_peak', 'min_P_Pa', 'fill_mach', 'n_steps'],
}
VALIDATION = {'chunc': ['peak', 'plateau', 'ratio'],
              'balls': ['peak', 'pk_tr', 'pk_mid'], 'haseg': ['peak']}


def report(out):
    rows = load(out)
    for motor in MOTORS:
        mr = {r['config']: r for r in rows if r['motor'] == motor}
        if 'base' not in mr:
            continue
        print(f"\n## {motor}\n")
        ms = METRICS[motor]
        print('| config | ' + ' | '.join(ms) + ' | minraw P | wall s | steps/s |')
        print('|' + '---|' * (len(ms) + 4))
        for c in CONFIGS:
            if c in mr:
                r = mr[c]
                print(f"| {c} | " + ' | '.join(f"{r[m]:.6g}" for m in ms)
                      + f" | {r['floor_min_raw_P']:.6g} | {r['wall_s']:.0f} | {r['steps_per_s']:.0f} |")
        print('\nrelative change vs base (%):\n')
        print('| config | ' + ' | '.join(ms) + ' |')
        print('|' + '---|' * (len(ms) + 1))
        flags = []
        for c in CONFIGS:
            if c in mr and c != 'base':
                cells = []
                for m in ms:
                    d = 100 * (mr[c][m] / mr['base'][m] - 1)
                    cells.append(f"{d:+.3f}")
                    if m in VALIDATION[motor] and abs(d) > 1.0:
                        flags.append(f"{motor} {c} {m} {d:+.3f}%")
                print(f"| {c} | " + ' | '.join(cells) + ' |')
        print('\nFLAGGED >1%: ' + ('; '.join(flags) if flags else 'none'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--motor')
    ap.add_argument('--config')
    ap.add_argument('--out', default=str(DEFAULT_OUT))
    ap.add_argument('--report', action='store_true')
    a = ap.parse_args()
    if a.report:
        return report(a.out)
    import srm_1d
    print(f"[provenance] srm_1d from {Path(srm_1d.__file__).parent}", flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    done = {(r['motor'], r['config']) for r in load(a.out)}
    for motor in ([a.motor] if a.motor else MOTORS):
        for c in ([a.config] if a.config else CONFIGS):
            if (motor, c) in done:
                print(f"skip {motor}/{c} (done)", flush=True)
                continue
            r = run_row(motor, c, a.out)
            with open(a.out, 'a') as f:
                f.write(json.dumps(r) + '\n')
            print(f"{motor}/{c}: peak {r['peak']:.5f} minP {r['min_P_Pa']:.1f} "
                  f"steps {r['n_steps']} wall {r['wall_s']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
