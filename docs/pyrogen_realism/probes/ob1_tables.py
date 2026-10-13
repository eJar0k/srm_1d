"""Build docs/pyrogen_realism/results/ob1r_tables.md from the post-rebase OB1 JSONL results (ob1r_*.jsonl)."""
import json
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parents[1] / 'results'


def load(name):
    return [json.loads(l) for l in open(R / name) if l.strip()]


def f(x, p=2):
    return '-' if x is None else f'{x:.{p}f}'


def ms(r, k):
    return None if r.get(k) is None else r[k] * 1e3


def lit(r):
    return r.get('frac_ignited', 0) >= 0.99


out = []

# T1
out.append('## T1 rate laws (mm/s)\n')
out.append('| P [MPa] | sim law 2e-5 P^0.5 | Nakano 71.1 P^0.589 (valid 0.005-0.2 MPa) | Zhu mean r0 (valid 2.8-50 MPa) | MIL-P-46994B (closed-bomb average) |')
out.append('|---|---|---|---|---|')
r0 = np.mean([6.848, 6.787, 5.389, 5.271]) * 1e-3
for P in (0.1, 0.2, 1.0, 2.8, 7.0, 20.0):
    sim = 2e-5 * (P * 1e6) ** 0.5 * 1e3
    nak = 71.1 * P ** 0.589
    zhu = r0 * min(P, 7.0) * 1e3 if P >= 2.8 else None
    out.append(f'| {P:g} | {sim:.1f} | {(f(nak,1) if P <= 0.2 else "(%.0f, out of range)" % nak)} | '
               f'{(f(zhu,1) if zhu else "(below range)")} | {"31.8-44.5" if P >= 2.8 else ""} |')
out.append('')

MODEL_DESC = {
    'pre_ob1': 'pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law',
    'lib_constant': 'I-A pellet, constant-area basket law (opt-in)',
    'lib': '**new default**: I-A pellet, exact cylinder area law',
    'lib_capZ': 'new default + rate capped at the Zhu plateau (opt-in)',
    'cube3_constant': '3 mm cube, constant-area law',
    'cube3': '3 mm cube, exact area law (L/D 1)',
    'cube3_capZ': '3 mm cube, exact area law + Zhu cap',
}
MODEL_ORDER = ['pre_ob1', 'lib_constant', 'lib', 'lib_capZ', 'cube3_constant', 'cube3', 'cube3_capZ']

# T2 realistic
re_ = {r['tag']: r for r in load('ob1r_realistic.jsonl')}
out.append('## T3 realistic rows (Chunc head basket, t_max 0.2 s, plateau = median 0.1-0.2 s; Peretz ignition h = main default)\n')
out.append('| pyrogen | g | DeMar | model | first ign ms | t90 ms | peak MPa | t_peak ms | plateau | peak/plateau | pyro t50 | t90 | duration ms |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for pyro in ('bpnv', 'mtv'):
    for g in (4.5, 6.0):
        for d in ('on', 'off'):
            for m in MODEL_ORDER:
                r = re_.get(f're_{pyro}_{g}_{d}_{m}')
                if r is None:
                    continue
                if 'error' in r:
                    out.append(f'| {pyro} | {g} | {d} | {MODEL_DESC[m]} | ERROR {r["error"][:50]} |')
                    continue
                if not lit(r):
                    out.append(f'| {pyro} | {g} | {d} | {MODEL_DESC[m]} | NOT LIT (peak {r["peak_MPa"]:.2f}) | | | | | | | | |')
                    continue
                out.append(f"| {pyro} | {g} | {d} | {MODEL_DESC[m]} | {f(ms(r,'t_first_ign_s'))} | {f(ms(r,'t90_s'))} | {r['peak_MPa']:.2f} | {f(ms(r,'t_peak_s'),1)} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} | {f(ms(r,'pyro_t50_s'),1)} | {f(ms(r,'pyro_t90_s'),1)} | {f(ms(r,'pyro_duration_s'),1)} |")
out.append('')

# T4 long
lg = {r['tag']: r for r in load('ob1r_long.jsonl')}
out.append('## T4 t_max 1.0 s, plateau = median 0.4-0.9 s (Chunc experiment: peak 8.89, plateau 8.66, ratio 1.027)\n')
out.append('| pyrogen | g | DeMar | model | peak MPa | t_peak ms | plateau | peak/plateau | pyro duration ms |')
out.append('|---|---|---|---|---|---|---|---|---|')
for pyro in ('bpnv', 'mtv'):
    for g in (4.5, 6.0):
        for d in ('on', 'off'):
            for m in ('pre_ob1', 'lib', 'cube3'):
                r = lg.get(f'lg_{pyro}_{g}_{d}_{m}')
                if r is None or 'error' in r:
                    continue
                if not lit(r):
                    out.append(f'| {pyro} | {g} | {d} | {MODEL_DESC[m]} | NOT LIT | | | | |')
                    continue
                out.append(f"| {pyro} | {g} | {d} | {MODEL_DESC[m]} | {r['peak_MPa']:.2f} | {f(ms(r,'t_peak_s'),1)} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} | {f(ms(r,'pyro_duration_s'),1)} |")
out.append('')

# T5 thresholds
th = {}
for r in load('ob1r_thresh.jsonl'):
    if r['tag'].startswith('thr_'):
        th[r['tag']] = r          # last summary row wins (finer plenum grid re-run)
out.append('## T5 ignition thresholds (min BPNV charge with >= 99 % of the 100 grain cells lit within 0.2 s; grid scan then bisection to 2 %)\n')
out.append('| topology | model | DeMar on [g] | DeMar off [g] |')
out.append('|---|---|---|---|')
for topo, models in (('head_basket', ('pre_ob1', 'lib_constant', 'lib', 'cube3')), ('forward_plenum', ('pre_ob1', 'lib', 'cube3'))):
    for m in models:
        def cell(d):
            r = th.get(f'thr_{topo}_{d}_{m}')
            if r is None:
                return '-'
            if r['threshold_g'] is None:
                return 'not lit (grid top)'
            return f"{r['threshold_g']:.3g}"
        out.append(f'| {topo} | {MODEL_DESC[m]} | {cell("on")} | {cell("off")} |')
out.append('')

# T6 mass sweep
mm = {r['tag']: r for r in load('ob1r_mass.jsonl')}
gs = (0.6, 0.9, 1.36, 1.8, 2.27, 2.93, 3.6, 4.54, 6.3, 9.0)
out.append('## T6 peak vs charge (BPNV head basket; peak MPa (peak/plateau); NL = not lit within 0.2 s)\n')
out.append('| DeMar | model | ' + ' | '.join(f'{g} g' for g in gs) + ' |')
out.append('|---|---|' + '---|' * len(gs))
for d in ('on', 'off'):
    for m in ('pre_ob1', 'lib', 'cube3'):
        cells = []
        for g in gs:
            r = mm.get(f'ms_{d}_{m}_{g}')
            cells.append('?' if r is None else (f"{r['peak_MPa']:.1f} ({r['peak_over_plateau']:.2f})" if lit(r) else 'NL'))
        out.append(f'| {d} | {m} | ' + ' | '.join(cells) + ' |')
out.append('')

# T7 canonical
vv = {r['tag']: r for r in load('ob1r_valid.jsonl')}
out.append('## T7 canonical motors, library BPNV, Sutton-mass forward plenum, Peretz ignition h (main default); old = pre-OB1 library pellet (3.2 mm, L/D 1), new = I-A pellet (3.175 mm, L/D 1.504)\n')
out.append('| motor | pellet | DeMar | rate law | lit | first ign ms | peak MPa | t_peak | plateau / trough / mid | ratios |')
out.append('|---|---|---|---|---|---|---|---|---|---|')
for m in ('chunc', 'balls', 'haseg'):
    for key in [('old32', 'on', 'lib'), ('lib', 'on', 'lib'), ('lib', 'on', 'capZ_lib'), ('old32', 'off', 'lib'), ('lib', 'off', 'lib')]:
        r = vv.get(f'{m}_{key[0]}_{key[1]}_{key[2]}')
        if r is None:
            continue
        pel = 'old 3.2 mm L/D 1' if key[0] == 'old32' else 'new I-A'
        law = 'library' if key[2] == 'lib' else 'capped (opt-in)'
        if not r['lit']:
            out.append(f"| {m} | {pel} | {key[1]} | {law} | NO | - | {r['peak']:.2f} | | | |")
            continue
        extra = ratios = ''
        tp = f"{r['t_peak']*1e3:.1f} ms" if m != 'haseg' else f"{r['t_peak']:.3f} s"
        if m == 'chunc':
            extra = f"{r['plateau']:.3f}"
            ratios = f"pk/pl {r['ratio']:.3f}"
        elif m == 'balls':
            extra = f"{r['trough']:.3f} / {r['mid']:.3f}"
            ratios = f"{r['pk_tr']:.2f} / {r['pk_mid']:.2f}"
        out.append(f"| {m} | {pel} | {key[1]} | {law} | yes | {r['first_ign_s']*1e3:.2f} | {r['peak']:.4f} | {tp} | {extra} | {ratios} |")
out.append('')

# T8 attribution / ablation on the new default
ex = {r['tag']: r for r in load('ob1r_extra.jsonl')}
out.append('## T8 attribution and ablation (BPNV 4.5 g head basket, DeMar on, t_max 0.2 s)\n')
out.append('| model | variant | first ign ms | peak MPa | plateau (0.1-0.2 s) | peak/plateau | pyro duration ms |')
out.append('|---|---|---|---|---|---|---|')
for m in ('pre_ob1', 'lib', 'lib_capZ'):
    base = re_[f're_bpnv_4.5_on_{m}']
    out.append(f"| {m} | Ma 2020 erosive (default) | {f(ms(base,'t_first_ign_s'))} | {base['peak_MPa']:.2f} | {base['plateau_MPa']:.2f} | {base['peak_over_plateau']:.2f} | {f(ms(base,'pyro_duration_s'),1)} |")
    for v in ('saint_robert', 'mukunda_paul'):
        r = ex[f'ex_{m}_{v}']
        out.append(f"| {m} | erosive = {v} | {f(ms(r,'t_first_ign_s'))} | {r['peak_MPa']:.2f} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} | {f(ms(r,'pyro_duration_s'),1)} |")
    r = ex[f'ex_{m}_pyro_only']
    out.append(f"| {m} | grain cannot ignite (pyrogen gas only) | - | {r['peak_MPa']:.2f} | - | - | - |")
for tag, desc in (('ex_lib_cap7_lib', 'rate flat above 7 MPa (library law)'), ('ex_lib_zhu_lo', 'Zhu law, r0 = 5.271e-3'),
                  ('ex_lib_zhu_mean', 'Zhu law, r0 = 6.074e-3 (mean)'), ('ex_lib_zhu_hi', 'Zhu law, r0 = 6.848e-3'),
                  ('ex_lib_bed3', 'diagnostic bed flame speed 3 m/s'), ('ex_lib_bed10', 'diagnostic bed flame speed 10 m/s'),
                  ('ex_lib_bed30', 'diagnostic bed flame speed 30 m/s')):
    r = ex[tag]
    if not lit(r):
        out.append(f'| lib | {desc} | NOT LIT | {r["peak_MPa"]:.2f} | | | |')
        continue
    out.append(f"| lib | {desc} | {f(ms(r,'t_first_ign_s'))} | {r['peak_MPa']:.2f} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} | {f(ms(r,'pyro_duration_s'),1)} |")
out.append('')

(R / 'ob1r_tables.md').write_text('\n'.join(out), encoding='utf-8')
print('wrote', R / 'ob1r_tables.md', len(out), 'lines')
