"""Build docs/pyrogen_realism/results/ob1_tables.md from the OB1 JSONL results. Run from the worktree root."""
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


OPT_ORDER = ['base', 'area23', 'cap7', 'capZ', 'zhu', 'zhu+area23', 'capZ+area23']
OPT_DESC = {
    'base': 'current pyrogen model (library law, constant area)',
    'area23': '(b) shrinking area A0 (m/m0)^(2/3), library law',
    'cap7': 'library law, rate flat above 7 MPa (IG1 probe)',
    'capZ': '(a) library law capped at Zhu plateau rate 42.5 mm/s (P_cap 4.5 MPa)',
    'zhu': '(a) Zhu piecewise law r = r0 P, flat above 7 MPa, constant area',
    'zhu+area23': '(a)+(b) Zhu law + shrinking area',
    'capZ+area23': '(a)+(b) capped library law + shrinking area',
}


def parse_tag(tag, prefix):
    rest = tag[len(prefix):]
    return rest


out = []

# T1 rate laws
out.append('## T1 rate laws (mm/s)\n')
out.append('| P [MPa] | sim law 2e-5 P^0.5 | Nakano 71.1 P^0.589 (valid 0.005-0.2 MPa) | YAML "candidate" 3.5e-5 P^0.589 | Zhu mean r0 (valid 2.8-50 MPa) | MIL-P-46994B (closed-bomb average) |')
out.append('|---|---|---|---|---|---|')
r0 = np.mean([6.848, 6.787, 5.389, 5.271]) * 1e-3
for P in (0.1, 0.2, 1.0, 2.8, 7.0, 20.0):
    sim = 2e-5 * (P * 1e6) ** 0.5 * 1e3
    nak = 71.1 * P ** 0.589
    yam = 3.5e-5 * (P * 1e6) ** 0.589 * 1e3
    zhu = r0 * min(P, 7.0) * 1e3 if P >= 2.8 else None
    nk = f(nak, 1) if P <= 0.2 else f'({nak:.0f}, out of range)'
    yk = f(yam, 1)
    zk = f(zhu, 1) if zhu else '(below range)'
    out.append(f'| {P:g} | {sim:.1f} | {nk} | {yk} | {zk} | {"31.8-44.5" if P >= 2.8 else ""} |')
out.append('')

# T2 ranking
rk = {r['tag']: r for r in load('ob1_rank.jsonl')}
out.append('## T2 ranking probes (BPNV 4.5 g, 3 mm cubes, head basket, t_max 0.2 s, plateau = median 0.1-0.2 s)\n')
for combo in ('current', 'forward'):
    out.append(f'### combo = {combo}\n')
    out.append('| option | first ign ms | t90 ms | peak MPa | d peak vs base | t_peak ms | peak/plateau | pyro t50 ms | pyro t90 ms | pyro duration ms |')
    out.append('|---|---|---|---|---|---|---|---|---|---|')
    base = rk[f'rk_{combo}_base']['peak_MPa']
    names = OPT_ORDER + ['area23+bed3', 'area23+bed10', 'area23+bed30']
    for n in names:
        r = rk[f'rk_{combo}_{n}']
        d = r['peak_MPa'] / base - 1
        lit = r['frac_ignited'] >= 0.99
        desc = OPT_DESC[n] if n in OPT_DESC else f'(b) + DIAGNOSTIC bed flame speed {n.split("bed")[1]} m/s (no sourced value)'
        if lit:
            out.append(f"| {n}: {desc} | {f(ms(r,'t_first_ign_s'))} | {f(ms(r,'t90_s'))} | {r['peak_MPa']:.2f} | {d*100:+.1f} % | {f(ms(r,'t_peak_s'),1)} | {r['peak_over_plateau']:.2f} | {f(ms(r,'pyro_t50_s'),1)} | {f(ms(r,'pyro_t90_s'),1)} | {f(ms(r,'pyro_duration_s'),1)} |")
        else:
            out.append(f"| {n}: {desc} | NOT LIT | - | {r['peak_MPa']:.2f} | - | - | - | - | - | - |")
    out.append('')

# T3 realistic
re_ = {r['tag']: r for r in load('ob1_realistic.jsonl')}
out.append('## T3 realistic rows (Chunc head basket, 3 mm cubes L/D 1, t_max 0.2 s, plateau = median 0.1-0.2 s)\n')
out.append('| pyrogen | g | combo | option | first ign ms | t90 ms | peak MPa | t_peak ms | plateau | peak/plateau | pyro t50 | t90 | duration ms |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for pyro in ('bpnv', 'mtv'):
    for g in (4.5, 6.0):
        for combo in ('current', 'forward'):
            for n in OPT_ORDER:
                r = re_.get(f're_{pyro}_{g}_{combo}_{n}')
                if r is None:
                    continue
                if r['frac_ignited'] < 0.99:
                    out.append(f'| {pyro} | {g} | {combo} | {n} | NOT LIT (peak {r["peak_MPa"]:.2f}) | | | | | | | | |')
                    continue
                out.append(f"| {pyro} | {g} | {combo} | {n} | {f(ms(r,'t_first_ign_s'))} | {f(ms(r,'t90_s'))} | {r['peak_MPa']:.2f} | {f(ms(r,'t_peak_s'),1)} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} | {f(ms(r,'pyro_t50_s'),1)} | {f(ms(r,'pyro_t90_s'),1)} | {f(ms(r,'pyro_duration_s'),1)} |")
out.append('')

# T4 long
lg = {r['tag']: r for r in load('ob1_long.jsonl')}
out.append('## T4 t_max 1.0 s, plateau = median 0.4-0.9 s (Chunc experiment: peak 8.89, plateau 8.66, ratio 1.027)\n')
out.append('| pyrogen | g | combo | option | peak MPa | t_peak ms | plateau | peak/plateau | pyro duration ms |')
out.append('|---|---|---|---|---|---|---|---|---|')
for pyro in ('bpnv', 'mtv'):
    for g in (4.5, 6.0):
        for combo in ('current', 'forward'):
            for n in ('base', 'area23', 'zhu+area23', 'capZ+area23'):
                r = lg.get(f'lg_{pyro}_{g}_{combo}_{n}')
                if r is None:
                    continue
                if r['frac_ignited'] < 0.99:
                    out.append(f'| {pyro} | {g} | {combo} | {n} | NOT LIT | | | | |')
                    continue
                out.append(f"| {pyro} | {g} | {combo} | {n} | {r['peak_MPa']:.2f} | {f(ms(r,'t_peak_s'),1)} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} | {f(ms(r,'pyro_duration_s'),1)} |")
out.append('')

# T5 thresholds
th = [r for r in load('ob1_thresh.jsonl') if r['tag'].startswith('thr_')]
out.append('## T5 head-basket ignition threshold (min BPNV charge with >= 99 % of 100 grain cells lit within 0.2 s; grid then bisection to 2 %)\n')
out.append('| option | current combo [g] | forward combo [g] |')
out.append('|---|---|---|')
tt = {(r['combo'], r['opt']): r for r in th}
for n in ('base', 'area23', 'capZ+area23', 'zhu+area23'):
    def cell(c):
        r = tt.get((c, n))
        return '-' if r is None or r['threshold_g'] is None else f"{r['threshold_g']:.3g}"
    out.append(f'| {n} | {cell("current")} | {cell("forward")} |')
out.append('')

# T6 mass sweep
ms_ = load('ob1_mass.jsonl')
out.append('## T6 peak vs charge (BPNV head basket; peak MPa (peak/plateau); NL = not lit within 0.2 s)\n')
gs = (0.6, 0.9, 1.36, 1.8, 2.27, 2.93, 3.6, 4.54, 6.3, 9.0)
out.append('| combo | option | ' + ' | '.join(f'{g} g' for g in gs) + ' |')
out.append('|---|---|' + '---|' * len(gs))
mm = {r['tag']: r for r in ms_}
for combo in ('current', 'forward'):
    for n in ('base', 'area23', 'capZ+area23'):
        cells = []
        for g in gs:
            r = mm[f'ms_{combo}_{n}_{g}']
            cells.append(f"{r['peak_MPa']:.1f} ({r['peak_over_plateau']:.2f})" if r['frac_ignited'] >= 0.99 else 'NL')
        out.append(f'| {combo} | {n} | ' + ' | '.join(cells) + ' |')
out.append('')

# T7 canonical
vv = {r['tag']: r for r in load('ob1_valid.jsonl')}
out.append('## T7 canonical motors, library BPNV, Sutton-mass forward plenum, current ignition combo (law = rate law; area law not applicable to the plenum)\n')
out.append('| motor | law | first ign ms | peak MPa | d peak vs lib | t_peak | plateau / trough / mid | ratios |')
out.append('|---|---|---|---|---|---|---|---|')
for m in ('chunc', 'balls', 'haseg'):
    b = vv[f'{m}_lib']['peak']
    for l in ('lib', 'cap7_lib', 'capZ_lib', 'zhu_mean'):
        r = vv[f'{m}_{l}']
        extra = ''
        ratios = ''
        if m == 'chunc':
            extra = f"{r['plateau']:.3f}"
            ratios = f"pk/pl {r['ratio']:.3f}"
        elif m == 'balls':
            extra = f"{r['trough']:.3f} / {r['mid']:.3f}"
            ratios = f"{r['pk_tr']:.2f} / {r['pk_mid']:.2f}"
        out.append(f"| {m} | {l} | {r['first_ign_s']*1e3:.2f} | {r['peak']:.4f} | {(r['peak']/b-1)*100:+.2f} % | {r['t_peak']*1e3 if m!='haseg' else r['t_peak']*1e3:.1f} ms | {extra} | {ratios} |")
out.append('')

# T8 attribution
ex = {r['tag']: r for r in load('ob1_extra.jsonl')}
out.append('## T8 attribution (BPNV 4.5 g head basket, current combo, t_max 0.2 s)\n')
out.append('| model | variant | peak MPa | plateau (0.1-0.2 s) | peak/plateau |')
out.append('|---|---|---|---|---|')
for n in ('base', 'capZ+area23'):
    b = rk[f'rk_current_{n}']
    out.append(f"| {n} | Ma 2020 erosive (default) | {b['peak_MPa']:.2f} | {b['plateau_MPa']:.2f} | {b['peak_over_plateau']:.2f} |")
    for v in ('saint_robert', 'mukunda_paul'):
        r = ex[f'ex_{n}_{v}']
        out.append(f"| {n} | erosive = {v} | {r['peak_MPa']:.2f} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} |")
    r = ex[f'ex_{n}_pyro_only']
    out.append(f"| {n} | grain cannot ignite (pyrogen gas only) | {r['peak_MPa']:.2f} | - | - |")
for law in ('zhu_lo', 'zhu_hi'):
    r = ex[f'ex_{law}+area23']
    out.append(f"| {law}+area23 | Ma 2020 (r0 = {'5.271e-3' if law=='zhu_lo' else '6.848e-3'} m/s/MPa) | {r['peak_MPa']:.2f} | {r['plateau_MPa']:.2f} | {r['peak_over_plateau']:.2f} |")
out.append('')

(R / 'ob1_tables.md').write_text('\n'.join(out), encoding='utf-8')
print('wrote', R / 'ob1_tables.md', len(out), 'lines')
