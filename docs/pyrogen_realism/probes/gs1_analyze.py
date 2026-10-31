"""GS1 analysis: tables, convergence orders, Richardson values, contact-front widths.
Usage: python docs/pyrogen_realism/probes/gs1_analyze.py
Reads results/gs1/{rows,profiles,full}.jsonl (+ rows.jsonl.guard.jsonl); writes results/gs1/tables.md."""
import json
import math
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parents[1] / 'results' / 'gs1'
NL = chr(10)


def load(name):
    p = R / name
    return [json.loads(l) for l in open(p) if l.strip()] if p.exists() else []


rows = load('rows.jsonl')
guard = {}
for g in load('rows.jsonl.guard.jsonl'):
    guard[g['tag']] = g


def peak_gb(tag):
    g = guard.get(tag)
    if not g or not g.get('guard'):
        return '-'
    return g['guard'].split('peak_job_gb=')[1].split()[0]


by = {}
for r in rows:
    tag = r['spec']['tag']
    base, _, rest = tag.partition('_N')
    by.setdefault(base, {})[int(rest)] = r

out = []


def f(v, d=2):
    return '-' if v is None or (isinstance(v, float) and math.isnan(v)) else f'{v:.{d}f}'


def order(vals, Ns):
    """Observed order from the last three grids (constant refinement ratio r); Richardson value from the finest two."""
    ks = [n for n in Ns if vals.get(n) is not None]
    if len(ks) < 3:
        return None, None
    n1, n2, n3 = ks[-3:]
    a, b, c = vals[n1], vals[n2], vals[n3]
    r = n3 / n2
    if abs(n2 / n1 - r) > 1e-9 or (b - a) == 0 or (c - b) == 0:
        return None, None
    ratio = (b - a) / (c - b)
    if ratio <= 1.0:
        return None, None
    p = math.log(ratio) / math.log(r)
    return p, c + (c - b) / (r ** p - 1.0)


metrics = [('t_first_ms', 't first [ms]'), ('t50_ms', 't50 [ms]'), ('t90_ms', 't90 [ms]'), ('front_speed_10_90_m_s', 'ignition-front v 10-90 [m/s]'),
           ('n_dark', 'dark cells'), ('dark_end_mm', 'dark-zone end [mm]'), ('peak_MPa', 'startup peak [MPa]'), ('t_peak_ms', 't_peak [ms]'),
           ('plateau_MPa', 'window median [MPa]'), ('peak_over_plateau', 'peak / window')]
for base, d in sorted(by.items()):
    Ns = sorted(d)
    out.append(NL + f'### {base}' + NL)
    out.append('| metric | ' + ' | '.join(f'N={n}' for n in Ns) + ' | order p | Richardson |' + NL + '|---|' + '---|' * (len(Ns) + 2))
    for key, lab in metrics:
        vals = {n: d[n].get(key) for n in Ns}
        p, rich = order(vals, Ns)
        dd = 3 if key in ('plateau_MPa', 'peak_over_plateau') else 2
        out.append(f'| {lab} | ' + ' | '.join(f(vals[n], dd) for n in Ns) + f' | {f(p, 2)} | {f(rich, dd)} |')
    out.append('| wall [s] / peak GB | ' + ' | '.join(f"{d[n]['wall_s']:.0f} / {peak_gb(d[n]['spec']['tag'])}" for n in Ns) + ' | | |')
    out.append('| steps | ' + ' | '.join(str(d[n]['steps']) for n in Ns) + ' | | |')

# ignition-front speed by quartile
out.append(NL + '## Ignition-front speed by quartile of the lit cells (m/s; later quartile = further downstream)' + NL)
out.append('| row | N | Q1 | Q2 | Q3 | Q4 |' + NL + '|---|---|---|---|---|---|')
for base, d in sorted(by.items()):
    for n in sorted(d):
        q = d[n].get('front_speed_quartiles_m_s') or [None] * 4
        out.append(f'| {base} | {n} | ' + ' | '.join(f(v, 0) for v in q) + ' |')


# contact front, from the pre-ignition profile rows
def crossing(x, th, level, start=0):
    """first x scanning downstream from `start` where th falls below `level` (linear interpolation)."""
    idx = np.flatnonzero(th[start:] < level)
    if idx.size == 0:
        return None
    i = start + int(idx[0])
    if i == 0:
        return float(x[0])
    a, b = th[i - 1], th[i]
    fr = (a - level) / (a - b) if a != b else 0.0
    return float(x[i - 1] + fr * (x[i] - x[i - 1]))


prof = load('profiles.jsonl')
if prof:
    out.append(NL + '## Contact front: gas-T rise and hot-gas fraction along the bore before first light (clean rows, profile rows)' + NL)
    out.append('x50 = position of the 50 % level of the rise (the front), w = distance between the 90 % and 10 % levels, lead = x(1 %) - x50 (how far the leading edge '
               'runs ahead of the front). T rise = (T - 293)/(T_max - 293); hot-gas fraction = 1 - Y_ambient (normalised). Positions in mm; dx in mm.' + NL)
    out.append('| row | t [ms] | N | dx | T: x50 | T: w | w/dx | T: lead | hot-gas: x50 | hot-gas: w | hot-gas: lead |' + NL + '|---|---|---|---|---|---|---|---|---|---|---|')
    summ = {}
    for r in sorted(prof, key=lambda r: (r['spec']['tag'].split('_N')[0], r['spec']['t_max'], r['spec']['N'])):
        pr = r['prof']
        x = np.array(pr['x_mm'])
        T = np.array(pr['T'])
        Y = np.array(pr['Y'])
        amb = [k for k, nme in enumerate(pr['species']) if 'ambient' in nme][0]
        dx = r['dx_mm']
        th = (T - 293.0) / (T.max() - 293.0)
        phi = 1.0 - Y[:, amb]
        phi = (phi - phi.min()) / (phi.max() - phi.min())
        res = {}
        for nm, arr in (('T', th), ('Y', phi)):
            st = int(np.argmax(arr))
            res[nm] = {lv: crossing(x, arr, lv, st) for lv in (0.9, 0.5, 0.1, 0.01)}

        def gap(nm, a, b):
            c = res[nm]
            return None if c[a] is None or c[b] is None else c[b] - c[a]
        tag = r['spec']['tag'].split('_N')[0]
        tms = r['spec']['t_max'] * 1e3
        wT = gap('T', 0.9, 0.1)
        out.append(f"| {tag} | {tms:.1f} | {r['spec']['N']} | {dx:.2f} | {f(res['T'][0.5], 0)} | {f(wT, 0)} | {f(None if wT is None else wT / dx, 1)} | "
                   f"{f(gap('T', 0.5, 0.01), 0)} | {f(res['Y'][0.5], 0)} | {f(gap('Y', 0.9, 0.1), 0)} | {f(gap('Y', 0.5, 0.01), 0)} |")
        summ[(tag, round(tms, 1), r['spec']['N'])] = (dx, res['T'][0.5], wT, gap('T', 0.5, 0.01))
    out.append(NL + '#### First-order numerical-diffusion comparison (T front)' + NL)
    out.append('Front speed u_f from x50 between 0.5 and 0.9 ms. First-order upwind diffusion D = u dx (1 - C)/2 with C = 0.3 (the CFL target) gives a 10-90 width of '
               '3.62 sqrt(D t), t = 0.9 ms (the age of the front). A numerically smeared front narrows as sqrt(dx) (width ratio 1.41 per halving of dx); a physical width stays '
               'constant (ratio 1.0).' + NL)
    out.append('| row | N | dx | u_f [m/s] | measured w at 0.9 ms [mm] | predicted numerical w [mm] | measured / predicted | w(coarser)/w(this) |' + NL + '|---|---|---|---|---|---|---|---|')
    for tag in sorted({k[0] for k in summ}):
        prev = None
        for N in (100, 200, 400, 800):
            a = summ.get((tag, 0.5, N))
            b = summ.get((tag, 0.9, N))
            if not a or not b or a[1] is None or b[1] is None or b[2] is None:
                continue
            uf = (b[1] - a[1]) / 0.4
            D = abs(uf) * (b[0] * 1e-3) * 0.7 / 2.0
            wp = 3.62 * math.sqrt(D * 0.9e-3) * 1e3
            out.append(f'| {tag} | {N} | {b[0]:.2f} | {uf:.0f} | {b[2]:.0f} | {wp:.0f} | {b[2] / wp:.1f} | {f(prev / b[2] if prev else None, 2)} |')
            prev = b[2]

full = list({r['spec']['tag']: r for r in load('full.jsonl')}.values())     # last row per tag
if full:
    out.append(NL + '## Full burns at the fixed (threshold) charges: startup, window median, late peak' + NL)
    out.append('| row | N | steps | wall [s] | first light [ms] | dark cells | startup peak [MPa] | t_peak [s] | window median [MPa] | peak/window | late peak (>0.5 s) [MPa] | t late [s] |' + NL + '|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in sorted(full, key=lambda r: r['spec']['tag']):
        out.append(f"| {r['spec']['tag']} | {r['spec']['N']} | {r['steps']} | {r['wall_s']:.0f} | {f(r.get('t_first_ms'), 2)} | {r['n_dark']} | {r['peak_MPa']:.3f} | "
                   f"{r['t_peak_ms'] / 1e3:.3f} | {f(r['plateau_MPa'], 3)} | {f(r.get('peak_over_plateau'), 3)} | {f(r.get('late_peak_MPa'), 3)} | {f(r.get('t_late_peak_s'), 2)} |")

(R / 'tables.md').write_text(NL.join(out) + NL, encoding='utf-8')
print(NL.join(out))


# ---- error vs the Richardson-extrapolated value (from the three finest grids), per metric and grid -----------------
def rich_val(vals):
    p, r = order(vals, [100, 200, 400, 800])
    return p, r


err_rows = []
summary = {}
keys = [('t_first_ms', 'first light'), ('t90_ms', 't90'), ('peak_MPa', 'startup peak'), ('t_peak_ms', 't_peak'), ('peak_over_plateau', 'peak/window'),
        ('front_speed_10_90_m_s', 'front speed')]
for key, lab in keys:
    for base, d in sorted(by.items()):
        vals = {n: d[n].get(key) for n in d}
        if len([v for v in vals.values() if v is not None]) < 4 or any(vals.get(n) is None for n in (200, 400, 800)):
            continue
        p, rv = rich_val(vals)
        if rv is None or p is None or p < 0.3:
            continue
        for n in (100, 200, 400, 800):
            if vals.get(n) is not None:
                summary.setdefault((lab, n), []).append((base, abs(vals[n] - rv) / abs(rv)))
out.append(NL + '## Error of each grid vs the Richardson-extrapolated value (rows with a clean 200-400-800 convergence, observed p >= 0.3)' + NL)
out.append('Worst case over the rows listed in the last column; the extrapolation assumes the asymptotic range, so treat the errors as indicative (observed p is 0.4-1.4, i.e. first-order-like).' + NL)
out.append('| metric | N=100 | N=200 | N=400 | N=800 | rows |' + NL + '|---|---|---|---|---|---|')
for key, lab in keys:
    cells = []
    nrows = 0
    for n in (100, 200, 400, 800):
        e = summary.get((lab, n), [])
        nrows = max(nrows, len(e))
        cells.append(f'{100 * max(v for _, v in e):.0f} % ({min(e, key=lambda t: -t[1])[0]})' if e else '-')
    out.append(f'| {lab} | ' + ' | '.join(cells) + f' | {nrows} |')
(R / 'tables.md').write_text(NL.join(out) + NL, encoding='utf-8')
print(NL.join(out[-12:]))
