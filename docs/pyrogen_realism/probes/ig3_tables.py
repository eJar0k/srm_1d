"""Build docs/pyrogen_realism/results/ig3_tables.md from the IG3 stage-1 JSONL results. Run from the worktree root."""
import json
import math
from collections import defaultdict
from pathlib import Path

R = Path('docs/pyrogen_realism/results')


def load(p):
    return [json.loads(l) for l in open(R / p) if l.strip()]


ms = lambda x: '-' if x is None or (isinstance(x, float) and math.isnan(x)) else f'{x * 1e3:.2f}'
out = []

# (c) head basket / forward plenum thresholds, DeMar off
b = load('ig3_sweep_bisect.jsonl')
out.append('## (c) Ignition threshold, Chunc, DeMar off, Peretz h (min charge with >=99 % of grain cells lit within 0.2 s; bisected to 2 %)\n')
out.append('| topology | Ap reading | min charge [g] |\n|---|---|---|')
d = defaultdict(list)
for r in b:
    d[(r['topo'], r['ap_ref'])].append((r['g'], r['frac_ignited']))
for k in sorted(d):
    ok = min(g for g, f in d[k] if f >= .99)
    out.append(f'| {k[0]} | {k[1]} | {ok:.3f} |')

# (c) per-motor forward-plenum thresholds
try:
    mt = load('ig3_motor_threshold.jsonl')
    out.append('\n## (c) Minimum forward-plenum BPNV charge that lights each canonical motor (>=99 % lit within 0.25 s, bisected to 3 %)\n')
    out.append('| motor | DeMar | Ap reading | min charge [g] |\n|---|---|---|---|')
    d = defaultdict(list)
    for r in mt:
        d[(r['motor'], r['demar'], r['ap_ref'])].append((r['g'], r['frac_ignited']))
    for k in sorted(d):
        ok = min(g for g, f in d[k] if f >= .99)
        out.append(f'| {k[0]} | {k[1]} | {k[2]} | {ok:.3f} |')
except FileNotFoundError:
    pass

# (b) realistic charges
rc = load('ig3_sweep_realistic.jsonl')
rc.sort(key=lambda r: (r['pyrogen'], r['g'], r['ap_ref']))
out.append('\n## (b) Realistic charges (Chunc head basket, 3 mm cubes L/D 1, DeMar off, Peretz h; t_max 0.2 s; plateau = median 0.1-0.2 s)\n')
out.append('| pyrogen | g | Ap reading | first ign [ms] | t90 [ms] | peak [MPa] | t_peak [ms] | plateau [MPa] | peak/plateau |\n|---|---|---|---|---|---|---|---|---|')
for r in rc:
    out.append(f"| {r['pyrogen']} | {r['g']} | {r['ap_ref']} | {ms(r['t_first_ign_s'])} | {ms(r['t90_s'])} | {r['peak_MPa']:.2f} | {ms(r['t_peak_s'])} | {r['plateau_MPa']:.2f} | {r['peak_MPa'] / r['plateau_MPa']:.2f} |")

# (a) canonical motors
v = load('ig3_valid.jsonl')
order = {'chunc': 0, 'balls': 1, 'haseg': 2}
v.sort(key=lambda r: (order[r['motor']], r['ap_ref']))
out.append('\n## (a) Canonical motors, library BPNV Sutton-mass forward-plenum igniter, DeMar step, Peretz h (B0 rows)\n')
out.append('| motor | Ap reading | igniter [g] | first ign [ms] | peak [MPa] | t_peak [ms] | plateau / trough / mid | ratios | steps |\n|---|---|---|---|---|---|---|---|---|')
for r in v:
    extra = {'chunc': f"plateau {r.get('plateau', 0):.3f}", 'balls': f"trough {r.get('trough', 0):.3f} / mid {r.get('mid', 0):.3f}", 'haseg': '-'}[r['motor']]
    rat = {'chunc': f"pk/pl {r.get('ratio', 0):.3f} (exp 1.027)", 'balls': f"pk/tr {r.get('pk_tr', 0):.2f}, pk/mid {r.get('pk_mid', 0):.2f} (exp ~1.50 / ~1.35)", 'haseg': 'peak 0.95x exp (6.44)'}[r['motor']]
    out.append(f"| {r['motor']} | {r['ap_ref']} | {r['pyro_mass_g']:.2f} | {ms(r['first_ign_s'])} | {r['peak']:.3f} | {r['t_peak'] * 1e3:.1f} | {extra} | {rat} | {r['steps']} |")
(R / 'ig3_tables.md').write_text('\n'.join(out) + '\n', encoding='utf-8')
print('\n'.join(out))
