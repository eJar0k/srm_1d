"""Build docs/pyrogen_realism/results/ig2_tables.md from the IG2 JSONL results. Run from the worktree root."""
import json, math, sys
from pathlib import Path
R = Path('docs/pyrogen_realism/results')
def load(p): return [json.loads(l) for l in open(R / p) if l.strip()]
ms = lambda x: '-' if x is None or (isinstance(x, float) and math.isnan(x)) else f'{x*1e3:.2f}'
out = []
# threshold
b = load('ig2_sweep_bisect.jsonl')
out.append('## Ignition threshold (min charge with >=99 % of grain cells lit within 0.2 s; bisected to 2 %)\n')
out.append('| topology | ignition h | DeMar | min charge [g] |\n|---|---|---|---|')
from collections import defaultdict
d = defaultdict(list)
for r in b: d[(r['topo'], r['ign_h'], r['demar'])].append((r['g'], r['frac_ignited']))
for k in sorted(d):
    ok = min(g for g, f in d[k] if f >= .99)
    out.append(f'| {k[0]} | {k[1]} | {k[2]} | {ok:.3f} |')
# realistic
rc = load('ig2_sweep_realistic.jsonl'); rc.sort(key=lambda r: (r['pyrogen'], r['g'], r['ign_h'], r['demar']))
out.append('\n## Realistic charges (head basket, 3 mm cubes L/D 1, Chunc; t_max 0.2 s; plateau = median 0.1-0.2 s)\n')
out.append('| pyrogen | g | ignition h | DeMar | first ign [ms] | t90 [ms] | peak [MPa] | t_peak [ms] | plateau [MPa] | peak/plateau |\n|---|---|---|---|---|---|---|---|---|---|')
for r in rc:
    out.append(f"| {r['pyrogen']} | {r['g']} | {r['ign_h']} | {r['demar']} | {ms(r['t_first_ign_s'])} | {ms(r['t90_s'])} | {r['peak_MPa']:.2f} | {ms(r['t_peak_s'])} | {r['plateau_MPa']:.2f} | {r['peak_MPa']/r['plateau_MPa']:.2f} |")
cr = load('ig2_chunc_real.jsonl'); cr.sort(key=lambda r: r['tag'])
out.append('\n## Realistic charges at t_max 1.0 s (plateau = median 0.4-0.9 s; Chunc exp: peak 8.89, plateau 8.66, ratio 1.027)\n')
out.append('| tag | peak [MPa] | t_peak [ms] | plateau | peak/plateau | first ign [ms] |\n|---|---|---|---|---|---|')
for r in cr:
    out.append(f"| {r['tag']} | {r['peak_MPa']:.2f} | {ms(r['t_peak_s'])} | {r['plateau_MPa']:.2f} | {r['peak_MPa']/r['plateau_MPa']:.2f} | {ms(r['t_first_ign_s'])} |")
v = load('ig2_valid.jsonl'); v.sort(key=lambda r: (r['motor'], r['ign_h'], r['demar']))
out.append('\n## Validation, canonical motors, library BPNV Sutton-mass forward-plenum igniter (B0 rows)\n')
out.append('| motor | ignition h | DeMar | igniter [g] | first ign [ms] | peak [MPa] | t_peak [s] | plateau / trough / mid | ratios | steps | wall [s] |\n|---|---|---|---|---|---|---|---|---|---|---|')
for r in v:
    extra = {'chunc': f"plateau {r.get('plateau', 0):.2f}", 'balls': f"trough {r.get('trough', 0):.2f} / mid {r.get('mid', 0):.2f}", 'haseg': '-'}[r['motor']]
    rat = {'chunc': f"pk/pl {r.get('ratio', 0):.3f}", 'balls': f"pk/tr {r.get('pk_tr', 0):.2f}, pk/mid {r.get('pk_mid', 0):.2f}", 'haseg': '-'}[r['motor']]
    out.append(f"| {r['motor']} | {r['ign_h']} | {r['demar']} | {r['pyro_mass_g']:.2f} | {ms(r['first_ign_s'])} | {r['peak']:.3f} | {r['t_peak']:.3f} | {extra} | {rat} | {r['steps']} | {r['wall_s']:.0f} |")
(R / 'ig2_tables.md').write_text('\n'.join(out) + '\n', encoding='utf-8')
print('\n'.join(out))
