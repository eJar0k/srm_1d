"""Print a JSONL results file as a table. Usage: ob1_table.py <jsonl> [tag-prefix]"""
import sys, json
rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
pre = sys.argv[2] if len(sys.argv) > 2 else ''
rows = [r for r in rows if r['tag'].startswith(pre)]
f = lambda x, p=2: '-' if x is None else (f'{x:.{p}f}' if isinstance(x, float) else str(x))
print('| tag | first ign ms | t90 ms | peak MPa | t_peak ms | plateau | pk/pl | pyro t50 ms | t90 ms | dur ms |')
print('|---|---|---|---|---|---|---|---|---|---|')
for r in sorted(rows, key=lambda r: r['tag']):
    if 'error' in r:
        print('|', r['tag'], '| ERROR', r['error'][:60]); continue
    ms = lambda k: None if r.get(k) is None else r[k] * 1e3
    print('|', r['tag'], '|', f(ms('t_first_ign_s')), '|', f(ms('t90_s')), '|', f(r['peak_MPa']), '|', f(ms('t_peak_s'), 1), '|',
          f(r['plateau_MPa']), '|', f(r['peak_over_plateau']), '|', f(ms('pyro_t50_s'), 1), '|', f(ms('pyro_t90_s'), 1), '|',
          f(ms('pyro_duration_s'), 1), '|')
