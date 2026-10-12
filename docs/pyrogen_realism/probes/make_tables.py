import json,sys
def load(f): return [json.loads(l) for l in open(f)]
def f(v): return '%.2f'%(1e3*v) if v else '-'
def tab(rows):
    print('| topology | lb (g) | ignites (frac) | t_first ms | t90 ms | peak MPa | t_peak ms | P_end(0.2s) | pyro dur ms |\n|---|---|---|---|---|---|---|---|---|')
    for d in sorted(rows,key=lambda d:(d['topo'],d['lb'])):
        print(f"| {d['topo']} | {d['lb']:g} ({d['lb']*453.59:.2f}) | {d['frac_ignited']:.2f} | {f(d['t_first_ign_s'])} | {f(d['t90_s'])} | {d['peak_MPa']:.2f} | {1e3*d['t_peak_s']:.1f} | {d['P_end_MPa']:.2f} | {1e3*d['pyro_duration_s']:.1f} |")
d='docs/pyrogen_realism/results/'
print('### baseline');tab(load(d+'sweep.jsonl')+load(d+'sweep_fine.jsonl'))
print('\n### hc x0.2');tab(load(d+'sweep_hc0.2.jsonl'))
