import json
f=lambda v: '%.2f'%(1e3*v) if v else '-'
print('| case | first ign ms | t90 ms | peak MPa | t_peak ms | pyro dur ms |\n|---|---|---|---|---|---|')
for fn in ['sens','cx','cy','cz']:
    for l in open(f'docs/pyrogen_realism/results/{fn}.jsonl'):
        d=json.loads(l)
        if d['tag']=='cz_warm': continue
        print(f"| {d['tag']} ({d['topo'][:4]},{d['lb']:g} lb) | {f(d['t_first_ign_s'])} | {f(d['t90_s'])} | {d['peak_MPa']:.2f} | {1e3*d['t_peak_s']:.1f} | {1e3*d['pyro_duration_s']:.1f} |")
