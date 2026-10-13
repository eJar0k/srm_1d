"""EN1 Sod ladder: {sensible,total} x {legacy,conservative momentum} on N cells, plus a CFL series.
    python docs/solver_verification/probes/en1_sod_table.py [out.json]
"""
import json, sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import en1_verification as v

NS = (100, 200, 400, 800, 1600)
rows = []
def one(N, total, cm, cfl=0.3, K=1):
    r = v.sod_case(N, total, cfl=cfl, cons_mom=cm, K=K)
    m = v.sod_metrics(r)
    d = v.sod_total_energy_drift(r)
    pl, pr = m['plateaus']['star-L'], m['plateaus']['star-R']
    return dict(N=N, total=total, cons_mom=cm, cfl=cfl, K=K, steps=r['steps'],
                xs=m['xs_num'], xs_ex=m['xs_exact'], S=m['S_num'], S_ex=m['S_exact'],
                l1=m['l1'], starL=pl, starR=pr,
                dE_rel=float(d['dE_rel']), dM_rel=float(d['dM_rel']))
for total in (False, True):
    for cm in (False, True):
        for N in NS:
            rows.append(one(N, total, cm))
for cfl in (0.3, 0.15, 0.075, 0.0375):
    rows.append(one(400, True, True, cfl=cfl))
for K in (1, 3):
    rows.append(one(400, True, True, K=K))
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(json.dumps(rows, default=float, indent=1))
def tag(r): return ('total' if r['total'] else 'sensible') + ('+cons' if r['cons_mom'] else '')
print("== ladder (CFL 0.3, K=1) ==")
prev = {}
for r in rows[:20]:
    t = tag(r)
    o = ''
    if t in prev:
        o = ' order ' + ' '.join(f"{k}:{math.log(prev[t]['l1'][k]/r['l1'][k])/math.log(2):.2f}" for k in ('rho','u','P','T'))
    prev[t] = r
    print(f"{t:15s} N={r['N']:5d} xs {r['xs']:.4f}/{r['xs_ex']:.4f} S {r['S']:.1f}/{r['S_ex']:.1f} L1 rho {r['l1']['rho']:.4f} u {r['l1']['u']:.4f} P {r['l1']['P']:.4f} T {r['l1']['T']:.4f}"
          f" | *L rho {r['starL']['rho'][0]:.4f}/{r['starL']['rho'][1]:.4f} u {r['starL']['u'][0]:.1f}/{r['starL']['u'][1]:.1f} P {r['starL']['P'][0]:.0f}/{r['starL']['P'][1]:.0f} T {r['starL']['T'][0]:.1f}/{r['starL']['T'][1]:.1f}"
          f" | *R rho {r['starR']['rho'][0]:.4f}/{r['starR']['rho'][1]:.4f} P {r['starR']['P'][0]:.0f}/{r['starR']['P'][1]:.0f} T {r['starR']['T'][0]:.1f}/{r['starR']['T'][1]:.1f}"
          f" | dE {r['dE_rel']:.2e} dM {r['dM_rel']:.2e}{o}")
print("== total+cons N=400 CFL / K series ==")
for r in rows[20:]:
    print(f"cfl {r['cfl']} K {r['K']}: S {r['S']:.1f}/{r['S_ex']:.1f} P*L {r['starL']['P'][0]:.0f}/{r['starL']['P'][1]:.0f} u*L {r['starL']['u'][0]:.1f}/{r['starL']['u'][1]:.1f} L1 {r['l1']} dE {r['dE_rel']:.2e} dM {r['dM_rel']:.2e}")
