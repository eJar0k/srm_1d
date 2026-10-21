"""EN1: one (kind, energy, momentum, N) duct case -> JSON row.  kind = shock | venturi.

    python docs/solver_verification/probes/en1_duct_ladder.py shock 1 0 200 out.json
args: kind total(0/1) conservative_momentum(0/1) N out.json [t_end]
(t_end default 0.06 s.  The marching loop is bounded: step cap, wall cap, dt floor; a
stalled or diverged case is returned with ``stall_reason`` set.)
"""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import en1_verification as v

kind, total, cm, N = sys.argv[1], sys.argv[2] == '1', sys.argv[3] == '1', int(sys.argv[4])
tend = float(sys.argv[6]) if len(sys.argv) > 6 else 0.06
kw = dict(shock_ar=1.5) if kind == 'shock' else dict(throat_M=0.8)
t0 = time.time()
r = v.duct_case(N, total, cons_mom=cm, t_end=tend, **kw)
row = dict(kind=kind, total=total, cons_mom=cm, N=N, wall=time.time() - t0,
           steps=r['steps'], diverged=bool(r['diverged']), stall=r['stall_reason'], t=r['tube'].t,
           resid=r['resid'][-1][1], P_back=r['P_back'])
if not r['diverged']:
    row['metrics'] = v.duct_metrics(r)
    row['x'] = r['x'].tolist(); row['P'] = r['P'].tolist(); row['T'] = r['T'].tolist()
    row['M'] = r['M'].tolist(); row['P_e'] = r['P_e'].tolist(); row['T_e'] = r['T_e'].tolist()
    row['M_e'] = r['M_e'].tolist()
Path(sys.argv[5]).write_text(json.dumps(row, default=float))
print(kind, total, cm, N, 'wall', round(row['wall'], 1), 'div', row['diverged'], flush=True)
