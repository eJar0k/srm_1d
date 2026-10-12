"""Print P_head / pyrogen mdot / nozzle massflow vs time. Usage: ob1_hist.py '<json spec>' [dt_ms]"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, ob1_case
spec = json.loads(sys.argv[1]); dtm = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
m, res = ob1_case.run_case(spec, return_res=True)
t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6
md = np.asarray(res['mdot_ig']); mf = np.asarray(res['massflow']); mp = np.asarray(res['m_pyrogen']) * 1e3
print('peak', m['peak_MPa'], 't_peak', m['t_peak_s'], 'pyro dur', m['pyro_duration_s'])
print('t_ms  P_head  mdot_ig[kg/s]  nozzle_mdot  m_pyro_left[g]')
for tt in np.arange(0, spec.get('t_max', 0.1) * 1e3, dtm):
    i = int(np.searchsorted(t, tt * 1e-3))
    if i >= t.size: break
    print(f'{tt:6.1f} {P[i]:7.2f} {md[i]:9.4f} {mf[i]:9.3f} {mp[i]:8.3f}')
