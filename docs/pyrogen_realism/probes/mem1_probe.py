"""MEM1: one run of a motor; prints est vs actual history steps. Run under the guard:
  python -m srm_1d.tools.guarded_run --mem-gb 4 --timeout-min 30 -- python docs/pyrogen_realism/probes/mem1_probe.py <ric> [t_max] [P_cutoff_Pa] [gui]
'gui' = the GUI path (srm1d_plugin.simulate with _solver_config_defaults, t_max 10)."""
import sys, json, time, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d.simulation as sim
from srm_1d.openmotor_adapter import run_from_ric, load_ric

ric = sys.argv[1]
t_max = float(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != '-' else None
pc = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != '-' else None
gui = 'gui' in sys.argv[4:]
cap = {}
import functools
orig = sim._run_simulation_once
@functools.wraps(orig)
def spy(*a, **k):
    r = orig(*a, **k); s = r['summary']
    cap.setdefault('runs', []).append(dict(est=s['history_estimated_steps'], cap=s['history_capacity'], steps=int(s['steps']), term=s['termination'], t_end=float(r['time'][-1])))
    return r
sim._run_simulation_once = spy
t0 = time.time()
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    if gui:
        from srm_1d.fmm_grain import _setup_openmotor_path
        _setup_openmotor_path()
        import srm_1d.srm1d_plugin as pl
        from motorlib import solvers, motor as om_motor
        d = load_ric(str(ROOT / ric))
        sr = solvers.get_solver('srm_1d-transient').simulate(om_motor.Motor(d), config=pl._solver_config_defaults())
    else:
        kw = dict(print_interval=0.0, verbose=False)
        if t_max: kw['t_max'] = t_max
        if pc: kw['P_cutoff'] = pc
        run_from_ric(str(ROOT / ric), **kw)
print(json.dumps(dict(ric=ric, gui=gui, runs=cap['runs'], wall_s=round(time.time() - t0, 1),
      warns=[str(x.message)[:120] for x in w if issubclass(x.category, RuntimeWarning)][:4])))
