"""IG4: pathfin54 on the GUI path (srm1d_plugin simulate() with the solver-config defaults).
Usage: ig4_pathfin.py [explicit_charge_g]   Run from the worktree root with PYTHONPATH=$PWD."""
import os, sys, time, json, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
from srm_1d.fmm_grain import _setup_openmotor_path
_setup_openmotor_path()
import srm_1d.srm1d_plugin as pl
from motorlib import solvers, motor as om_motor
from srm_1d.openmotor_adapter import load_ric

ric = ROOT / 'motors' / '2025.12.25 54-2800 ST2.0.ric'
d = load_ric(str(ric))
if len(sys.argv) > 1:
    d['igniter'] = dict(d['igniter']); d['igniter']['mass'] = float(sys.argv[1]) * 1e-3
motor = om_motor.Motor(d)
cfg = pl._solver_config_defaults()
solver = solvers.get_solver('srm_1d-transient')
# capture the raw run too
import srm_1d.simulation as sim
orig = sim.run_simulation
cap = {}
def spy(*a, **k):
    r = orig(*a, **k); cap['res'] = r; return r
import srm_1d.srm1d_plugin
sim.run_simulation = spy
t0 = time.time()
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    sr = solver.simulate(motor, config=cfg)
wall = time.time() - t0
s = cap['res']['summary']
t = np.asarray(cap['res']['time'])
print(json.dumps(dict(attempts=s['ignition_attempts'], masses_g=[round(m*1e3, 4) for m in s['ignition_attempt_masses_kg']],
      final_g=round(s['pyrogen_mass_final_kg']*1e3, 4), cells_lit=s['grain_cells_ignited'], grain_cells=s['grain_cells'],
      P_peak_MPa=round(s['P_peak']/1e6, 4), run_end_s=round(float(t[-1]), 4), termination=s['termination'],
      ignited=s['ignited'], wall_s=round(wall, 1), attempt_wall_s=[round(x, 1) for x in s['ignition_attempt_wall_s']],
      igniter_peak_mdot=s.get('igniter_peak_mdot_kg_s'), grain_mdot_after=s.get('grain_mdot_after_igniter_kg_s'),
      sr_success=bool(sr.success), alerts=[str(a.description) if hasattr(a, 'description') else str(a) for a in getattr(sr, 'alerts', [])][:3],
      warns=[str(x.message)[:200] for x in w if issubclass(x.category, RuntimeWarning)]), indent=1))
