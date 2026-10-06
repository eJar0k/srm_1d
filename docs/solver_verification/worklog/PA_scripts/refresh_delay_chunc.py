"""PA/P3: ignition-to-next-burn-refresh delay on canonical Chunc.

Run from the worktree root:  PYTHONPATH="$PWD" python <this file>
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.tools.ignition_diagnostics import ignition_refresh_delays

CANON = dict(roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None,
             T_ignition=756.0, print_interval=0.0, verbose=False)
res = run_from_ric(str(ROOT / 'motors' / 'machbusterNew.ric'), t_max=0.05,
                   P_cutoff=0.01e6, erosive_model='ma2020',
                   erosive_mp_cap=0.0, **CANON)[0]
out = ignition_refresh_delays(res)
s = res['summary']
print(f"steps={s['steps']} interval={out['interval']} "
      f"max_burn_refresh_gap_s={s['max_burn_refresh_gap_s']:.6e}")
print(f"n_ignitions={out['n_ignitions']} n_with_refresh={out['n_with_refresh']}")
print(f"ignition_rate_refresh_delay_max_s={out['max_s']:.6e}")
print(f"ignition_rate_refresh_delay_median_s={out['median_s']:.6e}")
