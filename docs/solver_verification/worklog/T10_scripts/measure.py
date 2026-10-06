"""T10 conservation measurement.

    PYTHONPATH=$PWD python docs/solver_verification/worklog/T10_scripts/measure.py chunc 1.0 out.json
Dumps: Sum of the geometry storage jump (J, kg) from T6's audit, Sum of the
T10 fill-debit channels (if present), summary mass_balance_error / mass_produced.
"""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.tools.ignition_diagnostics import sensible_enthalpy_audit

RICS = {"chunc": "machbusterNew.ric", "balls": "BALLSstick.ric"}
motor, t_max, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
print("[provenance]", srm_1d.__file__, flush=True)
res = run_from_ric(str(ROOT / "motors" / RICS[motor]), t_max=t_max, roughness=32e-6,
                   kappa=0.44, T_ignition=756.0, pyrogen="bpnv", pyrogen_mass=None,
                   print_interval=0.0, verbose=False)[0]
s = res["summary"]
aud = sensible_enthalpy_audit(res)["run"]
dt = np.asarray(res["dt"])
d = {"steps": int(s["steps"]), "P_peak": s["P_peak"], "t_peak": s["t_peak"],
     "mass_produced": s["mass_produced"], "mass_balance_error": s["mass_balance_error"],
     "mass_nozzle": s.get("mass_nozzle"),
     "geom_J_signed": aud["jumps_J"]["geometry"], "geom_kg": aud["jumps_kg"]["geometry"],
     "eos_J": aud["jumps_J"]["eos"], "mix_J": aud["jumps_J"]["mixture"],
     "total_J": aud["jumps_J"]["total"], "total_kg": aud["jumps_kg"]["total"],
     "budget": aud["budget_J"],
     "src_mismatch_max": aud["source_ledger_mismatch_max_abs_W"]}
if "regression_fill_debit_power" in res.keys():
    d["debit_J"] = float(np.sum(np.asarray(res["regression_fill_debit_power"]) * dt))
    d["debit_kg"] = float(np.sum(np.asarray(res["regression_fill_debit_massflow"]) * dt))
Path(out).write_text(json.dumps(d, indent=1, default=float))
print(json.dumps({k: v for k, v in d.items() if k != "budget"}, default=float))
