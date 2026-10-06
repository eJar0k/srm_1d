"""T6 measurement: run a motor and dump sensible_enthalpy_audit as JSON.

    PYTHONPATH=$PWD python docs/solver_verification/worklog/T6_scripts/audit_run.py chunc 1.0 out.json
    motors: chunc | balls   (committed .ric, canonical knobs, library BPNV)
Snapshots every 0.5 ms so the fill/spike phases are resolved (omitted-term
part D). Run from the worktree root.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
import numpy as np  # noqa: E402
import srm_1d  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402
from srm_1d.tools.ignition_diagnostics import sensible_enthalpy_audit  # noqa: E402

RICS = {"chunc": "machbusterNew.ric", "balls": "BALLSstick.ric"}


def main():
    motor, t_max, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    print("[provenance]", srm_1d.__file__, flush=True)
    res = run_from_ric(
        str(ROOT / "motors" / RICS[motor]), t_max=t_max,
        roughness=32e-6, kappa=0.44, T_ignition=756.0, pyrogen="bpnv",
        pyrogen_mass=None, print_interval=0.0, verbose=False,
        snapshot_interval=5e-4)[0]
    audit = sensible_enthalpy_audit(res)
    s = res["summary"]
    audit["meta"] = {k: s[k] for k in (
        "steps", "cells", "wall_time", "P_peak", "t_peak", "mass_produced",
        "pressure_floor_activations", "termination")}
    audit["meta"]["steps_per_s"] = s["steps"] / s["wall_time"]
    audit["meta"]["nozzle_mass_mismatch_max_abs_kg_s"] = float(
        np.max(np.abs(res["nozzle_mass_mismatch_kg_s"])))
    audit["meta"]["max_thermal_source_power_W"] = float(
        np.max(np.abs(res["thermal_source_power"])))
    audit["meta"]["max_nozzle_power_W"] = float(
        np.max(np.abs(res["nozzle_scalar_flux_power"])))
    Path(out).write_text(json.dumps(audit, indent=1, default=float))
    print(json.dumps(audit["meta"]))


if __name__ == "__main__":
    main()
