import sys, json, numpy as np
sys.path.insert(0, r"C:/tmp/sv/cx")
from scripts.run_fill_configuration_comparison import _variants, _resolve_options, CONFIG_PATH
from scripts.run_numerical_verification import _extract_metrics
from srm_1d.openmotor_adapter import run_from_ric
base = json.loads(CONFIG_PATH.read_text())["cases"]["chunc"]
V = _variants(base["run_options"])
motor = r"C:/tmp/sv/cx/" + base["motor_path"]
V["canon_bpnv"] = dict(V["current_bpnv"], roughness=32e-6, kappa=0.44, T_ignition=756.0, k_solid=0.271, P_cutoff=0.05e6)
for name in sys.argv[1:]:
    opts = dict(_resolve_options(V[name]), target_propellant_cells=100, cfl_target=0.3)
    r, perf, *_ = run_from_ric(motor, verbose=False, **opts)
    m = _extract_metrics(r, perf)
    keys = ["burn_update_interval_steps","ignition_rate_refresh_delay_count","ignition_rate_refresh_delay_max_s","ignition_rate_refresh_delay_median_s","max_abs_source_ledger_mismatch_w","max_abs_nozzle_boundary_mismatch_w","max_abs_between_step_storage_jump_j","sum_abs_between_step_storage_jump_j","max_abs_energy_residual_w","max_fill_mach","minimum_pressure_pa","limit_temperature_ceiling_total_activations","limit_temperature_floor_total_activations","limit_pressure_floor_total_activations","clipping_active_steps","integrated_abs_clipping_energy_j"]
    print(name, json.dumps({k:m.get(k) for k in keys}))
    gaps = r["summary"].get("max_burn_refresh_gap_s")
