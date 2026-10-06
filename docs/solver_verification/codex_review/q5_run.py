import sys, json, os, numpy as np, time
sys.path.insert(0, r"C:/tmp/sv/cx")
from scripts.run_fill_configuration_comparison import _variants, _resolve_options, CONFIG_PATH
from scripts.run_numerical_verification import _extract_metrics
from srm_1d.openmotor_adapter import run_from_ric
base = json.loads(CONFIG_PATH.read_text())["cases"]["chunc"]
V = _variants(base["run_options"])
motor = r"C:/tmp/sv/cx/" + base["motor_path"]
V["canon_bpnv"] = dict(V["current_bpnv"], roughness=32e-6, kappa=0.44, T_ignition=756.0, k_solid=0.271, P_cutoff=0.05e6)
cells=int(sys.argv[1]); cfl=float(sys.argv[2]); name=sys.argv[3]
opts = dict(_resolve_options(V[name]), target_propellant_cells=cells, cfl_target=cfl)
t0=time.time()
r, perf, *_ = run_from_ric(motor, verbose=False, **opts)
m = _extract_metrics(r, perf)
t=np.asarray(r["time"]); mm=np.asarray(r["max_mach"]); P=np.asarray(r["P_head"]); ip=int(np.argmax(P))
keys=["actual_cells","steps","peak_pressure_mpa","peak_time_s","max_fill_mach","minimum_pressure_pa","limit_pressure_floor_total_activations","limit_pressure_floor_active_cells","limit_pressure_floor_maximum_raw_deficit_pa","limit_pressure_floor_first_activation_time_s","limit_pressure_floor_last_activation_time_s","limit_pressure_floor_max_duration_x_m"]
out={k:m.get(k) for k in keys}
out.update(variant=name,cfl=cfl,mach_pre_peak=float(mm[t<t[ip]].max()),wall=round(time.time()-t0,1))
print("RES",json.dumps(out))
