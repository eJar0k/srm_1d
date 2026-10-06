import sys, json, time, numpy as np
sys.path.insert(0, r"C:/tmp/sv/cx")
import srm_1d, motorlib
print("srm_1d", srm_1d.__file__); print("motorlib", motorlib.__file__)
from pathlib import Path
from scripts.run_fill_configuration_comparison import _variants, _resolve_options, CONFIG_PATH
from srm_1d.openmotor_adapter import run_from_ric
base = json.loads(CONFIG_PATH.read_text())["cases"]["chunc"]
V = _variants(base["run_options"])
motor = r"C:/tmp/sv/cx/" + base["motor_path"]
cells = 100
which = sys.argv[1:] or ["current_bpnv"]
# extra: canonical knobs (v0.8 defaults) + bpnv
V["canon_bpnv"] = dict(V["current_bpnv"], roughness=32e-6, kappa=0.44, T_ignition=756.0, k_solid=0.271, P_cutoff=0.05e6)
for name in which:
    opts = dict(_resolve_options(V[name]), target_propellant_cells=cells, cfl_target=float(__import__("os").environ.get("CFL","0.3")))
    t0=time.time()
    r, perf, *_ = run_from_ric(motor, verbose=False, **opts)
    t = np.asarray(r["time"]); m = np.asarray(r["max_mach"]); P=np.asarray(r["P_head"])
    ip = int(np.argmax(P))
    w3 = t<=0.003; wp = t<t[ip]
    im = int(np.argmax(m))
    print(json.dumps(dict(variant=name, cells=cells, steps=len(t), wall=round(time.time()-t0,1),
        mach_le_3ms=float(m[w3].max()), mach_pre_peak=float(m[wp].max()), mach_all_50ms=float(m.max()),
        t_at_max_mach=float(t[im]), t_peak=float(t[ip]), P_peak_MPa=float(P[ip]/1e6),
        pyro=str(opts['pyrogen'] if isinstance(opts['pyrogen'],str) else opts['pyrogen'].name))))
