"""IG1 single-case runner. Usage: ig1_case.py '<json spec>' <out.jsonl> [npz_dir]

spec keys: topo (head_basket|forward_plenum), lb (pyrogen mass in lb), t_max,
tag, snap (snapshot interval s), plus optional probes: hc_scale, ramp_tau,
erosive_model, T_ign, pyro_over (dict of Pyrogen field overrides e.g. {"a":..,"n":..}).
Resumable: skips if tag already in the jsonl. Appends one JSON line of metrics.
Run from the worktree root with PYTHONPATH=$PWD (see docs/solver_verification/TASKS.md).
"""
import sys, json, time, os, dataclasses
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import srm_1d
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen

spec = json.loads(sys.argv[1]); out = Path(sys.argv[2])
npz_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else None
tag = spec['tag']
if out.exists():
    for ln in out.read_text().splitlines():
        if ln.strip() and json.loads(ln).get('tag') == tag:
            print('skip', tag); sys.exit(0)

LB = 0.45359237
ric = spec.get('ric', 'motors/machbusterNew.ric')
kw = dict(roughness=32e-6, kappa=0.44, T_ignition=spec.get('T_ign', 756.0),
          print_interval=0.0, verbose=False, t_max=spec.get('t_max', 0.2),
          snapshot_interval=spec.get('snap', 1e-3))
if 'hc_scale' in spec: kw['probe_hc_scale'] = spec['hc_scale']
if 'ramp_tau' in spec: kw['probe_ramp_tau'] = spec['ramp_tau']
if 'cut_t' in spec: kw['probe_cut_t'] = spec['cut_t']
if 'pcap' in spec: kw['probe_pcap'] = spec['pcap']
if 'area23' in spec: kw['probe_area23'] = spec['area23']
if 'erosive_model' in spec: kw['erosive_model'] = spec['erosive_model']
pyro = spec.get('pyrogen', 'bpnv')
if spec.get('pyro_over'):
    from srm_1d.openmotor_adapter import load_pyrogen
    p = load_pyrogen(pyro)
    pyro = dataclasses.replace(p, **spec['pyro_over'])
t0 = time.time()
res, perf, nozzle, geo, prop = run_from_ric(
    str(ROOT / ric), pyrogen=pyro, pyrogen_mass=spec['lb'] * LB,
    injection_topology=spec['topo'], **kw)
wall = time.time() - t0
t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6
s = res['summary']
snap0 = res['snapshots'][0]
isg = snap0['is_grain']
ig = np.asarray(res['ignition_time_by_cell'])[isg]
fin = ig[ig < 1e8]
frac = fin.size / max(ig.size, 1)
k = int(P.argmax())
m = dict(tag=tag, **{kk: spec[kk] for kk in spec if kk != 'tag'},
         peak_MPa=float(P[k]), t_peak_s=float(t[k]), P_end_MPa=float(P[-1]),
         frac_ignited=float(frac),
         t_first_ign_s=float(fin.min()) if fin.size else None,
         t50_s=float(np.percentile(ig, 50)) if frac >= .5 else None,
         t90_s=float(np.percentile(ig, 90)) if frac >= .9 else None,
         pyro_duration_s=float(s['pyrogen_duration']),
         pyro_peak_P_MPa=float(s['pyrogen_peak_P']) / 1e6,
         pyro_mass_burned_g=float(s['pyrogen_mass_burned']) * 1e3,
         steps=int(s['steps']), wall_s=wall, n_grain=int(ig.size))
with open(out, 'a') as f:
    f.write(json.dumps(m) + '\n')
if npz_dir is not None:
    npz_dir.mkdir(parents=True, exist_ok=True)
    keys = ['t', 'P', 'u', 'Mach', 'T', 'r_total', 'r_erosive', 'D_port', 'C_burn',
            'is_burning', 'is_grain', 'T_surf', 'mass_source', 'pyrogen_surface_heat_flux',
            'rho', 'probe_h_c', 'probe_q_conv', 'x']
    sn = {kk: np.array([sp[kk] for sp in res['snapshots']]) for kk in keys if kk != 'x'}
    sn['x'] = np.asarray(res['snapshots'][0]['x'])
    hist = {h: np.asarray(res[h]) for h in ['time', 'P_head', 'P_ig', 'mdot_ig', 'm_pyrogen',
            'n_ignited', 'pyrogen_surface_heat_power', 'max_mach', 'dt', 'convective_wall_sink_power']}
    np.savez_compressed(npz_dir / f'{tag}.npz', ignition_time=np.asarray(res['ignition_time_by_cell']),
                        **{'snap_' + a: b for a, b in sn.items()}, **{'h_' + a: b for a, b in hist.items()})
print(json.dumps(m))
