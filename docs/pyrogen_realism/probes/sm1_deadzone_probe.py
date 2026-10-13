# Reviewer probe (2026-10-12): BALLSstick dead zone, 4.8 g BPNV, DeMar off.
# Usage: python sm1_deadzone_probe.py <outdir> <charge_g>  (writes the committed BALLSstick.ric copy to <outdir>)
import sys, subprocess, dataclasses, numpy as np
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
out = sys.argv[1]
blob = subprocess.run(['git', 'show', 'HEAD:motors/BALLSstick.ric'], capture_output=True).stdout
ric = f'{out}/BALLSstick.ric'; open(ric, 'wb').write(blob)
pyro = dataclasses.replace(load_pyrogen('bpnv'), heat_delivery_mode='none')
res = run_from_ric(ric, pyrogen=pyro, pyrogen_mass=float(sys.argv[2])*1e-3, t_max=2.0,
                   P_cutoff=0.01e6, roughness=32e-6, kappa=0.44, T_ignition=756.0,
                   snapshot_interval=0.01, print_interval=0.0, verbose=False)[0]
ign = np.asarray(res['ignition_time_by_cell']); snaps = res['snapshots']
g = snaps[0]['is_grain']; gi = np.where(g)[0]
never = [i for i in gi if not (ign[i] >= 0)]
print('grain cells', gi[0], '..', gi[-1], ' never lit:', never)
x = snaps[0]['x']
show = [gi[0], gi[0]+3, gi[0]+6, gi[0]+9, gi[0]+12, gi[0]+20, gi[0]+50]
print('cell x[mm] D_port[mm] ign_t[ms]')
for i in show: print(f'{i:4d} {x[i]*1e3:7.1f} {snaps[0]["D_port"][i]*1e3:6.2f} {ign[i]*1e3:8.2f}')
for s in snaps:
    if s['t'] in () : pass
for s in snaps:
    t = s['t']
    if not (abs(t-0.02)<6e-3 or abs(t-0.05)<6e-3 or abs(t-0.1)<6e-3 or abs(t-0.3)<6e-3 or abs(t-1.0)<6e-3 or abs(t-1.8)<6e-3): continue
    print(f'\nt={t*1e3:.0f} ms  P_head={s["P"][0]/1e6:.2f} MPa')
    print(' cell  u[m/s]  T_gas[K]  T_surf[K]  burning')
    for i in [0, 1, 2] + show:
        print(f' {i:4d} {s["u"][i]:7.2f} {s["T"][i]:8.0f} {s["T_surf"][i]:8.0f}   {bool(s["is_burning"][i])}')
