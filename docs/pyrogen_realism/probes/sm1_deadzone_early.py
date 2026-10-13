# Reviewer probe (2026-10-12): first 12 ms of the BALLSstick dead zone, per cell u / T_gas / T_surf / Peretz h.
# Run sm1_deadzone_probe.py first (it writes <outdir>/BALLSstick.ric). Usage: python sm1_deadzone_early.py <outdir>
import sys, dataclasses, numpy as np
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
from srm_1d.simulation import _peretz_heat_transfer_coeff
out = sys.argv[1]; ric = f'{out}/BALLSstick.ric'
pyro = dataclasses.replace(load_pyrogen('bpnv'), heat_delivery_mode='none')
res = run_from_ric(ric, pyrogen=pyro, pyrogen_mass=4.8e-3, t_max=0.03, P_cutoff=0.01e6,
                   roughness=32e-6, kappa=0.44, T_ignition=756.0,
                   snapshot_interval=0.001, print_interval=0.0, verbose=False)[0]
snaps = res['snapshots']; x = snaps[0]['x']
At = np.pi/4*res['D_throat'][0]**2
cells = [3, 6, 9, 12, 15, 23, 53]
R, Cp, g = 8.314/0.0482, 2000.0, 1.17   # approximate mixture for an h estimate
print('h_peretz estimated with approximate R/Cp/gamma (relative values only)')
for s in snaps:
    t = s['t']
    if t > 0.012: break
    row = []
    for i in cells:
        A = np.pi/4*s['D_port'][i]**2
        h = _peretz_heat_transfer_coeff(s['P'][i], s['u'][i], s['T'][i], s['T_surf'][i], R, Cp, g,
                                        s['D_port'][i], x[i], A, At)
        row.append(f"{s['u'][i]:6.1f}/{s['T'][i]:5.0f}/{s['T_surf'][i]:4.0f}/{h:5.0f}")
    print(f"t={t*1e3:5.1f}ms P0={s['P'][0]/1e6:5.2f}  " + '  '.join(row))
print('cols: cells', cells, ' each = u[m/s]/T_gas[K]/T_surf[K]/h[W/m2K]')
print('Ap/At:', [round(np.pi/4*snaps[0]['D_port'][i]**2/At, 2) for i in cells])
