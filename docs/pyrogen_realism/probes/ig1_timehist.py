"""Analyse the ig1_case.py npz time histories (Measurement 2 and 3)."""
import sys, numpy as np
np.set_printoptions(linewidth=200, precision=4, suppress=True)
tag = sys.argv[1]; lb = float(sys.argv[2])
z = np.load(f'C:/tmp/sv_ig1/npz/{tag}.npz')
t = z['snap_t']; x = z['snap_x']; dx = x[1]-x[0]; N = len(x)
P = z['snap_P']; T = z['snap_T']; Ts = z['snap_T_surf']; u = z['snap_u']; rho = z['snap_rho']
h = z['snap_probe_h_c']; q = z['snap_probe_q_conv']; qpy = z['snap_pyrogen_surface_heat_flux']
rt = z['snap_r_total']; re = z['snap_r_erosive']; D = z['snap_D_port']; isg = z['snap_is_grain'][0] > .5
ign = z['ignition_time']; th = z['h_time']; mdot = z['h_mdot_ig']; mrem = z['h_m_pyrogen']
m0 = lb*0.45359237
rho_p = 1700.; Aport = np.pi/4*D[0]**2
Ap = Aport.mean()
Lc = m0/(0.6*rho_p*0.5*Ap); ncart = min(N, max(1, int(np.ceil(Lc/dx))))
print(f'{tag}: N={N} dx={dx*1e3:.1f}mm  basket cells 0..{ncart-1} (L_cart={Lc*1e3:.1f} mm)  m0={m0*1e3:.2f} g')
# --- pyrogen output
cum = m0 - mrem
k = np.argmax(mdot)
def tcross(frac):
    return th[np.argmax(cum >= frac*m0)] if (cum >= frac*m0).any() else np.nan
print(f'pyrogen: peak mdot={mdot.max():.4f} kg/s at {th[k]*1e3:.2f} ms; mean mdot (first-to-90%)={0.9*m0/(tcross(.9)-th[np.argmax(mdot>1e-9)]):.4f} kg/s')
print(f'  t10/t50/t90/t99 of gas-solid consumption [ms]: {tcross(.1)*1e3:.2f} {tcross(.5)*1e3:.2f} {tcross(.9)*1e3:.2f} {tcross(.99)*1e3:.2f}')
print(f'  bore P at basket (vol-avg): peak={np.max(P[:, :ncart].mean(1))/1e6:.3f} MPa @ {t[np.argmax(P[:, :ncart].mean(1))]*1e3:.2f} ms; head P peak {P[:,0].max()/1e6:.3f} MPa')
# --- ignition order
isg_idx = np.flatnonzero(isg)
order = isg_idx[np.argsort(ign[isg_idx])]
print('first-ignited cells (cell: t_ign ms):', ', '.join(f'{c}:{ign[c]*1e3:.3f}' for c in order[:8]))
print(f'ignition cascade: first {ign[order[0]]*1e3:.3f} ms, 50% {np.percentile(ign[isg_idx],50)*1e3:.3f}, 90% {np.percentile(ign[isg_idx],90)*1e3:.3f}, last {ign[order[-1]]*1e3:.3f}')
# --- per-cell flux channels over the 6 snapshots preceding ignition for first 4 cells
CAL=41840.
for c in order[:4]:
    ki = np.searchsorted(t, ign[c])
    print(f'\ncell {c} (x={x[c]*1e3:.0f} mm, ign at {ign[c]*1e3:.3f} ms): snapshots before ignition')
    print('  t[ms]   Tgas  Tsurf   u[m/s]  G[kg/m2s]  h_c[W/m2K]  q_conv[MW/m2]  q_pyro[MW/m2]  q_total/Keller_ceiling6.7')
    for kk in range(max(ki-6,0), min(ki+1, len(t))):
        G = rho[kk,c]*abs(u[kk,c])
        print(f'  {t[kk]*1e3:7.3f} {T[kk,c]:6.0f} {Ts[kk,c]:6.0f} {u[kk,c]:8.1f} {G:9.1f}  {h[kk,c]:11.0f} {q[kk,c]/1e6:10.2f} {qpy[kk,c]/1e6:10.2f} {(q[kk,c]+qpy[kk,c])/1e6/6.7:8.2f}')
# peak convective q over all cells before they ignite
qmax = np.array([np.max(np.where(np.arange(len(t))<=np.searchsorted(t, ign[c]), q[:,c], 0)) for c in isg_idx])
hmax = np.array([np.max(np.where(np.arange(len(t))<=np.searchsorted(t, ign[c]), h[:,c], 0)) for c in isg_idx])
print(f'\nover all grain cells (pre-ignition): peak q_conv median {np.median(qmax)/1e6:.1f} MW/m2, max {qmax.max()/1e6:.1f}; peak h_c median {np.median(hmax):.0f}, max {hmax.max():.0f}')
# --- hot gas front
front = np.array([x[np.flatnonzero(T[kk] > 350)].max() if (T[kk] > 350).any() else 0. for kk in range(len(t))])
print('\nhot-gas front (max x with T>350 K) vs t:')
for tt in [1e-4,3e-4,5e-4,1e-3,2e-3,3e-3,5e-3,1e-2,2e-2,3e-2,5e-2]:
    kk = np.searchsorted(t, tt)
    if kk < len(t): print(f'  t={t[kk]*1e3:7.2f} ms  front x={front[kk]*1e3:7.1f} mm  Tgas_head={T[kk,0]:.0f}  P_head={P[kk,0]/1e6:.2f} MPa  n_ign(cum)={(ign[isg_idx]<=t[kk]).sum()}')
# --- near-peak: aft cell G & erosive fraction
kp = np.argmax(P[:,0]); print(f'\nnear P_head peak (snapshot t={t[kp]*1e3:.2f} ms, P_head {P[kp,0]/1e6:.2f} MPa):')
for c in [0, N//4, N//2, 3*N//4, N-6, N-3]:
    G = rho[kp,c]*abs(u[kp,c])
    print(f'  cell {c:3d} x={x[c]*1e3:6.0f}mm  P={P[kp,c]/1e6:6.2f}  u={u[kp,c]:7.1f}  G={G:8.1f} kg/m2s  r_tot={rt[kp,c]*1e3:6.2f} mm/s  r_ero={re[kp,c]*1e3:6.2f}  ero_frac={re[kp,c]/max(rt[kp,c],1e-12):.2f}')
# plateau reference
kpl = np.searchsorted(t, t[-1]*0.95)
print(f'plateau ref (t={t[kpl]*1e3:.1f} ms): aft cell G={rho[kpl,N-6]*abs(u[kpl,N-6]):.0f}, ero_frac={re[kpl,N-6]/max(rt[kpl,N-6],1e-12):.2f}, P_head={P[kpl,0]/1e6:.2f}')
