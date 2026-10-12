"""H1 check vs Keller-Baer-Ryan 1966 (AIAA J. 4(8):1358-1364, DOI 10.2514/3.3680; docs/references/keller1966.pdf p.3 Eq. 11):
  h = 1.435e-4 * Tg^0.3 * G^0.905  cal/(cm^2 s K)   (air/O2/N2; Tg [K], G [g/cm^2 s]); average deviation 6%
Compare to the sim's bare Gnielinski h_c (probe channel) at the last pre-ignition snapshot of the first-ignited cells.
Keller validated flux band: 20-160 cal/cm2 s = 0.84-6.7 MW/m2 (Eq. 17 in the same paper)."""
import numpy as np, sys
CAL=41840.
def keller_h(Tg,G_kgm2s):
    return 1.435e-4*Tg**0.3*(G_kgm2s*0.1)**0.905*CAL   # kg/m2s -> g/cm2s = x0.1
for tag in ['th_hb_0.01','th_hb_0.001']:
    z=np.load(f'C:/tmp/sv_ig1/npz/{tag}.npz'); t=z['snap_t']; ign=z['ignition_time']; isg=z['snap_is_grain'][0]>.5
    idx=np.flatnonzero(isg); order=idx[np.argsort(ign[idx])][:8]
    print(f'\n{tag}: first 8 ignited cells, last snapshot with h_c>0 before ignition')
    print(' cell  t_ign[ms]  Tgas[K]  G[kg/m2s]  h_sim[W/m2K]  h_Keller  ratio  q_sim[MW/m2]  q_Keller=hK*dT[MW/m2]')
    rs=[]; qs=[]
    for c in order:
        ks=np.flatnonzero(z['snap_probe_h_c'][:,c]>0); k=ks[-1]
        Tg=z['snap_T'][k,c]; G=z['snap_rho'][k,c]*abs(z['snap_u'][k,c]); hs=z['snap_probe_h_c'][k,c]
        hk=keller_h(Tg,G); dT=Tg-z['snap_T_surf'][k,c]
        rs.append(hs/hk); qs.append(z['snap_probe_q_conv'][k,c]/1e6)
        print(f' {c:4d} {ign[c]*1e3:9.3f} {Tg:8.0f} {G:9.1f} {hs:12.0f} {hk:9.0f} {hs/hk:6.2f} {z["snap_probe_q_conv"][k,c]/1e6:9.2f} {hk*dT/1e6:12.2f}')
    print(f' median h ratio sim/Keller = {np.median(rs):.2f}; median q_sim = {np.median(qs):.2f} MW/m2 (Keller band ceiling 6.7)')
    # all grain cells, at their peak pre-ignition q
    rr=[]
    for c in idx:
        ks=np.flatnonzero((z['snap_probe_h_c'][:,c]>0)&(t<=ign[c]+1e-9))
        if ks.size==0: continue
        k=ks[np.argmax(z['snap_probe_q_conv'][ks,c])]
        rr.append(z['snap_probe_h_c'][k,c]/keller_h(z['snap_T'][k,c],z['snap_rho'][k,c]*abs(z['snap_u'][k,c])))
    print(f' all grain cells (at own peak pre-ignition q): h ratio median {np.median(rr):.2f} range {min(rr):.2f}-{max(rr):.2f}')
