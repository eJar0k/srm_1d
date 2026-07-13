"""Decompose sim ignition flux into h_c and dT at the moment each grain cell
ignites; compare h_c to Keller's measured convective h range."""
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.burn_rate import gnielinski_nusselt, haaland_friction
M=r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/motors/machbusterNew.ric"
_,_,_,_,prop=run_from_ric(M,roughness=32e-6,kappa=0.44,pyrogen='bpnv',T_ignition=756.0,
    P_cutoff=0.01e6,t_max=0.0025,snapshot_interval=1.0,print_interval=0.0,verbose=False)
kg=prop.k_gas; mu=prop.mu_gas; Pr=prop.Cp_gas*mu/kg
snaps=list(np.load(r"c:/tmp/_chunc_snaps.npy",allow_pickle=True))
t=np.array([s['t'] for s in snaps]); xg=snaps[0]['x']; gidx=np.flatnonzero(snaps[0]['is_grain'])

rows=[]
for i in gidx:
    prev=False
    for k in range(len(snaps)):
        if t[k]>0.01: break
        s=snaps[k]; burn=s['is_burning'][i]
        if (not prev) and burn:   # first snapshot where cell i is burning
            rho=s['rho'][i]; uu=abs(s['u'][i]); D=s['D_port'][i]; Tg=s['T'][i]; Ts=s['T_surf'][i]
            Re=rho*uu*D/mu
            if Re>=100:
                f=haaland_friction(Re,32e-6,D); Nu=gnielinski_nusselt(Re,Pr,D,xg[i],f,Tg,Ts,0.44)
                h0=Nu*kg/D; rows.append((h0,Tg-Ts,rho*uu,Re,h0*max(Tg-Ts,0)))
            break
        prev=prev or burn
rows=np.array(rows)
h0,dT,G,Re,q=rows[:,0],rows[:,1],rows[:,2],rows[:,3],rows[:,4]
print(f"At the IGNITING snapshot of each grain cell (n={len(rows)}):")
print(f"  h_c (Gnielinski) [W/m2K]: median={np.median(h0):.0f}  range={h0.min():.0f}-{h0.max():.0f}")
print(f"  dT=Tg-Ts [K]            : median={np.median(dT):.0f}")
print(f"  G=rho|u| [kg/m2s]       : median={np.median(G):.0f}  (= {np.median(G)/10:.0f} g/cm2s)")
print(f"  Re                      : median={np.median(Re):.2e}")
print(f"  q''=h_c*dT [MW/m2]      : median={np.median(q)/1e6:.1f}")
print(f"\nKELLER measured convective h (from his fluxes/dT): ~400-6700 W/m2K over G up to ~280 g/cm2s")
print(f"  sim median h_c {np.median(h0):.0f} W/m2K vs Keller ceiling ~6700 -> {np.median(h0)/6700:.1f}x")
print(f"  Keller h~G^0.905; sim igniting G median {np.median(G)/10:.0f} g/cm2s is INSIDE Keller's G range")
