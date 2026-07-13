"""Step 2 — sim convective ignition flux q'' vs Keller-Baer-Ryan 1966 regime.

Keller (convective AP/HTPB ignition, the motor-startup analog) validated the
surface-T criterion over Fs = 20-160 cal/cm2 s = 0.84-6.7 MW/m2, with
T_ign = 300 + 286.1*Fs^0.08 (K, Fs in cal/cm2 s).
Question: is the sim's Goodman/Gnielinski convective flux during the fill
INSIDE Keller's validated band, or far above it (=> steady Gnielinski Nu
over-predicting the transient developing-fill flux)?
"""
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.burn_rate import gnielinski_nusselt, haaland_friction

MOTOR = r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/motors/machbusterNew.ric"
# quick build to recover gas transport (numba cached -> fast)
_,_,_,_,prop = run_from_ric(MOTOR, roughness=32e-6, kappa=0.44, pyrogen='bpnv',
    T_ignition=756.0, P_cutoff=0.01e6, t_max=0.0025, snapshot_interval=1.0,
    print_interval=0.0, verbose=False)
k_gas=prop.k_gas; mu=prop.mu_gas; Cp=prop.Cp_gas; Pr=Cp*mu/k_gas
print(f"gas: k={k_gas:.4f} W/mK, mu={mu:.3e} Pa s, Cp={Cp:.1f}, Pr={Pr:.3f}")
CAL=41840.0  # W/m2 per cal/cm2 s

snaps=list(np.load(r"c:/tmp/_chunc_snaps.npy",allow_pickle=True))
t=np.array([s['t'] for s in snaps]); xg=snaps[0]['x']; isg=snaps[0]['is_grain']; gidx=np.flatnonzero(isg)
rough=32e-6; kappa=0.44

def qconv(s,i):
    rho=s['rho'][i]; uu=abs(s['u'][i]); D=s['D_port'][i]; Tg=s['T'][i]; Ts=s['T_surf'][i]
    Re=rho*uu*D/mu
    if Re<100: return 0.0,Re
    f=haaland_friction(Re,rough,D)
    Nu=gnielinski_nusselt(Re,Pr,D,xg[i],f,Tg,Ts,kappa)
    h0=Nu*k_gas/D
    return h0*max(Tg-Ts,0.0), Re

# peak convective flux each grain cell sees during 0-5 ms, and flux at its ignition
early=[k for k in range(len(snaps)) if t[k]<=0.006]
peakq=np.zeros(len(gidx)); ign_q=np.full(len(gidx),np.nan); ign_t=np.full(len(gidx),np.nan); Re_at=np.zeros(len(gidx))
for gi,i in enumerate(gidx):
    prev_burn=False
    for k in early:
        q,Re=qconv(snaps[k],i)
        if q>peakq[gi]: peakq[gi]=q; Re_at[gi]=Re
        if (not prev_burn) and snaps[k]['is_burning'][i]:
            ign_q[gi]=q; ign_t[gi]=t[k];
            prev_burn=True
        prev_burn=prev_burn or snaps[k]['is_burning'][i]

pk=peakq/1e6  # MW/m2
print(f"\n=== SIM convective ignition flux (grain cells, 0-6 ms) ===")
print(f"peak q'' per cell [MW/m2]: min={pk.min():.1f}  median={np.median(pk):.1f}  "
      f"mean={pk.mean():.1f}  max={pk.max():.1f}  p90={np.percentile(pk,90):.1f}")
print(f"peak q'' in cal/cm2 s    : median={np.median(peakq)/CAL:.0f}  max={peakq.max()/CAL:.0f}")
print(f"Re at peak flux          : median={np.median(Re_at):.2e}  max={Re_at.max():.2e}")
print(f"\nKELLER validated band    : 0.84-6.7 MW/m2 (20-160 cal/cm2 s), Re ~ shock-tube BL")
over=np.mean(pk>6.7)*100
print(f"--> {over:.0f}% of grain cells see peak q'' ABOVE Keller's 6.7 MW/m2 ceiling")
print(f"--> sim median {np.median(pk):.1f} MW/m2 = {np.median(pk)/6.7:.1f}x Keller's ceiling, "
      f"{np.median(pk)/0.84:.0f}x Keller's floor")

# Keller T_ign at the sim's flux level (Eq 17): does 756 K still hold?
Fs_cal=np.median(peakq)/CAL
Tign_keller=300+286.1*Fs_cal**0.08
print(f"\nKeller Eq.17 T_ign at sim median flux ({Fs_cal:.0f} cal/cm2 s) = {Tign_keller:.0f} K "
      f"(extrapolated; sim uses 756 K)")

# implied surface-T-threshold ignition time at sim flux vs Keller flux (inert semi-inf solid)
# t_ign = (pi/alpha)*[k_s*(Tign-T0)/(2*q)]^2 ; alpha from prop
ks=prop.k_solid; rho_p=prop.rho_propellant; Cps=prop.Cps; alpha=ks/(rho_p*Cps)
dT=756-293
def tign(q): return (np.pi*ks**2*dT**2)/(4*alpha*q**2)  # s
print(f"\nsolid: k_s={ks:.3f}, rho_p={rho_p:.0f}, Cps={Cps:.0f}, alpha={alpha:.2e} m2/s")
for q_MW in [1.0,3.0,6.7,np.median(pk),pk.max()]:
    print(f"  inert surface-T t_ign at q''={q_MW:5.1f} MW/m2 -> {tign(q_MW*1e6)*1e3:.2f} ms")
print(f"\nsim observed snap-on: whole grain lights in ~0.6-1.5 ms")
