"""Step-1 overlay plot: experimental vs sim aft-port G + pressure."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

snaps = list(np.load(r"c:/tmp/_chunc_snaps.npy", allow_pickle=True))
A_t = np.pi/4*0.017018034036068077**2; cstar = 1572.2730597254458
df = pd.read_excel(r"C:/Users/ejarocki/Downloads/ThomasMach5_edited.xlsx", sheet_name='thomasmach5')
te = df['time(ms)'].to_numpy(float)/1e3; Pe = df['Chamber Pressure - Measured (psi)'].to_numpy(float)*6894.76

t = np.array([s['t'] for s in snaps])
P = np.array([s['P'] for s in snaps]); rho=np.array([s['rho'] for s in snaps])
u=np.array([s['u'] for s in snaps]); Dp=np.array([s['D_port'] for s in snaps])
rer=np.array([s['r_erosive'] for s in snaps])
isg=snaps[0]['is_grain']; gidx=np.flatnonzero(isg); ia=gidx[-1]; xg=snaps[0]['x']
P_head=P[:,0]; Ap_aft=np.pi/4*Dp[:,ia]**2
Pe_on=np.interp(t,te,Pe,left=np.nan,right=np.nan)
G_exp=Pe_on*A_t/(cstar*Ap_aft); G_simqs=P_head*A_t/(cstar*Ap_aft)
ie=np.array([gidx[np.argmax(rer[k,gidx])] for k in range(len(snaps))])
G_act=np.array([rho[k,ie[k]]*abs(u[k,ie[k]]) for k in range(len(snaps))])

fig,(a1,a2)=plt.subplots(2,1,figsize=(9,8),sharex=True)
m=t<=0.1
a1.plot(t[m]*1e3,P_head[m]/1e6,'b-',label='sim P_head')
a1.plot(te*1e3,Pe/1e6,'k.--',ms=4,label='exp P (ThomasMach5)')
a1.set_ylabel('P (MPa)'); a1.legend(); a1.set_title('Chunc ignition transient — pressure & aft-port mass flux G'); a1.grid(alpha=.3); a1.set_xlim(0,100)
a2.plot(t[m]*1e3,G_act[m],'r-',lw=2,label='sim ACTUAL ρ|u| (erosion sees), peak-ero cell')
a2.plot(t[m]*1e3,G_simqs[m],'b--',label='sim quasi-steady P·A_t/(c*·A_port)')
a2.plot(t[m]*1e3,G_exp[m],'k.-',ms=4,label='EXP quasi-steady (from measured P)')
a2.set_ylabel('G (kg/m²·s)'); a2.set_xlabel('time (ms)'); a2.legend(); a2.grid(alpha=.3)
fig.tight_layout()
out=r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/artifacts/chunc_G_reconstruction.png"
import os; os.makedirs(os.path.dirname(out),exist_ok=True); fig.savefig(out,dpi=110)
print("saved",out)
