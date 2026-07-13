"""Analyze saved Chunc snapshots: reconstruct experimental vs sim aft-port G."""
import numpy as np, pandas as pd

snaps = list(np.load(r"c:/tmp/_chunc_snaps.npy", allow_pickle=True))
D_throat = 0.017018034036068077
A_t = np.pi/4*D_throat**2
cstar = 1572.2730597254458
print(f"A_t={A_t*1e6:.1f} mm^2, cstar={cstar:.1f} m/s")

# experimental high-res trace (MPa, s)
df = pd.read_excel(r"C:/Users/ejarocki/Downloads/ThomasMach5_edited.xlsx", sheet_name='thomasmach5')
te = df['time(ms)'].to_numpy(float)/1e3
Pe = df['Chamber Pressure - Measured (psi)'].to_numpy(float)*6894.76  # Pa

t = np.array([s['t'] for s in snaps])
def arr(key): return np.array([s[key] for s in snaps])  # (nsnap, ncell)
P = arr('P'); rho = arr('rho'); u = arr('u'); Dp = arr('D_port')
rer = arr('r_erosive'); rtot = arr('r_total'); Mach = arr('Mach')
isg = snaps[0]['is_grain']
gidx = np.flatnonzero(isg)
ia = gidx[-1]                       # aft-most grain cell
xg = snaps[0]['x']

P_head = P[:, 0]
A_port = np.pi/4*Dp**2              # (nsnap, ncell)
Ap_aft = A_port[:, ia]

# (1) experimental quasi-steady aft G  (interp exp P onto sim t; uses sim A_port_aft)
Pe_on_t = np.interp(t, te, Pe, left=np.nan, right=np.nan)
G_exp_qs = Pe_on_t * A_t / (cstar * Ap_aft)
# (2) sim quasi-steady aft G
G_sim_qs = P_head * A_t / (cstar * Ap_aft)
# (3) sim ACTUAL local rho*|u| — aft-most cell and peak-erosive cell
G_act_aft = rho[:, ia]*np.abs(u[:, ia])
# peak-erosive grain cell per snapshot
ie = np.array([gidx[np.argmax(rer[k, gidx])] for k in range(len(snaps))])
G_act_peakero = np.array([rho[k, ie[k]]*abs(u[k, ie[k]]) for k in range(len(snaps))])
x_peakero = xg[ie]
rer_peak = np.array([rer[k, ie[k]] for k in range(len(snaps))])
rtot_peak = np.array([rtot[k, ie[k]] for k in range(len(snaps))])
ef = rer_peak/np.maximum(rtot_peak, 1e-12)   # erosive fraction at peak-ero cell

def row(tt):
    k = int(np.argmin(np.abs(t-tt)))
    return dict(t_ms=t[k]*1e3, P_sim=P_head[k]/1e6, P_exp=Pe_on_t[k]/1e6,
               G_exp_qs=G_exp_qs[k], G_sim_qs=G_sim_qs[k],
               G_act_aft=G_act_aft[k], G_act_peakero=G_act_peakero[k],
               x_pe_mm=x_peakero[k]*1e3, Mach_aft=Mach[k, ia],
               ero_frac=ef[k])

print("\n  t_ms  P_sim P_exp | G_exp_qs G_sim_qs G_act_aft G_act_pkE | x_pkE  Mach_aft eroF")
for tt in [0.001,0.003,0.005,0.008,0.01155,0.015,0.02,0.03,0.05,0.1,0.2,0.4]:
    r=row(tt)
    print(f"  {r['t_ms']:6.1f} {r['P_sim']:5.2f} {r['P_exp'] if not np.isnan(r['P_exp']) else -1:5.2f} |"
          f" {r['G_exp_qs'] if not np.isnan(r['G_exp_qs']) else -1:7.1f} {r['G_sim_qs']:7.1f}"
          f" {r['G_act_aft']:8.1f} {r['G_act_peakero']:8.1f} | {r['x_pe_mm']:5.0f} {r['Mach_aft']:6.2f} {r['ero_frac']:5.2f}")

# peak of the actual erosion-driving G and where/when
kpk = np.argmax(G_act_peakero)
print(f"\nPEAK G_act_peakero = {G_act_peakero[kpk]:.1f} kg/m2s @ t={t[kpk]*1e3:.2f} ms, x={x_peakero[kpk]*1e3:.0f} mm")
print(f"  at that instant: P_sim={P_head[kpk]/1e6:.2f} MPa, Mach_aft={Mach[kpk,ia]:.2f}, ero_frac={ef[kpk]:.2f}")
# experimental aft-G ceiling (max over recorded plateau)
mfin = np.isfinite(G_exp_qs)
print(f"\nEXP aft-G ceiling (max over recorded trace) = {np.nanmax(G_exp_qs[mfin]):.1f} kg/m2s "
      f"@ t={t[mfin][np.argmax(G_exp_qs[mfin])]*1e3:.0f} ms")
# plateau sim values (300-400ms)
mp=(t>=0.3)&(t<=0.4)
print(f"SIM plateau (300-400ms): G_act_aft={G_act_aft[mp].mean():.1f}, G_sim_qs={G_sim_qs[mp].mean():.1f}, "
      f"G_exp_qs={np.nanmean(G_exp_qs[mp]):.1f} kg/m2s")
print(f"\nRATIOS at spike (t={t[kpk]*1e3:.1f}ms): G_act_peakero/G_exp_qs={G_act_peakero[kpk]/G_exp_qs[kpk]:.2f}x,"
      f"  G_act_peakero/G_sim_qs={G_act_peakero[kpk]/G_sim_qs[kpk]:.2f}x,"
      f"  G_sim_qs/G_exp_qs={G_sim_qs[kpk]/G_exp_qs[kpk]:.2f}x")
