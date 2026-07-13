"""PROBE sweep: transient h_c knock-down sensitivity (Chunc).
Scales the convective Gnielinski h_c (ignition gate + Ma erosive, shared Nu*k/D)
by `factor` while t < `window`. Measures spike vs plateau and snap-on time.
Target: experimental = NO spike (ratio ~1.0), plateau ~8.8 MPa.
"""
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric
M=r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/motors/machbusterNew.ric"

def run(factor, window):
    r,_,_,_,_=run_from_ric(M,roughness=32e-6,kappa=0.44,pyrogen='bpnv',T_ignition=756.0,
        P_cutoff=0.01e6,t_max=0.5,snapshot_interval=0.002,print_interval=0.0,verbose=False,
        h_establish_factor=factor, h_establish_window_s=window)
    s=r['summary']
    snaps=r['snapshots']; isg=snaps[0]['is_grain']; ng=isg.sum()
    tfull=np.nan
    for sn in snaps:
        if (sn['is_burning'][isg]).sum()==ng: tfull=sn['t']*1e3; break
    return s['P_peak']/1e6, s['t_peak']*1e3, s['P_mid']/1e6, tfull

print(f"{'factor':>6} {'window_ms':>9} | {'P_peak':>7} {'t_peak_ms':>9} {'P_mid':>6} {'ratio':>6} {'fullIgn_ms':>10}")
print("-"*70)
# baseline
for factor,window in [(1.0,0.0)]:
    Pp,tp,Pm,tf=run(factor,window)
    print(f"{factor:6.2f} {window*1e3:9.0f} | {Pp:7.2f} {tp:9.2f} {Pm:6.2f} {Pp/Pm:6.3f} {tf:10.1f}  <-- BASELINE")
for window in [0.005,0.010,0.020]:
    for factor in [0.5,0.3,0.15]:
        Pp,tp,Pm,tf=run(factor,window)
        print(f"{factor:6.2f} {window*1e3:9.0f} | {Pp:7.2f} {tp:9.2f} {Pm:6.2f} {Pp/Pm:6.3f} {tf:10.1f}")
print("\nExperimental target: ratio ~1.0 (NO spike), plateau ~8.8 MPa, pressurize ~10 ms")
