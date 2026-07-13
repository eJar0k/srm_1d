"""Survey Chunc igniter physics: sizing, impetus, gas generation vs literature."""
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
from srm_1d.propellant import R_UNIVERSAL
M=r"c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/motors/machbusterNew.ric"

# ---- 1. impetus check (gas-generation potency per unit mass) ----
pyro = load_pyrogen('bpnv')
Rs = R_UNIVERSAL/pyro.M
W_si = Rs*pyro.T_flame                       # J/kg  (impetus = R_specific*T_flame)
J_to_psi_in3 = 6894.76*1.6387064e-5          # 1 J = Pa*m3 -> psi*in3 per gram needs /1e-3
W_psi_in3_g = W_si*1e-3 / (6894.76*1.6387064e-5)  # J/g divided by (psi*in3 per J)... see below
# clean conversion: 1 J = 1/(6894.76*1.6387e-5) psi*in3 = 8.851 psi*in3 ; per gram: J/g * 8.851
W_psi_in3_g = (W_si/1000.0)*8.851
print("=== IMPETUS / FORCE CONSTANT (gas potency per gram) ===")
print(f"bpnv: M={pyro.M*1000:.1f} g/mol, T_flame={pyro.T_flame:.0f} K, gamma={pyro.gamma}")
print(f"  R_specific = {Rs:.1f} J/kgK ; W = R*T = {W_si/1000:.0f} J/g = {W_psi_in3_g:.0f} psi*in3/g")
print(f"  YAML cites DeMar literature impetus_W = {pyro.impetus_W:.0f} psi*in3/g")
print(f"  --> sim impetus is {W_psi_in3_g/pyro.impetus_W:.2f}x DeMar literature")
M_match = R_UNIVERSAL*pyro.T_flame/(pyro.impetus_W/8.851*1000)
print(f"  to match DeMar 5000 at T={pyro.T_flame:.0f}K: M would be {M_match*1000:.1f} g/mol (vs {pyro.M*1000:.1f})")

# ---- 2. igniter sizing for Chunc ----
r,perf,nozzle,geo,prop = run_from_ric(M,roughness=32e-6,kappa=0.44,pyrogen='bpnv',T_ignition=756.0,
    P_cutoff=0.01e6,t_max=0.5,snapshot_interval=0.2,print_interval=0.0,verbose=False)
s=r['summary']
print("\n=== igniter sizing (Sutton auto) ===")
print(f"  pyrogen_mass_initial = {s['pyrogen_mass_initial']*1000:.3f} g  (Sutton 0.12*V_F^0.7)")
print(f"  pyrogen_peak_P(plenum) = {s.get('pyrogen_peak_P',float('nan'))/1e6:.1f} MPa")
print(f"  pyrogen_duration = {s.get('pyrogen_duration',float('nan'))*1e3:.1f} ms")
print(f"  pyrogen_mass_burned = {s.get('pyrogen_mass_burned',0)*1000:.3f} g")

# ---- 3. igniter gas-generation history ----
t=np.array(r['time']); Ph=np.array(r['P_head'])
print("\n=== result keys with igniter data ===")
print([k for k in r.keys() if any(w in k.lower() for w in ('ig','pyro','mdot','plenum'))])
for key in ('mdot_ig','P_ig','m_pyrogen','T_ig'):
    if key in r:
        a=np.array(r[key]); print(f"  {key}: max={np.nanmax(a):.4g}, len={len(a)}")
