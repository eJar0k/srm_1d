import sys; sys.path.insert(0, __import__('os').getcwd())
import numpy as np
from srm_1d.solver import _piso_step_with_energy_diagnostics as f
N=20; dx=0.02; A=1e-3; R=300.; Cp=1200.; g=1.25
for Tsrc in (600.,3000.):
 for dt in (1e-6,1e-5):
  for S in (False,True):
   for K in (0,1,3):
    P=np.full(N,2e5); T=np.full(N,600.); rho=P/(R*T); u=np.zeros(N+1)
    src=np.zeros(N); src[0:3]=5.0
    th=src*Cp*Tsrc
    out=f(rho,u,P,T,np.full(N,A),np.full(N,0.035),src,th,np.zeros(N+1),np.full(N,0.02),dx,dt,np.full(N,g),np.full(N,R),np.full(N,Cp),np.full(N,3500.),1e-4,101325.,300.,N,0.0,K,1e-9,S)
    mb,ma,eos,mdot=out[16],out[17],out[18],out[20]
    srcm=np.sum(src)*dx
    carried=np.sum(out[0])*A*dx
    print(f"Tsrc {Tsrc:.0f} dt {dt:g} S {int(S)} K {K}: expected {srcm-mdot:.5f} ledger {(ma-mb)/dt:.5f} carried {(carried-mb)/dt:.5f} eos/dt {eos/dt:.5f} p2/dt {out[23]/dt:.5f} passes {out[21]} relT {out[22]:.1e}")
