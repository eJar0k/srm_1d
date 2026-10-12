"""H3 check (b): sim pyrogen law in DeMar's closed test vessel vs measured time-to-peak.
DeMar (docs/references/STARTING YOUR UP-GOERS_pub (1).pdf p.7,12): 75mm casing 5.75in x 2.75in ID
(V=34.1 in^3); BKNO3-Viton: peak 105 psi, time to peak 67 ms, impetus 5000 psi in3/g => M=PV/W=0.72 g.
Pure python ODE, no srm_1d import. Gas: sim bpnv (T 2719 K, M 0.0482, gas fraction 0.781), adiabatic."""
import numpy as np
V=np.pi/4*(2.75*.0254)**2*5.75*.0254; print(f'V={V*1e3:.3f} L')
Pamb=101325.; R=8.314/0.0482; T=2719.; fg=.781; rho=1700.
for m0 in [0.72e-3]:
  for label,d in [('3.2mm pellet cyl (sim default)',3.2e-3),('1mm grains',1e-3),('0.3mm powder',3e-4)]:
    A=m0*6/(rho*d)
    for a,n,lab in [(2e-5,.5,'sim a=2e-5 n=0.5')]:
        m=m0; Pv=Pamb; t=0.; dt=1e-5; P=[];ts=[]; mg=0
        # air mass heated? ignore: treat added gas as ideal at T; P = Pamb + mg*R*T/V (isothermal-ish upper bound)
        while m>1e-9 and t<2:
            mdot=rho*A*a*max(Pv,0)**n
            dm=min(mdot*dt,m); m-=dm; mg+=dm*fg
            Pv=Pamb+mg*R*T/V; t+=dt; P.append(Pv); ts.append(t)
        P=np.array(P); k=P.argmax()
        print(f'{label:32s} A={A*1e4:6.1f} cm2 {lab}: burn-out {ts[k]*1e3:7.1f} ms, P_peak {(P[k]-Pamb)/6894.76:6.1f} psi gauge (DeMar: 105 psi in 67 ms)')
