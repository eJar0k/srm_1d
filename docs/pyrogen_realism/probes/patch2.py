p='srm_1d/simulation.py'
s=open(p,encoding='utf-8').read()
def rep(old,new):
    global s
    assert s.count(old)==1,(old,s.count(old)); s=s.replace(old,new)
rep("""    probe_ramp_tau,               # uncontained pyrogen A_burn linear ramp time [s] (0 = off)
    probe_hq,                     # [2,N] scratch: recorded h_c, q_conv""","""    probe_ramp_tau,               # uncontained pyrogen A_burn linear ramp time [s] (0 = off)
    probe_cut_t,                  # uncontained pyrogen shut off at t > this [s] (0 = off)
    probe_pcap,                   # uncontained pyrogen burn-law pressure capped at this [Pa] (0 = off)
    probe_hq,                     # [3,N] scratch: rows 0,1 recorded h_c, q_conv; row 2 capped P""")
rep("""            new_remaining = _compute_uncontained_pyrogen_mdot(
                P, a_pyro, n_pyro, rho_pyro, A_burn_eff,""","""            if probe_cut_t > 0.0 and t > probe_cut_t:
                A_burn_eff = 0.0
            P_pyro = P
            if probe_pcap > 0.0:
                for i in range(N):
                    probe_hq[2, i] = P[i] if P[i] < probe_pcap else probe_pcap
                P_pyro = probe_hq[2]
            new_remaining = _compute_uncontained_pyrogen_mdot(
                P_pyro, a_pyro, n_pyro, rho_pyro, A_burn_eff,""")
rep("""        float(probe_hc_scale), float(probe_ramp_tau),
        np.zeros((2, N)),""","""        float(probe_hc_scale), float(probe_ramp_tau),
        float(probe_cut_t), float(probe_pcap),
        np.zeros((3, N)),""")
rep("""    probe_ramp_tau=0.0,
""","""    probe_ramp_tau=0.0,
    probe_cut_t=0.0,
    probe_pcap=0.0,
""")
open(p,'w',encoding='utf-8').write(s)
