p='srm_1d/simulation.py'
s=open(p,encoding='utf-8').read()
def rep(old,new,count=1):
    global s
    assert s.count(old)==count,(old,s.count(old))
    s=s.replace(old,new)
rep("N_SNAP_CHANNELS = 19","_SNAP_PROBE_H = 19         # IG1 probe: bare convective h_c [W/m^2K] of heated unignited cells\n_SNAP_PROBE_Q = 20         # IG1 probe: convective q'' = h_c*(T_gas-T_surf) [W/m^2]\nN_SNAP_CHANNELS = 21")
rep("""    rho_gas,                 # [N] cell gas density [kg/m^3]
    Cp_mix_gas,              # [N] cell mixture Cp [J/(kg K)]
):""","""    rho_gas,                 # [N] cell gas density [kg/m^3]
    Cp_mix_gas,              # [N] cell mixture Cp [J/(kg K)]
    probe_hc_scale,          # IG1 probe: multiplier on the IGNITION-heating h_c (1.0 = off)
    probe_hq,                # IG1 probe: [2,N] recorded (unscaled h_c, q_conv); diagnostic only
):""")
rep("""    for i in range(N):
        pyrogen_surface_heat_flux[i] = 0.0
        radiation_heat_flux[i] = 0.0
        radiation_sink_power[i] = 0.0
""","""    for i in range(N):
        probe_hq[0, i] = 0.0
        probe_hq[1, i] = 0.0
        pyrogen_surface_heat_flux[i] = 0.0
        radiation_heat_flux[i] = 0.0
        radiation_sink_power[i] = 0.0
""")
rep("""                if flame_spread_enabled:
                    h_c = h_c * flame_spread_augment[i]
                h_total = h_c""","""                if flame_spread_enabled:
                    h_c = h_c * flame_spread_augment[i]
                probe_hq[0, i] = h_c
                probe_hq[1, i] = h_c * (T[i] - T_surf[i])
                if probe_hc_scale != 1.0:
                    h_c = h_c * probe_hc_scale
                h_total = h_c""")
rep("""    progress_state,               # [2] float64: [0]=progress 0..1 (written),
                                  #              [1]=cancel flag (read; >0.5 aborts)
):""","""    progress_state,               # [2] float64: [0]=progress 0..1 (written),
                                  #              [1]=cancel flag (read; >0.5 aborts)
    # --- IG1 diagnosis probes (defaults: off) ---
    probe_hc_scale,               # multiplier on ignition-heating h_c (1.0 = off)
    probe_ramp_tau,               # uncontained pyrogen A_burn linear ramp time [s] (0 = off)
    probe_hq,                     # [2,N] scratch: recorded h_c, q_conv
):""")
rep("""            new_remaining = _compute_uncontained_pyrogen_mdot(
                P, a_pyro, n_pyro, rho_pyro, A_burn_per_cell,""","""            A_burn_eff = A_burn_per_cell
            if probe_ramp_tau > 0.0 and t < probe_ramp_tau:
                A_burn_eff = A_burn_per_cell * (t / probe_ramp_tau)
            new_remaining = _compute_uncontained_pyrogen_mdot(
                P, a_pyro, n_pyro, rho_pyro, A_burn_eff,""")
rep("""            rho, Cp_mix_arr,
        )""","""            rho, Cp_mix_arr,
            probe_hc_scale, probe_hq,
        )""")
rep("""                snap_data[snap_idx, _SNAP_RHO, i] = rho[i]
""","""                snap_data[snap_idx, _SNAP_RHO, i] = rho[i]
                snap_data[snap_idx, _SNAP_PROBE_H, i] = probe_hq[0, i]
                snap_data[snap_idx, _SNAP_PROBE_Q, i] = probe_hq[1, i]
""")
rep("""        progress_state,
    )

    wall_elapsed""","""        progress_state,
        float(probe_hc_scale), float(probe_ramp_tau),
        np.zeros((2, N)),
    )

    wall_elapsed""")
rep("""            'rho': snap_data[s, _SNAP_RHO, :].copy(),
        })""","""            'rho': snap_data[s, _SNAP_RHO, :].copy(),
            'probe_h_c': snap_data[s, _SNAP_PROBE_H, :].copy(),
            'probe_q_conv': snap_data[s, _SNAP_PROBE_Q, :].copy(),
        })""")
rep("""    progress_state=None,
    # Limiter""","""    progress_state=None,
    # IG1 diagnosis probes (probe-only, default off; not for production use)
    probe_hc_scale=1.0,
    probe_ramp_tau=0.0,
    # Limiter""")
open(p,'w',encoding='utf-8').write(s)
