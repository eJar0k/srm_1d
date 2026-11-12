"""TP1 motor table: thrust/performance from P_head (before) vs the nozzle-entrance
stagnation pressure P0_nozzle (after), from ONE run (the dynamics are identical; TP1 only
records an extra diagnostic and changes the post-processing).

  before = ``compute_thrust_series`` fed with P_head; after = ``compute_motor_performance``.
Also: the exit-cell mixture gamma over the first 50 ms, C*_sim with both pressures, the
throat history (erosion), and (pathfin54) the measured thrust curve comparison.

Run from the checkout root, one motor per process, through the guard:
    python -m srm_1d.tools.guarded_run --mem-gb 8 --timeout-min 60 -- \
        python docs/solver_verification/probes/tp1_motor.py chunc OUT.json
motor in {chunc, balls, haseg, pathfin}.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import srm_1d  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402
from srm_1d.nozzle import compute_thrust_series, _throat_aspect_loss  # noqa: E402,F401

_spec = importlib.util.spec_from_file_location(
    "mp_cap_crossmotor", ROOT / "docs" / "v0_7_4" / "probes" / "mp_cap_crossmotor.py")
_probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_probe)

RIC = {'chunc': 'machbusterNew.ric', 'balls': 'BALLSstick.ric', 'haseg': 'hasegawa_a.ric',
       'pathfin': '2025.12.25 54-2800 ST2.0.ric'}
T_MAX = {'chunc': 1.0, 'balls': 2.0, 'haseg': 10.0, 'pathfin': 6.0}
P_CUT = {'chunc': 0.01e6, 'balls': 0.01e6, 'haseg': 0.05e6, 'pathfin': 0.1e6}
PLATEAU_T = {'chunc': (0.4, 0.9), 'balls': (1.0, 1.1), 'haseg': (2.0, 3.0),
             'pathfin': (1.0, 2.0)}


def perf_from_pressure(res, perf, nozzle, prop, P):
    """Thrust / CF / Isp stats for an arbitrary pressure series (same D_throat history)."""
    from srm_1d.propellant import critical_flow_function, R_UNIVERSAL
    t = np.asarray(res['time'])
    D = np.asarray(res['D_throat'])
    rep = prop.representative_tab()
    c_star = np.sqrt(R_UNIVERSAL / rep.molecular_weight * rep.T_flame) / critical_flow_function(rep.gamma)
    th, cf, isp, _ = compute_thrust_series(
        np.asarray(P, float), D, len(t), rep.gamma, nozzle.D_exit,
        nozzle.divergence_losses(), nozzle.efficiency, nozzle.throat_length,
        nozzle.skin_losses(), float(res['P_ambient']), c_star, 9.80665)
    imp = float(np.sum(0.5 * (th[1:] + th[:-1]) * np.diff(t)))
    pk = float(th.max())
    burning = th > 0.05 * pk
    bi = np.where(burning)[0]
    bt = float(t[bi[-1]] - t[bi[0]])
    return dict(thrust=th, total_impulse=imp, peak_thrust=pk, burn_time=bt,
                average_thrust=imp / bt if bt > 0 else 0.0,
                average_Isp=float(np.mean(isp[burning])), average_C_F=float(np.mean(cf[burning])))


def cstar_stats(res, nozzle, P, win):
    t = np.asarray(res['time'])
    mdot = np.asarray(res['massflow'])
    D = np.asarray(res['D_throat'])
    At = 0.25 * np.pi * D * D
    sel = (t >= win[0]) & (t <= win[1]) & (mdot > 1e-9)
    return float(np.median(np.asarray(P)[sel] * At[sel] / mdot[sel]))


def main():
    motor, out = sys.argv[1], sys.argv[2]
    ric = _probe._committed_ric_path(RIC[motor])
    print(f"[provenance] srm_1d from {Path(srm_1d.__file__).parent}")
    kw = dict(t_max=T_MAX[motor], P_cutoff=P_CUT[motor], snapshot_interval=0.05)
    if motor == 'pathfin':
        # defaults: the motor's own igniter block (auto charge), canonical knobs by default
        kw.update(print_interval=0.0, verbose=False)
        res, perf, nozzle, geo, prop = run_from_ric(str(ric), **kw)
    else:
        res, perf, nozzle, geo, prop = run_from_ric(
            str(ric), erosive_model='ma2020', erosive_mp_cap=0.0, **kw, **_probe.CANON)
    t = np.asarray(res['time'])
    Ph = np.asarray(res['P_head'])
    P0 = np.asarray(res['P0_nozzle'])
    T0 = np.asarray(res['T0_nozzle'])
    gam = np.asarray(res['gamma_nozzle'])
    Rg = np.asarray(res['R_nozzle'])
    mdot = np.asarray(res['massflow'])
    summary = res['summary']
    win = PLATEAU_T[motor]
    sel = (t >= win[0]) & (t <= win[1])

    before = perf_from_pressure(res, perf, nozzle, prop, Ph)
    after = dict(thrust=perf['thrust'], total_impulse=perf['total_impulse'],
                 peak_thrust=perf['peak_thrust'], burn_time=perf['burn_time'],
                 average_thrust=perf['average_thrust'], average_Isp=perf['average_Isp'],
                 average_C_F=perf['average_C_F'])
    rep = prop.representative_tab()
    row = dict(motor=motor, n_steps=int(summary['steps']), termination=summary.get('termination'),
               auto_charge_g=float(summary.get('pyrogen_mass_initial', float('nan'))) * 1e3,
               plateau_window=win, c_star_prop=float(perf['c_star']),
               gamma_prop_rep=float(rep.gamma),
               D_throat_initial_mm=float(nozzle.D_throat) * 1e3,
               D_throat_final_mm=float(np.asarray(res['D_throat'])[-1]) * 1e3,
               erosion_coeff=float(nozzle.erosion_coeff), slag_coeff=float(nozzle.slag_coeff),
               P0_le_Phead_all=bool(np.all(P0 <= Ph * (1 + 1e-12))),
               P0_over_Phead_min=float(np.min(P0[Ph > 1e5] / Ph[Ph > 1e5])),
               Phead_over_P0_plateau=float(np.median(Ph[sel] / P0[sel])),
               Phead_over_P0_first50ms_max=float(np.max(Ph[(t < 0.05) & (Ph > 2e5)] / P0[(t < 0.05) & (Ph > 2e5)])) if np.any((t < 0.05) & (Ph > 2e5)) else None,
               Phead_over_P0_max_overall=float(np.max(Ph[Ph > 1e6] / P0[Ph > 1e6])),
               T0_plateau_K=float(np.median(T0[sel])), T_flame_prop=float(rep.T_flame),
               peakP_head_MPa=float(Ph.max()) / 1e6,
               before={k: v for k, v in before.items() if k != 'thrust'},
               after={k: v for k, v in after.items() if k != 'thrust'})
    # independent check of the recorded P0/P_exit: the NB1 probe's reconstruction from the
    # snapshot nearest the window centre (exit-cell T, D_port, gamma, R) and the run's mdot
    snaps = res['snapshots']
    snap = min(snaps, key=lambda s: abs(s['t'] - 0.5 * (win[0] + win[1])))
    A_e = 0.25 * np.pi * float(snap['D_port'][-1]) ** 2
    gm, Rm = float(res['gamma_mix_final'][-1]), float(res['R_mix_final'][-1])
    Pe_h = np.asarray(res['P_exit'])
    Me = float(np.median(mdot[sel] / (A_e * Pe_h[sel] * np.sqrt(gm / (Rm * float(snap['T'][-1]))))))
    row['nb1_style_check'] = dict(M_exit=Me, P0_over_P=float((1 + 0.5 * (gm - 1) * Me ** 2) ** (gm / (gm - 1))),
                                  recorded_P0_over_P_exit=float(np.median(P0[sel] / Pe_h[sel])),
                                  A_e_over_A_t=float(A_e / (0.25 * np.pi * float(np.median(np.asarray(res['D_throat'])[sel])) ** 2)))
    # C*: sim vs prop, with both pressures (plateau median)
    row['cstar_sim_P0_plateau'] = cstar_stats(res, nozzle, P0, win)
    row['cstar_sim_Phead_plateau'] = cstar_stats(res, nozzle, Ph, win)
    row['cstar_sim_P0_avg_burn'] = float(perf['c_star_sim_avg'])
    # mixture c* from the exit cell (gamma, R, T0): c* = sqrt(R T0)/Gamma(gamma)
    from srm_1d.propellant import critical_flow_function
    Gam = np.array([critical_flow_function(float(g)) for g in gam[sel][::max(1, sel.sum() // 200)]])
    cm = np.sqrt(Rg[sel][::max(1, sel.sum() // 200)] * T0[sel][::max(1, sel.sum() // 200)]) / Gam
    row['cstar_mixture_exitcell_plateau'] = float(np.median(cm))
    # Gamma in the first 50 ms (and the rep-tab value the CF uses)
    for lbl, m in (('first50ms', t < 0.05), ('plateau', sel)):
        g = gam[m]
        row['gamma_' + lbl] = dict(min=float(g.min()), median=float(np.median(g)), max=float(g.max()))
    # CF with the mixture gamma vs scalar (first 50 ms pressure > 0.2 MPa), diagnostic only
    # thrust over pressure effective CF*At at plateau
    At = 0.25 * np.pi * np.asarray(res['D_throat']) ** 2
    row['F_over_Phead_At_plateau'] = {
        'before': float(np.median(before['thrust'][sel] / (Ph[sel] * At[sel]))),
        'after_over_Phead': float(np.median(after['thrust'][sel] / (Ph[sel] * At[sel]))),
        'after_over_P0': float(np.median(after['thrust'][sel] / (P0[sel] * At[sel])))}
    row['F_over_Phead_plateau_N_per_MPa'] = {
        'before': float(np.median(before['thrust'][sel] / (Ph[sel] / 1e6))),
        'after': float(np.median(after['thrust'][sel] / (Ph[sel] / 1e6)))}

    if motor == 'pathfin':
        d = np.loadtxt(ROOT / 'static_fire_data' / 'pathfin54.csv', delimiter=',', skiprows=1)
        te, Fe, Pe = d[:, 0], d[:, 1], d[:, 2]
        ie = float(np.sum(0.5 * (Fe[1:] + Fe[:-1]) * np.diff(te)))
        pke = float(Fe.max())
        bi = np.where(Fe > 0.05 * pke)[0]
        bte = float(te[bi[-1]] - te[bi[0]])
        row['measured'] = dict(total_impulse=ie, peak_thrust=pke, burn_time=bte,
                               average_thrust=ie / bte, t_last=float(te[-1]),
                               peak_P_MPa=float(Pe.max()) / 1e6, F_end=float(Fe[-1]), P_end_MPa=float(Pe[-1]) / 1e6)
        # the measured series to a common window: sim thrust over [t0 + te[0], t0 + te[-1]],
        # t0 = sim time at which thrust first exceeds 5 % of peak minus the measured onset
        row['sim_t_first_5pct'] = float(t[np.where(after['thrust'] > 0.05 * after['peak_thrust'])[0][0]])
        row['meas_t_first_5pct'] = float(te[bi[0]])
        # measured thrust/pressure at the plateau: median over F>0.8*Fmax region
        m80 = Fe > 0.8 * pke
        row['measured']['F_over_P_N_per_MPa_plateau'] = float(np.median(Fe[m80] / (Pe[m80] / 1e6)))
        m80s = after['thrust'] > 0.8 * after['peak_thrust']
        row['sim_F_over_Phead_N_per_MPa_plateau_after'] = float(np.median(after['thrust'][m80s] / (Ph[m80s] / 1e6)))
        row['sim_F_over_Phead_N_per_MPa_plateau_before'] = float(np.median(before['thrust'][m80s] / (Ph[m80s] / 1e6)))
        row['F_over_P_At_measured'] = float(np.median(Fe[m80] / Pe[m80]) / float(At[0]))
        # measured gauge vs absolute: first sample pressure (Pa) as a baseline check
        row['measured']['P_first_Pa'] = float(Pe[0])
        row['measured_trace'] = dict(t=te.tolist(), F=Fe.tolist(), P=(Pe / 1e6).tolist())

    np.savez_compressed(out + '.npz', t=t, P_head=Ph, P_exit=np.asarray(res['P_exit']), P0=P0, T0=T0,
                        gamma=gam, R=Rg, mdot=mdot, D_throat=np.asarray(res['D_throat']),
                        F_before=before['thrust'], F_after=after['thrust'])
    idx = np.linspace(0, len(t) - 1, min(2000, len(t))).astype(int)
    row['trace'] = dict(t=t[idx].tolist(), Phead_MPa=(Ph[idx] / 1e6).tolist(),
                        P0_MPa=(P0[idx] / 1e6).tolist(),
                        F_before=before['thrust'][idx].tolist(), F_after=after['thrust'][idx].tolist())
    with open(out, 'w') as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: v for k, v in row.items() if k not in ('trace', 'measured_trace')}, indent=1))


if __name__ == '__main__':
    main()
