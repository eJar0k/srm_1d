"""EN1 Part 3 -- effect of the energy / momentum formulation on the motors.

One full transient run per (motor, case) writes a JSON row with everything the
EN1 worklog tables need:

  * B0 rows (Chunc / BALLSstick / Hasegawa A, committed .ric, canonical knobs,
    library BPNV, the shipped auto-sized charge) -- peak, t_peak, plateau, ratios,
    n_steps, charge actually used;
  * ignition timing (first light, t90 of the cells that light) and the cells never lit;
  * the Chunc fill Mach map from 1 ms snapshots (max Mach, supersonic extent, the
    downstream end of the supersonic region = shock position, bore-gas T at the first
    lit cell);
  * the T6 energy ledger (EOS reconciliation jump and the other storage jumps, fill phase);
  * wall time and steps/s of the run.

Cases (``--case``):
    sensible   energy_formulation='sensible', momentum_formulation='legacy'   (shipped)
    total      energy_formulation='total',    momentum_formulation='legacy'
    total_cons energy_formulation='total',    momentum_formulation='conservative'
    sens_cons  energy_formulation='sensible', momentum_formulation='conservative'

``--mass`` selects the igniter charge: ``auto`` (the shipped auto-size + upsize
fallback) or a charge in grams (explicit: never resized, so every formulation
gets the identical igniter).

Run from the worktree root with the pyenv 3.10.5 interpreter, e.g.
    python docs/solver_verification/probes/en1_motor_tables.py chunc total --out r.json
"""
import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import srm_1d  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402
from srm_1d.tools.ignition_diagnostics import sensible_enthalpy_audit  # noqa: E402

CANON = dict(roughness=32e-6, kappa=0.44, T_ignition=756.0, print_interval=0.0,
             verbose=False, pyrogen='bpnv')
CASES = {
    'sensible': ('sensible', 'legacy'),
    'total': ('total', 'legacy'),
    'total_cons': ('total', 'conservative'),
    'sens_cons': ('sensible', 'conservative'),
}
# hist = preallocated per-step history rows (default 5e6 rows ~ 2.3 GB, over the 4 GB
# coexist budget when two jobs overlap); sized ~1.4x the step count of each motor.
MOTORS = {
    'chunc': dict(ric='machbusterNew.ric', t_max=1.0, P_cutoff=0.01e6, hist=800_000),
    'balls': dict(ric='BALLSstick.ric', t_max=2.0, P_cutoff=0.01e6, hist=1_300_000),
    'haseg': dict(ric='hasegawa_a.ric', t_max=10.0, P_cutoff=0.05e6, hist=3_300_000),
}


def fill_map(res, t_end=0.010):
    """Per-snapshot Mach statistics of the bore for t <= t_end."""
    rows = []
    for s in res['snapshots']:
        if s['t'] > t_end:
            break
        M = np.abs(np.asarray(s['Mach']))
        x = np.asarray(s['x'])
        dx = float(x[1] - x[0])
        sup = np.flatnonzero(M > 1.0)
        rows.append(dict(
            t_ms=float(s['t'] * 1e3),
            M_max=float(M.max()), x_Mmax_mm=float(x[int(M.argmax())] * 1e3),
            n_supersonic=int(sup.size),
            supersonic_extent_mm=float(sup.size * dx * 1e3),
            x_shock_mm=(float((x[sup.max()] + 0.5 * dx) * 1e3) if sup.size else None),
            T_max=float(np.max(s['T'])), T_min=float(np.min(s['T'])),
            P_max_MPa=float(np.max(s['P']) / 1e6)))
    return rows


def ignition_stats(res):
    it = np.asarray(res['ignition_time_by_cell'], float)
    lit = it < 1.0e9
    out = dict(n_lit=int(lit.sum()), n_cells=int(it.size))
    if lit.any():
        tl = it[lit]
        out.update(t_first_ms=float(tl.min() * 1e3),
                   t90_ms=float(np.percentile(tl, 90) * 1e3),
                   first_cell=int(np.flatnonzero(lit)[int(tl.argmin())]))
    s = res['summary']
    out.update(grain_cells_never_lit=int(s['grain_cells_never_lit']),
               grain_never_lit_cell_range=s.get('grain_never_lit_cell_range'),
               grain_fraction_lit=float(s['grain_fraction_lit']),
               ignited=bool(s.get('ignited', True)))
    return out


def first_cell_temperature(res, ign):
    """Bore-gas T at the first-lit cell at the snapshots around its ignition time."""
    if 'first_cell' not in ign:
        return None
    c, t0 = ign['first_cell'], ign['t_first_ms'] * 1e-3
    snaps = res['snapshots']
    before = [s for s in snaps if s['t'] <= t0]
    after = [s for s in snaps if s['t'] >= t0]
    out = {}
    if before:
        out['T_before_K'] = (float(before[-1]['t'] * 1e3), float(before[-1]['T'][c]))
    if after:
        out['T_after_K'] = (float(after[0]['t'] * 1e3), float(after[0]['T'][c]))
    # max gas T over lit cells at the first snapshot after the first light
    return out


def ledger(res):
    a = sensible_enthalpy_audit(res)
    ph = a['phases']
    out = {'phase_bounds': a['phase_bounds']}
    for name in ('fill', 'spike', 'plateau'):
        b = ph[name]
        out[name] = {k: b[k] for k in b if not isinstance(b[k], (dict, list))}
        out[name].update({k: b[k] for k in b if isinstance(b[k], dict)})
    out['run'] = {k: a['run'][k] for k in a['run']}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('motor', choices=sorted(MOTORS))
    ap.add_argument('case', choices=sorted(CASES))
    ap.add_argument('--mass', default='auto',
                    help="'auto' or the igniter charge in grams")
    ap.add_argument('--t-max', type=float, default=None)
    ap.add_argument('--snap', type=float, default=1.0e-3)
    ap.add_argument('--cells', type=int, default=100)
    ap.add_argument('--out', default=None)
    ap.add_argument('--no-ledger', action='store_true')
    ap.add_argument('--driver', default='static', choices=('static', 'recovery'),
                    help='RT1 ignition driver temperature')
    ap.add_argument('--solid', default=None, choices=(None, 'goodman', 'peretz', 'conduction'),
                    help='solid_heating_model (None = the base default)')
    ap.add_argument('--jet', default='off', choices=('off', 'max', 'replace'),
                    help='HE1 head_end_jet_h')
    ap.add_argument('--hist', type=int, default=None, help='history_capacity override')
    ap.add_argument('--aft', default='lump', help="G1 aft_end: cone|cone_nolump|lump|none")
    ap.add_argument('--head-gap-mm', type=float, default=None, help='G1 head gap override [mm]')
    a = ap.parse_args()
    cfg = MOTORS[a.motor]
    ef, mf = CASES[a.case]
    t_max = a.t_max or cfg['t_max']
    mass = None if a.mass == 'auto' else float(a.mass) * 1e-3
    print(f"[provenance] srm_1d from {Path(srm_1d.__file__).parent}", flush=True)
    t0 = time.time()
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        res = run_from_ric(str(ROOT / 'motors' / cfg['ric']), t_max=t_max,
                           P_cutoff=cfg['P_cutoff'], pyrogen_mass=mass,
                           energy_formulation=ef, momentum_formulation=mf,
                           snapshot_interval=a.snap,
                           **({} if a.hist is None else dict(history_capacity=a.hist)),   # default: MEM1 auto-sized history
                           ignition_driver_temperature=a.driver,
                           head_end_jet_h=a.jet,
                           **({} if a.solid is None else dict(solid_heating_model=a.solid)),
                           aft_end=a.aft, head_gap=(None if a.head_gap_mm is None else a.head_gap_mm * 1e-3),
                           target_propellant_cells=a.cells, **CANON)[0]
    wall = time.time() - t0
    s = res['summary']
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    row = dict(motor=a.motor, case=a.case, mass=a.mass, t_max=t_max,
               driver=a.driver, jet=a.jet, solid=a.solid or s.get('solid_heating_model', 'default'),
               energy=ef, momentum=mf, wall_total_s=wall,
               wall_time_s=float(s['wall_time']), steps=int(s['steps']),
               termination=s['termination'],
               pyrogen_mass_final_g=float(s.get('pyrogen_mass_final_kg', float('nan'))) * 1e3,
               ignition_attempts=int(s.get('ignition_attempts', 1)),
               pyrogen_mass_auto=bool(s.get('pyrogen_mass_auto', True)),
               pyrogen_mass_initial_g=float(s.get('pyrogen_mass_initial', float('nan'))) * 1e3,
               warnings=[str(w.message)[:200] for w in wl
                         if 'never ignited' in str(w.message)
                         or 'auto-sized' in str(w.message)])
    k = int(P.argmax())
    row['peak'] = float(P[k])
    row['t_peak'] = float(t[k])
    if a.motor == 'chunc':
        pl = float(np.median(P[(t >= 0.4) & (t <= 0.9)]))
        row.update(plateau=pl, ratio=float(P[k]) / pl)
    elif a.motor == 'balls':
        ms = t < 0.3
        kp = int(np.argmax(P[ms]))
        Ppk, tpk = float(P[ms][kp]), float(t[ms][kp])
        after = (t > tpk) & (t < 1.2)
        Ptr = float(np.min(P[after]))
        Pmid = float(P[int(np.argmin(np.abs(t - 1.03)))])
        row.update(peak=Ppk, t_peak=tpk, trough=Ptr, mid=Pmid, pk_tr=Ppk / Ptr,
                   pk_mid=Ppk / Pmid)
    if mass is not None:
        # run_from_ric on c5f4af9 ignores an explicit charge when the .ric block says
        # mass -1 unless pyrogen='bpnv' is passed (CANON does); verify it was honoured.
        assert not row['pyrogen_mass_auto'] and abs(row['pyrogen_mass_initial_g'] - float(a.mass)) < 1e-6, row
    if a.motor == 'haseg' and (t > 0.5).any():
        row['late_peak'] = float(P[t > 0.5].max())   # the progressive late peak (exp ~6.44 MPa)
        row['t_late_peak'] = float(t[t > 0.5][int(P[t > 0.5].argmax())])
    mm = np.asarray(res['max_mach'])
    pre = t <= (row['t_peak'] if a.motor != 'haseg' else 0.05)
    row['fill_mach_max'] = float(mm[t <= row['t_peak']].max()) if (t <= row['t_peak']).any() else None
    row['mach_max_all'] = float(mm.max())
    row['floor_activations'] = int(s.get('pressure_floor_activations', 0))
    row['min_P_Pa_first_50ms'] = float(np.min(np.asarray(res['min_pressure'])[t <= 0.05])) \
        if 'min_pressure' in res else None
    ign = ignition_stats(res)
    row['ignition'] = ign
    row['ign_time_ms_by_cell'] = [None if v >= 1.0e9 else float(v) * 1e3
                                  for v in np.asarray(res['ignition_time_by_cell'])]
    d14 = np.asarray(res['ignition_diag'])[:, 14]
    row['jet_h_at_ignition_max'] = float(np.nanmax(d14)) if np.isfinite(d14).any() else None
    row['first_cell_T'] = first_cell_temperature(res, ign)
    row['fill_map'] = fill_map(res)
    # EOS / storage-jump channels (whole run and fill phase) for the T6 table
    if not a.no_ledger:
        try:
            row['ledger'] = ledger(res)
        except Exception as exc:   # keep the row, report the failure
            row['ledger_error'] = repr(exc)
    print(json.dumps({k: v for k, v in row.items() if k not in ('fill_map', 'ledger')},
                     default=float), flush=True)
    if a.out:
        Path(a.out).write_text(json.dumps(row, default=float, indent=1))


if __name__ == '__main__':
    main()
