"""P5 probe: energy-coupled PISO iteration (thread F discriminator).

Rows are keyed (part, motor, target_cells, cfl, tmax, K, S) with
K = piso_energy_coupling_iters, S = piso_corr2_storage (0/1 ablation: the
omitted storage term of the second pressure corrector, see worklog/P5.md).
Each finished row is appended to a JSONL accumulator; present keys are skipped
(resumable). Rows run ONE AT A TIME in this process.

Parts (all canonical knobs: roughness 32e-6, kappa 0.44, T_ignition 756,
library BPNV, ma2020, committed .ric):
  A  Chunc 100 cells, CFL 0.3, t_max 0.05  : EOS jump vs K, fill closure
  B  Chunc grids x CFL, t_max 0.05          : fill Mach / min P / floor
  C  Chunc 100/400 cells, CFL 0.3, t_max 1  : peak, t_peak, ignition, plateau
  D  BALLSstick 100 cells, CFL 0.3, t_max 2 : peak, pk/trough, pk/mid, fill Mach
Usage (repo root, PYTHONPATH=$PWD, SRM1D_OPENMOTOR_PATH exported):
  python docs/solver_verification/probes/p5_energy_coupling.py run --part A
  python .../p5_energy_coupling.py run --part B --grid 100,400 --cfl 0.3 --K 0,3 --S 0,1
  python .../p5_energy_coupling.py report
"""
import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np  # noqa: E402

OUT = Path("C:/tmp/sv_p5/results.jsonl")
CANON = dict(roughness=32e-6, kappa=0.44, T_ignition=756.0,
             print_interval=0.0, verbose=False, pyrogen='bpnv', pyrogen_mass=None)
MOTORS = {'chunc': ('machbusterNew.ric', {}),
          'balls': ('BALLSstick.ric', {'P_cutoff': 0.01e6})}


def key_of(r):
    return (r['part'], r['motor'], int(r['target_cells']), float(r['cfl']),
            float(r['tmax']), int(r['K']), int(r['S']))


def load(path=OUT):
    rows = []
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _fin(x):
    x = float(x)
    return x if np.isfinite(x) else None


def run_row(part, motor, cells, cfl, tmax, K, S):
    from srm_1d.openmotor_adapter import run_from_ric
    from srm_1d.tools.ignition_diagnostics import sensible_enthalpy_audit
    ric, extra = MOTORS[motor]
    kw = dict(CANON, **extra)
    t0 = time.time()
    res = run_from_ric(str(ROOT / 'motors' / ric), target_propellant_cells=cells,
                       cfl_target=cfl, t_max=tmax, piso_energy_coupling_iters=K,
                       piso_corr2_storage=bool(S), **kw)[0]
    wall = time.time() - t0
    s = res['summary']
    t = np.asarray(res['time']); P = np.asarray(res['P_head']) / 1e6
    mm = np.asarray(res['max_mach']); mp = np.asarray(res['min_pressure'])
    n = len(t)
    k = int(P.argmax())
    row = dict(part=part, motor=motor, target_cells=cells, cells=int(s['cells']),
               cfl=cfl, tmax=tmax, K=K, S=S, n_steps=int(s['steps']), wall_s=wall,
               steps_per_s=int(s['steps']) / wall, termination=s.get('termination'),
               peak=float(P[k]), t_peak=float(t[k]), t_end=float(t[-1]))
    pre = t < t[k]
    row['fill_mach_prepeak'] = float(mm[pre].max()) if pre.any() else None
    row['t_fill_mach'] = float(t[pre][mm[pre].argmax()]) if pre.any() else None
    m3 = t <= 0.003
    row['fill_mach_3ms'] = float(mm[m3].max()) if m3.any() else None
    j = int(mp.argmin())
    row.update(min_P_Pa=float(mp[j]), t_min_P=float(t[j]),
               floor_min_raw_P=_fin(s['pressure_floor_min_raw_P']),
               floor_activations=int(s['pressure_floor_activations']),
               floor_steps=int(s['pressure_floor_steps']),
               floor_mass_kg=float(s['pressure_floor_mass_added_kg']))
    # coupling cost / convergence
    cp = np.asarray(res['piso_coupling_passes']); rd = np.asarray(res['piso_coupling_rel_dT'])
    row.update(mean_passes=float(cp.mean()), max_passes=int(cp.max()),
               frac_steps_unconverged=float(np.mean(rd >= 1e-6)),
               mean_rel_dT=float(rd.mean()), max_rel_dT=float(rd.max()))
    fm = t <= 0.003
    row['mean_passes_fill3ms'] = float(cp[fm].mean()) if fm.any() else None
    # EOS storage jump channels (last slot incomplete -> excluded)
    gJ = np.asarray(res['storage_jump_eos_J'])[:n - 1]
    gk = np.asarray(res['mass_jump_eos_kg'])[:n - 1]
    p2 = np.asarray(res['storage_jump_eos_p2_kg'])[:n - 1]
    row.update(eos_J_abs=float(np.abs(gJ).sum()), eos_J_signed=float(gJ.sum()),
               eos_J_max=float(np.abs(gJ).max()),
               eos_kg_abs=float(np.abs(gk).sum()), eos_kg_signed=float(gk.sum()),
               eos_p2_kg_signed=float(p2.sum()), eos_p2_kg_abs=float(np.abs(p2).sum()))
    # fill-phase closure (T6 audit): jumps / gross source and no-jump closure error
    try:
        au = sensible_enthalpy_audit(res)
        fb = au['phases']['fill']['budget_J']
        row.update(fill_bounds=au['phase_bounds'],
                   fill_src=fb['sources_in_gross'], fill_eos_J=fb['eos'],
                   fill_jumps_total_J=fb['jumps_total'],
                   fill_closure_nojump_frac=fb['closure_error_without_jumps_frac_of_source'],
                   fill_closure_jump_frac=fb['closure_error_with_jumps_frac_of_source'],
                   run_closure_nojump_frac=au['run']['budget_J']['closure_error_without_jumps_frac_of_source'])
    except Exception as e:  # noqa
        row['audit_error'] = repr(e)
    # ignition timing (cells that ignite by end of run)
    ig = np.asarray(res['ignition_time_by_cell'], dtype=float)
    ign = ig[ig < 1e9]
    row['n_ignited'] = int(ign.size)
    row['n_cells_total'] = int(ig.size)
    if ign.size:
        row['t_first_ign'] = float(ign.min())
        row['t_ign_50'] = float(np.percentile(ign, 50))
        row['t_ign_90'] = float(np.percentile(ign, 90))
        row['t_ign_all'] = float(ign.max())
    if motor == 'chunc' and tmax >= 0.9:
        w = (t >= 0.4) & (t <= 0.9)
        row['plateau'] = float(np.median(P[w]))
    if motor == 'balls':
        ms = t < 0.3
        kp = int(np.argmax(P[ms]))
        Ppk, tpk = float(P[ms][kp]), float(t[ms][kp])
        after = (t > tpk) & (t < 1.2)
        Ptr = float(np.min(P[after]))
        Pmid = float(P[int(np.argmin(np.abs(t - 1.03)))])
        row.update(peak=Ppk, t_peak=tpk, trough=Ptr, mid=Pmid, pk_tr=Ppk / Ptr, pk_mid=Ppk / Pmid)
        pre = t < tpk
        row['fill_mach_prepeak'] = float(mm[pre].max()) if pre.any() else None
    return row


def ints(s):
    return [int(x) for x in s.split(',')]


def floats(s):
    return [float(x) for x in s.split(',')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run', 'report'])
    ap.add_argument('--part', default='A')
    ap.add_argument('--motor', default=None)
    ap.add_argument('--grid', default=None)
    ap.add_argument('--cfl', default=None)
    ap.add_argument('--K', default='0,1,3')
    ap.add_argument('--S', default='0,1')
    ap.add_argument('--tmax', type=float, default=None)
    ap.add_argument('--out', default=str(OUT))
    a = ap.parse_args()
    if a.cmd == 'report':
        return report(a.out)
    import srm_1d
    print(f"[provenance] srm_1d from {Path(srm_1d.__file__).parent}", flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    defaults = {'A': ('chunc', '100', '0.3', 0.05), 'B': ('chunc', '100,400,800', '0.3,0.15', 0.05),
                'C': ('chunc', '100,400', '0.3', 1.0), 'D': ('balls', '100', '0.3', 2.0)}
    dm, dg, dc, dt = defaults[a.part]
    motor = a.motor or dm
    grid = ints(a.grid or dg)
    cfls = floats(a.cfl or dc)
    tmax = a.tmax or dt
    done = {key_of(r) for r in load(a.out)}
    for cells in grid:
        for cfl in cfls:
            for S in ints(a.S):
                for K in ints(a.K):
                    r0 = dict(part=a.part, motor=motor, target_cells=cells, cfl=cfl,
                              tmax=tmax, K=K, S=S)
                    if key_of(r0) in done:
                        print('skip', key_of(r0), flush=True)
                        continue
                    r = run_row(a.part, motor, cells, cfl, tmax, K, S)
                    with open(a.out, 'a') as f:
                        f.write(json.dumps(r) + '\n')
                    print(f"{key_of(r0)}: peak {r['peak']:.4f} tpk {r['t_peak']*1e3:.2f}ms "
                          f"Mach {r['fill_mach_prepeak']} minP {r['min_P_Pa']:.0f} "
                          f"floor {r['floor_activations']} steps {r['n_steps']} "
                          f"wall {r['wall_s']:.0f}s passes {r['mean_passes']:.2f}", flush=True)


def report(out):
    rows = load(out)

    def g(v, f='{:.5g}'):
        return '-' if v is None else f.format(v)
    for part in 'ABCD':
        pr = [r for r in rows if r['part'] == part]
        if not pr:
            continue
        print(f"\n### Part {part}")
        pr.sort(key=lambda r: (r['target_cells'], -r['cfl'], r['S'], r['K']))
        print("| cells | CFL | S | K | peak MPa | t_peak ms | plateau | fill Mach | min P Pa | raw floor P | floor act | "
              "first ign ms | t90 ms | EOS J abs | EOS J signed | EOS kg abs | p2 kg | fill closure(no jumps) | passes | steps | steps/s |")
        print("|" + "---|" * 21)
        for r in pr:
            print(f"| {r['cells']} | {r['cfl']} | {r['S']} | {r['K']} | {r['peak']:.4f} | {r['t_peak']*1e3:.2f} | "
                  f"{g(r.get('plateau'), '{:.4f}')} | {g(r['fill_mach_prepeak'], '{:.3f}')} | {r['min_P_Pa']:.0f} | "
                  f"{g(r['floor_min_raw_P'], '{:.0f}')} | {r['floor_activations']} | "
                  f"{g(None if r.get('t_first_ign') is None else r['t_first_ign']*1e3, '{:.3f}')} | "
                  f"{g(None if r.get('t_ign_90') is None else r['t_ign_90']*1e3, '{:.2f}')} | "
                  f"{r['eos_J_abs']:.4g} | {r['eos_J_signed']:.4g} | {r['eos_kg_abs']:.3g} | {r['eos_p2_kg_signed']:.3g} | "
                  f"{g(r.get('fill_closure_nojump_frac'), '{:.4f}')} | {r['mean_passes']:.2f} | {r['n_steps']} | {r['steps_per_s']:.0f} |")
        if part == 'D':
            for r in pr:
                print(f"BALLS S{r['S']} K{r['K']}: peak {r['peak']:.4f} pk/trough {r['pk_tr']:.4f} pk/mid {r['pk_mid']:.4f} "
                      f"fill Mach {r['fill_mach_prepeak']}")


if __name__ == '__main__':
    main()
