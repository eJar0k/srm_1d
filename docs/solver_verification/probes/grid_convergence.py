"""T8 probe: Chunc grid convergence (part A) + fill-Mach attribution (part B).

ADAPTED (not cherry-picked) from the friend's Codex-fork scripts
``scripts/run_spatial_extension.py`` (8b0619c: ladder, resume, successive
change, factor-2 observed order + Richardson) and
``scripts/run_fill_configuration_comparison.py`` (cab26d0: legacy-BPNV variant
construction), switched to the CANONICAL config: committed ``.ric``, knobs
roughness 32e-6 / kappa 0.44 / T_ignition 756 / default k_solid, ma2020,
``pyrogen='bpnv', pyrogen_mass=None`` (library BPNV) unless the igniter variant
says otherwise.

Fill Mach is recorded in several windows so definitions can be compared:
  fill_mach_prepeak : max(max_mach[t <  t_peak])        (T7 / REOPENED-6 style)
  fill_mach_3ms     : max(max_mach[t <= 3 ms])          (Codex's metric)
  fill_mach_20ms    : max(max_mach[t <= 20 ms])         (plugin-path triage)
each with the time of the maximum. ``max_mach`` is the max over CELLS of
|Mach| (cell-centred), the same history array in both trees.

Every row is appended to a JSONL accumulator; a (part, igniter, cells, cfl,
tmax) key already present is skipped (resumable).

Usage (from the worktree root, pyenv 3.10.5, SRM1D_OPENMOTOR_PATH exported):
  python docs/solver_verification/probes/grid_convergence.py A --grid 50,100 --cfl 0.3
  python docs/solver_verification/probes/grid_convergence.py B --grid 50,100,200
  python docs/solver_verification/probes/grid_convergence.py report
Options: --tmax (A default 1.0, B default 0.05), --out JSONL path.
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

DEFAULT_OUT = Path("C:/tmp/sv_t8/results.jsonl")
RIC = ROOT / "motors" / "machbusterNew.ric"
CANON = dict(roughness=32e-6, kappa=0.44, T_ignition=756.0,
             print_interval=0.0, verbose=False)


def _pyrogen_for(igniter):
    """Return the kwargs selecting the igniter variant."""
    from srm_1d.openmotor_adapter import (SUPERSEDED_PYROGENS,
                                          _pyrogen_from_mapping)
    if igniter == 'lib_bpnv':
        return dict(pyrogen='bpnv', pyrogen_mass=None)
    if igniter == 'lib_mtv':
        return dict(pyrogen='mtv', pyrogen_mass=None)
    if igniter == 'legacy_bpnv':
        # pre-bfc2f3f BPNV snapshot (== what pre-fix .ric files embedded)
        snap = SUPERSEDED_PYROGENS['BPNV'][0]
        return dict(pyrogen=_pyrogen_from_mapping(snap), pyrogen_mass=None)
    if igniter == 'embedded':   # pyrogen=None -> the .ric's own igniter block
        return dict(pyrogen=None, pyrogen_mass=None)
    raise ValueError(igniter)


def _win_max(t, M, mask):
    if not mask.any():
        return None, None
    idx = np.flatnonzero(mask)
    k = idx[int(np.argmax(M[idx]))]
    return float(M[k]), float(t[k])


def _finite(x):
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def run_row(part, igniter, cells, cfl, tmax, knobs='canon'):
    from srm_1d.openmotor_adapter import run_from_ric
    kw = dict(CANON)
    if knobs == 'codex':   # Codex baseline_configs Chunc knobs (stale; for attribution only)
        kw.update(roughness=35e-6, kappa=0.45, T_ignition=900.0, k_solid=0.4,
                  P_cutoff=50000.0)
    kw.update(_pyrogen_for(igniter))
    t0 = time.time()
    res = run_from_ric(str(RIC), target_propellant_cells=cells, cfl_target=cfl,
                       t_max=tmax, **kw)[0]
    wall = time.time() - t0
    s = res['summary']
    t = np.asarray(res['time'])
    P = np.asarray(res['P_head']) / 1e6
    M = np.asarray(res['max_mach'])
    mp = np.asarray(res['min_pressure'])
    k = int(P.argmax())
    row = dict(part=part, igniter=igniter, knobs=knobs, target_cells=cells, cfl=cfl,
               tmax=tmax, cells=int(s['cells']),
               peak=float(P[k]), t_peak=float(t[k]),
               n_steps=int(s['steps']), wall_s=wall,
               termination=s.get('termination'),
               history_cap_reached=bool(s.get('history_cap_reached')),
               history_capacity=s.get('history_capacity'),
               history_estimated_steps=s.get('history_estimated_steps'),
               pyrogen_name=s.get('pyrogen_name'),
               pyrogen_mass_initial=_finite(s.get('pyrogen_mass_initial')),
               min_P_Pa=float(mp.min()), t_minP=float(t[int(mp.argmin())]),
               floor_activations=int(s['pressure_floor_activations']),
               floor_steps=int(s['pressure_floor_steps']),
               floor_min_raw_P=_finite(s['pressure_floor_min_raw_P']),
               floor_first_time_s=_finite(s['pressure_floor_first_time_s']),
               t_end=float(t[-1]),
               plateau=(float(np.median(P[(t >= 0.4) & (t <= 0.9)]))
                        if tmax >= 0.9 and (t >= 0.4).any() else None))
    # run_simulation's history array can fill before t_max on fine grids
    # (est_steps assumes a=1000 m/s); then the plateau window is partial.
    row['truncated'] = bool(float(t[-1]) < 0.999 * tmax)
    row['plateau_partial'] = bool(row['plateau'] is not None and float(t[-1]) < 0.9)
    for name, mask in (('prepeak', t < t[k]), ('3ms', t <= 0.003),
                       ('20ms', t <= 0.020)):
        m, tm = _win_max(t, M, mask)
        row['fill_mach_' + name] = m
        row['t_fill_mach_' + name] = tm
    row['dt_min'] = _finite(s['dt_min'])
    return row


def load_rows(path):
    rows = []
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if line.strip():
                rows.append(json.loads(line))
    return rows


def key_of(r):
    return (r['part'], r['igniter'], r.get('knobs', 'canon'), int(r['target_cells']),
            float(r['cfl']), float(r['tmax']))


# ----------------------------------------------------------------- report ---
def observed_order(f3, f2, f1):
    """f1 finest, f3 coarsest, refinement ratio 2 (target cells)."""
    d1, d2 = f2 - f1, f3 - f2
    if d1 == 0 or d2 == 0 or d1 * d2 <= 0:
        return None, None, 'non-monotone' if d1 * d2 < 0 else 'degenerate'
    p = math.log(d2 / d1) / math.log(2.0)
    if p <= 0:
        return p, None, 'p<=0 (diverging)'
    return p, f1 + (f1 - f2) / (2.0 ** p - 1.0), 'ok'


def _g(v, f='{:.5g}'):
    return '-' if v is None else f.format(v)


def report(path):
    allrows = load_rows(path)
    rows = [r for r in allrows if r['part'] == 'A' and r.get('knobs', 'canon') == 'canon']
    for cfl in sorted({r['cfl'] for r in rows}, reverse=True):
        sel = [r for r in rows if r['cfl'] == cfl]
        by = {}
        for r in sel:      # prefer the full-length row
            c = r['target_cells']
            if (c not in by or r['tmax'] > by[c]['tmax']
                    or (r['tmax'] == by[c]['tmax'] and 't_end' in r)):
                by[c] = r
        cs = sorted(by)
        print(f"\n### Part A, CFL {cfl}  (* = t_max<0.9: plateau n/a; ^ = history array full before t_max, plateau partial, t_end in last column)")
        print("| target | cells | tmax | peak MPa | t_peak ms | plateau MPa | "
              "Mach pre-peak (t ms) | Mach 3ms (t ms) | Mach 20ms (t ms) | "
              "min P Pa | min raw P Pa | floor act (steps) | n_steps | wall s | t_end s |")
        print("|" + "---|" * 15)
        for c in cs:
            r = by[c]
            star = ('*' if r['tmax'] < 0.9 else '') + (
                '^' if r.get('plateau_partial') else '')
            print(f"| {c}{star} | {r['cells']} | {r['tmax']:g} | {r['peak']:.5f} | "
                  f"{r['t_peak']*1e3:.4f} | {_g(r['plateau'], '{:.5f}')} | "
                  f"{r['fill_mach_prepeak']:.4f} ({r['t_fill_mach_prepeak']*1e3:.3f}) | "
                  f"{r['fill_mach_3ms']:.4f} ({r['t_fill_mach_3ms']*1e3:.3f}) | "
                  f"{r['fill_mach_20ms']:.4f} ({r['t_fill_mach_20ms']*1e3:.3f}) | "
                  f"{r['min_P_Pa']:.1f} | {_g(r['floor_min_raw_P'], '{:.1f}')} | "
                  f"{r['floor_activations']} ({r['floor_steps']}) | "
                  f"{r['n_steps']} | {r['wall_s']:.0f} | {_g(r.get('t_end'), '{:.4f}')} |")
        metrics = ('peak', 't_peak', 'plateau', 'fill_mach_prepeak',
                   'fill_mach_3ms', 'fill_mach_20ms', 'min_P_Pa')
        if len(cs) >= 3:
            c3, c2, c1 = cs[-3:]
            print(f"\nObserved order / Richardson from target cells "
                  f"{c3}/{c2}/{c1} (f1 = {c1} finest):")
            print("| metric | f3 | f2 | f1 | p | f_inf | status |")
            print("|---|---|---|---|---|---|---|")
            for m in metrics:
                v = [by[c][m] for c in (c3, c2, c1)]
                if m == 'plateau' and any(by[c].get('plateau_partial') for c in (c3, c2, c1)):
                    print(f"| {m} | excluded (partial window) | | | | | |")
                    continue
                if m == 'min_P_Pa' and any(by[c]['floor_activations'] > 0 for c in (c3, c2, c1)):
                    print(f"| {m} | {v[0]:.6g} | {v[1]:.6g} | {v[2]:.6g} | - | - | floor-censored |")
                    continue
                if any(x is None for x in v):
                    print(f"| {m} | n/a (needs full-length rows) | | | | | |")
                    continue
                p, fi, st = observed_order(*v)
                print(f"| {m} | {v[0]:.6g} | {v[1]:.6g} | {v[2]:.6g} | "
                      f"{_g(p, '{:.3f}')} | {_g(fi, '{:.6g}')} | {st} |")
    rows = [r for r in allrows if r['part'] == 'B']
    if rows:
        print("\n### Part B (t_max 0.05)")
        print("| igniter | knobs | target | cells | cfl | peak MPa | t_peak ms | "
              "Mach 3ms (t ms) | Mach 20ms (t ms) | Mach pre-peak (t ms) | "
              "min P Pa | floor act | n_steps | wall s |")
        print("|" + "---|" * 14)
        for r in sorted(rows, key=lambda r: (r.get('knobs', 'canon'), r['igniter'], r['cfl'],
                                             r['target_cells'])):
            print(f"| {r['igniter']} | {r.get('knobs', 'canon')} | {r['target_cells']} | {r['cells']} | "
                  f"{r['cfl']} | {r['peak']:.5f} | {r['t_peak']*1e3:.3f} | "
                  f"{r['fill_mach_3ms']:.4f} ({r['t_fill_mach_3ms']*1e3:.3f}) | "
                  f"{r['fill_mach_20ms']:.4f} ({r['t_fill_mach_20ms']*1e3:.3f}) | "
                  f"{r['fill_mach_prepeak']:.4f} ({r['t_fill_mach_prepeak']*1e3:.3f}) | "
                  f"{r['min_P_Pa']:.1f} | {r['floor_activations']} | "
                  f"{r['n_steps']} | {r['wall_s']:.0f} |")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['A', 'B', 'report'])
    ap.add_argument('--grid', default=None, help='comma list of target cells')
    ap.add_argument('--cfl', default='0.3', help='comma list')
    ap.add_argument('--igniters', default=None,
                    help='B: comma list of lib_bpnv,legacy_bpnv,lib_mtv,embedded')
    ap.add_argument('--tmax', type=float, default=None)
    ap.add_argument('--knobs', default='canon', choices=['canon', 'codex'])
    ap.add_argument('--out', default=str(DEFAULT_OUT))
    a = ap.parse_args()
    if a.mode == 'report':
        report(a.out)
        return
    import srm_1d
    print(f"[provenance] srm_1d from {Path(srm_1d.__file__).parent}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    done = {key_of(r) for r in load_rows(a.out)}
    default_grid = '50,100,200,400,800,1600' if a.mode == 'A' else '50,100,200'
    grid = [int(x) for x in (a.grid or default_grid).split(',')]
    cfls = [float(x) for x in a.cfl.split(',')]
    if a.mode == 'A':
        igs, tmax = ['lib_bpnv'], (a.tmax if a.tmax is not None else 1.0)
    else:
        igs = (a.igniters or 'lib_bpnv,legacy_bpnv,lib_mtv').split(',')
        tmax = a.tmax if a.tmax is not None else 0.05
    for cfl in cfls:
        for ig in igs:
            for c in grid:
                k = (a.mode, ig, a.knobs, c, cfl, tmax)
                if k in done:
                    print(f"[skip] {k}")
                    continue
                print(f"[run] {k}", flush=True)
                row = run_row(a.mode, ig, c, cfl, tmax, a.knobs)
                with open(a.out, 'a') as f:
                    f.write(json.dumps(row) + '\n')
                print(f"  cells {row['cells']} peak {row['peak']:.8f} "
                      f"t_peak {row['t_peak']:.6f} n_steps {row['n_steps']} "
                      f"wall {row['wall_s']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
