"""
cross_motor_lhs.py — v0.7.5 cross-motor re-LHS (shared-knob calibration).
=========================================================================

Models the v0.7.4 Phase C.2 2x2 sweep
(``hasegawa_a_lhs_mode_transport_2x2.py``) but the cells are MOTORS instead
of mode x transport. Sweeps the shared Ma + Goodman-ignition physics knobs
(roughness, kappa, T_ignition, k_solid) across the fired motors that ship
experimental traces — Hasegawa A, Zerox, Chunc (machbusterNew) — and finds
the single knob set that best fits all of them at once.

Post Phase 3/4 each motor is self-describing: transport comes from its
`.ric` per-tab block at the new **frozen** default, and the igniter comes
from its `.ric` `data.igniter` block (so ``pyrogen=None`` — no per-motor
pyrogen kwargs here). Pyrogen sizing is NOT swept; this isolates the
cross-motor physics calibration.

Cross-motor combine
-------------------
``run_lhs`` draws a fixed-seed Latin Hypercube, so the SAME knob sets are
evaluated for every motor. Rows are keyed by the knob tuple, each motor's
fitness is normalized by its median (equal weighting across motors), and
summed to a combined score; the lowest combined score is the cross-motor
optimum.

Conventions inherited from the 2x2 driver:
    SRM_LHS_SAMPLES   samples per motor (default 1000; <=4 = smoke test)
    SRM_LHS_WORKERS   worker count (default 16; 'auto' = cpu_count)
    SRM_LHS_FITNESS   'mse' (default, robust cross-motor) | 'segmented'
    SRM_LHS_PROGRESS  'brief' (default) | 'verbose' | 'none'
    SRM_LHS_MOTORS    comma subset of {hasegawa_a,zerox,chunc}

Launch (overnight, large pool):
    python -m examples.cross_motor_lhs
Pre-flight smoke test (3 runs, ~minutes):
    SRM_LHS_SAMPLES=1 python -m examples.cross_motor_lhs

Outputs (artifacts/cross_motor_lhs/):
    <motor>/lhs.csv            per-motor LHS rows
    <motor>/best_diagnostics.png   per-motor combined-best trace vs exp
    cross_motor_combined.csv   knob sets ranked by combined score
    comparison.png             best shared knobs, one panel per motor
    rank1_knobs.md             combined top-5 knob table
"""

import os
import multiprocessing
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.tools.sensitivity import (
    DEFAULT_PRESSURE_SEGMENTS, run_lhs, mse_fitness,
    segmented_pressure_fitness,
)
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import (
    HASEGAWA_MOTOR_A_EXPERIMENTAL, ZEROX_EXPERIMENTAL, CHUNC_EXPERIMENTAL,
)

MOTORS_DIR = Path(__file__).resolve().parents[1] / 'motors'
OUTPUT_ROOT = Path('artifacts') / 'cross_motor_lhs'
SEED = 42

# --- worker / sample config (2x2 conventions) ---
DEFAULT_WORKERS = 16
N_SAMPLES = int(os.environ.get('SRM_LHS_SAMPLES', '1000'))
_WORKERS_ENV = os.environ.get('SRM_LHS_WORKERS')
if _WORKERS_ENV is None:
    MAX_WORKERS = DEFAULT_WORKERS
elif _WORKERS_ENV.lower() == 'auto':
    MAX_WORKERS = multiprocessing.cpu_count()
else:
    MAX_WORKERS = int(_WORKERS_ENV)
if N_SAMPLES <= 4:           # smoke test → single worker, clean streams
    MAX_WORKERS = 1

FITNESS_MODE = os.environ.get('SRM_LHS_FITNESS', 'mse')
PROGRESS_MODE = os.environ.get('SRM_LHS_PROGRESS', 'brief')
SIM_VERBOSE = os.environ.get('SRM_SIM_VERBOSE', '0').lower() in {'1', 'true', 'yes'}

# Shared physics knobs (NOT pyrogen sizing). Bounds match the v0.7.4
# Phase C.2 2x2 lit-bounded band: roughness 15-100 um, kappa near the 0.45
# center, T_ignition AP/HTPB-typical, k_solid the AP/HTPB+Al defensible
# band (0.26-0.32, reference_ap_htpb_k_solid_literature).
SHARED_BOUNDS = {
    'roughness':  (15e-6, 100e-6),
    'kappa':      (0.40, 0.50),
    'T_ignition': (750.0, 950.0),
    'k_solid':    (0.26, 0.32),
}
_KNOB_KEYS = tuple(SHARED_BOUNDS.keys())

SEGMENT_WEIGHTS = {
    'mse_spike': 0.25, 'mse_post_spike': 0.35,
    'mse_plateau': 0.20, 'mse_taildown': 0.20,
}

# Per-motor fit config: .ric, experimental dict, ignition-peak alignment
# window [s], sim horizon. time_offset aligns experimental ignition to t=0.
MOTOR_CONFIGS = {
    'hasegawa_a': {
        'ric': 'hasegawa_a.ric', 'exp': HASEGAWA_MOTOR_A_EXPERIMENTAL,
        'time_offset': 0.0, 'peak_align_window': (0.02, 0.18), 't_max': 6.0,
    },
    'zerox': {
        'ric': 'zerox.ric', 'exp': ZEROX_EXPERIMENTAL,
        'time_offset': ZEROX_EXPERIMENTAL.get('time_offset', -0.3),
        'peak_align_window': (0.0, 0.6), 't_max': 8.5,
    },
    'chunc': {
        'ric': 'machbusterNew.ric', 'exp': CHUNC_EXPERIMENTAL,
        'time_offset': 0.0, 'peak_align_window': (0.0, 0.20), 't_max': 2.6,
    },
}
MOTOR_COLORS = {'hasegawa_a': '#1f77b4', 'zerox': '#ff7f0e', 'chunc': '#2ca02c'}


def _knob_key(row):
    return tuple(round(float(row[k]), 12) for k in _KNOB_KEYS)


def _fitness_fn(cfg, t_exp, p_exp):
    if FITNESS_MODE == 'segmented':
        return segmented_pressure_fitness(
            t_exp, p_exp, t_min=0.01, segments=DEFAULT_PRESSURE_SEGMENTS,
            weights=SEGMENT_WEIGHTS, peak_align_window=cfg['peak_align_window'])
    if FITNESS_MODE == 'mse':
        return mse_fitness(t_exp, p_exp, t_min=0.01,
                           peak_align_window=cfg['peak_align_window'])
    raise ValueError("SRM_LHS_FITNESS must be 'mse' or 'segmented'")


def _run_motor(name, cfg):
    cell_dir = OUTPUT_ROOT / name
    cell_dir.mkdir(parents=True, exist_ok=True)
    exp = cfg['exp']
    t_exp = np.asarray(exp['time']) + cfg['time_offset']
    p_exp = np.asarray(exp['pressure'])
    print(f"\n[motor {name}] {cfg['ric']}  N={N_SAMPLES} workers={MAX_WORKERS} "
          f"fitness={FITNESS_MODE}")
    rows = run_lhs(
        motor_path=str(MOTORS_DIR / cfg['ric']),
        bounds=SHARED_BOUNDS,
        n_samples=N_SAMPLES,
        fitness_fn=_fitness_fn(cfg, t_exp, p_exp),
        n_workers=MAX_WORKERS,
        seed=SEED,
        csv_path=str(cell_dir / 'lhs.csv'),
        progress_mode=PROGRESS_MODE,
        sim_verbose=SIM_VERBOSE,
        # Each motor uses its own .ric igniter block + frozen transport.
        pyrogen=None,
        t_max=cfg['t_max'], P_cutoff=0.05e6,
        snapshot_interval=cfg['t_max'], print_interval=cfg['t_max'] * 10,
    )
    return rows


def _median_norm(values):
    finite = values[np.isfinite(values)]
    med = np.median(finite) if finite.size else 1.0
    return values / med if med > 0 else values


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    requested = os.environ.get('SRM_LHS_MOTORS')
    names = ([n.strip() for n in requested.split(',')] if requested
             else list(MOTOR_CONFIGS.keys()))

    print('=' * 60)
    print('v0.7.5 cross-motor re-LHS (shared physics knobs, frozen default)')
    print(f'  motors    = {names}')
    print(f'  N/motor   = {N_SAMPLES}   workers = {MAX_WORKERS}')
    print(f'  fitness   = {FITNESS_MODE}')
    print(f'  bounds    = {SHARED_BOUNDS}')
    print('=' * 60)

    per_motor = {name: _run_motor(name, MOTOR_CONFIGS[name]) for name in names}

    # --- combine by shared knob set ---
    fitness_by_motor = {}
    for name, rows in per_motor.items():
        fitness_by_motor[name] = {
            _knob_key(r): float(r['fitness'])
            for r in rows if not r.get('error')
        }
    common = set.intersection(*[set(d) for d in fitness_by_motor.values()])
    if not common:
        print("\nNo knob set succeeded across all motors — check bounds/sims.")
        return

    keys = sorted(common)
    norm = {name: _median_norm(np.array([fitness_by_motor[name][k] for k in keys]))
            for name in names}
    combined = np.sum([norm[name] for name in names], axis=0)
    order = np.argsort(combined)

    rows_out = []
    for idx in order:
        k = keys[idx]
        row = dict(zip(_KNOB_KEYS, k))
        row['combined'] = float(combined[idx])
        for name in names:
            row[f'fitness_{name}'] = fitness_by_motor[name][k]
        rows_out.append(row)

    cols = list(_KNOB_KEYS) + ['combined'] + [f'fitness_{n}' for n in names]
    combined_csv = OUTPUT_ROOT / 'cross_motor_combined.csv'
    with open(combined_csv, 'w', encoding='utf-8') as f:
        f.write(','.join(cols) + '\n')
        for row in rows_out:
            f.write(','.join(f"{row[c]:.6g}" for c in cols) + '\n')
    print(f"\nSaved {combined_csv} ({len(rows_out)} shared-knob sets)")

    # rank-1 knobs markdown (combined top-5)
    md = ['# Cross-motor re-LHS — combined top-5 shared knobs', '',
          f'**Samples per motor**: {N_SAMPLES}  **Workers**: {MAX_WORKERS}  '
          f'**Fitness**: {FITNESS_MODE}', '',
          '| Rank | combined | roughness [um] | kappa | T_ign [K] | k_solid | '
          + ' | '.join(f'fit_{n}' for n in names) + ' |',
          '|------|----------|----------------|-------|-----------|---------|'
          + '|'.join('-------' for _ in names) + '|']
    for rank, row in enumerate(rows_out[:5], start=1):
        md.append(
            f"| {rank} | {row['combined']:.3f} | {row['roughness']*1e6:.1f} | "
            f"{row['kappa']:.3f} | {row['T_ignition']:.0f} | {row['k_solid']:.3f} | "
            + ' | '.join(f"{row[f'fitness_{n}']:.3f}" for n in names) + ' |')
    (OUTPUT_ROOT / 'rank1_knobs.md').write_text('\n'.join(md) + '\n', encoding='utf-8')

    print("\n" + "=" * 60)
    print("--- TOP 5 CROSS-MOTOR KNOB SETS (lower combined = better) ---")
    print("=" * 60)
    for rank, row in enumerate(rows_out[:5], start=1):
        print(f"Rank {rank}: combined={row['combined']:.3f}  "
              f"roughness={row['roughness']*1e6:.1f}um kappa={row['kappa']:.3f} "
              f"T_ign={row['T_ignition']:.0f}K k_solid={row['k_solid']:.3f}")
        print("  per-motor fitness: " +
              "  ".join(f"{n}={row[f'fitness_{n}']:.3f}" for n in names))

    # --- best shared knobs: per-motor trace overlay ---
    best = rows_out[0]
    best_params = {k: best[k] for k in _KNOB_KEYS}
    fig, axes = plt.subplots(1, len(names), figsize=(6 * len(names), 5),
                             squeeze=False)
    for ax, name in zip(axes[0], names):
        cfg = MOTOR_CONFIGS[name]
        exp = cfg['exp']
        result, *_ = run_from_ric(
            str(MOTORS_DIR / cfg['ric']), pyrogen=None,
            t_max=cfg['t_max'], P_cutoff=0.05e6,
            snapshot_interval=cfg['t_max'], print_interval=cfg['t_max'] * 10,
            verbose=False, **best_params)
        ax.plot(np.asarray(exp['time']) + cfg['time_offset'], exp['pressure'],
                'k.-', linewidth=2, label='experimental')
        ax.plot(result['time'], result['P_head'] / 1e6,
                color=MOTOR_COLORS.get(name, 'b'), linewidth=1.5,
                label='1D PISO (best shared)')
        ax.set_title(f"{name}  fit={best[f'fitness_{name}']:.3f}")
        ax.set_xlabel('Time [s]'); ax.set_ylabel('P_head [MPa]')
        ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle(
        f"Cross-motor best shared knobs (N={N_SAMPLES})  "
        f"roughness={best['roughness']*1e6:.0f}um kappa={best['kappa']:.2f} "
        f"T_ign={best['T_ignition']:.0f}K k_solid={best['k_solid']:.2f}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUTPUT_ROOT / 'comparison.png', dpi=200)
    plt.close(fig)
    print(f"Saved {OUTPUT_ROOT / 'comparison.png'}")
    print(f"Saved {OUTPUT_ROOT / 'rank1_knobs.md'}")


if __name__ == '__main__':
    main()
