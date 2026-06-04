"""
cross_motor_lhs_v075.py — v0.7.5 cross-motor re-LHS (the v0.7.4 closing task).
=============================================================================

Runs on the v0.7.x base (branch v0.7.0-phase4 / tag v0.7.4) — NOT the
v0.8.0 frontend branch. This is the calibration line: validated physics,
sidecar transport, explicit pyrogen, no channel/API churn.

Implements the v0.7.4 TASKS.md "Cross-motor regression + re-LHS" step:
sweep the SHARED Ma + Goodman-ignition knobs across the fired motors that
ship experimental traces (Hasegawa A, Zerox, Chunc/machbusterNew) and find
the single knob set that best fits all three — with the v0.7.4 spike fix
(F+Z) ENABLED and FROZEN transport (the v0.7.3.3 frozen-wins finding).

Locked per the documented task:
- Transport = FROZEN (each motor's `<motor>.frozen.transport.yaml`).
- Igniter   = explicit `pyrogen='bpnv'` (v0.7.x has no .ric igniter block).
- Spike fix = `flame_front_enabled=True`, `zn_enabled=True`,
              `kappa_zn=1.0` FIXED — Z-N relaxation strength is NOT a free
              LHS knob (feedback / v0.7.4 TASKS).
- Physical bounds = roughness >= 15 um, kappa near the 0.45 center
              (feedback_roughness_kappa_physical_bounds).

Cross-motor combine: a fixed-seed LHS evaluates the SAME knob sets for
every motor; rows are keyed by the knob tuple, each motor's fitness is
median-normalized (equal weighting), and summed — lowest combined score is
the cross-motor optimum.

Conventions (2x2 driver):
    SRM_LHS_SAMPLES  samples per motor (default 1000; <=4 = smoke test)
    SRM_LHS_WORKERS  workers (default 16; 'auto' = cpu_count)
    SRM_LHS_FITNESS  'mse' (default) | 'segmented'
    SRM_LHS_PROGRESS 'brief' | 'verbose' | 'none'
    SRM_LHS_MOTORS   comma subset of {hasegawa_a,zerox,chunc}

Launch (overnight): python -m srm_1d.examples.cross_motor_lhs_v075
Smoke test:         SRM_LHS_SAMPLES=1 python -m srm_1d.examples.cross_motor_lhs_v075

Outputs: artifacts/cross_motor_lhs_v075/{<motor>/lhs.csv, cross_motor_combined.csv,
         comparison.png, rank1_knobs.md}
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
OUTPUT_ROOT = Path('artifacts') / 'cross_motor_lhs_v075'
SEED = 42

DEFAULT_WORKERS = 16
# 8h budget sizing (measured 2026-06-03 on this 16-core box, single-thread,
# post-warmup): hasegawa_a ~30.9 s/sim, zerox ~76.3 s/sim, chunc ~26.8 s/sim
# => ~133.9 s per knob-tuple (all 3 motors). At ideal 16x parallel, 8h fits
# ~3440 tuples; default carries a ~0.85 margin (~3000) for real parallel
# inefficiency (HT/memory contention) + P_cutoff per-sim spread. Motors run
# sequentially and the CSV checkpoint is crash-safe but NOT resumable, so size
# to FINISH within budget (overrun loses the last motor's combine).
N_SAMPLES = int(os.environ.get('SRM_LHS_SAMPLES', '3000'))
_WORKERS_ENV = os.environ.get('SRM_LHS_WORKERS')
if _WORKERS_ENV is None:
    MAX_WORKERS = DEFAULT_WORKERS
elif _WORKERS_ENV.lower() == 'auto':
    MAX_WORKERS = multiprocessing.cpu_count()
else:
    MAX_WORKERS = int(_WORKERS_ENV)
if N_SAMPLES <= 4:
    MAX_WORKERS = 1

FITNESS_MODE = os.environ.get('SRM_LHS_FITNESS', 'mse')
PROGRESS_MODE = os.environ.get('SRM_LHS_PROGRESS', 'brief')
SIM_VERBOSE = os.environ.get('SRM_SIM_VERBOSE', '0').lower() in {'1', 'true', 'yes'}

# Shared physics knobs only. Physical bounds enforced (roughness >= 15 um,
# kappa near 0.45). k_solid in the AP/HTPB+Al defensible band.
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

# v0.7.4 spike-fix knobs — LOCKED across the sweep (not LHS variables).
SPIKE_FIX = dict(
    flame_front_enabled=True,
    zn_enabled=True,
    kappa_zn=1.0,   # FIXED — do not let the LHS treat as a free knob.
)

MOTOR_CONFIGS = {
    'hasegawa_a': {
        'ric': 'hasegawa_a.ric', 'frozen': 'hasegawa_a.frozen.transport.yaml',
        'exp': HASEGAWA_MOTOR_A_EXPERIMENTAL, 'time_offset': 0.0,
        'peak_align_window': (0.02, 0.18), 't_max': 6.0,
    },
    'zerox': {
        'ric': 'zerox.ric', 'frozen': 'zerox.frozen.transport.yaml',
        'exp': ZEROX_EXPERIMENTAL,
        'time_offset': ZEROX_EXPERIMENTAL.get('time_offset', -0.3),
        'peak_align_window': (0.0, 0.6), 't_max': 8.5,
    },
    'chunc': {
        'ric': 'machbusterNew.ric', 'frozen': 'machbusterNew.frozen.transport.yaml',
        'exp': CHUNC_EXPERIMENTAL, 'time_offset': 0.0,
        'peak_align_window': (0.0, 0.20), 't_max': 2.6,
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
    print(f"\n[motor {name}] {cfg['ric']} (frozen)  N={N_SAMPLES} "
          f"workers={MAX_WORKERS} fitness={FITNESS_MODE}  F+Z on")
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
        # Locked: explicit pyrogen, FROZEN transport, F+Z spike fix.
        pyrogen='bpnv',
        transport_path=str(MOTORS_DIR / cfg['frozen']),
        t_max=cfg['t_max'], P_cutoff=0.05e6,
        snapshot_interval=cfg['t_max'], print_interval=cfg['t_max'] * 10,
        **SPIKE_FIX,
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
    print('v0.7.5 cross-motor re-LHS  (FROZEN transport, F+Z on, kappa_zn=1)')
    print(f'  motors  = {names}')
    print(f'  N/motor = {N_SAMPLES}  workers = {MAX_WORKERS}  fitness = {FITNESS_MODE}')
    print(f'  bounds  = {SHARED_BOUNDS}')
    print('=' * 60)

    per_motor = {name: _run_motor(name, MOTOR_CONFIGS[name]) for name in names}

    fitness_by_motor = {
        name: {_knob_key(r): float(r['fitness'])
               for r in rows if not r.get('error')}
        for name, rows in per_motor.items()
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
    with open(OUTPUT_ROOT / 'cross_motor_combined.csv', 'w', encoding='utf-8') as f:
        f.write(','.join(cols) + '\n')
        for row in rows_out:
            f.write(','.join(f"{row[c]:.6g}" for c in cols) + '\n')

    md = ['# v0.7.5 cross-motor re-LHS — combined top-5 (FROZEN, F+Z on)', '',
          f'**Samples/motor**: {N_SAMPLES}  **Workers**: {MAX_WORKERS}  '
          f'**Fitness**: {FITNESS_MODE}  **kappa_zn**: 1.0 (fixed)', '',
          '| Rank | combined | roughness [um] | kappa | T_ign [K] | k_solid | '
          + ' | '.join(f'fit_{n}' for n in names) + ' |',
          '|------|----------|----------------|-------|-----------|---------|'
          + '|'.join('----' for _ in names) + '|']
    for rank, row in enumerate(rows_out[:5], start=1):
        md.append(
            f"| {rank} | {row['combined']:.3f} | {row['roughness']*1e6:.1f} | "
            f"{row['kappa']:.3f} | {row['T_ignition']:.0f} | {row['k_solid']:.3f} | "
            + ' | '.join(f"{row[f'fitness_{n}']:.3f}" for n in names) + ' |')
    (OUTPUT_ROOT / 'rank1_knobs.md').write_text('\n'.join(md) + '\n', encoding='utf-8')

    print("\n--- TOP 5 CROSS-MOTOR KNOB SETS (lower combined = better) ---")
    for rank, row in enumerate(rows_out[:5], start=1):
        print(f"Rank {rank}: combined={row['combined']:.3f}  "
              f"roughness={row['roughness']*1e6:.1f}um kappa={row['kappa']:.3f} "
              f"T_ign={row['T_ignition']:.0f}K k_solid={row['k_solid']:.3f}  | "
              + "  ".join(f"{n}={row[f'fitness_{n}']:.3f}" for n in names))

    best = rows_out[0]
    best_params = {k: best[k] for k in _KNOB_KEYS}
    fig, axes = plt.subplots(1, len(names), figsize=(6 * len(names), 5),
                             squeeze=False)
    for ax, name in zip(axes[0], names):
        cfg = MOTOR_CONFIGS[name]
        exp = cfg['exp']
        result, *_ = run_from_ric(
            str(MOTORS_DIR / cfg['ric']), pyrogen='bpnv',
            transport_path=str(MOTORS_DIR / cfg['frozen']),
            t_max=cfg['t_max'], P_cutoff=0.05e6,
            snapshot_interval=cfg['t_max'], print_interval=cfg['t_max'] * 10,
            verbose=False, **SPIKE_FIX, **best_params)
        ax.plot(np.asarray(exp['time']) + cfg['time_offset'], exp['pressure'],
                'k.-', linewidth=2, label='experimental')
        ax.plot(result['time'], result['P_head'] / 1e6,
                color=MOTOR_COLORS.get(name, 'b'), linewidth=1.5,
                label='1D PISO (best shared)')
        ax.set_title(f"{name}  fit={best[f'fitness_{name}']:.3f}")
        ax.set_xlabel('Time [s]'); ax.set_ylabel('P_head [MPa]')
        ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle(
        f"v0.7.5 cross-motor best shared knobs (FROZEN, F+Z; N={N_SAMPLES})  "
        f"roughness={best['roughness']*1e6:.0f}um kappa={best['kappa']:.2f} "
        f"T_ign={best['T_ignition']:.0f}K k_solid={best['k_solid']:.2f}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUTPUT_ROOT / 'comparison.png', dpi=200)
    plt.close(fig)
    print(f"\nSaved {OUTPUT_ROOT}/cross_motor_combined.csv, rank1_knobs.md, comparison.png")


if __name__ == '__main__':
    main()