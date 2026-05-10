"""
zerox_lhs.py — 6-variable Latin Hypercube calibration of Zerox
================================================================

Searches a 6-dimensional parameter space for the best fit to the
Zerox static-fire pressure trace. Includes "pinned" variant sweeps
where one parameter is held at its current (Hasegawa-A inherited)
value while the other five are LHS-sampled — useful for assessing
whether each knob is essential or optional for a good fit.

History: v0.7.0 replaced the v0.6.0 ignition knobs with pyrogen mass,
pyrogen throat area, and surface ignition temperature.

Bounds are informed by the prior `zerox_a_sweep` (a near 1.0 once
throat erosion is corrected) and `zerox_erosion_sweep` (sweet spot
near erosion × 2.65).

Persistent artifacts (so plots can be regenerated without re-running):

  zerox_lhs.csv         — one row per sample: params + fitness + tag
  zerox_lhs.npz         — params, fitness, resampled pressure traces,
                          tags, time grid, param names. Load with
                          `data = np.load('zerox_lhs.npz', allow_pickle=True)`
  zerox_lhs_raw.pkl     — pickled list[dict], full result records
  zerox_lhs_top.png     — top-10 best traces overlaid on experimental
  zerox_lhs_pinned.png  — 8-panel: main best + 7 pinned-variant bests
  zerox_lhs_ensemble.png — alpha-faded ensemble of top-50

Usage:
    "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" \\
        -m srm_1d.examples.zerox_lhs

Wall-time estimate (auto-printed at start). With 6 cores at ~100 s/sim
and 392 samples, ~110 min. With more cores, proportionally less.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import csv
import os
import pickle
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc

from srm_1d.openmotor_adapter import (
    build_pyrogen_chamber,
    load_pyrogen,
    load_ric,
    load_transport,
    ric_to_sim_args,
)
from srm_1d.simulation import run_simulation
from srm_1d.plotting import ZEROX_EXPERIMENTAL


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox.ric'
TRANSPORT_PATH = MOTOR_PATH.with_name('zerox.transport.yaml')

# Current ("Hasegawa-A inherited") values used for pinned variants.
CURRENT = {
    'erosion_coeff_scale': 1.0,
    'a_scale': 1.0,
    'pyrogen_mass': 0.0024,
    'pyrogen_throat_area': 2.0e-5,
    'T_ignition': 850.0,
    'kappa': 0.45,
}

# Bounds for the pyrogen-based LHS.
BOUNDS = {
    'erosion_coeff_scale': (1.5, 3.5),
    'a_scale':             (0.85, 1.10),
    'pyrogen_mass':        (1e-4, 5e-3),
    'pyrogen_throat_area': (1e-6, 5e-5),
    'T_ignition':          (700.0, 950.0),
    'kappa':               (0.30, 0.60),
}

N_MAIN_SAMPLES = 280
N_PINNED_PER_PARAM = 16   # 16 × 6 = 96 pinned samples
SEED = 42

# Locked simulation kwargs (not swept).
LOCKED_SIM_KWARGS = dict(
    roughness=20e-6,
    cfl_target=0.5,
    dt_max=1e-4,
    t_max=8.0,
    P_cutoff=0.01e6,
    snapshot_interval=10.0,   # effectively no snapshots — saves memory in 392 sims
    print_interval=10.0,
)

# Resample every sim's pressure trace onto this fixed grid before
# storing — keeps the npz small and replots fast.
RESAMPLE_T = np.linspace(0.0, 8.0, 500)


# ============================================================
# Worker (top-level for ProcessPoolExecutor pickling)
# ============================================================

def _run_one_lhs(args):
    """Run one LHS sample and return resampled trace + fitness."""
    idx, params, tag, motor_path, transport_path, locked_kwargs = args
    try:
        motor = load_ric(motor_path)
        gas_props = load_transport(transport_path)

        sim_kwargs = dict(locked_kwargs)
        sim_kwargs['kappa'] = params['kappa']
        sim_kwargs['T_ignition'] = params['T_ignition']

        sim_args = ric_to_sim_args(motor, gas_props=gas_props, **sim_kwargs)
        # erosion_coeff and a are mutated on the loaded objects.
        sim_args['nozzle'].erosion_coeff *= params['erosion_coeff_scale']
        for tab in sim_args['propellant'].tabs:
            tab.a *= params['a_scale']
        sim_args['pyrogen_chamber'] = build_pyrogen_chamber(
            load_pyrogen('bpnv'), sim_args['geo'], sim_args['nozzle'],
            pyrogen_mass=params['pyrogen_mass'],
            pyrogen_throat_area=params['pyrogen_throat_area'],
        )

        geo = sim_args.pop('geo')
        propellant = sim_args.pop('propellant')
        result = run_simulation(geo, propellant, **sim_args)

        # Resample full P_head onto fixed grid for storage.
        P_resampled = np.interp(RESAMPLE_T, result['time'],
                                result['P_head'] / 1e6, left=0.0, right=0.0)

        # Compute MSE fitness vs experimental (after time_offset alignment).
        t_exp = np.asarray(ZEROX_EXPERIMENTAL['time']) \
            + ZEROX_EXPERIMENTAL.get('time_offset', 0.0)
        p_exp = np.asarray(ZEROX_EXPERIMENTAL['pressure'])
        p_sim_at_exp = np.interp(t_exp, result['time'],
                                  result['P_head'] / 1e6, left=0.0, right=0.0)
        mask = t_exp >= 0.0
        fitness = float(np.mean((p_sim_at_exp[mask] - p_exp[mask]) ** 2))
        if not np.isfinite(fitness):
            fitness = 1e6

        return idx, params, tag, P_resampled, fitness, None
    except Exception as e:
        return idx, params, tag, np.zeros_like(RESAMPLE_T), 1e6, str(e)


# ============================================================
# Sample generators
# ============================================================

def build_main_samples(n, seed):
    keys = list(BOUNDS.keys())
    sampler = qmc.LatinHypercube(d=len(keys), seed=seed)
    raw = sampler.random(n=n)
    scaled = qmc.scale(raw,
                       [BOUNDS[k][0] for k in keys],
                       [BOUNDS[k][1] for k in keys])
    return [(dict(zip(keys, row.tolist())), 'main') for row in scaled]


def build_pinned_samples(n_per, seed):
    out = []
    keys = list(BOUNDS.keys())
    for i, pin_key in enumerate(keys):
        free_keys = [k for k in keys if k != pin_key]
        sampler = qmc.LatinHypercube(d=len(free_keys), seed=seed + 1000 * (i + 1))
        raw = sampler.random(n=n_per)
        scaled = qmc.scale(raw,
                           [BOUNDS[k][0] for k in free_keys],
                           [BOUNDS[k][1] for k in free_keys])
        for row in scaled:
            params = dict(zip(free_keys, row.tolist()))
            params[pin_key] = CURRENT[pin_key]
            out.append((params, f'pinned_{pin_key}'))
    return out


# ============================================================
# Save artifacts (belt-and-suspenders order)
# ============================================================

def save_artifacts(results):
    keys = list(BOUNDS.keys())

    # 1) Pickle raw first — minimal failure surface.
    try:
        with open('zerox_lhs_raw.pkl', 'wb') as f:
            pickle.dump(results, f)
        print("Saved zerox_lhs_raw.pkl")
    except Exception as e:
        print(f"WARN: pickle failed: {e}")

    # 2) CSV (params + fitness + tag).
    try:
        with open('zerox_lhs.csv', 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['idx', 'tag', 'fitness', 'error'] + keys)
            for r in results:
                row = [r['idx'], r['tag'], r['fitness'], r['error'] or '']
                row.extend(r['params'][k] for k in keys)
                w.writerow(row)
        print("Saved zerox_lhs.csv")
    except Exception as e:
        print(f"WARN: csv failed: {e}")

    # 3) NPZ with resampled traces.
    try:
        params_arr = np.array([[r['params'][k] for k in keys] for r in results])
        fitness_arr = np.array([r['fitness'] for r in results])
        P_arr = np.stack([r['P_resampled'] for r in results])
        tags_arr = np.array([r['tag'] for r in results])
        np.savez_compressed(
            'zerox_lhs.npz',
            params=params_arr,
            fitness=fitness_arr,
            pressures=P_arr,
            tags=tags_arr,
            time_grid=RESAMPLE_T,
            param_names=np.array(keys),
        )
        print("Saved zerox_lhs.npz (compressed)")
    except Exception as e:
        print(f"WARN: npz failed: {e}")


# ============================================================
# Plot artifacts
# ============================================================

def plot_results(results):
    keys = list(BOUNDS.keys())
    t_exp = np.asarray(ZEROX_EXPERIMENTAL['time']) \
        + ZEROX_EXPERIMENTAL.get('time_offset', 0.0)
    p_exp = np.asarray(ZEROX_EXPERIMENTAL['pressure'])

    valid = [r for r in results if r['error'] is None]
    valid.sort(key=lambda r: r['fitness'])

    # ---------- Plot 1: top-10 overall ----------
    try:
        top_K = 10
        fig, ax = plt.subplots(figsize=(11, 7))
        cmap = plt.cm.viridis(np.linspace(0.10, 0.85, top_K))
        for color, r in zip(cmap, valid[:top_K]):
            ax.plot(RESAMPLE_T, r['P_resampled'], '-', color=color,
                    linewidth=1.6, alpha=0.85,
                    label=f"{r['tag']} MSE={r['fitness']:.3f}")
        ax.plot(t_exp, p_exp, 'k-', linewidth=2.5, marker='o', markersize=4,
                markevery=3, label='Experimental')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Head-End Pressure [MPa]')
        ax.set_title(f'Zerox LHS — Top {top_K} curves overall '
                     f'(MSE in MPa², lower is better)')
        ax.legend(loc='upper right', fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3); ax.set_xlim(0, 8.5); ax.set_ylim(bottom=0)
        plt.tight_layout()
        fig.savefig('zerox_lhs_top.png', dpi=150)
        plt.close(fig)
        print("Saved zerox_lhs_top.png")
    except Exception as e:
        print(f"WARN: top plot failed: {e}")

    # ---------- Plot 2: best per pinned variant ----------
    try:
        fig, axes = plt.subplots(4, 2, figsize=(14, 16))
        axes = axes.flatten()

        # Panel 0 — main 7-var best 5
        main_results = [r for r in valid if r['tag'] == 'main']
        ax = axes[0]
        for r in main_results[:5]:
            ax.plot(RESAMPLE_T, r['P_resampled'], '-', alpha=0.7, linewidth=1.4)
        ax.plot(t_exp, p_exp, 'k-', linewidth=2.5, label='Experimental')
        if main_results:
            ax.set_title(f'Main 7-var (best 5)\nbest MSE={main_results[0]["fitness"]:.3f}',
                         fontsize=11)
        ax.set_xlabel('Time [s]'); ax.set_ylabel('P [MPa]')
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
        ax.set_xlim(0, 8.5); ax.set_ylim(bottom=0)

        for i, key in enumerate(keys):
            ax = axes[i + 1]
            pinned = [r for r in valid if r['tag'] == f'pinned_{key}']
            for r in pinned[:5]:
                ax.plot(RESAMPLE_T, r['P_resampled'], '-', alpha=0.7, linewidth=1.4)
            ax.plot(t_exp, p_exp, 'k-', linewidth=2.5, label='Experimental')
            if pinned:
                ax.set_title(f'{key} pinned at {CURRENT[key]:.4g} (best 5)\n'
                             f'best MSE={pinned[0]["fitness"]:.3f}',
                             fontsize=10)
            else:
                ax.set_title(f'{key} pinned (no samples)', fontsize=10)
            ax.set_xlabel('Time [s]'); ax.set_ylabel('P [MPa]')
            ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
            ax.set_xlim(0, 8.5); ax.set_ylim(bottom=0)

        plt.tight_layout()
        fig.savefig('zerox_lhs_pinned.png', dpi=150)
        plt.close(fig)
        print("Saved zerox_lhs_pinned.png")
    except Exception as e:
        print(f"WARN: pinned plot failed: {e}")

    # ---------- Plot 3: ensemble of top-50 ----------
    try:
        ensemble_K = min(50, len(valid))
        fig, ax = plt.subplots(figsize=(11, 7))
        for r in valid[:ensemble_K]:
            ax.plot(RESAMPLE_T, r['P_resampled'], '-', color='steelblue',
                    linewidth=0.7, alpha=0.30)
        # Highlight best 1
        if valid:
            ax.plot(RESAMPLE_T, valid[0]['P_resampled'], '-', color='red',
                    linewidth=2.0, label=f"Best  MSE={valid[0]['fitness']:.3f}")
        ax.plot(t_exp, p_exp, 'k-', linewidth=2.5, marker='o', markersize=4,
                markevery=3, label='Experimental')
        ax.set_xlabel('Time [s]'); ax.set_ylabel('P [MPa]')
        ax.set_title(f'Zerox LHS — Top {ensemble_K} ensemble '
                     f'(faded blue) + best (red) + experimental (black)')
        ax.legend(loc='upper right'); ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 8.5); ax.set_ylim(bottom=0)
        plt.tight_layout()
        fig.savefig('zerox_lhs_ensemble.png', dpi=150)
        plt.close(fig)
        print("Saved zerox_lhs_ensemble.png")
    except Exception as e:
        print(f"WARN: ensemble plot failed: {e}")

    # ---------- Print best 5 to console ----------
    print(f"\nBest 5 parameter sets (overall):")
    print(f"  {'rank':<5}{'tag':<24}{'MSE':<10}", end="")
    for k in keys:
        print(f"{k[:14]:<16}", end="")
    print()
    for rank, r in enumerate(valid[:5]):
        print(f"  {rank+1:<5}{r['tag']:<24}{r['fitness']:<10.4f}", end="")
        for k in keys:
            print(f"{r['params'][k]:<16.4g}", end="")
        print()


# ============================================================
# Driver
# ============================================================

def main():
    n_workers = os.cpu_count() or 6
    main_samples = build_main_samples(N_MAIN_SAMPLES, SEED)
    pinned_samples = build_pinned_samples(N_PINNED_PER_PARAM, SEED)
    all_samples = main_samples + pinned_samples

    print(f"==== Zerox LHS run ====")
    print(f"  Workers: {n_workers}")
    print(f"  Main 7-var samples:   {len(main_samples)}")
    print(f"  Pinned variants:      {len(BOUNDS)} × {N_PINNED_PER_PARAM} = "
          f"{len(pinned_samples)}")
    print(f"  Total:                {len(all_samples)}")
    est_min = len(all_samples) * 110.0 / n_workers / 60.0
    print(f"  Est. wall time @110s/sim: {est_min:.1f} min")
    print()

    work = [
        (i, params, tag, str(MOTOR_PATH), str(TRANSPORT_PATH), LOCKED_SIM_KWARGS)
        for i, (params, tag) in enumerate(all_samples)
    ]

    t0 = time.time()
    results = []
    completed = 0
    with ProcessPoolExecutor(max_workers=n_workers) as exe:
        futures = [exe.submit(_run_one_lhs, w) for w in work]
        for fut in as_completed(futures):
            idx, params, tag, P_resampled, fitness, err = fut.result()
            results.append({
                'idx': idx, 'params': params, 'tag': tag,
                'P_resampled': P_resampled, 'fitness': fitness, 'error': err,
            })
            completed += 1
            if completed % 20 == 0 or completed == len(work):
                best = min((r['fitness'] for r in results), default=1e6)
                elapsed = time.time() - t0
                eta = elapsed / completed * (len(work) - completed)
                print(f"  {completed}/{len(work)}  best MSE={best:.4f}  "
                      f"elapsed={elapsed/60:.1f}min  eta={eta/60:.1f}min")

    results.sort(key=lambda r: r['idx'])
    print(f"\nAll {len(results)} samples done. Total wall: "
          f"{(time.time()-t0)/60:.1f} min")

    save_artifacts(results)
    plot_results(results)


if __name__ == '__main__':
    main()
