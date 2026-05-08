"""
hasegawa_a_lhs.py — Hasegawa Motor A Latin Hypercube optimization.
==================================================================

5-variable LHS sweep over the Ma erosive-burning + ignition parameters,
fitting head-end pressure trace MSE against Hasegawa et al. (2006)
experimental data. The Rank-1 result of an N=500 run of this script
is the v0.6.0 calibration baseline (see DEVNOTES "Calibration State").

Usage:
    python -m srm_1d.examples.hasegawa_a_lhs

By default runs N=500 samples across all CPU cores. Override
``N_SAMPLES`` or ``MAX_WORKERS`` below for shorter test runs.

Output:
    hasegawa_a_lhs.csv — one row per sample, all params + fitness
    hasegawa_a_lhs_top5.png — top-5 trace overlay
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.tools.sensitivity import run_lhs, mse_fitness
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import HASEGAWA_MOTOR_A_EXPERIMENTAL


MOTOR_PATH = str(Path(__file__).resolve().parents[1] / 'motors' / 'hasegawa_a.ric')
EXPERIMENTAL_TIME_OFFSET = 0.02  # align experimental ignition with sim t=0
N_SAMPLES = 500
MAX_WORKERS = None  # default: os.cpu_count()


def main():
    t_exp = HASEGAWA_MOTOR_A_EXPERIMENTAL['time'] + EXPERIMENTAL_TIME_OFFSET
    p_exp = HASEGAWA_MOTOR_A_EXPERIMENTAL['pressure']

    bounds = {
        'roughness':         (5e-6, 50e-6),
        'igniter_mass':      (0.001, 0.050),
        'ignition_ramp_tau': (0.001, 0.10),
        'P_ignition':        (0.005e6, 0.1e6),
        'igniter_tau':       (0.001, 0.20),
    }

    rows = run_lhs(
        motor_path=MOTOR_PATH,
        bounds=bounds,
        n_samples=N_SAMPLES,
        fitness_fn=mse_fitness(t_exp, p_exp, t_min=0.01),
        n_workers=MAX_WORKERS,
        seed=42,
        csv_path='hasegawa_a_lhs.csv',
        # Locked sim kwargs
        kappa=0.45, t_max=6.0, P_cutoff=0.05e6,
        snapshot_interval=2.0, print_interval=20.0,
    )

    # ============================================================
    # Top-5 summary + plot
    # ============================================================
    sorted_rows = sorted(rows, key=lambda r: r['fitness'])
    print()
    print("=" * 50)
    print("--- TOP 5 BEST FITS ---")
    print("=" * 50)
    for rank, r in enumerate(sorted_rows[:5], start=1):
        print(f"Rank {rank} (MSE: {r['fitness']:.4f}):")
        print(f"  Roughness    = {r['roughness']*1e6:.1f} μm")
        print(f"  Ign Mass     = {r['igniter_mass']*1000:.1f} g")
        print(f"  Ign Ramp Tau = {r['ignition_ramp_tau']*1000:.1f} ms")
        print(f"  P_ignition   = {r['P_ignition']/1e6:.3f} MPa")
        print(f"  Ign Tau      = {r['igniter_tau']*1000:.1f} ms")
        print("-" * 30)

    # Re-run the top-5 to recapture full traces (results aren't stored
    # to keep memory bounded during the sweep)
    plt.figure(figsize=(14, 9))
    plt.plot(t_exp, p_exp, 'k.-', linewidth=2.5, zorder=10,
             label='Experimental (Hasegawa)')
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for rank, r in enumerate(sorted_rows[:5]):
        params = {k: r[k] for k in bounds.keys()}
        result, *_ = run_from_ric(
            MOTOR_PATH,
            kappa=0.45, t_max=6.0, P_cutoff=0.05e6,
            snapshot_interval=2.0, print_interval=20.0,
            **params,
        )
        plt.plot(result['time'], result['P_head'] / 1e6,
                 color=colors[rank], linewidth=1.5, alpha=0.9,
                 label=f"Rank {rank+1} | MSE={r['fitness']:.3f}")

    plt.title(f"Hasegawa A — {len(bounds)}-Variable LHS (Top 5, N={N_SAMPLES})",
              fontsize=16)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("Head-End Pressure [MPa]", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 5.5)
    plt.legend(loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.savefig("hasegawa_a_lhs_top5.png", dpi=300)
    plt.close()
    print("\nSaved hasegawa_a_lhs_top5.png")


if __name__ == "__main__":
    main()
