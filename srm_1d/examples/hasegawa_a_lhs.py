"""
hasegawa_a_lhs.py — Hasegawa Motor A Latin Hypercube optimization.
==================================================================

LHS sweep over Ma erosive-burning and v0.7.0 pyrogen ignition
parameters, fitting the Hasegawa et al. (2006) head-end pressure trace.
The default fitness is segmented so the spike, post-spike shoulder,
plateau, and taildown are visible as separate diagnostics.

Usage:
    python -m srm_1d.examples.hasegawa_a_lhs

By default runs N=500 samples across all CPU cores. Override
``SRM_HASEGAWA_LHS_SAMPLES`` or ``SRM_HASEGAWA_LHS_WORKERS`` for shorter
test runs. Set ``SRM_LHS_PROGRESS=brief|verbose|none`` and
``SRM_SIM_VERBOSE=1`` to restore per-run solver summary blocks.

Output:
    artifacts/hasegawa_a_lhs/hasegawa_a_lhs.csv
    artifacts/hasegawa_a_lhs/hasegawa_a_lhs_top5.png
    artifacts/hasegawa_a_lhs/hasegawa_a_lhs_metrics.png
    artifacts/hasegawa_a_lhs/hasegawa_a_lhs_diagnostics.png
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.tools.sensitivity import (
    DEFAULT_PRESSURE_SEGMENTS,
    run_lhs,
    mse_fitness,
    pressure_trace_metrics,
    segmented_pressure_fitness,
)
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import HASEGAWA_MOTOR_A_EXPERIMENTAL


MOTOR_PATH = str(Path(__file__).resolve().parents[1] / 'motors' / 'hasegawa_a.ric')
EXPERIMENTAL_TIME_OFFSET = 0.0  # align experimental ignition with sim t=0
N_SAMPLES = int(os.environ.get('SRM_HASEGAWA_LHS_SAMPLES', '500'))
MAX_WORKERS = os.environ.get('SRM_HASEGAWA_LHS_WORKERS')
MAX_WORKERS = None if MAX_WORKERS in (None, '', 'auto') else int(MAX_WORKERS)
DEFAULT_OUTPUT_PREFIX = str(
    Path('artifacts') / 'hasegawa_a_lhs' / 'hasegawa_a_lhs'
)
OUTPUT_PREFIX = os.environ.get('SRM_HASEGAWA_LHS_PREFIX', DEFAULT_OUTPUT_PREFIX)
FITNESS_MODE = os.environ.get('SRM_HASEGAWA_LHS_FITNESS', 'segmented')
PROGRESS_MODE = os.environ.get('SRM_LHS_PROGRESS', 'brief')
SIM_VERBOSE = os.environ.get('SRM_SIM_VERBOSE', '0').lower() in {'1', 'true', 'yes'}

SEGMENT_WEIGHTS = {
    'mse_spike': 0.25,
    'mse_post_spike': 0.35,
    'mse_plateau': 0.20,
    'mse_taildown': 0.20,
}


def _plot_metric_tradeoffs(rows, prefix):
    valid = [r for r in rows if not r.get('error')]
    if not valid:
        return

    def arr(key):
        return np.array([float(r.get(key, np.nan)) for r in valid])

    def set_robust_limits(ax, x_values, y_values, lo=0.0):
        for setter, values in ((ax.set_xlim, x_values), (ax.set_ylim, y_values)):
            finite = values[np.isfinite(values)]
            if finite.size == 0:
                continue
            hi = np.percentile(finite, 90.0) if finite.size > 8 else np.max(finite)
            if np.isfinite(hi) and hi > lo:
                setter(lo, hi * 1.15)

    fitness = arr('fitness')
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    post = arr('mse_post_spike')
    plateau = arr('mse_plateau')
    sc = axes[0, 0].scatter(
        post, plateau, c=fitness,
        cmap='viridis_r', s=24, alpha=0.85,
    )
    axes[0, 0].set_xlabel('Post-spike MSE [MPa^2]')
    axes[0, 0].set_ylabel('Plateau MSE [MPa^2]')
    axes[0, 0].set_title('Post-spike vs plateau')
    set_robust_limits(axes[0, 0], post, plateau)
    fig.colorbar(sc, ax=axes[0, 0], label='fitness')

    tail = arr('mse_taildown')
    axes[0, 1].scatter(tail, plateau,
                       c=fitness, cmap='viridis_r', s=24, alpha=0.85)
    axes[0, 1].set_xlabel('Taildown MSE [MPa^2]')
    axes[0, 1].set_ylabel('Plateau MSE [MPa^2]')
    axes[0, 1].set_title('Taildown vs plateau')
    set_robust_limits(axes[0, 1], tail, plateau)

    duration = arr('pyrogen_duration_ms')
    axes[1, 0].scatter(duration, post,
                       c=arr('pyrogen_peak_P_MPa'), cmap='plasma',
                       s=24, alpha=0.85)
    axes[1, 0].set_xlabel('Pyrogen duration [ms]')
    axes[1, 0].set_ylabel('Post-spike MSE [MPa^2]')
    axes[1, 0].set_title('Igniter duration vs shoulder error')
    set_robust_limits(axes[1, 0], duration, post)

    axes[1, 1].scatter(arr('peak_error_pct'), arr('trough_error_pct'),
                       c=fitness, cmap='viridis_r', s=24, alpha=0.85)
    axes[1, 1].axhline(0.0, color='0.5', linewidth=0.8)
    axes[1, 1].axvline(0.0, color='0.5', linewidth=0.8)
    axes[1, 1].set_xlabel('Spike peak error [%]')
    axes[1, 1].set_ylabel('Post-spike trough error [%]')
    axes[1, 1].set_title('Spike amplitude vs shoulder floor')

    for ax in axes.flat:
        ax.grid(True, alpha=0.3)
    fig.suptitle('Hasegawa A segmented metric tradeoffs')
    fig.tight_layout()
    path = f'{prefix}_metrics.png'
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f'Saved {path}')


def _plot_best_diagnostics(result, t_exp, p_exp, prefix, t_offset=0.0):
    """Best-run diagnostic plot.

    v0.7.1 Phase 5: takes ``t_offset`` (seconds) so the sim trace and
    its residual are drawn on the SAME aligned time axis the fitness
    function scored against. Without this, the plotted sim trace would
    show the un-aligned (raw) ignition timing and the residual would
    bake in the ignition-phasing error — misleading because the
    optimizer judged the run by the aligned residual instead.
    """
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=False)

    t_raw = result['time']
    t_aligned = t_raw + t_offset
    p_sim = result['P_head'] / 1e6
    p_at_exp = np.interp(t_exp, t_aligned, p_sim)

    axes[0].plot(t_aligned, p_sim, 'b-', linewidth=2,
                 label=f'simulation (aligned, t_offset={t_offset*1000:+.1f} ms)')
    axes[0].plot(t_raw, p_sim, color='0.7', linewidth=0.8, alpha=0.7,
                 label='simulation (raw, pre-alignment)')
    axes[0].plot(t_exp, p_exp, 'ko-', linewidth=1.5, markersize=3,
                 label='experimental')
    axes[0].set_ylabel('P_head [MPa]')
    axes[0].legend(loc='best', fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t_exp, p_at_exp - p_exp, 'r.-', linewidth=1.2,
                 label='aligned residual (what the optimizer scored)')
    axes[1].axhline(0.0, color='0.4', linewidth=0.8)
    axes[1].set_ylabel('Residual [MPa]')
    axes[1].legend(loc='best', fontsize=9)
    axes[1].grid(True, alpha=0.3)

    ax_pig = axes[2]
    ax_mdot = ax_pig.twinx()
    ax_pig.plot(t_aligned, result['P_ig'] / 1e6, color='tab:purple',
                linewidth=1.6, label='P_ig')
    ax_mdot.plot(t_aligned, result['mdot_ig'] * 1000.0, color='tab:orange',
                 linewidth=1.2, label='mdot_ig')
    ax_pig.set_ylabel('P_ig [MPa]')
    ax_mdot.set_ylabel('mdot_ig [g/s]')
    ax_pig.grid(True, alpha=0.3)
    lines, labels = ax_pig.get_legend_handles_labels()
    more_lines, more_labels = ax_mdot.get_legend_handles_labels()
    ax_pig.legend(lines + more_lines, labels + more_labels, loc='best')

    snapshots = result.get('snapshots', [])
    if snapshots:
        snap_t = np.array([s['t'] for s in snapshots])
        ign_frac = np.array([
            np.mean(s['is_burning'][s['is_grain']]) if np.any(s['is_grain']) else 0.0
            for s in snapshots
        ])
        max_tsurf = np.array([np.max(s['T_surf']) for s in snapshots])
        ax_frac = axes[3]
        ax_tsurf = ax_frac.twinx()
        ax_frac.plot(snap_t, ign_frac, 'g-', linewidth=1.8,
                     label='burning grain fraction')
        ax_tsurf.plot(snap_t, max_tsurf, color='tab:red',
                      linewidth=1.2, label='max T_surf')
        ax_frac.set_ylabel('Burning fraction [-]')
        ax_tsurf.set_ylabel('max T_surf [K]')
        lines, labels = ax_frac.get_legend_handles_labels()
        more_lines, more_labels = ax_tsurf.get_legend_handles_labels()
        ax_frac.legend(lines + more_lines, labels + more_labels, loc='best')
    axes[3].set_xlabel('Time [s]')
    axes[3].grid(True, alpha=0.3)

    fig.suptitle('Hasegawa A best-run diagnostics')
    fig.tight_layout()
    path = f'{prefix}_diagnostics.png'
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f'Saved {path}')


def main():
    Path(OUTPUT_PREFIX).parent.mkdir(parents=True, exist_ok=True)

    t_exp = HASEGAWA_MOTOR_A_EXPERIMENTAL['time'] + EXPERIMENTAL_TIME_OFFSET
    p_exp = HASEGAWA_MOTOR_A_EXPERIMENTAL['pressure']

    # v0.7.1 Phase 5 bounds (k_solid-tightened diagnostic sweep,
    # 2026-05-22). The prior literature-bounded re-sweep produced a
    # top-5 ensemble with fitted k_solid spanning 0.175-0.776 W/(m·K),
    # a 4.4x range for a single propellant composition. Subsequent
    # AP/HTPB k_solid literature dive (memory:
    # reference_ap_htpb_k_solid_literature) established that published
    # cold-measurement bounds for AP/HTPB+Al cluster in 0.20-0.40,
    # centered ~0.25-0.30. No source supports k > 0.6 except as a
    # free fit parameter compensating for absent physics.
    #
    # This sweep tightens k_solid to the full literature-defensible
    # band (0.20-0.40) as a DIAGNOSTIC experiment:
    #
    #   - If best fitness lands near the prior 0.0644 (within ~20%),
    #     the model is honest and previous calibration is recoverable.
    #   - If best fitness rises sharply (>~0.1), confirmed model gap.
    #     The right response is NOT to re-widen k_solid; it's to
    #     slate Al2O3 two-phase lag (Pardue 1992) and/or Z-N dynamic
    #     burn rate as v0.7.2 priority work.
    #
    # All other bounds inherited from the prior sweep (commit af49a68):
    #   roughness 5-100 um, heat_flux 30-200 cal/cm^2/s, etc.
    bounds = {
        # Ma erosive-burning knobs
        'roughness':           (5e-6, 100e-6),
        'kappa':               (0.30, 0.60),
        # Pyrogen plenum sizing
        'pyrogen_mass':        (0.001, 0.050),
        'pyrogen_throat_area': (1e-6, 5e-5),
        'pyrogen_volume':      (1e-6, 1.5e-5),
        'pyrogen_heat_flux_cal_cm2_s': (30.0, 200.0),
        # Goodman ignition / surface-conduction
        'T_ignition':          (800.0, 1100.0),
        'k_solid':             (0.20, 0.40),
    }

    # v0.7.1 Phase 5: per-sample ignition-timing alignment. The fitness
    # finds the sim and experimental peak times inside this window
    # before computing MSE — keeps the optimizer focused on trace
    # SHAPE rather than ignition PHASING. Hasegawa A's global peak is
    # the ignition spike; the (0.02, 0.18) window scoops it cleanly.
    PEAK_ALIGN_WINDOW = (0.02, 0.18)

    metrics_fn = pressure_trace_metrics(
        t_exp, p_exp, t_min=0.01, segments=DEFAULT_PRESSURE_SEGMENTS,
        peak_align_window=PEAK_ALIGN_WINDOW,
    )
    if FITNESS_MODE == 'mse':
        fitness_fn = mse_fitness(
            t_exp, p_exp, t_min=0.01,
            peak_align_window=PEAK_ALIGN_WINDOW,
        )
    elif FITNESS_MODE == 'segmented':
        fitness_fn = segmented_pressure_fitness(
            t_exp, p_exp, t_min=0.01,
            segments=DEFAULT_PRESSURE_SEGMENTS,
            weights=SEGMENT_WEIGHTS,
            peak_align_window=PEAK_ALIGN_WINDOW,
        )
    else:
        raise ValueError("SRM_HASEGAWA_LHS_FITNESS must be 'segmented' or 'mse'")

    rows = run_lhs(
        motor_path=MOTOR_PATH,
        bounds=bounds,
        n_samples=N_SAMPLES,
        fitness_fn=fitness_fn,
        metrics_fn=metrics_fn,
        n_workers=MAX_WORKERS,
        seed=42,
        csv_path=f'{OUTPUT_PREFIX}.csv',
        progress_mode=PROGRESS_MODE,
        sim_verbose=SIM_VERBOSE,
        # Locked sim kwargs (kappa now in the LHS bounds dict above)
        pyrogen='bpnv',
        t_max=6.0, P_cutoff=0.05e6,
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
        print(f"Rank {rank} (fitness: {r['fitness']:.4f}):")
        print(f"  Roughness    = {r['roughness']*1e6:.1f} um")
        print(f"  Kappa        = {r['kappa']:.3f}")
        print(f"  Pyro Mass    = {r['pyrogen_mass']*1000:.1f} g")
        print(f"  Pyro Throat  = {r['pyrogen_throat_area']*1e6:.2f} mm^2")
        print(f"  Pyro Volume  = {r['pyrogen_volume']*1e6:.1f} cm^3")
        print(f"  Pyro HeatFlx = {r['pyrogen_heat_flux_cal_cm2_s']:.1f} cal/cm^2/s")
        print(f"  T_ignition   = {r['T_ignition']:.0f} K")
        print(f"  k_solid      = {r['k_solid']:.3f} W/(m.K)")
        print(f"  MSE segments = spike {r.get('mse_spike', np.nan):.3f}, "
              f"post {r.get('mse_post_spike', np.nan):.3f}, "
              f"plateau {r.get('mse_plateau', np.nan):.3f}, "
              f"tail {r.get('mse_taildown', np.nan):.3f}")
        print(f"  Peak/trough  = {r.get('peak_error_pct', np.nan):+.1f}% / "
              f"{r.get('trough_error_pct', np.nan):+.1f}%")
        t_off = r.get('t_offset_applied_s', np.nan)
        if np.isfinite(t_off):
            print(f"  t_offset     = {t_off*1000.0:+.1f} ms (sim peak shifted to match exp)")
        print("-" * 30)

    _plot_metric_tradeoffs(rows, OUTPUT_PREFIX)

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
            t_max=6.0, P_cutoff=0.05e6,
            snapshot_interval=2.0, print_interval=20.0,
            pyrogen='bpnv',
            verbose=SIM_VERBOSE,
            **params,
        )
        # v0.7.1 Phase 5: draw each rank on its aligned time axis (the
        # axis the fitness function judged) so visual comparison
        # matches the reported MSE. The offset was captured in the
        # row's metrics dict during the LHS sweep.
        t_offset_rank = float(r.get('t_offset_applied_s', 0.0) or 0.0)
        plt.plot(result['time'] + t_offset_rank, result['P_head'] / 1e6,
                 color=colors[rank], linewidth=1.5, alpha=0.9,
                 label=f"Rank {rank+1} | fitness={r['fitness']:.3f} | "
                       f"t_off={t_offset_rank*1000:+.0f} ms")

    plt.title(f"Hasegawa A — {len(bounds)}-Variable LHS (Top 5, N={N_SAMPLES})",
              fontsize=16)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("Head-End Pressure [MPa]", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 5.5)
    plt.legend(loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_PREFIX}_top5.png", dpi=300)
    plt.close()
    print(f"\nSaved {OUTPUT_PREFIX}_top5.png")

    if sorted_rows:
        best_params = {k: sorted_rows[0][k] for k in bounds.keys()}
        best_result, *_ = run_from_ric(
            MOTOR_PATH,
            t_max=6.0, P_cutoff=0.05e6,
            snapshot_interval=0.02, print_interval=20.0,
            pyrogen='bpnv',
            verbose=SIM_VERBOSE,
            **best_params,
        )
        best_t_offset = float(sorted_rows[0].get('t_offset_applied_s', 0.0) or 0.0)
        _plot_best_diagnostics(
            best_result, t_exp, p_exp, OUTPUT_PREFIX,
            t_offset=best_t_offset,
        )


if __name__ == "__main__":
    main()
