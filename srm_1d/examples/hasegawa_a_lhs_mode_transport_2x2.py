"""
hasegawa_a_lhs_mode_transport_2x2.py — v0.7.4 Phase C.2 LHS sweep

Full 2x2 (heat_delivery_mode x transport-choice):
    cell DF: demar      x frozen
    cell DE: demar      x effective
    cell RF: radiation  x frozen
    cell RE: radiation  x effective

Active cells: all four. The frozen-collapse blocker discovered
during the 2026-05-27 smoke test was resolved later that day by the
v0.7.3.2 architecture fix (Kn-based pyrogen throat sizing +
tighter CFL defaults). See
srm_1d/docs/v0_7_4/references/frozen_collapse_investigation.md
for the resolution and the regression-prevention pytest gate at
srm_1d/tests/test_canonical_examples.py.

Each cell runs an LHS sweep over (roughness, kappa, T_ignition,
k_solid). Saint-Robert (a, n) held at seed (v0.7.4 Phase C.3 will
revisit with Mizushima 2016 BPNV values).

Default: PRE-FLIGHT SMOKE TEST (SRM_LHS_SAMPLES=1, 4 runs total,
~5 minutes). The driver exposes env vars so this scales to the
full 200/cell sweep with a one-variable change.

Pyrogen-form geometry uses Mizushima 2016 industry-standard BKNO3
pellets (3.2 mm cylinder, L/D=1.0).

User-confirmed scope decisions (2026-05-27):
- Cu/Al thermite EXCLUDED (mass-balance hack in YAML would be a
  confounding variable).
- injection_topology = 'forward_plenum' (calibration-canonical for
  Hasegawa A; other topologies need their own per-motor sweeps).
- roughness upper bound extended 50um -> 100um (capture rough
  propellants); smoke test must verify no numerical-collapse.

Env vars:
    SRM_LHS_SAMPLES   samples per cell (default 1 = smoke test)
    SRM_LHS_WORKERS   worker count (default 4 = half of 8 phys cores)
                      Set to 'auto' for multiprocessing.cpu_count().
    SRM_LHS_FITNESS   'segmented' (default) | 'mse'
    SRM_LHS_PROGRESS  'brief' (default) | 'verbose' | 'none'

Usage:
    # Smoke test (default)
    python -m srm_1d.examples.hasegawa_a_lhs_mode_transport_2x2

    # Full sweep at half workers
    SRM_LHS_SAMPLES=200 \\
        python -m srm_1d.examples.hasegawa_a_lhs_mode_transport_2x2

    # Full sweep at 6 workers (75% load)
    SRM_LHS_SAMPLES=200 SRM_LHS_WORKERS=6 \\
        python -m srm_1d.examples.hasegawa_a_lhs_mode_transport_2x2

Outputs (under artifacts/hasegawa_a_lhs_2x2/):
    <cell>/lhs.csv             — all LHS samples for the cell
    <cell>/best_diagnostics.png — top-fit trace + experimental
    comparison.png             — 4-cell top-fit overlay
    rank1_knobs.md             — markdown table of rank-1 knobs
    sanity_check.md            — smoke-test sanity verification log
"""

import os
import multiprocessing
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
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
from srm_1d.plotting import HASEGAWA_MOTOR_A_EXPERIMENTAL


# ============================================================
# Configuration
# ============================================================

MOTORS_DIR = Path(__file__).resolve().parents[1] / 'motors'
MOTOR_PATH = str(MOTORS_DIR / 'hasegawa_a.ric')
FROZEN_TRANSPORT = str(MOTORS_DIR / 'hasegawa_a.frozen.transport.yaml')
EFFECTIVE_TRANSPORT = str(MOTORS_DIR / 'hasegawa_a.transport.yaml')

EXPERIMENTAL_TIME_OFFSET = 0.0

# 8 physical cores; 16 logical CPUs with hyperthreading.
# DEFAULT_WORKERS = 16 saturates all logical CPUs for the long sweep
# (user-confirmed 2026-05-27: "use all workers available"). Drop to
# 8 (= physical core count) for normal foreground use; lower
# (e.g., 4) for smoke tests when keeping the machine responsive.
DEFAULT_WORKERS = 16

# Default samples per cell tuned for ~6 hr wall time at 16 workers
# with the 6-cell grid: per-run ~30-60s wall single-thread, so
# 16 workers ≈ ~1500-3000 samples/hr → 6 hr × 16 cores ≈ 6000-9000
# samples total ≈ 1000-1500 per cell. We target 1000 per cell;
# override via SRM_LHS_SAMPLES for smoke tests or shorter runs.
N_SAMPLES = int(os.environ.get('SRM_LHS_SAMPLES', '1000'))
_WORKERS_ENV = os.environ.get('SRM_LHS_WORKERS')
if _WORKERS_ENV is None:
    MAX_WORKERS = DEFAULT_WORKERS
elif _WORKERS_ENV.lower() == 'auto':
    MAX_WORKERS = multiprocessing.cpu_count()
else:
    MAX_WORKERS = int(_WORKERS_ENV)

# For N=1 smoke test, single-worker keeps the output streams clean
# (no async noise during pre-flight verification).
if N_SAMPLES <= 4:
    MAX_WORKERS = 1

FITNESS_MODE = os.environ.get('SRM_LHS_FITNESS', 'segmented')
PROGRESS_MODE = os.environ.get('SRM_LHS_PROGRESS', 'brief')
SIM_VERBOSE = os.environ.get('SRM_SIM_VERBOSE', '0').lower() in {'1', 'true', 'yes'}

OUTPUT_ROOT = Path('artifacts') / 'hasegawa_a_lhs_2x2'

# Peak-time alignment window — same as the v0.7.1 Phase 5 effective sweep
PEAK_ALIGN_WINDOW = (0.02, 0.18)

# LHS BOUNDS (4-D, user-confirmed Phase C.2 scope)
LHS_BOUNDS = {
    'roughness':  (15e-6, 100e-6),    # 15-100 um (user-extended upper bound)
    'kappa':      (0.40, 0.50),       # narrow band around 0.45 lit center
    'T_ignition': (750.0, 950.0),     # AP/HTPB typical
    'k_solid':    (0.26, 0.32),       # lit-tightened band
}

SEGMENT_WEIGHTS = {
    'mse_spike':      0.25,
    'mse_post_spike': 0.35,
    'mse_plateau':    0.20,
    'mse_taildown':   0.20,
}

# Cell definitions — name, topology, mode, transport_path, label, color.
# 4 forward_plenum cells (the 2x2 mode x transport grid) +
# 2 head_basket cells (frozen vs effective; mode='demar' for both since
#   that's the BPNV YAML default — gives a clean direct comparison
#   against the forward_plenum demar cells).
# v0.7.3.2 architecture fix (Kn-throat + tighter CFL) unblocked
# frozen cells; v0.7.3.3 made mode-axis active for forward_plenum.
CELLS = [
    # cell_name             topology          mode         transport            label                          color
    ('demar_fwd_frozen',    'forward_plenum', 'demar',     FROZEN_TRANSPORT,    'fwd_plenum DeMar x Frozen',     '#1f77b4'),
    ('demar_fwd_effective', 'forward_plenum', 'demar',     EFFECTIVE_TRANSPORT, 'fwd_plenum DeMar x Effective',  '#ff7f0e'),
    ('rad_fwd_frozen',      'forward_plenum', 'radiation', FROZEN_TRANSPORT,    'fwd_plenum Rad x Frozen',       '#2ca02c'),
    ('rad_fwd_effective',   'forward_plenum', 'radiation', EFFECTIVE_TRANSPORT, 'fwd_plenum Rad x Effective',    '#d62728'),
    ('demar_hb_frozen',     'head_basket',    'demar',     FROZEN_TRANSPORT,    'head_basket DeMar x Frozen',    '#9467bd'),
    ('demar_hb_effective',  'head_basket',    'demar',     EFFECTIVE_TRANSPORT, 'head_basket DeMar x Effective', '#8c564b'),
]


# ============================================================
# Per-cell sweep runner
# ============================================================

def run_cell(cell_name, topology, mode, transport_path, label,
             fitness_fn, metrics_fn):
    """Run LHS sweep for one cell of the 6-cell grid."""
    cell_dir = OUTPUT_ROOT / cell_name
    cell_dir.mkdir(parents=True, exist_ok=True)
    csv_path = str(cell_dir / 'lhs.csv')

    # Build the pyrogen object with the cell's heat_delivery_mode
    # override. Sent through to run_from_ric per LHS draw — and
    # importantly, the Pyrogen dataclass pickles cleanly for the
    # ProcessPoolExecutor workers.
    pyrogen_obj = load_pyrogen('bpnv')
    pyrogen_obj.heat_delivery_mode = mode

    transport_basename = Path(transport_path).name
    print(f"\n[cell {cell_name}] topology='{topology}', mode='{mode}', "
          f"transport='{transport_basename}'")
    print(f"  LHS samples = {N_SAMPLES}, workers = {MAX_WORKERS}")
    print(f"  Pyrogen.heat_delivery_mode = {pyrogen_obj.heat_delivery_mode}")
    print(f"  Pyrogen.particle_diameter_m = {pyrogen_obj.particle_diameter_m}")
    print(f"  Output dir = {cell_dir}")

    rows = run_lhs(
        motor_path=MOTOR_PATH,
        bounds=LHS_BOUNDS,
        n_samples=N_SAMPLES,
        fitness_fn=fitness_fn,
        metrics_fn=metrics_fn,
        n_workers=MAX_WORKERS,
        seed=42,
        csv_path=csv_path,
        progress_mode=PROGRESS_MODE,
        sim_verbose=SIM_VERBOSE,
        # Locked sim kwargs
        pyrogen=pyrogen_obj,
        transport_path=transport_path,
        injection_topology=topology,
        t_max=6.0, P_cutoff=0.05e6,
        snapshot_interval=2.0, print_interval=20.0,
    )

    sorted_rows = sorted(rows, key=lambda r: r['fitness'])
    return sorted_rows, cell_dir


def render_best_diagnostics(rank1, t_exp, p_exp, topology, mode,
                            transport_path, cell_dir, label):
    """Render the top-fit pressure trace for this cell."""
    pyrogen_obj = load_pyrogen('bpnv')
    pyrogen_obj.heat_delivery_mode = mode

    params = {k: rank1[k] for k in LHS_BOUNDS.keys()}
    result, perf, *_ = run_from_ric(
        MOTOR_PATH,
        transport_path=transport_path,
        injection_topology=topology,
        t_max=6.0, P_cutoff=0.05e6,
        snapshot_interval=2.0, print_interval=20.0,
        pyrogen=pyrogen_obj,
        verbose=False,
        **params,
    )

    t_sim = result['time']
    p_sim = result['P_head'] / 1e6
    t_off = rank1.get('t_offset_applied_s', 0.0)
    if not np.isfinite(t_off):
        t_off = 0.0
    t_aligned = t_sim + t_off

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(t_exp, p_exp, 'k.-', linewidth=1.8, markersize=3,
            label='experimental')
    ax.plot(t_aligned, p_sim, '-', linewidth=2,
            label=f'sim ({label}, t_off={t_off*1000:+.1f} ms)')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Head pressure [MPa]')
    ax.set_title(f'Hasegawa A — {label} — rank-1 fit (fitness {rank1["fitness"]:.4f})')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = cell_dir / 'best_diagnostics.png'
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path), result


def render_comparison(cell_results, t_exp, p_exp):
    """Render 4-cell top-fit overlay vs experimental."""
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.plot(t_exp, p_exp, 'k.-', linewidth=2.5, markersize=4,
            label='Experimental (Hasegawa)', zorder=10)
    for cell_name, _topology, mode, _tp, label, color in CELLS:
        rank1 = cell_results[cell_name]['rank1']
        result = cell_results[cell_name]['result']
        t_sim = result['time']
        p_sim = result['P_head'] / 1e6
        t_off = rank1.get('t_offset_applied_s', 0.0) or 0.0
        ax.plot(t_sim + t_off, p_sim, '-', linewidth=1.6, color=color,
                label=f'{label} (fit {rank1["fitness"]:.3f})')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Head pressure [MPa]')
    ax.set_title('Hasegawa A — 2x2 mode x transport rank-1 comparison')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = OUTPUT_ROOT / 'comparison.png'
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def render_rank1_knobs_md(cell_results):
    """Write a markdown table of rank-1 knobs per cell."""
    lines = [
        '# v0.7.4 Phase C.2 — rank-1 knobs per cell',
        '',
        f'**Samples per cell**: {N_SAMPLES}',
        f'**Workers**: {MAX_WORKERS}',
        '',
        '| Cell | Mode | Transport | fitness | roughness [um] | kappa | T_ignition [K] | k_solid [W/m/K] | t_offset [ms] |',
        '|------|------|-----------|---------|----------------|-------|----------------|-----------------|----------------|',
    ]
    for cell_name, topology, mode, transport_path, label, _color in CELLS:
        r = cell_results[cell_name]['rank1']
        t_off = r.get('t_offset_applied_s', 0.0) or 0.0
        lines.append(
            f"| {cell_name} | {mode} | {Path(transport_path).name} | "
            f"{r['fitness']:.4f} | {r['roughness']*1e6:.1f} | "
            f"{r['kappa']:.3f} | {r['T_ignition']:.0f} | "
            f"{r['k_solid']:.3f} | {t_off*1000:+.1f} |"
        )
    path = OUTPUT_ROOT / 'rank1_knobs.md'
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return str(path)


def render_sanity_check(cell_results):
    """Write a sanity-check log for the smoke test (per plan).

    Verifies (1) modes actually toggle (peak P differs),
    (2) peak-time alignment computed, (3) transport YAMLs loaded
    correctly, (4) no numerical collapse.
    """
    lines = [
        '# v0.7.4 Phase C.2 Smoke-Test Sanity Check',
        '',
        f'N_SAMPLES per cell: {N_SAMPLES}  (smoke-test mode if <= 4)',
        f'Workers: {MAX_WORKERS}',
        '',
        '## (1) Mode toggle — peak P should differ between modes',
        '',
        '| Cell | mode | transport | P_peak [MPa] | t_peak [s] | termination |',
        '|------|------|-----------|--------------|------------|-------------|',
    ]
    peak_p_per_cell = {}
    for cell_name, topology, mode, transport_path, label, _color in CELLS:
        result = cell_results[cell_name]['result']
        summary = result['summary']
        peak_p = summary['P_peak'] / 1e6
        t_peak = summary['t_peak']
        term_code = summary.get('termination_code', '?')
        peak_p_per_cell[cell_name] = peak_p
        lines.append(
            f"| {cell_name} | {mode} | {Path(transport_path).name} | "
            f"{peak_p:.3f} | {t_peak:.3f} | code={term_code} |"
        )

    # Mode toggle verification — works for any non-empty set of cells.
    demar_peaks = [v for k, v in peak_p_per_cell.items() if k.startswith('demar')]
    rad_peaks   = [v for k, v in peak_p_per_cell.items() if k.startswith('radiation')]
    if demar_peaks and rad_peaks:
        demar_avg = float(np.mean(demar_peaks))
        rad_avg = float(np.mean(rad_peaks))
        demar_rad_diff = abs(demar_avg - rad_avg) / max(demar_avg, 1e-6)
    else:
        demar_rad_diff = float('nan')
    lines += [
        '',
        f"Mode-toggle |Demar avg - Radiation avg| / Demar avg = {demar_rad_diff:.2%}",
        '',
        '✓ HEALTHY if > 5%; ✗ SUSPICIOUS if < 1% (toggle may be a no-op).',
        '',
        '## (2) Peak-time alignment — t_offset should be nonzero',
        '',
        '| Cell | t_offset_applied_s [ms] |',
        '|------|--------------------------|',
    ]
    nonzero_offsets = 0
    for cell_name, *_ in CELLS:
        rank1 = cell_results[cell_name]['rank1']
        t_off = rank1.get('t_offset_applied_s', 0.0) or 0.0
        if abs(t_off) > 1e-9:
            nonzero_offsets += 1
        lines.append(f"| {cell_name} | {t_off*1000:+.2f} |")

    lines += [
        '',
        f"Nonzero t_offsets: {nonzero_offsets} / {len(CELLS)} cells.",
        '',
        '✓ HEALTHY if at least 1 cell has nonzero t_offset.',
        '✗ SUSPICIOUS if all are 0 — alignment may be silently falling back.',
        '',
        '## (3) Transport YAML resolution',
        '',
        'Each cell printed the resolved transport path in its run header. ',
        'Confirm the _frozen cells used `hasegawa_a.frozen.transport.yaml` and ',
        'the _effective cells used `hasegawa_a.transport.yaml`.',
        '',
        '## (4) Numerical-collapse check',
        '',
        '| Cell | termination_code | OK? |',
        '|------|------------------|-----|',
    ]
    n_collapse = 0
    for cell_name, *_ in CELLS:
        summary = cell_results[cell_name]['result']['summary']
        term = summary.get('termination_code', '?')
        ok = '✓' if term != 4 else '✗ COLLAPSE'
        if term == 4:
            n_collapse += 1
        lines.append(f"| {cell_name} | {term} | {ok} |")
    lines += [
        '',
        f"Numerical collapses: {n_collapse} / {len(CELLS)} cells.",
        '',
        '✓ HEALTHY if all cells terminated cleanly (code != 4).',
        '✗ If collapses occurred, narrow `roughness` upper bound back to 50e-6 m.',
        '',
        '## Roughness bound observation',
        '',
        f'LHS roughness bounds = ({LHS_BOUNDS["roughness"][0]*1e6:.0f}, '
        f'{LHS_BOUNDS["roughness"][1]*1e6:.0f}) um.',
        '',
    ]
    path = OUTPUT_ROOT / 'sanity_check.md'
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return str(path)


# ============================================================
# Main
# ============================================================

def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    t_exp = HASEGAWA_MOTOR_A_EXPERIMENTAL['time'] + EXPERIMENTAL_TIME_OFFSET
    p_exp = HASEGAWA_MOTOR_A_EXPERIMENTAL['pressure']

    metrics_fn = pressure_trace_metrics(
        t_exp, p_exp, t_min=0.01,
        segments=DEFAULT_PRESSURE_SEGMENTS,
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
        raise ValueError("SRM_LHS_FITNESS must be 'segmented' or 'mse'")

    print('=' * 60)
    print('v0.7.4 Phase C.2 — Hasegawa A 2x2 LHS sweep')
    print(f'  Mode x Transport = DeMar/Radiation x Frozen/Effective')
    print(f'  N_SAMPLES per cell = {N_SAMPLES}')
    print(f'  Workers = {MAX_WORKERS}')
    print(f'  Fitness = {FITNESS_MODE}')
    print(f'  Bounds = {LHS_BOUNDS}')
    print('=' * 60)

    cell_results = {}
    for cell_name, topology, mode, transport_path, label, _color in CELLS:
        sorted_rows, cell_dir = run_cell(
            cell_name, topology, mode, transport_path, label,
            fitness_fn, metrics_fn,
        )
        rank1 = sorted_rows[0]
        # Render the top-fit pressure trace + capture the result
        # dict for the comparison overlay.
        _diag_path, result = render_best_diagnostics(
            rank1, t_exp, p_exp, topology, mode, transport_path,
            cell_dir, label,
        )
        cell_results[cell_name] = {
            'rank1': rank1,
            'result': result,
            'all_rows': sorted_rows,
        }
        print(f"  [cell {cell_name}] rank-1 fitness = {rank1['fitness']:.4f}, "
              f"P_peak = {result['summary']['P_peak']/1e6:.2f} MPa, "
              f"t_peak = {result['summary']['t_peak']:.3f} s")

    # Comparison overlay
    cmp_path = render_comparison(cell_results, t_exp, p_exp)
    print(f"\nComparison plot: {cmp_path}")

    knobs_path = render_rank1_knobs_md(cell_results)
    print(f"Rank-1 knobs table: {knobs_path}")

    if N_SAMPLES <= 4:
        sanity_path = render_sanity_check(cell_results)
        print(f"Smoke-test sanity check: {sanity_path}")
        print()
        print(f"Pre-flight smoke test complete. Review the {len(CELLS)} traces in")
        print(f"  {OUTPUT_ROOT}/comparison.png")
        print("and the sanity-check log before scheduling the full sweep.")
    else:
        print()
        print(f"Full sweep complete ({N_SAMPLES} samples per cell, "
              f"{N_SAMPLES * len(CELLS)} total runs).")


if __name__ == '__main__':
    main()
