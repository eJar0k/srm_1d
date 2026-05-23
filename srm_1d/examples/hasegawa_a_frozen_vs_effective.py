"""
hasegawa_a_frozen_vs_effective.py — v0.7.1 Phase 5 Task 1 A/B test
====================================================================

Runs Hasegawa Motor A twice — once with the FROZEN RPA transport
(`hasegawa_a.frozen.transport.yaml`, k=0.3685, Cp=2060) and once with
the EFFECTIVE RPA pair (`hasegawa_a.transport.yaml`, k=0.6517, Cp=2764).
Both runs use the same v0.7.1 Phase 5 full3_kbound LHS rank-1
calibration parameters.

Originally written to test whether the LHS optimizer's habit of pegging
k_solid against its lower bound was the model compensating for
under-heat-transfer through a frozen-k_gas film. v0.7.1 Phase 5
close-out: confirmed YES, switching to effective shifts k_solid to the
literature center but exposes a structural 11% ignition-spike
under-prediction (v0.7.2 work).

NOTE: as of v0.7.1, EFFECTIVE is the default. The auto-resolved sibling
of hasegawa_a.ric is `hasegawa_a.transport.yaml` which now contains
effective values. The frozen sibling `hasegawa_a.frozen.transport.yaml`
is preserved for this diagnostic.

Outputs:
    artifacts/hasegawa_a_freeff/<stamp>/
        comparison.png   — both sim traces overlaid on experimental
        metrics.txt      — P_peak / t_peak / mse for each run
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import HASEGAWA_MOTOR_A_EXPERIMENTAL
from srm_1d.run_artifacts import artifact_dir


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'hasegawa_a.ric'
FROZEN_YAML = Path(__file__).resolve().parents[1] / 'motors' / 'hasegawa_a.frozen.transport.yaml'
EFFECTIVE_YAML = Path(__file__).resolve().parents[1] / 'motors' / 'hasegawa_a.transport.yaml'

# Full3_kbound LHS rank-1 (idx 110, fitness 0.0635) — literature-bounded
# k_solid sweep, the best of the six v0.7.1 Phase 5 sweeps.
LHS_PARAMS = dict(
    roughness=3.0e-5,                      # 30 µm
    kappa=0.4464,
    pyrogen='bpnv',
    pyrogen_mass=9.455e-3,                 # 9.46 g
    pyrogen_throat_area=1.4e-5,            # 14 mm²
    pyrogen_volume=1.1e-5,                 # 11 cm³
    pyrogen_heat_flux_cal_cm2_s=39.155,
    T_ignition=837.65,
    k_solid=0.2058,
)

EXPERIMENTAL_TIME_OFFSET = 1.1  # matches hasegawa_motor_a.py


def run_with_transport(label, transport_path):
    print(f"\n=== Running {label} ({transport_path.name}) ===")
    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        transport_path=str(transport_path),
        P_cutoff=0.05e6,
        snapshot_interval=0.5,
        print_interval=0.5,
        verbose=False,
        **LHS_PARAMS,
    )
    return result, perf


def summarize(label, result):
    t = result['time']
    P = result['P_head'] / 1e6  # MPa
    P_peak = float(P.max())
    t_peak = float(t[int(P.argmax())])
    return dict(label=label, P_peak=P_peak, t_peak=t_peak, t=t, P=P)


def main():
    frozen = summarize('frozen', run_with_transport('frozen', FROZEN_YAML)[0])
    effective = summarize('effective', run_with_transport('effective', EFFECTIVE_YAML)[0])

    out = artifact_dir('hasegawa_a_freeff')
    print(f"\nArtifacts: {out}")

    # ----- Comparison plot --------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 7))

    exp = HASEGAWA_MOTOR_A_EXPERIMENTAL
    ax.plot(exp['time'] + EXPERIMENTAL_TIME_OFFSET, exp['pressure'],
            'k-', linewidth=1.5, marker='o', markersize=3,
            markevery=max(1, len(exp['time']) // 25),
            label=exp.get('label', 'Hasegawa 2006 (digitized)'))

    ax.plot(frozen['t'], frozen['P'], 'b-', linewidth=2,
            label=f"Frozen k=0.3685, Cp=2060 (P_peak={frozen['P_peak']:.2f} MPa)")
    ax.plot(effective['t'], effective['P'], 'r-', linewidth=2,
            label=f"Effective k=0.6517, Cp=2764 (P_peak={effective['P_peak']:.2f} MPa)")

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Head-End Pressure [MPa]', fontsize=12)
    ax.set_title(
        "Hasegawa A — Frozen vs Effective Transport (v0.7.1 Phase 5 Task 1)\n"
        "Both runs: full3_kbound rank-1 params (roughness=30µm, kappa=0.446, "
        "k_solid=0.206, 9.5g BPNV)",
        fontsize=11,
    )
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=10)

    comparison_path = out / 'comparison.png'
    fig.savefig(comparison_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {comparison_path}")

    # ----- Metrics text -----------------------------------------------------
    lines = [
        "Hasegawa A — Frozen vs Effective Transport (v0.7.1 Phase 5 Task 1)",
        "=" * 72,
        "",
        f"Params (full3_kbound rank-1, fitness=0.0635 against frozen YAML):",
        f"  roughness            = {LHS_PARAMS['roughness']*1e6:.1f} µm",
        f"  kappa                = {LHS_PARAMS['kappa']}",
        f"  pyrogen_mass         = {LHS_PARAMS['pyrogen_mass']*1000:.2f} g",
        f"  pyrogen_throat_area  = {LHS_PARAMS['pyrogen_throat_area']*1e6:.2f} mm²",
        f"  pyrogen_volume       = {LHS_PARAMS['pyrogen_volume']*1e6:.2f} cm³",
        f"  pyrogen_heat_flux    = {LHS_PARAMS['pyrogen_heat_flux_cal_cm2_s']:.1f} cal/cm²/s",
        f"  T_ignition           = {LHS_PARAMS['T_ignition']} K",
        f"  k_solid              = {LHS_PARAMS['k_solid']}",
        "",
        f"{'Case':<12} {'k_gas':>10} {'Cp_gas':>10} {'P_peak [MPa]':>14} {'t_peak [s]':>12}",
        "-" * 60,
        f"{'Frozen':<12} {'0.3685':>10} {'2060':>10} {frozen['P_peak']:>14.3f} {frozen['t_peak']:>12.3f}",
        f"{'Effective':<12} {'0.6517':>10} {'2764':>10} {effective['P_peak']:>14.3f} {effective['t_peak']:>12.3f}",
        "",
        f"Experimental P_peak (from CSV): {max(exp['pressure']):.3f} MPa",
        f"Experimental t_peak (CSV, no offset): {exp['time'][int(np.argmax(exp['pressure']))]:.3f} s",
    ]
    metrics_path = out / 'metrics.txt'
    metrics_path.write_text('\n'.join(lines))
    print(f"  Saved: {metrics_path}")
    print()
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
