"""
ISP_Super_Loki.py -- run the 'ISP_Super_Loki.ric' through the v0.7.3
Phase A head_basket topology.

Super Loki's physical igniter is a head-end BKNO3 pellet charge in a
consumable moisture cup — NO defined orifice or pressure-containing
aft cap. Modeling it as ``forward_plenum`` is wrong physics; the
v0.7.3 Phase A ``head_basket`` topology (uncontained pyrogen, each
pellet burning at local bore pressure) is the appropriate fit. See
the ``PyrogenChamber`` docstring at
``srm_1d/igniter_plenum.py`` L52-L120 for the architectural split.

Usage:
    python -m srm_1d.examples.ISP_Super_Loki
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
    plot_flow_snapshots, plot_field_heatmap,
    ISP_SUPER_LOKI_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir


CASE_NAME = 'ISP_Super_Loki'
MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'ISP_Super_Loki.ric'


def main():
    OUTPUT_DIR = artifact_dir(CASE_NAME)

    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        roughness=30e-6,
        kappa=0.45,
        pyrogen='mtv',
        pyrogen_mass=None,
        T_ignition=850.0,
        P_cutoff=0.01e6,
        # v0.7.3 Phase A: physical-topology fit. head_basket distributes
        # pyrogen pellets uniformly across the head-end cartridge cells;
        # each pellet burns at its host cell's local bore pressure.
        injection_topology='head_basket',
        cartridge_length_m=-1.0,  # derive from pyrogen mass + bore geometry
        snapshot_interval=0.01,   # dense snapshots for ignition diagnostics
        print_interval=0.2,
        verbose=True,
    )

    summary = result['summary']
    print(
        f"ISP_Super_Loki [head_basket]: "
        f"P_peak={summary['P_peak'] / 1e6:.2f} MPa, "
        f"t_peak={summary['t_peak']:.3f} s, "
        f"t_burn={summary['t_burn']:.3f} s, "
        f"impulse={perf['total_impulse']:.1f} N*s, "
        f"designation={perf['motor_designation']}"
    )

    plot_pressure(
        result,
        title="ISP_Super_Loki Pressure Trace (head_basket)",
        experimental=ISP_SUPER_LOKI_EXPERIMENTAL,
        save_path=str(OUTPUT_DIR / "pressure.png"),
    )

    plot_flow_snapshot(
        result,
        t_target=0.2,
        save_path=str(OUTPUT_DIR / "flow_t0p2.png"),
    )

    plot_flow_snapshots(
        result,
        t_targets=[0.01, 0.05, 0.20, 0.50, 1.00],
        fields=('P', 'u_cell', 'T', 'is_burning'),
        title="ISP_Super_Loki — flow evolution",
        save_path=str(OUTPUT_DIR / "flow_multi.png"),
    )

    plot_field_heatmap(
        result,
        fields=('P', 'u_cell', 'T', 'is_burning'),
        title="ISP_Super_Loki — x-t heatmap (full burn)",
        save_path=str(OUTPUT_DIR / "heatmap_full.png"),
    )

    plot_field_heatmap(
        result,
        fields=('P', 'u_cell', 'T', 'is_burning'),
        t_max=0.5,
        title="ISP_Super_Loki — x-t heatmap (ignition transient, t<=0.5s)",
        save_path=str(OUTPUT_DIR / "heatmap_ignition.png"),
    )

    plot_summary(
        result,
        performance=perf,
        title="ISP_Super_Loki Simulation Summary (head_basket)",
        experimental=ISP_SUPER_LOKI_EXPERIMENTAL,
        save_path=str(OUTPUT_DIR / "summary.png"),
    )

    plt.close('all')
    print(f"Plots saved under {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
