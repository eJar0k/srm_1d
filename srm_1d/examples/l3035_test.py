"""
l3035_test.py -- run the L3035 .ric through the v0.7.0 workflow.

This mirrors `BALLSstick.py` and is intentionally simple: load the .ric,
auto-discover the sibling transport YAML, use the built-in BPNV pyrogen,
and write plots under `artifacts/l3035/`.

Usage:
    python -m srm_1d.examples.l3035_test
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import plot_pressure, plot_flow_snapshot, plot_summary


CASE_NAME = 'l3035'
MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'L3035.ric'
OUTPUT_DIR = Path('artifacts') / CASE_NAME


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        roughness=50e-6,
        kappa=0.45,
        pyrogen='bpnv',
        pyrogen_mass=None,
        T_ignition=850.0,
        P_cutoff=0.01e6,
        snapshot_interval=0.2,
        print_interval=0.2,
        verbose=False,
    )

    summary = result['summary']
    print(
        f"L3035: P_peak={summary['P_peak'] / 1e6:.2f} MPa, "
        f"t_burn={summary['t_burn']:.3f} s, "
        f"impulse={perf['total_impulse']:.1f} N*s, "
        f"designation={perf['motor_designation']}"
    )

    plot_pressure(
        result,
        title="L3035 Pressure Trace",
        save_path=OUTPUT_DIR / "l3035_pressure.png",
    )

    plot_flow_snapshot(
        result,
        t_target=0.2,
        save_path=OUTPUT_DIR / "l3035_flow.png",
    )

    plot_summary(
        result,
        performance=perf,
        title="L3035 Simulation Summary",
        save_path=OUTPUT_DIR / "l3035_summary.png",
    )

    plt.close('all')
    print(f"Plots saved under {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
