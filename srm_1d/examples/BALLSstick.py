"""
BALLSstick.py -- run the BALLSstick .ric through the v0.7.0 workflow.

This is the quick pattern for any openMotor .ric file:

1. Put `<motor>.ric` in `srm_1d/motors/`.
2. Add sibling `<motor>.transport.yaml` with `mu`, `k`, and `Cp`.
3. Call `run_from_ric(..., pyrogen='bpnv')`.
4. Save pressure, flow, and summary plots under `artifacts/<motor>/`.

Usage:
    python -m srm_1d.examples.BALLSstick

This is an exploratory bounded run, not a calibrated prediction. Increase
or remove `t_max` after checking the pressure trace and geometry.
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import plot_pressure, plot_flow_snapshot, plot_summary


CASE_NAME = 'BALLSstick'
MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / f'{CASE_NAME}.ric'
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
        t_max=3.0,
        snapshot_interval=0.5,
        print_interval=1.0,
        verbose=False,
    )

    summary = result['summary']
    print(
        f"{CASE_NAME}: P_peak={summary['P_peak'] / 1e6:.2f} MPa, "
        f"t_burn={summary['t_burn']:.3f} s, "
        f"impulse={perf['total_impulse']:.1f} N*s"
    )

    plot_pressure(
        result,
        title=f"{CASE_NAME} Pressure Trace",
        save_path=OUTPUT_DIR / f"{CASE_NAME}_pressure.png",
    )

    plot_flow_snapshot(
        result,
        t_target=0.5,
        save_path=OUTPUT_DIR / f"{CASE_NAME}_flow.png",
    )

    plot_summary(
        result,
        performance=perf,
        title=f"{CASE_NAME} Simulation Summary",
        save_path=OUTPUT_DIR / f"{CASE_NAME}_summary.png",
    )

    plt.close('all')
    print(f"Plots saved under {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
