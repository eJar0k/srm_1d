"""
bates_4seg.py — 4-Segment BATES Motor Simulation
===================================================

Runs the example 4-segment BATES motor (`srm_1d/motors/example_bates.ric`)
and demonstrates gap-flow physics.

Usage:
    python -m srm_1d.examples.bates_4seg

Outputs:
    bates_pressure.png  — pressure trace
    bates_flow.png      — flow field snapshot
    bates_thrust.png    — thrust and Isp
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import plot_pressure, plot_thrust, plot_flow_snapshot


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'example_bates.ric'


def main():
    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        roughness=20e-6,
        kappa=0.45,
        pyrogen='bpnv',
        pyrogen_mass=None,
        T_ignition=850.0,
        P_cutoff=0.05e6,
        snapshot_interval=0.1,
        print_interval=0.1,
    )

    plot_pressure(result, title="4-Segment BATES",
                  save_path="bates_pressure.png")

    plot_flow_snapshot(result, t_target=0.3,
                       title="BATES — Flow at t ≈ 0.3s",
                       save_path="bates_flow.png")

    plot_thrust(result, perf, title="4-Segment BATES",
                save_path="bates_thrust.png")

    plt.close('all')
    print("\nAll plots saved.")


if __name__ == '__main__':
    main()
