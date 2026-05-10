"""
Zerox_test.py — Single-shot Zerox static-fire trace
====================================================

Quick sim of the Zerox motor (forward Finocyl + aft BATES, ~1.45 kg
"Risky Batman V3") with the current borrowed Hasegawa A calibration.
Used as a baseline for the Zerox tuning workflow — see
zerox_diagnostic.py for the spike/throat/FMM diagnostic plots.

The igniter parameters here are inherited from Hasegawa A's v0.6.0
LHS rank-1 fit and are known-wrong for Zerox. They will be replaced
by an independent LHS fit at the end of the diagnostic phase.

Usage:
    python -m srm_1d.examples.Zerox_test
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric, save_csv
from srm_1d.plotting import (
    plot_pressure, plot_thrust, plot_flow_snapshot,
    ZEROX_EXPERIMENTAL,
)


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox.ric'


def main():
    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        # Borrowed Hasegawa A calibration — known-wrong for Zerox.
        roughness=20e-6,
        kappa=0.45,
        pyrogen='bpnv',
        pyrogen_mass=None,
        T_ignition=850.0,
        cfl_target=0.5,
        dt_max=1e-4,
        t_max=8.0,
        P_cutoff=0.01e6,
        snapshot_interval=0.5,
        print_interval=0.2,
    )

    plot_pressure(
        result, title="Zerox Test",
        experimental=ZEROX_EXPERIMENTAL,
        time_offset=ZEROX_EXPERIMENTAL.get('time_offset', 0.0),
        save_path="Zerox_pressure.png",
    )
    plot_flow_snapshot(
        result, t_target=0.3,
        title="Zerox — Flow at t ≈ 0.3s",
        save_path="Zerox_flow.png",
    )
    plot_thrust(
        result, perf, title="Zerox Test",
        save_path="Zerox_thrust.png",
    )
    plt.close('all')
    print("\nAll plots saved.")

    save_csv("output.csv", result, perf, geo=geo, propellant=prop)
    print("CSV saved.")


if __name__ == '__main__':
    main()
