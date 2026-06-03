"""
hasegawa_motor_a.py — Validate against Hasegawa et al. (2006) Motor A
=======================================================================

Loads the Hasegawa Motor A definition from `srm_1d/motors/hasegawa_a.ric`
(geometry + propellant + nozzle) plus its sibling
`hasegawa_a.transport.yaml` (RPA effective gas transport), runs the
1D PISO simulation, and compares the head-end pressure trace against
digitized experimental data from Ma et al. (2020) Fig. 10.

The roughness / igniter parameters below are the Rank-1 fit from the
v0.6.0 Latin Hypercube sweep (see srm_1d/tools/sensitivity.py and
the hasegawa_a_lhs example). They are the canonical calibration
target — see DEVNOTES "Calibration State" and the
`project_hasegawa_calibration_state` memory.

Usage:
    python -m examples.hasegawa_motor_a

Outputs:
    artifacts/hasegawa_a/<timestamp>_<sha>[-dirty]/
        pressure.png  — pressure trace vs experimental
        flow.png      — flow field snapshot at t ≈ 2s
        summary.png   — 4-panel summary

Each run lands in its own stamped subdirectory so re-runs don't
overwrite earlier traces (important during Phase 5 LHS where
calibration runs accumulate quickly).
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
    HASEGAWA_MOTOR_A_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir, verify_run_health


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'hasegawa_a.ric'
EXPERIMENTAL_TIME_OFFSET = 1.1  # align experimental ignition with sim t=0


def main():
    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        roughness=37.1e-6,
        kappa=0.45,
        pyrogen='bpnv',
        pyrogen_mass=None,
        T_ignition=850.0,
        P_cutoff=0.05e6,
        snapshot_interval=0.2,
        print_interval=0.2,
    )

    out = artifact_dir('hasegawa_a')

    plot_pressure(
        result,
        title="Motor A — 1D PISO vs Experimental (Hasegawa 2006)",
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=EXPERIMENTAL_TIME_OFFSET,
        save_path=str(out / 'pressure.png'),
    )

    plot_flow_snapshot(
        result, t_target=2.0,
        save_path=str(out / 'flow.png'),
    )

    plot_summary(
        result, performance=perf,
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=EXPERIMENTAL_TIME_OFFSET,
        title="Hasegawa Motor A — Simulation Summary",
        save_path=str(out / 'summary.png'),
    )

    plt.close('all')
    print(f"\nAll plots saved to {out}")
    # v0.7.3.2 run-health gate: surfaces collapsed runs that would
    # otherwise save misleading plots and trick the user into reading
    # a stale-looking artifact dir as a successful calibration. Don't
    # raise — the artifacts are still useful for debugging — but
    # banner clearly when something is wrong.
    verify_run_health(result, motor_name='Hasegawa A (canonical)')


if __name__ == '__main__':
    main()
