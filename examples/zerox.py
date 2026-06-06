"""
zerox.py — Canonical Zerox static-fire reproduction
=====================================================

Loads the LHS-calibrated Zerox motor from
`srm_1d/motors/zerox_LHS.ric` (calibrated `erosionCoeff` and
propellant `a` baked in) plus its sibling `.transport.yaml` (RPA
frozen gas transport), runs the 1D PISO simulation with the LHS
rank-1 igniter parameters, and compares the head-end pressure trace
against digitized experimental data from the user's static fire.

The igniter / kappa parameters below are the rank-1 fit from the
v0.6.0 7-variable LHS sweep (see srm_1d/examples/zerox_lhs.py and
DEVNOTES "Calibration State"). MSE = 0.071 MPa² vs experimental
(RMS error ≈ 0.27 MPa, ~7% of peak).

Usage:
    python -m examples.zerox

Outputs (under artifacts/zerox/<timestamp>_<sha>[-dirty]/):
    pressure.png  — pressure trace vs experimental
    flow.png      — flow field snapshot at t ≈ 0.3s
    summary.png   — 4-panel summary
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
    ZEROX_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox_LHS.ric'


def main():
    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        # v0.7.5 cross-motor re-LHS optimum (shared fired-motor fit;
        # docs/v0_7_5/RESULT.md). Supersedes the v0.6.0 per-motor fit
        # (roughness 20um / kappa 0.3286 -- the latter below the 0.40
        # physical floor).
        roughness=32e-6,
        kappa=0.44,
        pyrogen='bpnv',
        pyrogen_mass=None,
        T_ignition=756.0,
        P_cutoff=0.01e6,
        snapshot_interval=0.2,
        print_interval=0.2,
    )

    # Per-run artifact dir matches the run_template.py convention.
    out = artifact_dir('zerox')

    plot_pressure(
        result,
        title="Zerox — 1D PISO vs Experimental (LHS-calibrated)",
        experimental=ZEROX_EXPERIMENTAL,
        time_offset=ZEROX_EXPERIMENTAL.get('time_offset', 0.0),
        save_path=str(out / "pressure.png"),
    )

    plot_flow_snapshot(
        result, t_target=0.3,
        save_path=str(out / "flow.png"),
    )

    plot_summary(
        result, performance=perf,
        experimental=ZEROX_EXPERIMENTAL,
        time_offset=ZEROX_EXPERIMENTAL.get('time_offset', 0.0),
        title="Zerox — Simulation Summary",
        save_path=str(out / "summary.png"),
    )

    plt.close('all')
    print(f"\nAll plots saved to {out}")


if __name__ == '__main__':
    main()
