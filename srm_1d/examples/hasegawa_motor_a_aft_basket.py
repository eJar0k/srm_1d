"""
hasegawa_motor_a_aft_basket.py — v0.7.3 Phase A diagnostic variant

Reruns Hasegawa Motor A with ``injection_topology='aft_basket'``: the
pyrogen pellets are placed in the aft end of the bore (cells
[N - n_cart, N - 1]) instead of upstream. Each pellet burns at its
host cell's local bore pressure; no momentum injection (uncontained
model — see PyrogenChamber docstring at
srm_1d/igniter_plenum.py L52-L120).

**Diagnostic question** per
``srm_1d/docs/v0_7_2/candidates_post_phaseA.md`` §3: if the
simultaneous-ignition pressure-spike artifact persists under reversed
mass-injection topology (mass entering near the nozzle instead of
the head-end), the artifact lives in the per-cell Goodman ignition
kernel, not the pyrogen source model — so the next v0.7.3+ candidate
should target per-cell coupling (candidate 4 / Z-N burn-rate lag /
solid-phase axial conduction). If the spike disappears or shifts
qualitatively, the artifact was driven by head-end mass concentration
and the unified-igniter refactor (candidate 6) gets prioritized.

Same knobs as ``hasegawa_motor_a.py`` (roughness=37.1um, kappa=0.45,
T_ignition=850, BPNV pyrogen). Only the topology differs.

Usage:
    python -m srm_1d.examples.hasegawa_motor_a_aft_basket
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
    plot_flow_snapshots, plot_field_heatmap,
    HASEGAWA_MOTOR_A_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir


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
        # v0.7.3 Phase A diagnostic — reverse the mass-injection topology
        injection_topology='aft_basket',
        cartridge_length_m=-1.0,   # derive from pyrogen mass + bore geometry
        snapshot_interval=0.005,   # dense coverage of the back->front cascade
        print_interval=0.2,
    )

    out = artifact_dir('hasegawa_a_aft_basket')

    plot_pressure(
        result,
        title="Hasegawa Motor A — aft_basket diagnostic (v0.7.3 Phase A)",
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=EXPERIMENTAL_TIME_OFFSET,
        save_path=str(out / 'pressure.png'),
    )

    plot_flow_snapshot(
        result, t_target=0.05,
        save_path=str(out / 'flow_t0p05.png'),
    )

    plot_flow_snapshot(
        result, t_target=2.0,
        save_path=str(out / 'flow_t2.png'),
    )

    plot_flow_snapshots(
        result,
        t_targets=[0.005, 0.020, 0.050, 0.100, 0.500],
        fields=('P', 'u_cell', 'T', 'is_burning'),
        title="Hasegawa A aft_basket — back->front cascade",
        save_path=str(out / 'flow_multi.png'),
    )

    plot_field_heatmap(
        result,
        fields=('P', 'u_cell', 'T', 'is_burning'),
        title="Hasegawa A aft_basket — x-t heatmap (full burn)",
        save_path=str(out / 'heatmap_full.png'),
    )

    plot_field_heatmap(
        result,
        fields=('P', 'u_cell', 'T', 'is_burning'),
        t_max=0.5,
        title="Hasegawa A aft_basket — x-t heatmap (ignition, t<=0.5s)",
        save_path=str(out / 'heatmap_ignition.png'),
    )

    plot_summary(
        result, performance=perf,
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=EXPERIMENTAL_TIME_OFFSET,
        title="Hasegawa A aft_basket — Simulation Summary",
        save_path=str(out / 'summary.png'),
    )

    summary = result['summary']
    print(
        f"hasegawa_motor_a_aft_basket: "
        f"P_peak={summary['P_peak'] / 1e6:.2f} MPa @ t={summary['t_peak']:.3f} s, "
        f"t_burn={summary['t_burn']:.3f} s, "
        f"impulse={perf['total_impulse']:.1f} N*s"
    )
    print(f"Plots saved under {out}")
    plt.close('all')


if __name__ == '__main__':
    main()
