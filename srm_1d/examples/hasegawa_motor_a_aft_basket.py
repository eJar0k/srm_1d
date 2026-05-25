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

from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
    plot_flow_snapshots, plot_field_heatmap,
    HASEGAWA_MOTOR_A_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'hasegawa_a.ric'
EXPERIMENTAL_TIME_OFFSET = 1.1  # align experimental ignition with sim t=0


def _run_one(mode, out):
    """Run Hasegawa A aft_basket with one heat-delivery mode."""
    pyrogen_obj = load_pyrogen('bpnv')
    pyrogen_obj.heat_delivery_mode = mode

    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        roughness=37.1e-6,
        kappa=0.45,
        pyrogen=pyrogen_obj,
        pyrogen_mass=None,
        T_ignition=850.0,
        P_cutoff=0.05e6,
        injection_topology='aft_basket',
        cartridge_length_m=-1.0,
        snapshot_interval=0.005,
        print_interval=0.2,
        verbose=False,
    )

    summary = result['summary']
    print(
        f"Hasegawa A [aft_basket / {mode}]: "
        f"P_peak={summary['P_peak']/1e6:.2f} MPa @ t={summary['t_peak']:.3f}s, "
        f"impulse={perf['total_impulse']:.1f} N*s"
    )

    plot_pressure(
        result,
        title=f"Hasegawa A — aft_basket / {mode}",
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=EXPERIMENTAL_TIME_OFFSET,
        save_path=str(out / f'pressure_{mode}.png'),
    )
    plot_flow_snapshots(
        result,
        t_targets=[0.005, 0.020, 0.050, 0.100, 0.500],
        fields=('P', 'u', 'T', 'is_burning'),
        title=f"Hasegawa A aft_basket — back->front cascade ({mode})",
        save_path=str(out / f'flow_multi_{mode}.png'),
    )
    plot_field_heatmap(
        result,
        fields=('P', 'u', 'T', 'T_surf', 'is_burning'),
        t_max=0.5,
        title=f"Hasegawa A aft_basket — x-t heatmap (ignition, {mode})",
        save_path=str(out / f'heatmap_ignition_{mode}.png'),
    )
    plot_summary(
        result, performance=perf,
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=EXPERIMENTAL_TIME_OFFSET,
        title=f"Hasegawa A aft_basket Summary ({mode})",
        save_path=str(out / f'summary_{mode}.png'),
    )
    plt.close('all')
    return summary


def main():
    out = artifact_dir('hasegawa_a_aft_basket')
    print("v0.7.3 Phase B.6 — Hasegawa A aft_basket diagnostic A/B/control")
    print("  Diagnostic Q: does the simultaneous-ignition artifact "
          "persist under reversed (aft) mass-injection topology?")
    print("  Forward_plenum baseline: P_peak ~ 6.20 MPa @ t ~ 0.03 s "
          "(v0.7.0 calibrated; pre-Phase-B.0)")
    print()

    s_none = _run_one('none', out)
    s_demar = _run_one('demar', out)
    s_radiation = _run_one('radiation', out)

    print()
    print(f"Plots saved under {out}")
    print()
    print("A/B summary:")
    print(f"  none:      P_peak = {s_none['P_peak']/1e6:7.3f} MPa "
          f"@ t={s_none['t_peak']:.3f} s")
    print(f"  demar:     P_peak = {s_demar['P_peak']/1e6:7.3f} MPa "
          f"@ t={s_demar['t_peak']:.3f} s")
    print(f"  radiation: P_peak = {s_radiation['P_peak']/1e6:7.3f} MPa "
          f"@ t={s_radiation['t_peak']:.3f} s")
    print("  baseline:  P_peak ~ 6.20 MPa (forward_plenum, v0.7.0)")


if __name__ == '__main__':
    main()
