"""
test_plotting_channels.py — v0.8.0 Phase 2 channel-native plotting tests.

Verifies the generic ``plot_channels`` engine and that the bespoke figures
(pressure / thrust / summary / flow-snapshot / heatmap) are reproduced
through the channel model — accepting either a results dict or a
``SimulationChannels``, producing identical plotted data, with units
carried on the channels and converted at draw time.

Uses the non-interactive Agg backend; one real short simulation built once
for the module.
"""

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from srm_1d.channels import build_channels  # noqa: E402
from srm_1d import plotting  # noqa: E402


@pytest.fixture(scope="module")
def real_run():
    pytest.importorskip("numba")
    pytest.importorskip("skfmm")
    from srm_1d import run_simulation
    from srm_1d.nozzle import compute_motor_performance
    from srm_1d.tests.test_simulation_phase3 import _small_motor, _test_chamber
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0, t_max=0.010, P_cutoff=1.0,
        dt_max=2.0e-5, burn_update_interval=1, snapshot_interval=0.002,
        cfl_target=0.5, verbose=False,
    )
    perf = compute_motor_performance(result, nozzle, prop)
    return geo, prop, nozzle, result, perf


def teardown_function(_):
    plt.close('all')


# ---- plot_channels generic engine -----------------------------------

def test_plot_channels_unit_conversion(real_run):
    *_, result, _ = real_run
    sc = build_channels(result)
    fig, axes = plotting.plot_channels(
        sc, 'P_head', display_units={'P_head': 'MPa'}, overlay=True)
    line = axes[0].lines[0]
    # Plotted y-data is Pa -> MPa (carried on the channel), x is time [s].
    np.testing.assert_allclose(line.get_ydata(), result['P_head'] / 1e6,
                               rtol=1e-9)
    np.testing.assert_array_equal(line.get_xdata(), result['time'])
    plt.close(fig)


def test_plot_channels_stacked_multi(real_run):
    *_, result, _ = real_run
    fig, axes = plotting.plot_channels(result, ['P_head', 'Kn', 'massflow'])
    assert len(axes) == 3
    np.testing.assert_array_equal(axes[0].lines[0].get_ydata(),
                                  result['P_head'])  # default unit = Pa
    plt.close(fig)


def test_plot_channels_accepts_dict_and_channels(real_run):
    *_, result, _ = real_run
    f1, a1 = plotting.plot_channels(result, 'Kn', overlay=True)
    f2, a2 = plotting.plot_channels(build_channels(result), 'Kn', overlay=True)
    np.testing.assert_array_equal(a1[0].lines[0].get_ydata(),
                                  a2[0].lines[0].get_ydata())
    plt.close(f1)
    plt.close(f2)


# ---- bespoke figures, dict vs channels parity -----------------------

def test_plot_pressure_dict_vs_channels(real_run):
    *_, result, _ = real_run
    f1, a1 = plotting.plot_pressure(result)
    f2, a2 = plotting.plot_pressure(build_channels(result))
    y1 = a1.lines[0].get_ydata()
    y2 = a2.lines[0].get_ydata()
    np.testing.assert_array_equal(y1, y2)
    # Primary trace is MPa via the generic path.
    np.testing.assert_allclose(y1, result['P_head'] / 1e6, rtol=1e-9)
    plt.close(f1)
    plt.close(f2)


def test_plot_thrust_via_generic_path(real_run):
    *_, result, perf = real_run
    fig, axes = plotting.plot_thrust(result, perf)
    # 2-panel: thrust then Isp, both sourced from the perf channels.
    np.testing.assert_array_equal(axes[0].lines[0].get_ydata(), perf['thrust'])
    np.testing.assert_array_equal(axes[1].lines[0].get_ydata(), perf['Isp'])
    plt.close(fig)


def test_plot_summary_runs_dict_and_channels(real_run):
    *_, result, perf = real_run
    f1, _ = plotting.plot_summary(result, performance=perf)
    f2, _ = plotting.plot_summary(build_channels(result), performance=perf)
    assert f1 is not None and f2 is not None
    plt.close(f1)
    plt.close(f2)


def test_plot_flow_snapshot_channel_sourced(real_run):
    *_, result, _ = real_run
    f1, _ = plotting.plot_flow_snapshot(result, snap_index=1)
    f2, _ = plotting.plot_flow_snapshot(build_channels(result), snap_index=1)
    assert f1 is not None and f2 is not None
    plt.close(f1)
    plt.close(f2)


def test_plot_field_heatmap_channel_sourced(real_run):
    *_, result, _ = real_run
    fig, axes = plotting.plot_field_heatmap(result, fields=('P', 'T'))
    assert fig is not None
    plt.close(fig)


def test_plot_flow_snapshots_channel_sourced(real_run):
    *_, result, _ = real_run
    snaps = result['snapshots']
    t_targets = [snaps[0]['t'], snaps[-1]['t']]
    fig, axes = plotting.plot_flow_snapshots(result, t_targets,
                                             fields=('P', 'T'))
    assert axes.shape == (2, 2)
    plt.close(fig)
