import copy
from types import SimpleNamespace

import numpy as np
import pytest

from srm_1d import run_simulation
from srm_1d.igniter_plenum import PyrogenChamber
from srm_1d.nozzle import Nozzle
from srm_1d.propellant import Pyrogen
from srm_1d.tests._motor_fixtures import (
    example_bates_geo,
    hasegawa_propellant_1,
    single_cylinder_geo,
)
from srm_1d.tools.ignition_diagnostics import (
    analyze_ignition_spike,
    ignition_spread_metrics,
    source_timeseries,
)


def _snap(t, burning, r_total=None, r_erosive=None, endface=None):
    burning = np.asarray(burning, dtype=bool)
    n = len(burning)
    if r_total is None:
        r_total = np.zeros(n)
    if r_erosive is None:
        r_erosive = np.zeros(n)
    if endface is None:
        endface = np.zeros(n)
    return {
        "t": float(t),
        "x": np.arange(n, dtype=float),
        "P": np.full(n, 101325.0),
        "T": np.full(n, 300.0),
        "T_surf": np.full(n, 293.0),
        "is_burning": burning,
        "is_grain": np.ones(n, dtype=bool),
        "C_burn": np.ones(n),
        "r_total": np.asarray(r_total, dtype=float),
        "r_erosive": np.asarray(r_erosive, dtype=float),
        "endface_msource": np.asarray(endface, dtype=float),
    }


def _synthetic_result(time, pressure, mdot, snapshots=()):
    time = np.asarray(time, dtype=float)
    return {
        "time": time,
        "P_head": np.asarray(pressure, dtype=float),
        "P_exit": np.asarray(pressure, dtype=float),
        "mdot_ig": np.asarray(mdot, dtype=float),
        "P_ig": np.asarray(mdot, dtype=float) * 1.0e6 + 101325.0,
        "m_pyrogen": np.maximum(0.0, 0.001 - np.cumsum(mdot) * 0.01),
        "snapshots": list(snapshots),
        "P_ambient": 101325.0,
        "summary": {},
    }


def _test_chamber():
    pyro = Pyrogen(
        name="diagnostic-test",
        a=3.0e-5,
        n=0.5,
        rho=1700.0,
        T_flame=2800.0,
        M=0.030,
        gamma=1.25,
        impetus_W=5000.0,
    )
    return PyrogenChamber(
        pyrogen=pyro,
        m_pyrogen_initial=0.003,
        A_burn_initial=5.0e-4,
        A_throat=2.0e-5,
        V_plenum=3.0e-6,
        burn_law="end_burning",
    )


def _small_motor(inhibited=True):
    if inhibited:
        geo = single_cylinder_geo(
            D_bore=0.030, D_outer=0.060, length=0.120,
            target_propellant_cells=12,
        )
    else:
        geo = example_bates_geo(target_propellant_cells=32)
    prop = hasegawa_propellant_1()
    nozzle = Nozzle(D_throat=0.010, D_exit=0.020, efficiency=0.95)
    return geo, prop, nozzle


def _run_small(**kwargs):
    geo, prop, nozzle = _small_motor(kwargs.pop("inhibited", True))
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=kwargs.pop("T_ignition", 294.0),
        t_max=kwargs.pop("t_max", 0.006),
        P_cutoff=kwargs.pop("P_cutoff", 1.0),
        dt_max=kwargs.pop("dt_max", 2.0e-5),
        burn_update_interval=kwargs.pop("burn_update_interval", 1),
        snapshot_interval=kwargs.pop("snapshot_interval", 0.001),
        cfl_target=kwargs.pop("cfl_target", 0.5),
        verbose=False,
        **kwargs,
    )
    return result, geo, prop


def test_classifies_pyrogen_driven_spike_when_peak_overlaps_mdot():
    result = _synthetic_result(
        [0.0, 0.1, 0.2],
        [101325.0, 3.0e6, 1.5e6],
        [0.0, 0.02, 0.01],
    )
    diagnostics = analyze_ignition_spike(result)
    assert diagnostics["classification"]["primary_driver"] == "pyrogen_combustion"
    assert diagnostics["classification"]["pyrogen_combustion"]


def test_classifies_erosive_snap_on_after_pyrogen_burnout():
    result = _synthetic_result(
        [0.0, 0.1, 0.2, 0.3],
        [101325.0, 0.4e6, 1.0e6, 3.0e6],
        [0.02, 0.0, 0.0, 0.0],
        snapshots=[
            _snap(0.1, [True, True], r_total=[1.0, 1.0], r_erosive=[0.0, 0.0]),
            _snap(0.3, [True, True], r_total=[10.0, 10.0], r_erosive=[9.0, 9.0]),
        ],
    )
    diagnostics = analyze_ignition_spike(
        result, propellant=SimpleNamespace(rho_propellant=1.0)
    )
    assert diagnostics["classification"]["primary_driver"] == "erosive_snap_on"
    assert diagnostics["classification"]["erosive_snap_on"]


def test_detects_instant_ignition_collapse():
    result = _synthetic_result(
        [0.0, 0.1, 0.2],
        [101325.0, 1.0e6, 0.8e6],
        [0.0, 0.0, 0.0],
        snapshots=[
            _snap(0.0, [False, False, False, False]),
            _snap(0.1, [True, True, True, True]),
            _snap(0.2, [True, True, True, True]),
        ],
    )
    spread = ignition_spread_metrics(result)
    assert spread["instant_ignition_collapse"]
    assert spread["spread_10_90_s"] == pytest.approx(0.0)


def test_source_timeseries_estimates_sources_without_mutating_result():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
        snapshots=[
            _snap(0.1, [True, True], r_total=[3.0, 3.0],
                  r_erosive=[1.0, 1.0], endface=[0.5, 0.5]),
        ],
    )
    original = copy.deepcopy(result)
    sources = source_timeseries(
        result, propellant=SimpleNamespace(rho_propellant=2.0)
    )
    assert sources["normal_sidewall_kg_s"][0] == pytest.approx(8.0)
    assert sources["erosive_sidewall_kg_s"][0] == pytest.approx(4.0)
    assert sources["endface_kg_s"][0] == pytest.approx(1.0)
    assert sources["pyrogen_kg_s"][0] == pytest.approx(1.0)
    np.testing.assert_array_equal(result["P_head"], original["P_head"])
    np.testing.assert_array_equal(
        result["snapshots"][0]["r_total"],
        original["snapshots"][0]["r_total"],
    )


def test_default_diagnostic_controls_preserve_behavior():
    base, _geo, _prop = _run_small()
    explicit, _geo, _prop = _run_small(
        initial_gas_temperature=None,
        diagnostic_disable_erosive=False,
        diagnostic_disable_endfaces=False,
    )
    np.testing.assert_allclose(base["P_head"], explicit["P_head"])
    np.testing.assert_allclose(base["mdot_ig"], explicit["mdot_ig"])


def test_initial_gas_temperature_can_delay_goodman_heating():
    hot, _geo, _prop = _run_small(T_ignition=500.0)
    ambient, _geo, prop = _run_small(
        T_ignition=500.0,
        initial_gas_temperature=propellant_initial_temperature(),
    )
    hot_spread = ignition_spread_metrics(hot)
    ambient_spread = ignition_spread_metrics(ambient)
    assert ambient["summary"]["initial_gas_temperature"] == pytest.approx(293.0)
    assert (
        np.isnan(ambient_spread["t50_s"])
        or ambient_spread["t50_s"] >= hot_spread["t50_s"]
    )


def propellant_initial_temperature():
    return hasegawa_propellant_1().T_initial


def test_disabling_erosive_preserves_normal_burn_and_zeroes_erosive_increment():
    baseline, _geo, _prop = _run_small(roughness=100e-6, T_ignition=294.0)
    disabled, _geo, _prop = _run_small(
        roughness=100e-6,
        T_ignition=294.0,
        diagnostic_disable_erosive=True,
    )
    baseline_r_erosive = max(np.max(s["r_erosive"]) for s in baseline["snapshots"])
    disabled_r_erosive = max(np.max(s["r_erosive"]) for s in disabled["snapshots"])
    disabled_r_total = max(np.max(s["r_total"]) for s in disabled["snapshots"])
    assert baseline_r_erosive > 0.0
    assert disabled_r_erosive == pytest.approx(0.0)
    assert disabled_r_total > 0.0


def test_disabling_endfaces_removes_endface_source_without_changing_sidewall_geometry():
    baseline, _geo, _prop = _run_small(inhibited=False, T_ignition=294.0)
    disabled, _geo, _prop = _run_small(
        inhibited=False,
        T_ignition=294.0,
        diagnostic_disable_endfaces=True,
    )
    assert max(np.max(s["endface_msource"]) for s in baseline["snapshots"]) > 0.0
    assert max(np.max(s["endface_msource"]) for s in disabled["snapshots"]) == pytest.approx(0.0)
    np.testing.assert_allclose(
        baseline["snapshots"][0]["C_burn"],
        disabled["snapshots"][0]["C_burn"],
        atol=1.0e-6,
    )
