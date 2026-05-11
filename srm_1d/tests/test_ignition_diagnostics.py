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
    energy_momentum_timeseries,
    ignition_spread_metrics,
    pressure_landmarks,
    source_timeseries,
)


def _snap(t, burning, r_total=None, r_erosive=None, endface=None,
          pyrogen_surface_heat_flux=None, radiation_heat_flux=None):
    burning = np.asarray(burning, dtype=bool)
    n = len(burning)
    if r_total is None:
        r_total = np.zeros(n)
    if r_erosive is None:
        r_erosive = np.zeros(n)
    if endface is None:
        endface = np.zeros(n)
    if pyrogen_surface_heat_flux is None:
        pyrogen_surface_heat_flux = np.zeros(n)
    if radiation_heat_flux is None:
        radiation_heat_flux = np.zeros(n)
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
        "pyrogen_surface_heat_flux": np.asarray(pyrogen_surface_heat_flux, dtype=float),
        "radiation_heat_flux": np.asarray(radiation_heat_flux, dtype=float),
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


def _test_chamber(**pyrogen_overrides):
    pyro = Pyrogen(
        name="diagnostic-test",
        a=3.0e-5,
        n=0.5,
        rho=1700.0,
        T_flame=2800.0,
        M=0.030,
        gamma=1.25,
        impetus_W=5000.0,
        heat_flux_cal_cm2_s=69.4,
    )
    for key, value in pyrogen_overrides.items():
        setattr(pyro, key, value)
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
    assert "energy" in diagnostics


def test_classifies_erosive_snap_on_after_pyrogen_burnout():
    result = _synthetic_result(
        [0.0, 0.1, 0.2],
        [101325.0, 3.0e6, 1.0e6],
        [0.02, 0.0, 0.0],
        snapshots=[
            _snap(0.1, [True, True], r_total=[10.0, 10.0], r_erosive=[9.0, 9.0]),
            _snap(0.2, [True, True], r_total=[1.0, 1.0], r_erosive=[0.0, 0.0]),
        ],
    )
    diagnostics = analyze_ignition_spike(
        result, propellant=SimpleNamespace(rho_propellant=1.0)
    )
    assert diagnostics["classification"]["primary_driver"] == "erosive_snap_on"
    assert diagnostics["classification"]["erosive_snap_on"]


def test_pressure_landmarks_label_startup_peak_separately_from_global_peak():
    result = _synthetic_result(
        [0.0, 0.05, 0.2, 2.0],
        [101325.0, 1.0e6, 0.8e6, 5.0e6],
        [0.02, 0.0, 0.0, 0.0],
    )
    pressure = pressure_landmarks(result, startup_margin_s=0.10)
    assert pressure["startup_window_peak_time_s"] == pytest.approx(0.05)
    assert pressure["startup_window_peak_pressure_mpa"] == pytest.approx(1.0)
    assert pressure["global_peak_time_s"] == pytest.approx(2.0)
    assert pressure["global_peak_pressure_mpa"] == pytest.approx(5.0)
    assert pressure["peak_time_s"] == pressure["startup_window_peak_time_s"]


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
    assert sources["pyrogen_surface_heat_power_w"][0] == pytest.approx(0.0)
    assert sources["pyrogen_surface_heat_flux_w_m2"][0] == pytest.approx(0.0)
    assert sources["radiation_heat_power_w"][0] == pytest.approx(0.0)
    assert sources["radiation_heat_flux_w_m2"][0] == pytest.approx(0.0)
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
        diagnostic_disable_momentum=False,
        diagnostic_disable_pyrogen_surface_heating=False,
        diagnostic_disable_adjacent_radiation=False,
        igniter_axial_momentum_fraction=1.0,
    )
    np.testing.assert_allclose(base["P_head"], explicit["P_head"])
    np.testing.assert_allclose(base["mdot_ig"], explicit["mdot_ig"])


def test_source_timeseries_reports_pyrogen_surface_heating_power():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
        snapshots=[
            _snap(
                0.1, [False, False],
                pyrogen_surface_heat_flux=[1000.0, 500.0],
            ),
        ],
    )
    sources = source_timeseries(result)
    assert sources["pyrogen_surface_heat_power_w"][0] == pytest.approx(1500.0)
    assert sources["pyrogen_surface_heat_flux_w_m2"][0] == pytest.approx(1000.0)


def test_source_timeseries_reports_adjacent_radiation_heating_power():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
        snapshots=[
            _snap(
                0.1, [False, False],
                radiation_heat_flux=[2000.0, 1000.0],
            ),
        ],
    )
    sources = source_timeseries(result)
    assert sources["radiation_heat_power_w"][0] == pytest.approx(3000.0)
    assert sources["radiation_heat_flux_w_m2"][0] == pytest.approx(2000.0)


def test_energy_momentum_timeseries_preserves_result_ledgers():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
    )
    result["pyrogen_enthalpy_power"] = np.array([10.0, 20.0])
    result["pyrogen_momentum_residual"] = np.array([0.0, 1.0e-12])
    audit = energy_momentum_timeseries(result)
    np.testing.assert_allclose(audit["pyrogen_enthalpy_power"], [10.0, 20.0])
    np.testing.assert_allclose(audit["pyrogen_momentum_residual"], [0.0, 1.0e-12])
    np.testing.assert_allclose(audit["energy_residual"], [0.0, 0.0])


def test_ambient_initial_gas_ignites_from_pyrogen_surface_heating():
    ambient, _geo, prop = _run_small(
        T_ignition=500.0,
        initial_gas_temperature=propellant_initial_temperature(),
    )
    disabled, _geo, _prop = _run_small(
        T_ignition=500.0,
        initial_gas_temperature=propellant_initial_temperature(),
        diagnostic_disable_pyrogen_surface_heating=True,
    )
    ambient_spread = ignition_spread_metrics(ambient)
    disabled_spread = ignition_spread_metrics(disabled)
    assert ambient["summary"]["initial_gas_temperature"] == pytest.approx(293.0)
    assert np.isfinite(ambient_spread["first_ignition_time_s"])
    assert max(
        np.max(s["pyrogen_surface_heat_flux"]) for s in ambient["snapshots"]
    ) > 0.0
    assert max(
        np.max(s["pyrogen_surface_heat_flux"]) for s in disabled["snapshots"]
    ) == pytest.approx(0.0)
    heated_cells = set()
    for snap in ambient["snapshots"]:
        heated = np.flatnonzero(snap["pyrogen_surface_heat_flux"] > 0.0)
        if heated.size:
            heated_cells.update(int(i) for i in heated)
    assert len(heated_cells) == 1
    assert (
        np.isnan(disabled_spread["first_ignition_time_s"])
        or disabled_spread["first_ignition_time_s"] >= ambient_spread["first_ignition_time_s"]
    )


def propellant_initial_temperature():
    return hasegawa_propellant_1().T_initial


def test_custom_pyrogen_missing_heat_flux_requires_explicit_disable():
    geo, prop, nozzle = _small_motor()
    chamber = _test_chamber(heat_flux_cal_cm2_s=None)
    with pytest.raises(ValueError, match="heat_flux_cal_cm2_s"):
        run_simulation(
            geo, prop, nozzle, chamber,
            T_ignition=500.0,
            initial_gas_temperature=prop.T_initial,
            t_max=0.001,
            P_cutoff=1.0,
            dt_max=2.0e-5,
            burn_update_interval=1,
            snapshot_interval=0.001,
            verbose=False,
        )

    result = run_simulation(
        geo, prop, nozzle, chamber,
        T_ignition=500.0,
        initial_gas_temperature=prop.T_initial,
        diagnostic_disable_pyrogen_surface_heating=True,
        t_max=0.001,
        P_cutoff=1.0,
        dt_max=2.0e-5,
        burn_update_interval=1,
        snapshot_interval=0.001,
        verbose=False,
    )
    assert result["summary"]["diagnostic_disable_pyrogen_surface_heating"]


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


def test_disabling_momentum_removes_igniter_momentum_source():
    baseline, _geo, _prop = _run_small(T_ignition=294.0)
    disabled, _geo, _prop = _run_small(
        T_ignition=294.0,
        diagnostic_disable_momentum=True,
    )
    baseline_momentum = max(np.max(s["momentum_source"]) for s in baseline["snapshots"])
    disabled_momentum = max(np.max(s["momentum_source"]) for s in disabled["snapshots"])
    assert baseline_momentum > 0.0
    assert disabled_momentum == pytest.approx(0.0)
