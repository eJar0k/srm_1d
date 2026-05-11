import numpy as np
import pytest

from srm_1d import run_simulation
from srm_1d.igniter_plenum import PyrogenChamber
from srm_1d.nozzle import Nozzle
from srm_1d.propellant import Pyrogen
from srm_1d.simulation import (
    _goodman_ignition_sources_and_mass,
    _cal_cm2_s_to_w_m2,
    _pyrogen_surface_heat_power,
    _pyrogen_surface_thermal_sink,
    _thermal_source_power,
)
from srm_1d.tests._motor_fixtures import (
    hasegawa_propellant_1,
    single_cylinder_geo,
)


def _test_chamber():
    pyro = Pyrogen(
        name="phase3-test",
        a=3.0e-5,
        n=0.5,
        rho=1700.0,
        T_flame=2800.0,
        M=0.030,
        gamma=1.25,
        impetus_W=5000.0,
        heat_flux_cal_cm2_s=69.4,
    )
    return PyrogenChamber(
        pyrogen=pyro,
        m_pyrogen_initial=0.003,
        A_burn_initial=5.0e-4,
        A_throat=2.0e-5,
        V_plenum=3.0e-6,
        burn_law="end_burning",
    )


def _small_motor():
    geo = single_cylinder_geo(
        D_bore=0.030, D_outer=0.060, length=0.120,
        target_propellant_cells=12,
    )
    prop = hasegawa_propellant_1()
    nozzle = Nozzle(D_throat=0.010, D_exit=0.020, efficiency=0.95)
    return geo, prop, nozzle


def test_legacy_igniter_kwargs_are_not_accepted():
    geo, prop, nozzle = _small_motor()
    with pytest.raises(TypeError):
        run_simulation(
            geo, prop, nozzle, _test_chamber(),
            igniter_tau=0.1,
            t_max=0.001,
        )


def test_pyrogen_heat_flux_unit_conversion():
    assert _cal_cm2_s_to_w_m2(1.0) == pytest.approx(41840.0)
    assert _cal_cm2_s_to_w_m2(69.4) == pytest.approx(69.4 * 41840.0)


def test_pyrogen_surface_heat_power_uses_measured_flux_when_uncapped():
    power, flux = _pyrogen_surface_heat_power(
        0.1, 2800.0, 300.0, 0.1, 0.01, 2000.0, 1.0e6,
    )
    assert power == pytest.approx(1000.0)
    assert flux == pytest.approx(1.0e6)


def test_pyrogen_surface_heat_power_is_sensible_enthalpy_capped():
    power, flux = _pyrogen_surface_heat_power(
        0.001, 800.0, 300.0, 0.1, 0.01, 2000.0, 1.0e9,
    )
    assert power == pytest.approx(1000.0)
    assert flux == pytest.approx(1.0e6)


def test_pyrogen_surface_thermal_sink_conserves_energy_units():
    power_w = 1200.0
    sink = _pyrogen_surface_thermal_sink(
        power_w, Cp_gas=2000.0, dx=0.02, pyrogen_thermal_source=100.0,
    )
    assert sink * 2000.0 * 0.02 == pytest.approx(power_w)


def test_thermal_source_power_matches_solver_units():
    thermal_source = np.array([10.0, -2.5])
    assert _thermal_source_power(thermal_source, 2000.0, 0.01, 2) == pytest.approx(150.0)


def test_adjacent_radiation_heats_only_neighbors_and_conserves_sink():
    N = 3
    dx = 0.01
    P = np.full(N, 101325.0)
    T = np.full(N, 300.0)
    T_surf = np.full(N, 293.0)
    delta = np.full(N, 1.0e-6)
    has_ignited = np.array([False, True, False])
    is_burning = np.array([False, True, False])
    is_grain = np.array([True, True, True])
    ignition_time = np.full(N, 1.0e10)
    r_total = np.array([0.0, 0.01, 0.0])
    r_erosive = np.zeros(N)
    mass_source = np.zeros(N)
    thermal_source = np.zeros(N)
    C_burn = np.ones(N)
    endface = np.zeros(N)
    pyrogen_flux = np.zeros(N)
    radiation_flux = np.zeros(N)
    radiation_sink_power = np.zeros(N)
    radiation_emitter = np.zeros(N, dtype=np.bool_)
    x = np.array([0.005, 0.015, 0.025])
    Re = np.zeros(N)
    D_hyd = np.full(N, 0.03)
    f = np.zeros(N)

    out = _goodman_ignition_sources_and_mass(
        P, T, T_surf, delta, has_ignited, is_burning, is_grain,
        ignition_time, r_total, r_erosive,
        mass_source, thermal_source,
        C_burn, endface, pyrogen_flux,
        radiation_flux, radiation_sink_power, radiation_emitter,
        x, Re, D_hyd, f,
        0.0, 1.0e-6, 1700.0, 3041.0, 293.0,
        0.5, 0.3685, 50e-6, 0.45, 1.0e-7, 0.3,
        10000.0, N, dx, 0.0, 300.0, 2060.0,
        0.0, 0.45,
    )

    radiation_heat_power = out[4]
    assert radiation_flux[0] > 0.0
    assert radiation_flux[1] == pytest.approx(0.0)
    assert radiation_flux[2] > 0.0
    assert radiation_sink_power[1] == pytest.approx(radiation_heat_power)
    assert radiation_sink_power[0] == pytest.approx(0.0)
    assert radiation_sink_power[2] == pytest.approx(0.0)
    assert thermal_source[1] < 1700.0 * r_total[1] * C_burn[1] * 3041.0


def test_pyrogen_driven_run_reports_ignition_and_pyrogen_state():
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0,
        t_max=0.010,
        P_cutoff=1.0,
        dt_max=2.0e-5,
        burn_update_interval=1,
        snapshot_interval=0.002,
        cfl_target=0.5,
    )

    summary = result['summary']
    assert summary['pyrogen_mass_burned'] > 0.0
    assert summary['pyrogen_peak_P'] > 101325.0
    assert np.max(result['P_head']) > 101325.0
    assert np.max(result['mdot_ig']) > 0.0
    assert len(result['P_ig']) == len(result['time'])

    assert result['snapshots']
    last = result['snapshots'][-1]
    assert 'T_surf' in last
    assert 'pyrogen_surface_heat_flux' in last
    assert 'radiation_heat_flux' in last
    assert np.any(last['is_burning'])
    assert np.max(last['T_surf']) >= 294.0
    assert np.max(last['pyrogen_surface_heat_flux']) >= 0.0
    assert result['summary']['radiation_emissivity'] == pytest.approx(prop.radiation_emissivity)
    assert 'pyrogen_enthalpy_power' in result
    assert 'gas_surface_heat_sink_power' in result
    assert 'energy_residual' in result
    assert 'pyrogen_momentum_residual' in result
