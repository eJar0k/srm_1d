import numpy as np
import pytest

from srm_1d import run_simulation
from srm_1d.igniter_plenum import PyrogenChamber
from srm_1d.nozzle import Nozzle
from srm_1d.propellant import Pyrogen
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
    assert np.any(last['is_burning'])
    assert np.max(last['T_surf']) >= 294.0
