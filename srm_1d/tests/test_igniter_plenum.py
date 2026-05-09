import math

import numpy as np
import pytest

from srm_1d import Pyrogen
from srm_1d.igniter_plenum import (
    PyrogenChamber,
    _choked_orifice_mdot,
    _step_plenum_ode,
    chamber_params,
    initial_plenum_state,
    pyrogen_params,
    step_plenum,
    sutton_pyrogen_mass,
)
from srm_1d.propellant import R_UNIVERSAL


def _test_pyrogen(**overrides):
    values = {
        "name": "test",
        "a": 2.0e-6,
        "n": 0.5,
        "rho": 1700.0,
        "T_flame": 2400.0,
        "M": 0.030,
        "gamma": 1.25,
        "impetus_W": 5000.0,
    }
    values.update(overrides)
    return Pyrogen(**values)


def _analytic_choked_steady_pressure(pyro, A_burn, A_throat):
    Gamma = math.sqrt(
        pyro.gamma
        * (2.0 / (pyro.gamma + 1.0))
        ** ((pyro.gamma + 1.0) / (pyro.gamma - 1.0))
    )
    R_specific = R_UNIVERSAL / pyro.M
    base = pyro.rho * pyro.a * A_burn * math.sqrt(
        R_specific * pyro.T_flame
    ) / (A_throat * Gamma)
    return base ** (1.0 / (1.0 - pyro.n))


def test_pyrogen_dataclass_importable():
    pyro = _test_pyrogen()
    assert pyro.name == "test"
    assert pyro.impetus_W == pytest.approx(5000.0)


def test_closed_bomb_conserves_mass_and_pressure_rises():
    pyro = _test_pyrogen()
    chamber = PyrogenChamber(
        pyrogen=pyro,
        m_pyrogen_initial=0.010,
        A_burn_initial=1.0e-4,
        A_throat=0.0,
        V_plenum=1.0e-5,
        burn_law="end_burning",
    )
    state = initial_plenum_state(chamber, P_initial=101325.0, T_initial=300.0)
    total_initial = state[0] + state[1]
    pressures = []

    for _ in range(200):
        state, mdot_out, _mdot_gen, P_ig = step_plenum(
            chamber, state, dt=2.5e-4, P_main=101325.0
        )
        assert mdot_out == pytest.approx(0.0)
        assert np.all(state >= 0.0)
        pressures.append(P_ig)

    assert state[0] + state[1] == pytest.approx(total_initial, rel=1e-6)
    assert pressures[-1] > pressures[0]
    assert np.all(np.diff(np.array(pressures)) >= -1e-6)


def test_constant_area_choked_outflow_matches_steady_pressure():
    pyro = _test_pyrogen(a=5.0e-6)
    A_burn = 1.0e-4
    A_throat = 5.0e-7
    chamber = PyrogenChamber(
        pyrogen=pyro,
        m_pyrogen_initial=0.100,
        A_burn_initial=A_burn,
        A_throat=A_throat,
        V_plenum=1.0e-6,
        burn_law="end_burning",
    )
    state = initial_plenum_state(chamber, P_initial=101325.0, T_initial=300.0)
    p_params = pyrogen_params(pyro)
    c_params = chamber_params(chamber)

    P_ig = 0.0
    for _ in range(8000):
        state, _mdot_out, _mdot_gen, P_ig = _step_plenum_ode(
            state, p_params, c_params, 1.0e-5, 101325.0
        )

    expected = _analytic_choked_steady_pressure(pyro, A_burn, A_throat)
    assert P_ig == pytest.approx(expected, rel=0.05)
    assert state[2] == pytest.approx(pyro.T_flame, rel=0.05)


def test_burnout_drives_pyrogen_and_outflow_to_zero_smoothly():
    pyro = _test_pyrogen(a=2.0e-5)
    chamber = PyrogenChamber(
        pyrogen=pyro,
        m_pyrogen_initial=2.0e-5,
        A_burn_initial=5.0e-4,
        A_throat=2.0e-6,
        V_plenum=2.0e-6,
        burn_law="end_burning",
    )
    state = initial_plenum_state(chamber, P_initial=101325.0, T_initial=300.0)
    p_params = pyrogen_params(pyro)
    c_params = chamber_params(chamber)
    peak_mdot = 0.0
    mdot_out = 0.0

    for _ in range(10000):
        state, mdot_out, _mdot_gen, _P_ig = _step_plenum_ode(
            state, p_params, c_params, 2.0e-5, 101325.0
        )
        assert np.all(state >= 0.0)
        peak_mdot = max(peak_mdot, mdot_out)

    assert state[0] == pytest.approx(0.0, abs=1e-10)
    assert peak_mdot > 0.0
    assert mdot_out < peak_mdot * 1e-2


def test_subsonic_orifice_fallback_matches_formula():
    P_ig = 1.0e6
    P_main = 0.90e6
    T_ig = 2200.0
    A_t = 1.0e-6
    gamma = 1.25
    M = 0.030

    pressure_ratio = P_main / P_ig
    term = pressure_ratio ** (2.0 / gamma) - pressure_ratio ** (
        (gamma + 1.0) / gamma
    )
    expected = P_ig * A_t * math.sqrt(
        (2.0 * gamma / (gamma - 1.0)) * M / (R_UNIVERSAL * T_ig)
    ) * math.sqrt(term)

    mdot = _choked_orifice_mdot(
        P_ig, T_ig, A_t, gamma, R_UNIVERSAL, M, P_main
    )
    assert mdot == pytest.approx(expected, rel=1e-12)
    assert _choked_orifice_mdot(
        P_ig, T_ig, A_t, gamma, R_UNIVERSAL, M, P_ig
    ) == pytest.approx(0.0)


def test_sutton_pyrogen_mass_returns_kg():
    expected_grams = 0.12 * 100.0 ** 0.7
    assert sutton_pyrogen_mass(100.0) == pytest.approx(
        expected_grams * 1e-3
    )


def test_unsupported_burn_law_rejected_at_python_boundary():
    with pytest.raises(ValueError, match="burn_law"):
        PyrogenChamber(
            pyrogen=_test_pyrogen(),
            m_pyrogen_initial=0.001,
            A_burn_initial=1.0e-4,
            A_throat=1.0e-6,
            V_plenum=1.0e-6,
            burn_law="cylindrical",
        )
