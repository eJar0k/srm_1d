"""
test_rt1_recovery.py -- RT1 opt-in ignition driver temperature (adiabatic-wall T_aw).

``ignition_driver_temperature='static'`` (default) is the shipped ignition heating;
``'recovery'`` drives it with T_aw = T + r u^2/(2 Cp), r = Pr^(1/3) turbulent (bulk
Re > 2300) or Pr^(1/2) laminar.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d import simulation as sim  # noqa: E402


class TestDriverTemperature:
    def test_static_is_identity(self):
        T, r, turb = sim._driver_temperature(0, 1, 2500.0, 300.0, 2000.0, 0.7, 1.2, 1e6)
        assert T == 2500.0 and r == 0.0 and turb == 0.0

    def test_turbulent_peretz_uses_cube_root_and_local_gamma(self):
        g = 1.2
        Pr = g / (1.77 * g - 0.45)
        T, r, turb = sim._driver_temperature(1, 1, 2500.0, 300.0, 2000.0, 0.7, g, 1e5)
        assert turb == 1.0
        assert r == pytest.approx(Pr ** (1.0 / 3.0))
        assert T == pytest.approx(2500.0 + r * 300.0 ** 2 / (2 * 2000.0))

    def test_laminar_uses_square_root_and_run_pr_for_gnielinski(self):
        T, r, turb = sim._driver_temperature(1, 0, 1000.0, 100.0, 1500.0, 0.64, 1.2, 1500.0)
        assert turb == 0.0
        assert r == pytest.approx(0.64 ** 0.5)
        assert T == pytest.approx(1000.0 + r * 100.0 ** 2 / 3000.0)

    def test_transition_is_2300(self):
        _, _, below = sim._driver_temperature(1, 1, 1000.0, 10.0, 1500.0, 0.7, 1.2, 2299.0)
        _, _, above = sim._driver_temperature(1, 1, 1000.0, 10.0, 1500.0, 0.7, 1.2, 2301.0)
        assert below == 0.0 and above == 1.0

    def test_recovery_factor_is_below_one_and_T_aw_below_stagnation(self):
        T, r, _ = sim._driver_temperature(1, 1, 2000.0, 500.0, 2000.0, 0.7, 1.2, 1e6)
        assert 0.0 < r < 1.0
        assert 2000.0 < T < 2000.0 + 500.0 ** 2 / (2 * 2000.0)


def test_unknown_driver_mode_rejected():
    with pytest.raises(ValueError):
        sim._ignition_driver_code('bogus')


@pytest.fixture(scope="module")
def short_runs():
    from srm_1d.openmotor_adapter import run_from_ric
    ric = os.path.join(os.path.dirname(__file__), "..", "motors", "machbusterNew.ric")
    out = {}
    for d in ('static', 'recovery'):
        out[d] = run_from_ric(
            ric, target_propellant_cells=60, t_max=0.012, roughness=32e-6, kappa=0.44,
            T_ignition=756.0, print_interval=0.0, verbose=False, pyrogen='bpnv',
            pyrogen_mass=2.623e-3, ignition_driver_temperature=d,
            history_capacity=100000)[0]
    return out


class TestMotorPlumbing:
    def test_stamp(self, short_runs):
        for d, r in short_runs.items():
            assert r['summary']['ignition_driver_temperature'] == d

    def test_recovery_changes_the_ignition(self, short_runs):
        a = np.asarray(short_runs['static']['ignition_time_by_cell'])
        b = np.asarray(short_runs['recovery']['ignition_time_by_cell'])
        assert not np.array_equal(a, b)
        assert b.min() < a.min()          # the hotter driver lights the first cell earlier

    def test_diag_columns_record_driver_and_regime(self, short_runs):
        d = np.asarray(short_runs['recovery']['ignition_diag'])
        lit = np.isfinite(d[:, 0])
        assert lit.any()
        T_gas, T_drv, r = d[lit, 4], d[lit, 10], d[lit, 11]
        assert np.all(T_drv >= T_gas) and np.all(r > 0.0) and np.all(r < 1.0)
        assert set(np.unique(d[lit, 13])) <= {0.0, 1.0}
        s = np.asarray(short_runs['static']['ignition_diag'])
        ls = np.isfinite(s[:, 0])
        assert np.allclose(s[ls, 10], s[ls, 4]) and np.all(s[ls, 11] == 0.0)

    def test_bad_mode_raises(self):
        from srm_1d.openmotor_adapter import run_from_ric
        ric = os.path.join(os.path.dirname(__file__), "..", "motors", "machbusterNew.ric")
        with pytest.raises(ValueError):
            run_from_ric(ric, t_max=0.001, verbose=False, ignition_driver_temperature='bogus')
