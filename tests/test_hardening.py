"""
test_hardening.py -- solver-verification PB / P1 hardening (ports of Codex
e46562b / e749f71 / 86c587f, adapted): ``pressure_floor_pa`` kwarg (default
1e3 bit-identical, validated finite > 0), ``burn_update_interval >= 1``
validation, ``D_hyd <= 0 -> (r0, 0)`` in ``burn_rate_cell``.
"""
import numpy as np
import pytest

pytest.importorskip("numba")

from srm_1d import run_simulation  # noqa: E402
from srm_1d.burn_rate import burn_rate_cell  # noqa: E402
from srm_1d.solver import _piso_step_with_energy_diagnostics  # noqa: E402
from tests.test_simulation_phase3 import _small_motor, _test_chamber  # noqa: E402

R_GAS = 300.0
CP_GAS = 2000.0
GAMMA = 1.2


def _floor_step(pressure, floor_pa=None):
    """One venting cell whose raw post-correction pressure undershoots."""
    T = 300.0
    args = (
        np.array([pressure / (R_GAS * T)]), np.zeros(2), np.array([pressure]),
        np.array([T]), np.array([1.0e-3]), np.array([0.035]),
        np.zeros(1), np.zeros(1), np.zeros(2), np.zeros(1),
        0.01, 1.0e-5, np.array([GAMMA]), np.array([R_GAS]),
        np.array([CP_GAS]), np.array([3030.0]),
        0.0, 101325.0, T, 1,
    )
    if floor_pa is None:
        return _piso_step_with_energy_diagnostics(*args)
    return _piso_step_with_energy_diagnostics(*args, 0.0, floor_pa)


class TestPressureFloorKwarg:
    def test_default_floor_is_the_historical_1000_pa_bit_identical(self):
        default = _floor_step(500.0)
        explicit = _floor_step(500.0, floor_pa=1.0e3)
        assert default[2][0] == 1000.0
        for x, y in zip(default, explicit):
            assert np.array_equal(np.asarray(x), np.asarray(y))

    def test_override_changes_the_clamp_value(self):
        out = _floor_step(500.0, floor_pa=750.0)
        assert out[2][0] == pytest.approx(750.0)
        assert out[12] == 1  # n_floor_cells

    @pytest.mark.parametrize("p", [0.0, -1.0, np.nan, np.inf])
    def test_pressure_floor_must_be_finite_and_positive(self, p):
        geo, prop, nozzle = _small_motor()
        with pytest.raises(
                ValueError, match="pressure_floor_pa must be finite and positive"):
            run_simulation(geo, prop, nozzle, _test_chamber(),
                           pressure_floor_pa=p, verbose=False)

    def test_summary_reports_the_floor(self):
        geo, prop, nozzle = _small_motor()
        kw = dict(T_ignition=294.0, t_max=0.001, P_cutoff=1.0,
                  snapshot_interval=0.001, verbose=False)
        r = run_simulation(geo, prop, nozzle, _test_chamber(), **kw)
        assert r['summary']['pressure_floor_pa'] == 1.0e3
        r = run_simulation(geo, prop, nozzle, _test_chamber(),
                           pressure_floor_pa=5.0e3, **kw)
        assert r['summary']['pressure_floor_pa'] == 5.0e3


class TestOtherHardening:
    def test_nonpositive_hydraulic_diameter_disables_erosive_term(self):
        A_SR = 4.9e-3 / (4.9e6 ** 0.3)
        tabs = (np.array([0.0]), np.array([2e7]), np.array([A_SR]),
                np.array([0.3]), 1)
        for D in (0.0, -0.01):
            r_total, r_erosive = burn_rate_cell(
                5.0e6, 5.0e5, D, 0.5, 20e-6,
                0.4943, 0.3685, 2060.0, 3041.0, 1000.0,
                1700.0, 1500.0, 293.0, *tabs, 0.45, 9.0e-5, 0.0, 0)
            assert r_total == pytest.approx(A_SR * 5.0e6 ** 0.3)
            assert r_erosive == 0.0

    @pytest.mark.parametrize("k", [0, -3, 1.5])
    def test_burn_update_interval_must_be_positive_integer(self, k):
        geo, prop, nozzle = _small_motor()
        with pytest.raises(ValueError, match="burn_update_interval"):
            run_simulation(geo, prop, nozzle, _test_chamber(),
                           burn_update_interval=k, verbose=False)
