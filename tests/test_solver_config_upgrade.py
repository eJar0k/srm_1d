"""test_solver_config_upgrade.py -- docs/igniter_provider I4.

The plugin's GUI run-config defaults are derived from ``run_simulation`` (they
once drifted to the pre-v0.7.5 50 um / 0.45 / 850 K), and a saved config that
still holds those untouched superseded defaults is upgraded.
"""
import copy
import inspect

import pytest

pytest.importorskip("skfmm")
pytest.importorskip("numba")

from srm_1d import srm1d_plugin  # noqa: E402
from srm_1d.fmm_grain import _setup_openmotor_path  # noqa: E402
from srm_1d.simulation import run_simulation  # noqa: E402

# The user's saved preferences (solverConfigs.srm_1d-transient), copied here.
USER_SAVED = {'P_cutoff': 344750.0, 'T_ignition': 850.0, 'cfl_target': 0.3,
              'kappa': 0.45, 'roughness': 50.0, 't_max': 10.0}


@pytest.fixture(scope="module")
def solver():
    try:
        _setup_openmotor_path()
    except ImportError as exc:       # pragma: no cover
        pytest.skip(f"openMotor checkout unavailable: {exc}")
    from motorlib import solvers  # type: ignore
    return solvers.get_solver(srm1d_plugin.SOLVER_NAME)


def test_schema_defaults_equal_run_simulation_defaults(solver):
    schema = solver.get_config_schema()
    params = inspect.signature(run_simulation).parameters
    got = schema.getProperties()
    assert set(got) == {'t_max', 'P_cutoff', 'cfl_target', 'roughness',
                        'kappa', 'T_ignition', 'ignition_heat_transfer'}
    # IG3: the ignition-h closure is a choice field whose default is read from
    # run_simulation and whose choices are exactly srm_1d's validated modes.
    from srm_1d import simulation as _sim
    assert got['ignition_heat_transfer'] == params['ignition_heat_transfer'].default == 'peretz'
    assert set(schema.props['ignition_heat_transfer'].values) == set(_sim._IGN_H_MODES)
    for key, value in got.items():
        if key == 'ignition_heat_transfer':
            continue
        want = params[key].default
        if key == 'roughness':
            want = want * 1.0e6                      # m -> um
        assert value == pytest.approx(want, rel=1e-12), key
    # the canonical v0.7.5 knobs, spelled out
    assert got['roughness'] == pytest.approx(32.0)
    assert (got['kappa'], got['T_ignition']) == (0.44, 756.0)


def test_superseded_values_are_not_current_defaults():
    current = srm1d_plugin._solver_config_defaults()
    for key, olds in srm1d_plugin.SUPERSEDED_SOLVER_CONFIG.items():
        for old in olds:
            assert old != pytest.approx(current[key], rel=1e-9), key


def test_users_saved_config_upgrades_exactly_the_three_stale_keys(solver):
    saved = copy.deepcopy(USER_SAVED)
    saved['P_cutoff'] = 344750.0
    before = copy.deepcopy(saved)
    new, upgraded = solver.upgrade_solver_config(saved)
    assert saved == before                                   # pure
    assert sorted(upgraded) == ['T_ignition', 'kappa', 'roughness']
    cur = srm1d_plugin._solver_config_defaults()
    for key in upgraded:
        assert new[key] == cur[key]
    assert new['roughness'] == pytest.approx(32.0)
    assert (new['kappa'], new['T_ignition']) == (0.44, 756.0)
    for key in ('P_cutoff', 'cfl_target', 't_max'):          # untouched
        assert new[key] == saved[key]
    assert new['P_cutoff'] == 344750.0


def test_customised_values_are_untouched(solver):
    vals = dict(USER_SAVED, roughness=40.0, kappa=0.5, T_ignition=900.0)
    new, upgraded = solver.upgrade_solver_config(vals)
    assert upgraded == []
    assert new == vals


def test_idempotent_and_tolerant(solver):
    once, _ = solver.upgrade_solver_config(USER_SAVED)
    twice, up2 = solver.upgrade_solver_config(once)
    assert up2 == [] and twice == once
    # GUI float noise still recognised; missing keys / empty dict fine.
    noisy = {'roughness': 50.00000000000001}
    assert solver.upgrade_solver_config(noisy)[1] == ['roughness']
    assert solver.upgrade_solver_config({}) == ({}, [])
    assert solver.upgrade_solver_config(None) == ({}, [])
