"""Permit rule of srm_1d.tools.refresh_ric_igniters (docs/igniter_provider I3b).

An ADDED ``/config/<key>`` is permitted only when its value equals the fork's
``MotorConfig`` default for that key (read from the fork, not hardcoded);
anything else outside the igniter material / openMotor-added allow-list is a
STOP."""
import pytest

pytest.importorskip("skfmm")
pytest.importorskip("numba")

from srm_1d.fmm_grain import _setup_openmotor_path  # noqa: E402
from srm_1d.tools.refresh_ric_igniters import (  # noqa: E402
    _is_permitted, diff_paths, motor_config_defaults)


@pytest.fixture(scope="module")
def defaults():
    try:
        _setup_openmotor_path()
    except ImportError as exc:       # pragma: no cover
        pytest.skip(f"openMotor checkout unavailable: {exc}")
    import motorlib.motor
    return motor_config_defaults(motorlib.motor)


def test_defaults_come_from_motorconfig(defaults):
    assert {'maxMachNumber', 'timestep', 'ambPressure', 'taperSlices'} <= set(defaults)


def test_added_config_key_at_default_is_permitted(defaults):
    after = {'config': {'maxMachNumber': defaults['maxMachNumber']}}
    assert _is_permitted('/config/maxMachNumber', 'added', False, after, defaults)


def test_added_config_key_not_at_default_is_refused(defaults):
    after = {'config': {'maxMachNumber': defaults['maxMachNumber'] + 1.0}}
    assert not _is_permitted('/config/maxMachNumber', 'added', False, after, defaults)


def test_config_change_or_removal_is_refused(defaults):
    after = {'config': {'maxMachNumber': defaults['maxMachNumber']}}
    for kind in ('changed', 'removed'):
        assert not _is_permitted('/config/maxMachNumber', kind, False, after, defaults)


def test_unknown_config_key_and_other_sections_are_refused(defaults):
    after = {'config': {'notAKey': 0.0}, 'nozzle': {'throat': 0.0}}
    assert not _is_permitted('/config/notAKey', 'added', False, after, defaults)
    assert not _is_permitted('/nozzle/throat', 'added', False, after, defaults)
    assert not _is_permitted('/config/a/b', 'added', False, after, defaults)


def test_existing_rules_unchanged(defaults):
    assert _is_permitted('/igniter/pyrogen/M', 'changed', True, {}, defaults)
    assert not _is_permitted('/igniter/pyrogen/M', 'changed', False, {}, defaults)
    assert _is_permitted('/solverConfigs', 'added', False, {}, defaults)
    assert not _is_permitted('/solverConfigs', 'changed', False, {}, defaults)


def test_diff_paths_reports_added_removed_changed():
    d = dict(diff_paths({'a': 1, 'b': 2, 'c': [1]}, {'a': 1, 'b': 3, 'd': 0, 'c': [2]}))
    assert d == {'/b': 'changed', '/c[0]': 'changed', '/d': 'added'}
