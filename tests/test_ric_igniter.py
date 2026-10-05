"""
test_ric_igniter.py — v0.8.0 Phase 4 igniter-as-data (D3).

Verifies the motor ``data.igniter`` block: default seeding, round-trip
through the Pyrogen material + chamber sizing, that all repo motors carry a
parseable igniter block after migration, and that run_from_ric reads the
igniter from the motor file when no pyrogen kwarg is given.
"""

import os
import glob
import dataclasses
import warnings

import pytest
import yaml

from srm_1d.openmotor_adapter import (
    default_igniter_block, load_igniter, _seed_motor_igniter,
    _pyrogen_to_block, _block_to_pyrogen, load_ric, build_pyrogen_chamber,
    StaleIgniterMaterialWarning,
)
from srm_1d.propellant import Pyrogen

MOTORS = os.path.join(os.path.dirname(__file__), '..', 'motors')


def test_default_block_has_material_and_sizing():
    block = default_igniter_block('bpnv')
    assert block['pyrogen']['name'] == 'BPNV'
    assert block['injection_topology'] == 'forward_plenum'
    assert block['mass'] == -1.0  # auto sentinel
    for key in ('throat_area', 'volume', 'burn_area', 'cartridge_length_m'):
        assert block[key] == -1.0


def test_pyrogen_block_round_trip():
    block = default_igniter_block('bpnv')
    pyro = _block_to_pyrogen(block['pyrogen'])
    assert isinstance(pyro, Pyrogen)
    re = _pyrogen_to_block(pyro)
    # Round-trip preserves every material field.
    for key, val in block['pyrogen'].items():
        assert re[key] == val


def test_load_igniter_maps_sentinels_to_auto():
    block = default_igniter_block('bpnv')
    motor = {'igniter': block}
    pyro, sizing = load_igniter(motor)
    assert pyro.name == 'BPNV'
    # -1.0 sentinels become None (auto) for build_pyrogen_chamber.
    assert sizing['pyrogen_mass'] is None
    assert sizing['pyrogen_throat_area'] is None
    assert sizing['pyrogen_burn_area'] is None
    assert sizing['injection_topology'] == 'forward_plenum'
    # cartridge_length_m keeps the -1.0 sentinel (build derives it).
    assert sizing['cartridge_length_m'] == -1.0


def test_load_igniter_none_when_absent():
    assert load_igniter({'propellant': {}}) is None


def test_seed_motor_igniter_idempotent():
    data = {'propellant': {}}
    _seed_motor_igniter(data)
    assert 'igniter' in data
    first = data['igniter']
    _seed_motor_igniter(data)
    assert data['igniter'] is first  # not overwritten


def test_explicit_sizing_survives_load():
    block = default_igniter_block('bpnv')
    block['mass'] = 0.006
    block['injection_topology'] = 'head_basket'
    pyro, sizing = load_igniter({'igniter': block})
    assert sizing['pyrogen_mass'] == pytest.approx(0.006)
    assert sizing['injection_topology'] == 'head_basket'


def test_all_repo_motors_have_parseable_igniter():
    rics = glob.glob(os.path.join(MOTORS, '*.ric'))
    assert rics
    for ric in rics:
        motor = load_ric(ric)
        ig = load_igniter(motor)
        assert ig is not None, f"{os.path.basename(ric)} has no igniter block"
        pyro, sizing = ig
        assert pyro.name and pyro.a > 0 and pyro.gamma > 1.0


def test_build_chamber_from_motor_file_igniter():
    """The loaded igniter builds a valid PyrogenChamber (self-describing)."""
    from tests.test_simulation_phase3 import _small_motor
    geo, _, nozzle = _small_motor()
    motor = load_ric(os.path.join(MOTORS, 'hasegawa_a.ric'))
    pyro, sizing = load_igniter(motor)
    chamber = build_pyrogen_chamber(pyro, geo, nozzle, **sizing)
    assert chamber.m_pyrogen_initial > 0
    assert chamber.injection_topology == 'forward_plenum'


# ---------------------------------------------------------------------------
# sv/T1 -- lossless igniter block round-trip + stale-material warning
# ---------------------------------------------------------------------------

_INIT_FIELDS = [f for f in dataclasses.fields(Pyrogen) if f.init]


def _every_field_nondefault():
    """A Pyrogen with EVERY init field set to a value != its default.

    Generated from ``dataclasses.fields`` so a field added to ``Pyrogen``
    later is covered automatically (and fails the round-trip test if the
    adapter drops it, as ``gas_mass_fraction`` once was)."""
    overrides = {'form': 'powder', 'heat_delivery_mode': 'radiation'}
    kwargs = {}
    for i, f in enumerate(_INIT_FIELDS):
        if f.name == 'name':
            val = 'ROUNDTRIP'
        elif f.name in overrides:
            val = overrides[f.name]
        elif f.default is dataclasses.MISSING:
            val = 1.0 + 0.37 * (i + 1)            # required float field
        elif f.default is None:
            val = 1234.5 + i                       # Optional[float]
        elif isinstance(f.default, str):
            val = f.default + '_alt'
        elif isinstance(f.default, float):
            val = f.default * 1.37 + 0.0123 * (i + 1)
        else:
            raise AssertionError(
                f"test helper cannot synthesize a value for Pyrogen.{f.name} "
                f"(default {f.default!r}); extend _every_field_nondefault")
        kwargs[f.name] = val
    p = Pyrogen(**kwargs)
    for f in _INIT_FIELDS:                         # really non-default
        if f.default is not dataclasses.MISSING:
            assert getattr(p, f.name) != f.default, f.name
    return p


def test_pyrogen_block_round_trip_every_field():
    p = _every_field_nondefault()
    block = _pyrogen_to_block(p)
    # direct round-trip
    assert dataclasses.asdict(_block_to_pyrogen(block)) == dataclasses.asdict(p)
    # and through the YAML text a .ric is actually saved as
    reloaded = yaml.safe_load(yaml.safe_dump(block))
    q = _block_to_pyrogen(reloaded)
    for f in _INIT_FIELDS:
        assert getattr(q, f.name) == getattr(p, f.name), f.name
    assert q.gas_mass_fraction == p.gas_mass_fraction != 1.0


def test_pyrogen_block_keys_equal_init_fields():
    names = {f.name for f in _INIT_FIELDS}
    assert set(_pyrogen_to_block(_every_field_nondefault())) == names
    block = default_igniter_block('bpnv')['pyrogen']
    assert set(block) == names
    assert block['gas_mass_fraction'] == pytest.approx(0.781)


def _stale_motor(**pyrogen_edits):
    block = default_igniter_block('bpnv')
    block['pyrogen'].update(pyrogen_edits)
    return {'igniter': block}


def test_current_library_block_emits_no_stale_warning():
    with warnings.catch_warnings():
        warnings.simplefilter('error', StaleIgniterMaterialWarning)
        pyro, _ = load_igniter({'igniter': default_igniter_block('bpnv')})
    assert pyro.gas_mass_fraction == pytest.approx(0.781)


def test_prefix_bfc2f3f_block_emits_one_stale_warning_and_keeps_values():
    motor = _stale_motor(M=0.03)
    del motor['igniter']['pyrogen']['gas_mass_fraction']   # key absent -> 1.0
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        pyro, _ = load_igniter(motor)
    stale = [w for w in rec if issubclass(w.category, StaleIgniterMaterialWarning)]
    assert len(stale) == 1
    msg = str(stale[0].message)
    assert 'BPNV' in msg and 'bpnv.yaml' in msg
    assert 'M: embedded 0.03 vs library 0.0482' in msg
    assert 'gas_mass_fraction: embedded 1.0 vs library 0.781' in msg
    assert msg.endswith(
        "The embedded values are used. Re-select the igniter material or "
        "re-seed the motor's igniter block to pick up library corrections.")
    # warn only -- the embedded values are what is used
    assert pyro.M == 0.03 and pyro.gas_mass_fraction == 1.0


def test_stale_warning_match_rules():
    # name is matched case-insensitively and itself is not a compared field
    with pytest.warns(StaleIgniterMaterialWarning, match='kappa_jet'):
        load_igniter(_stale_motor(name='bpnv', kappa_jet=3.0))
    # float tolerance: 1e-12 relative is "equal"
    with warnings.catch_warnings():
        warnings.simplefilter('error', StaleIgniterMaterialWarning)
        load_igniter(_stale_motor(M=0.0482 * (1.0 + 1e-12)))
    # None vs a value is a difference (library Cp_gas is absent -> None)
    with pytest.warns(StaleIgniterMaterialWarning, match='Cp_gas'):
        load_igniter(_stale_motor(Cp_gas=2000.0))
    # a string field
    with pytest.warns(StaleIgniterMaterialWarning, match='heat_delivery_mode'):
        load_igniter(_stale_motor(heat_delivery_mode='radiation'))
    # no library entry of that name -> no warning
    with warnings.catch_warnings():
        warnings.simplefilter('error', StaleIgniterMaterialWarning)
        load_igniter(_stale_motor(name='MY_CUSTOM_IGNITER', M=0.03))
