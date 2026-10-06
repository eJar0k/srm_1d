"""
test_igniter_provider.py -- plugin-provided igniter (docs/igniter_provider).

srm_1d is the single source of truth for the igniter model: the ``Pyrogen``
fields, the chamber sizing defaults, the built-in material library, the GUI
metadata, and the history of SUPERSEDED library snapshots. This module pins

* the GUI metadata <-> ``Pyrogen`` / ``_IGNITER_SIZING_DEFAULTS`` coverage and
  that every library value fits its schema range (``FloatProperty.setValue``
  silently DROPS an out-of-range value, so a bad range would corrupt a GUI
  round-trip without any error),
* the five optional ``get_igniter_*`` hooks on the registered solver,
* ``upgrade_igniter_block`` (exact-match, pure, idempotent), and
* ``load_igniter``'s upgrade-then-stale-check behaviour.

Needs an openMotor checkout (``motorlib``) for the schema/hook tests; skipped
cleanly otherwise, exactly like tests/test_srm1d_plugin.py.
"""

import copy
import dataclasses
import warnings

import pytest

pytest.importorskip("skfmm")
pytest.importorskip("numba")

from srm_1d.fmm_grain import _setup_openmotor_path  # noqa: E402
from srm_1d import srm1d_plugin  # noqa: E402
from srm_1d.igniter_plenum import _VALID_TOPOLOGIES, _burn_law_code  # noqa: E402
from srm_1d.simulation import _heat_delivery_code  # noqa: E402
from srm_1d.propellant import Pyrogen  # noqa: E402
from srm_1d.openmotor_adapter import (  # noqa: E402
    SUPERSEDED_PYROGENS, IgniterMaterialUpgradedWarning,
    StaleIgniterMaterialWarning, _IGNITER_SIZING_DEFAULTS, _pyrogen_to_block,
    default_igniter_block, library_materials, load_igniter, load_pyrogen,
    upgrade_igniter_block,
)

CHAMBER_KEYS = tuple(_IGNITER_SIZING_DEFAULTS)


@pytest.fixture(scope="module")
def solver():
    """The registered srm_1d solver (needs motorlib; skip if unavailable)."""
    try:
        _setup_openmotor_path()
    except ImportError as exc:       # pragma: no cover - checkout absent
        pytest.skip(f"openMotor checkout unavailable: {exc}")
    import srm_1d.srm1d_plugin  # noqa: F401  (registers on import)
    from motorlib import solvers  # type: ignore
    s = solvers.get_solver(srm1d_plugin.SOLVER_NAME)
    assert s is not None
    return s


def _legacy_block(**chamber):
    """A motor ``igniter`` block whose material is the pre-bfc2f3f BPNV."""
    block = {'pyrogen': copy.deepcopy(SUPERSEDED_PYROGENS['BPNV'][0])}
    block.update(copy.deepcopy(_IGNITER_SIZING_DEFAULTS))
    block.update(chamber)
    return block


def _current_bpnv():
    return _pyrogen_to_block(load_pyrogen('bpnv'))


# ----------------------------------------------------------------
# (a) metadata coverage
# ----------------------------------------------------------------

def test_material_meta_plus_passthrough_covers_pyrogen_fields():
    init_fields = {f.name for f in dataclasses.fields(Pyrogen) if f.init}
    meta = set(srm1d_plugin._PYROGEN_FIELD_META)
    passthrough = set(srm1d_plugin._PYROGEN_GUI_PASSTHROUGH)
    assert meta | passthrough == init_fields
    assert not (meta & passthrough)
    # Every passthrough field carries a stated reason.
    assert all(srm1d_plugin._PYROGEN_GUI_PASSTHROUGH.values())
    # The field added by bfc2f3f is editable, not silently passed through.
    assert 'gas_mass_fraction' in meta


def test_chamber_meta_matches_sizing_defaults():
    assert set(srm1d_plugin._IGNITER_CHAMBER_META) == set(_IGNITER_SIZING_DEFAULTS)


def test_enum_choices_are_what_srm1d_accepts():
    meta = srm1d_plugin._IGNITER_CHAMBER_META
    assert tuple(meta['injection_topology']['choices']) == tuple(_VALID_TOPOLOGIES)
    for law in meta['burn_law']['choices']:
        _burn_law_code(law)                      # raises if srm_1d rejects it
    with pytest.raises(ValueError):
        _burn_law_code('not-a-law')
    codes = []
    for mode in srm1d_plugin._PYROGEN_FIELD_META['heat_delivery_mode']['choices']:
        codes.append(_heat_delivery_code(mode))  # raises if srm_1d rejects it
    assert len(set(codes)) == len(codes)
    with pytest.raises(ValueError):
        _heat_delivery_code('not-a-mode')


def test_schemas_follow_meta(solver):
    from motorlib.properties import (  # type: ignore
        FloatProperty, EnumProperty, StringProperty)
    for schema, meta in (
            (solver.get_igniter_material_schema(),
             srm1d_plugin._PYROGEN_FIELD_META),
            (solver.get_igniter_chamber_schema(),
             srm1d_plugin._IGNITER_CHAMBER_META)):
        assert list(schema.props) == list(meta)
        for key, m in meta.items():
            prop = schema.props[key]
            assert prop.dispName == m['label']
            if 'choices' in m:
                assert isinstance(prop, EnumProperty)
                assert list(prop.values) == list(m['choices'])
            elif 'min' in m:
                assert isinstance(prop, FloatProperty)
                assert (prop.unit, prop.min, prop.max) == (
                    m['unit'], m['min'], m['max'])
            else:
                assert isinstance(prop, StringProperty)


def test_schema_defaults_are_srm1d_defaults(solver):
    chamber = solver.get_igniter_chamber_schema().getProperties()
    assert chamber == _IGNITER_SIZING_DEFAULTS
    material = solver.get_igniter_material_schema().getProperties()
    for f in dataclasses.fields(Pyrogen):
        if f.init and f.default is not dataclasses.MISSING and f.default is not None:
            assert material[f.name] == f.default, f.name


# ----------------------------------------------------------------
# (b) range fit: nothing in the library may be dropped by FloatProperty
# ----------------------------------------------------------------

def test_every_library_material_fits_the_schema_exactly(solver):
    mats = library_materials()
    assert mats
    for mat in mats:
        schema = solver.get_igniter_material_schema()
        values = {k: v for k, v in mat.items() if v is not None}
        # Anything the schema does not carry must be a declared passthrough.
        assert set(values) - set(schema.props) <= set(
            srm1d_plugin._PYROGEN_GUI_PASSTHROUGH), mat['name']
        schema.setProperties(values)
        got = schema.getProperties()
        for key in values:
            if key in schema.props:
                assert got[key] == values[key], (mat['name'], key)


def test_default_igniter_chamber_values_fit_the_schema_exactly(solver):
    block = default_igniter_block('bpnv')
    chamber = {k: v for k, v in block.items() if k != 'pyrogen'}
    assert set(chamber) == set(CHAMBER_KEYS)
    schema = solver.get_igniter_chamber_schema()
    schema.setProperties(chamber)
    assert schema.getProperties() == chamber


# ----------------------------------------------------------------
# (c) hooks
# ----------------------------------------------------------------

def test_hooks_return_the_library_and_default(solver):
    assert solver.capabilities.get('igniter')
    mats = solver.get_igniter_materials()
    assert mats == library_materials()
    names = [m['name'] for m in mats]
    assert len(names) == len(set(names))
    assert {'BPNV', 'MTV'} <= set(names)
    # Fresh dicts every call: mutating one result cannot leak into the next.
    mats[0]['name'] = 'MUTATED'
    assert solver.get_igniter_materials()[0]['name'] != 'MUTATED'
    assert solver.get_default_igniter() == default_igniter_block('bpnv')


def test_schema_hooks_return_fresh_objects(solver):
    for hook in (solver.get_igniter_material_schema,
                 solver.get_igniter_chamber_schema):
        a, b = hook(), hook()
        assert a is not b
        key = next(iter(a.props))
        assert a.props[key] is not b.props[key]


def test_upgrade_hook_delegates_to_the_adapter(solver):
    assert solver.upgrade_igniter_block(_legacy_block()) == \
        upgrade_igniter_block(_legacy_block())


# ----------------------------------------------------------------
# (d) upgrade_igniter_block
# ----------------------------------------------------------------

def test_legacy_bpnv_block_is_upgraded_and_idempotent():
    chamber = {'mass': 0.0123, 'injection_topology': 'head_basket',
               'burn_law': 'end_burning', 'cartridge_length_m': 0.1}
    block = _legacy_block(**chamber)
    new, changed = upgrade_igniter_block(block)
    assert changed is True
    assert new['pyrogen'] == _current_bpnv()
    assert new['pyrogen']['gas_mass_fraction'] == 0.781
    for key in CHAMBER_KEYS:                       # sizing never touched
        assert new[key] == block[key]
    assert set(new) == set(block)
    again, changed2 = upgrade_igniter_block(new)   # idempotent
    assert changed2 is False
    assert again == new


def test_legacy_with_cp_gas_none_still_upgrades():
    block = _legacy_block()
    block['pyrogen']['Cp_gas'] = None
    new, changed = upgrade_igniter_block(block)
    assert changed is True
    assert new['pyrogen'] == _current_bpnv()


def test_name_match_is_case_insensitive():
    block = _legacy_block()
    block['pyrogen']['name'] = 'bpnv'
    new, changed = upgrade_igniter_block(block)
    assert changed is True
    assert new['pyrogen'] == _current_bpnv()


def test_float_tolerance_is_relative_1e9():
    near = _legacy_block()
    near['pyrogen']['M'] = 0.030 * (1.0 + 1e-12)
    assert upgrade_igniter_block(near)[1] is True
    far = _legacy_block()
    far['pyrogen']['M'] = 0.030 * (1.0 + 1e-6)
    assert upgrade_igniter_block(far)[1] is False


# ``name`` is the lookup key (case-insensitive); a different name is simply a
# different material, so it is excluded here.
@pytest.mark.parametrize(
    'key', sorted(k for k in SUPERSEDED_PYROGENS['BPNV'][0] if k != 'name'))
def test_any_altered_legacy_field_blocks_the_upgrade(key):
    block = _legacy_block()
    v = block['pyrogen'][key]
    block['pyrogen'][key] = v + '_x' if isinstance(v, str) else v * 1.01
    snapshot = copy.deepcopy(block)
    new, changed = upgrade_igniter_block(block)
    assert changed is False
    assert new == snapshot                         # untouched (customised)


def test_missing_legacy_field_blocks_the_upgrade():
    block = _legacy_block()
    del block['pyrogen']['kappa_jet']
    assert upgrade_igniter_block(block)[1] is False


def test_extra_non_none_key_blocks_the_upgrade():
    # gas_mass_fraction=1.0 is the Pyrogen default but is NOT absent/None:
    # a block that carries it explicitly is not the exact legacy snapshot.
    block = _legacy_block()
    block['pyrogen']['gas_mass_fraction'] = 1.0
    assert upgrade_igniter_block(block)[1] is False
    block = _legacy_block()
    block['pyrogen']['my_custom_key'] = 'x'
    assert upgrade_igniter_block(block)[1] is False


def test_other_materials_are_untouched():
    for name in ('mtv', 'mtv_fast', 'thermite_cual', 'bpnv'):
        block = default_igniter_block(name)
        new, changed = upgrade_igniter_block(block)
        assert changed is False, name
        assert new == block


def test_upgrade_does_not_mutate_its_input():
    block = _legacy_block()
    snapshot = copy.deepcopy(block)
    new, _ = upgrade_igniter_block(block)
    assert block == snapshot
    new['pyrogen']['name'] = 'MUTATED'
    new['mass'] = 99.0
    assert block == snapshot


def test_degenerate_blocks_pass_through():
    assert upgrade_igniter_block(None) == (None, False)
    assert upgrade_igniter_block({'mass': 1.0}) == ({'mass': 1.0}, False)
    assert upgrade_igniter_block({'pyrogen': None}) == ({'pyrogen': None}, False)
    assert upgrade_igniter_block({'pyrogen': {'a': 1}}) == (
        {'pyrogen': {'a': 1}}, False)


def test_superseded_registry_is_well_formed():
    lib = {m['name'].lower(): m for m in library_materials()}
    for name, snapshots in SUPERSEDED_PYROGENS.items():
        assert name.lower() in lib, f"{name} has no current library entry"
        assert snapshots
        for snap in snapshots:
            assert snap['name'].lower() == name.lower()
            # A snapshot equal to the CURRENT entry would "upgrade" forever.
            assert upgrade_igniter_block({'pyrogen': lib[name.lower()]})[1] is False
            assert snap != lib[name.lower()]
    for mat in library_materials():                # current entries are stable
        assert upgrade_igniter_block({'pyrogen': mat})[1] is False


# ----------------------------------------------------------------
# (e) load_igniter
# ----------------------------------------------------------------

def _record(block):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        out = load_igniter({'igniter': block})
    return out, caught


def test_load_igniter_upgrades_legacy_block_with_one_warning():
    block = _legacy_block(mass=0.004)
    snapshot = copy.deepcopy(block)
    (pyro, sizing), caught = _record(block)
    upgraded = [w for w in caught
                if issubclass(w.category, IgniterMaterialUpgradedWarning)]
    stale = [w for w in caught
             if issubclass(w.category, StaleIgniterMaterialWarning)]
    assert len(upgraded) == 1
    assert not stale                      # stays silent for an upgraded block
    msg = str(upgraded[0].message)
    assert 'BPNV' in msg and 'superseded' in msg and 'persists' in msg
    assert pyro.gas_mass_fraction == 0.781
    assert pyro.M == 0.0482
    assert sizing['pyrogen_mass'] == 0.004      # chamber sizing carried through
    assert block == snapshot                    # caller's dict not modified


def test_load_igniter_current_block_emits_neither_warning():
    _, caught = _record(default_igniter_block('bpnv'))
    assert not [w for w in caught
                if issubclass(w.category, (IgniterMaterialUpgradedWarning,
                                           StaleIgniterMaterialWarning))]


def test_load_igniter_customised_legacy_keeps_stale_warning_only():
    block = _legacy_block()
    block['pyrogen']['M'] = 0.031            # user tweaked the old snapshot
    (pyro, _), caught = _record(block)
    assert pyro.M == 0.031                   # embedded value still authoritative
    assert not [w for w in caught
                if issubclass(w.category, IgniterMaterialUpgradedWarning)]
    assert len([w for w in caught
                if issubclass(w.category, StaleIgniterMaterialWarning)]) == 1


def test_gui_unit_conversion_noise_is_still_recognised_and_upgraded():
    """The user's live GUI library (pyrogens.yaml) holds the legacy BPNV with
    float noise from the GUI's display-unit round trip (values copied here;
    the user's file is never read by code). rel tol 1e-9 must absorb it."""
    noisy = copy.deepcopy(SUPERSEDED_PYROGENS['BPNV'][0])
    noisy['rho'] = 1699.9999999999998
    noisy['particle_diameter_m'] = 0.0032000000000000006
    assert noisy != SUPERSEDED_PYROGENS['BPNV'][0]     # the noise is real
    new, changed = upgrade_igniter_block({'pyrogen': noisy})
    assert changed is True
    assert new['pyrogen'] == _current_bpnv()
