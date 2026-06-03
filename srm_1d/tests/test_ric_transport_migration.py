"""
test_ric_transport_migration.py — v0.8.0 Phase 3 .ric transport migration.

Verifies the version-locked per-tab transport migration: folding a
.transport.yaml sidecar (+ frozen sibling) into a .ric, sharing transport
across motors with the same propellant, preserving openMotor's YAML tags,
the frozen-default variant, and the D7 hard-fault on missing transport.
"""

import os
import shutil

import pytest

yaml = pytest.importorskip("yaml")

from srm_1d.openmotor_adapter import (  # noqa: E402
    migrate_ric_transport, build_transport_library, load_ric,
    convert_propellant, _ric_raw_loader, RIC_FORMAT_VERSION,
    _TRANSPORT_TAB_KEYS,
)

MOTORS = os.path.join(os.path.dirname(__file__), '..', 'motors')


def _write_old_ric(path, pname='Test Prop'):
    """A minimal pre-0.7.0 .ric (no transport keys) with the python tags."""
    text = (
        "data:\n"
        "  config: {ambPressure: 101325.0}\n"
        "  grains:\n"
        "  - properties: {coreDiameter: 0.02, diameter: 0.05, "
        "inhibitedEnds: Both, length: 0.3}\n"
        "    type: BATES\n"
        "  nozzle: {throat: 0.01, exit: 0.02, efficiency: 0.95, "
        "divAngle: 15.0, convAngle: 60.0, slagCoeff: 0.0, "
        "erosionCoeff: 0.0, throatLength: 0.0}\n"
        "  propellant:\n"
        f"    density: 1700.0\n    name: {pname}\n"
        "    tabs:\n"
        "    - {a: 1.0e-05, n: 0.4, k: 1.2, t: 2800.0, m: 25.0, "
        "minPressure: 0.0, maxPressure: 10000000.0}\n"
        "type: !!python/object/apply:uilib.fileIO.fileTypes\n- 3\n"
        "version: !!python/tuple\n- 0\n- 6\n- 1\n"
    )
    with open(path, 'w', encoding='utf-8') as f:
        f.write(text)


def _write_sidecar(path, mu, k, cp):
    with open(path, 'w', encoding='utf-8') as f:
        yaml.safe_dump({'mu': mu, 'k': k, 'Cp': cp}, f)


def test_migrate_folds_sidecars_and_preserves_tags(tmp_path):
    ric = str(tmp_path / 'm.ric')
    _write_old_ric(ric)
    _write_sidecar(str(tmp_path / 'm.transport.yaml'), 9e-5, 0.65, 2700)   # effective
    _write_sidecar(str(tmp_path / 'm.frozen.transport.yaml'), 9e-5, 0.37, 2050)  # frozen

    migrate_ric_transport(ric, variant='frozen', write=True)

    # Tags preserved (file stays openMotor-loadable) + version bumped.
    text = open(ric, encoding='utf-8').read()
    assert '!!python/object/apply:uilib.fileIO.fileTypes' in text
    assert '!!python/tuple' in text
    raw = yaml.load(open(ric, encoding='utf-8'), Loader=_ric_raw_loader())
    assert tuple(raw['version']) == RIC_FORMAT_VERSION

    motor = load_ric(ric)
    tab = motor['propellant']['tabs'][0]
    assert motor['propellant']['transportVariant'] == 'frozen'
    assert tab['mu'] == pytest.approx(9e-5)
    assert tab['kThermalFrozen'] == pytest.approx(0.37)
    assert tab['cpFrozen'] == pytest.approx(2050)
    assert tab['kThermalEffective'] == pytest.approx(0.65)
    assert tab['cpEffective'] == pytest.approx(2700)

    # Frozen-default solver inputs come from the frozen slot.
    prop = convert_propellant(motor['propellant'])
    assert prop.k_gas == pytest.approx(0.37)
    assert prop.Cp_gas == pytest.approx(2050)


def test_lone_primary_fills_both_slots(tmp_path):
    """A motor with only a primary sidecar: both slots get that value, so
    the frozen default still runs (behavior-preserving)."""
    ric = str(tmp_path / 'm.ric')
    _write_old_ric(ric)
    _write_sidecar(str(tmp_path / 'm.transport.yaml'), 8e-5, 0.6, 2600)

    migrate_ric_transport(ric, variant='frozen', write=True)
    tab = load_ric(ric)['propellant']['tabs'][0]
    assert tab['kThermalFrozen'] == pytest.approx(0.6)
    assert tab['kThermalEffective'] == pytest.approx(0.6)


def test_no_sidecar_keeps_sentinel_and_hard_faults(tmp_path):
    ric = str(tmp_path / 'm.ric')
    _write_old_ric(ric)
    migrate_ric_transport(ric, transport=None, variant='frozen', write=True)
    tab = load_ric(ric)['propellant']['tabs'][0]
    for key in _TRANSPORT_TAB_KEYS:
        assert tab[key] == 0.0  # sentinel
    with pytest.raises(ValueError, match="transport"):
        convert_propellant(load_ric(ric)['propellant'])


def test_shared_propellant_fills_sidecarless_motor(tmp_path):
    """Two motors with the same propellant; only one has a sidecar. The
    library shares it so the other gets transport too."""
    a = str(tmp_path / 'a.ric')
    b = str(tmp_path / 'b.ric')
    _write_old_ric(a, pname='Shared Prop')
    _write_old_ric(b, pname='Shared Prop')
    _write_sidecar(str(tmp_path / 'a.transport.yaml'), 9e-5, 0.65, 2700)
    _write_sidecar(str(tmp_path / 'a.frozen.transport.yaml'), 9e-5, 0.37, 2050)

    lib = build_transport_library(str(tmp_path))
    assert 'Shared Prop' in lib and lib['Shared Prop'] is not None
    # b has no sidecar but inherits via the library
    migrate_ric_transport(b, transport=lib['Shared Prop'], variant='frozen')
    tab = load_ric(b)['propellant']['tabs'][0]
    assert tab['kThermalFrozen'] == pytest.approx(0.37)
    assert tab['kThermalEffective'] == pytest.approx(0.65)


def test_repo_motors_all_migrated_and_loadable():
    """The real repo motors were bulk-migrated: every one is at the current
    format version, has non-sentinel transport, and loads without fault."""
    import glob
    rics = glob.glob(os.path.join(MOTORS, '*.ric'))
    assert rics
    for ric in rics:
        motor = load_ric(ric)
        prop = motor['propellant']
        assert prop.get('transportVariant') in ('frozen', 'effective')
        # Reads transport from the .ric (no sidecar) without hard-faulting.
        converted = convert_propellant(prop)
        assert converted.mu_gas > 0 and converted.k_gas > 0 and converted.Cp_gas > 0
