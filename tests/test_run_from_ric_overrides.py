"""run_from_ric igniter-sizing precedence: explicit kwarg > .ric block > defaults.

Regression for a silent bug (found 2026-10-13): with a motor whose .ric
igniter block says ``mass: -1`` (auto), ``run_from_ric(..., pyrogen_mass=X)``
ignored X and auto-sized the charge. Every explicit igniter-sizing kwarg was
dropped whenever the motor file carried an igniter block.
"""
import os
import warnings

import pytest

from srm_1d.openmotor_adapter import _merge_igniter_sizing, run_from_ric

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHUNC = os.path.join(ROOT, "motors", "machbusterNew.ric")


def test_merge_explicit_beats_block_beats_defaults():
    block = {'pyrogen_mass': None, 'injection_topology': 'forward_plenum',
             'basket_fill_fraction': 0.5}
    merged = _merge_igniter_sizing(block, pyrogen_mass=2.0e-3,
                                   injection_topology=None,
                                   basket_area_law='constant')
    assert merged['pyrogen_mass'] == 2.0e-3                 # explicit wins
    assert merged['injection_topology'] == 'forward_plenum'  # block kept
    assert merged['basket_fill_fraction'] == 0.5
    assert merged['basket_area_law'] == 'constant'          # explicit added
    assert block['pyrogen_mass'] is None                    # input not mutated


def test_merge_without_block_passes_only_explicit_values():
    merged = _merge_igniter_sizing(None, pyrogen_mass=None, cartridge_length_m=None,
                                   injection_topology='head_basket')
    assert merged == {'injection_topology': 'head_basket'}


def test_explicit_charge_overrides_the_ric_igniter_block():
    # Whatever the motor file's block says (auto in the committed file, an
    # explicit charge in a user's edited copy), the kwarg must win.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_from_ric(CHUNC, pyrogen_mass=2.0e-3, t_max=0.002,
                           print_interval=0.0, verbose=False)[0]
    s = res['summary']
    assert not s['pyrogen_mass_auto']
    assert s['ignition_attempts'] == 1
    assert s['pyrogen_mass_initial'] == pytest.approx(2.0e-3)
