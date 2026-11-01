"""
test_d1_defaults.py -- D1 (2026-10-31): the EN1 / RT1 / HE1 options are the run_simulation defaults.

The pre-D1 behaviour stays reachable (and bit-for-bit) through the alternatives:
``energy_formulation='sensible'``, ``momentum_formulation='legacy'``,
``ignition_driver_temperature='static'``, ``head_end_jet_h='off'``.
"""
import inspect
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d.simulation import run_simulation  # noqa: E402

RIC = os.path.join(os.path.dirname(__file__), "..", "motors", "machbusterNew.ric")
PRE_D1 = dict(energy_formulation='sensible', momentum_formulation='legacy',
              ignition_driver_temperature='static', head_end_jet_h='off')


def test_signature_defaults():
    p = inspect.signature(run_simulation).parameters
    assert p['energy_formulation'].default == 'total'
    assert p['momentum_formulation'].default == 'conservative'
    assert p['ignition_driver_temperature'].default == 'recovery'
    assert p['head_end_jet_h'].default == 'max'


@pytest.fixture(scope="module")
def runs():
    from srm_1d.openmotor_adapter import run_from_ric
    kw = dict(target_propellant_cells=60, t_max=0.012, roughness=32e-6, kappa=0.44,
              T_ignition=756.0, print_interval=0.0, verbose=False, pyrogen='bpnv',
              pyrogen_mass=2.623e-3, history_capacity=100000)
    return {'default': run_from_ric(RIC, **kw)[0],
            'explicit_d1': run_from_ric(
                RIC, energy_formulation='total', momentum_formulation='conservative',
                ignition_driver_temperature='recovery', head_end_jet_h='max', **kw)[0],
            'pre_d1': run_from_ric(RIC, **PRE_D1, **kw)[0]}


def test_default_is_exactly_the_explicit_d1_options(runs):
    a, b = runs['default'], runs['explicit_d1']
    assert np.array_equal(np.asarray(a['P_head']), np.asarray(b['P_head']))
    assert np.array_equal(np.asarray(a['ignition_time_by_cell']),
                          np.asarray(b['ignition_time_by_cell']))


def test_pre_d1_options_differ_and_are_stamped(runs):
    assert not np.array_equal(np.asarray(runs['default']['P_head']),
                              np.asarray(runs['pre_d1']['P_head']))
    s = runs['pre_d1']['summary']
    assert (s['energy_formulation'], s['momentum_formulation'],
            s['ignition_driver_temperature'], s['head_end_jet_h']) == (
        'sensible', 'legacy', 'static', 'off')
    s = runs['default']['summary']
    assert (s['energy_formulation'], s['momentum_formulation'],
            s['ignition_driver_temperature'], s['head_end_jet_h']) == (
        'total', 'conservative', 'recovery', 'max')


def test_energy_convention_follows_the_formulation(runs):
    assert runs['default']['summary']['energy_convention'] == 'total_energy_CvT_plus_KE'
    assert runs['pre_d1']['summary']['energy_convention'] == 'sensible_enthalpy_CpT'
    assert 'kinetic_energy' not in runs['default']['summary']['energy_omitted_terms']
    assert 'kinetic_energy' in runs['pre_d1']['summary']['energy_omitted_terms']
