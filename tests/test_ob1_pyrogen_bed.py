"""
tests/test_ob1_pyrogen_bed.py -- OB1 opt-in pyrogen burn-model options
======================================================================

* ``Pyrogen.P_plateau`` : the burn RATE is evaluated at min(P, P_plateau)
  (Zhu et al. 2008 B/KNO3 closed-bomb law: r = r0*P up to 7 MPa, flat above).
* ``PyrogenChamber.basket_area_law = 'shrinking'`` : per-cell mass bookkeeping,
  burning area A = A0 (m/m0)^(2/3) (exact for cubes / spheres / L/D = 1
  cylinders under uniform face regression).
* ``PyrogenChamber.bed_flame_speed_m_s`` : diagnostic sequential lighting.

Defaults must be unchanged: every option is OFF unless asked for.
Hand-checked formulas, default identity, and pyrogen mass conservation.
"""
import dataclasses
import math

import numpy as np
import pytest

from srm_1d.igniter_plenum import (
    PyrogenChamber, _plenum_rates, chamber_params, pyrogen_params,
)
from srm_1d.openmotor_adapter import load_pyrogen
from srm_1d.propellant import Pyrogen
from srm_1d.simulation import (
    _compute_uncontained_pyrogen_mdot, _compute_uncontained_pyrogen_mdot_cells,
)

A_SR, N_SR, RHO = 2.0e-5, 0.5, 1700.0
# Zhu et al. 2008 Table 2 (m/s/MPa): 6.848e-3, 6.787e-3, 5.389e-3, 5.271e-3
ZHU_R0_MEAN = (6.848e-3 + 6.787e-3 + 5.389e-3 + 5.271e-3) / 4.0


def _pyro(**kw):
    d = dict(name='t', a=A_SR, n=N_SR, rho=RHO, T_flame=2700.0, M=0.048, gamma=1.17)
    d.update(kw)
    return Pyrogen(**d)


# ---------------------------------------------------------------- defaults
def test_defaults_are_off():
    assert _pyro().P_plateau == 0.0
    for name in ('bpnv', 'mtv'):
        assert load_pyrogen(name).P_plateau == 0.0
    ch = PyrogenChamber(pyrogen=_pyro(), m_pyrogen_initial=1e-3, A_burn_initial=1e-3,
                        A_throat=1e-5, V_plenum=1e-5)
    assert ch.basket_area_law == 'constant'
    assert ch.bed_flame_speed_m_s == 0.0
    pp, cp = pyrogen_params(ch.pyrogen), chamber_params(ch)
    assert pp.shape == (8,) and pp[7] == 0.0
    assert cp.shape == (7,) and cp[5] == 0.0 and cp[6] == 0.0


def test_validation():
    base = dict(pyrogen=_pyro(), m_pyrogen_initial=1e-3, A_burn_initial=1e-3, A_throat=1e-5, V_plenum=1e-5)
    with pytest.raises(ValueError):
        PyrogenChamber(**base, basket_area_law='cubic')
    with pytest.raises(ValueError):
        PyrogenChamber(**base, bed_flame_speed_m_s=-1.0)
    with pytest.raises(ValueError):                      # sequential lighting needs the per-cell path
        PyrogenChamber(**base, bed_flame_speed_m_s=5.0)
    PyrogenChamber(**base, basket_area_law='shrinking', bed_flame_speed_m_s=5.0)
    with pytest.raises(ValueError):
        PyrogenChamber(**{**base, 'pyrogen': _pyro(P_plateau=-1.0)})


# ---------------------------------------------------------------- plateau
def test_plateau_off_is_identical_and_on_caps_the_rate():
    N = 3
    P = np.array([1.0e6, 7.0e6, 20.0e6])
    A = 1.0e-3
    m_off, m_zero, m_on = np.zeros(N), np.zeros(N), np.zeros(N)
    _compute_uncontained_pyrogen_mdot(P, A_SR, N_SR, RHO, A, 1.0, 1e-6, 0, 2, N, m_off)
    _compute_uncontained_pyrogen_mdot(P, A_SR, N_SR, RHO, A, 1.0, 1e-6, 0, 2, N, m_zero, 1.0, 0.0)
    assert np.array_equal(m_off, m_zero)                       # P_plateau = 0 -> unchanged
    _compute_uncontained_pyrogen_mdot(P, A_SR, N_SR, RHO, A, 1.0, 1e-6, 0, 2, N, m_on, 1.0, 7.0e6)
    # hand check: rate = a * min(P, 7e6)^0.5 ; mdot = rho * rate * A
    expect = RHO * A * A_SR * np.sqrt(np.minimum(P, 7.0e6))
    np.testing.assert_allclose(m_on, expect, rtol=1e-12)
    assert m_on[2] == pytest.approx(m_on[1])                    # flat above the plateau
    assert m_on[0] == pytest.approx(m_off[0])                   # below it: unchanged


def test_zhu_law_in_the_a_n_plateau_parameterisation():
    """r = r0*P (2.8-7 MPa) and 7*r0 (7-50 MPa), r0 in m/s/MPa  ==  a = r0/1e6, n = 1, P_plateau = 7e6."""
    a, n, Pc = ZHU_R0_MEAN / 1e6, 1.0, 7.0e6
    N = 3
    P = np.array([2.8e6, 7.0e6, 50.0e6])
    mdot = np.zeros(N)
    _compute_uncontained_pyrogen_mdot(P, a, n, RHO, 1.0, 1.0, 1e-6, 0, 2, N, mdot, 1.0, Pc)
    rates = mdot / RHO                                           # A = 1 m^2
    np.testing.assert_allclose(rates, [ZHU_R0_MEAN * 2.8, ZHU_R0_MEAN * 7.0, ZHU_R0_MEAN * 7.0], rtol=1e-12)
    assert rates[1] == pytest.approx(0.0425, rel=2e-3)          # 42.5 mm/s plateau (mean of Zhu Table 2)


def test_plenum_plateau_hand_check():
    pyro = _pyro(P_plateau=7.0e6)
    ch = PyrogenChamber(pyrogen=pyro, m_pyrogen_initial=1e-3, A_burn_initial=2e-4, A_throat=1e-6, V_plenum=1e-5)
    # state where P_ig = 30 MPa: m_gas = P V M / (R T)
    from srm_1d.propellant import R_UNIVERSAL
    T = 2700.0
    m_gas = 30.0e6 * ch.V_plenum * pyro.M / (R_UNIVERSAL * T)
    args = (ch.m_pyrogen_initial, m_gas, T, pyro.a, pyro.n, pyro.rho, pyro.T_flame, pyro.M, pyro.gamma,
            ch.m_pyrogen_initial, ch.A_burn_initial, ch.A_throat, ch.V_plenum, 0.0, 101325.0)
    gen_plateau, _, P_ig = _plenum_rates(*args, 7.0e6)
    gen_off, _, _ = _plenum_rates(*args)
    assert P_ig == pytest.approx(30.0e6, rel=1e-9)
    assert gen_plateau == pytest.approx(RHO * ch.A_burn_initial * A_SR * math.sqrt(7.0e6), rel=1e-12)
    assert gen_off == pytest.approx(RHO * ch.A_burn_initial * A_SR * math.sqrt(30.0e6), rel=1e-12)


# ---------------------------------------------------------------- per-cell area law
def _cells_step(P, m_cell, m0, t, dt, lo, hi, N, mdot, gmf=1.0, Pc=0.0, t_light=None, A0=1.0e-3):
    if t_light is None:
        t_light = np.zeros(N)
    return _compute_uncontained_pyrogen_mdot_cells(
        P, A_SR, N_SR, RHO, A0, m_cell, m0, t, dt, lo, hi, N, mdot, gmf, Pc, t_light)


def test_area_follows_two_thirds_power_hand_check():
    N, m0 = 4, 1.0e-3
    P = np.full(N, 4.0e6)
    for frac in (1.0, 0.5, 0.125):
        m_cell = np.zeros(N)
        m_cell[:] = frac * m0
        mdot = np.zeros(N)
        _cells_step(P, m_cell, m0, 0.0, 1.0e-9, 0, 3, N, mdot)
        expect = RHO * A_SR * math.sqrt(4.0e6) * 1.0e-3 * frac ** (2.0 / 3.0)
        np.testing.assert_allclose(mdot, expect, rtol=1e-12)
    # frac = 0.125 -> (1/8)^(2/3) = 1/4 exactly
    assert mdot[0] == pytest.approx(0.25 * RHO * A_SR * math.sqrt(4.0e6) * 1.0e-3, rel=1e-12)


def test_cell_depletes_exactly_and_conserves_mass_with_gas_fraction():
    N, m0, gmf, dt = 3, 2.0e-5, 0.781, 2.0e-5
    P = np.array([1.0e6, 3.0e6, 9.0e6])        # cells burn out at different times (each on its own host-cell P)
    m_cell = np.full(N, m0)
    mdot = np.zeros(N)
    burned_gas = 0.0
    t = 0.0
    done_time = np.full(N, -1.0)
    for step in range(400000):
        rem = _cells_step(P, m_cell, m0, t, dt, 0, 2, N, mdot, gmf, A0=1.0e-4)
        burned_gas += mdot.sum() * dt
        t += dt
        for i in range(N):
            if m_cell[i] == 0.0 and done_time[i] < 0:
                done_time[i] = t
        if rem <= 0.0:
            break
    assert rem == 0.0 and np.all(m_cell == 0.0)
    # gas delivered = gas_mass_fraction * solid consumed (3 cells * m0), to round-off
    assert burned_gas == pytest.approx(gmf * N * m0, rel=1e-12)
    assert done_time[2] < done_time[1] < done_time[0]          # higher host-cell P burns out first
    # after depletion nothing is produced
    _cells_step(P, m_cell, m0, t, dt, 0, 2, N, mdot, gmf, A0=1.0e-4)
    assert np.all(mdot == 0.0)


def test_shrinking_cube_burn_time_matches_the_analytic_web_time():
    """One cube (A = 6 s^2, m = rho s^3) at constant pressure regresses uniformly: web time = (s/2)/r."""
    s0 = 3.0e-3
    m0 = RHO * s0 ** 3
    A0 = 6.0 * s0 ** 2
    P0 = 5.0e6
    r = A_SR * math.sqrt(P0)
    N, dt = 1, 2.0e-6
    m_cell = np.array([m0])
    mdot = np.zeros(N)
    P = np.array([P0])
    t = 0.0
    while m_cell[0] > 0.0:
        _cells_step(P, m_cell, m0, t, dt, 0, 0, N, mdot, 1.0, A0=A0)
        t += dt
        assert t < 1.0
    assert t == pytest.approx((s0 / 2.0) / r, rel=2e-3)      # Euler step error << 0.2 %


def test_bed_flame_speed_delays_unlit_cells():
    N, m0 = 4, 1.0e-3
    P = np.full(N, 2.0e6)
    t_light = np.array([0.0, 1.0e-3, 2.0e-3, 3.0e-3])
    m_cell = np.full(N, m0)
    mdot = np.zeros(N)
    _cells_step(P, m_cell, m0, 1.5e-3, 1e-6, 0, 3, N, mdot, t_light=t_light)
    assert mdot[0] > 0 and mdot[1] > 0 and mdot[2] == 0.0 and mdot[3] == 0.0
    assert m_cell[2] == m0 and m_cell[3] == m0                 # untouched until lit


# ---------------------------------------------------------------- run level
@pytest.fixture(scope='module')
def _runs():
    from srm_1d.openmotor_adapter import run_from_ric
    from pathlib import Path
    ric = str(Path(__file__).resolve().parents[1] / 'motors' / 'machbusterNew.ric')
    pyro = dataclasses.replace(load_pyrogen('bpnv'), particle_diameter_m=3e-3, particle_LD_ratio=1.0)
    kw = dict(pyrogen=pyro, pyrogen_mass=4.5e-3, injection_topology='head_basket', roughness=32e-6, kappa=0.44,
              T_ignition=756.0, print_interval=0.0, verbose=False, t_max=0.03, snapshot_interval=1e-3)
    out = {}
    out['default'] = run_from_ric(ric, **kw)[0]
    out['explicit_const'] = run_from_ric(ric, basket_area_law='constant', **kw)[0]
    out['shrinking'] = run_from_ric(ric, basket_area_law='shrinking', **kw)[0]
    return out


def test_run_default_identity(_runs):
    a, b = _runs['default'], _runs['explicit_const']
    assert np.array_equal(np.asarray(a['P_head']), np.asarray(b['P_head']))
    assert a['summary']['steps'] == b['summary']['steps']


def test_run_shrinking_conserves_pyrogen_mass_and_changes_the_burn(_runs):
    s = _runs['shrinking']['summary']
    assert s['pyrogen_mass_burned'] == pytest.approx(4.5e-3, rel=1e-9)
    mp = np.asarray(_runs['shrinking']['m_pyrogen'])
    assert mp[0] > 0.99 * 4.5e-3 and mp[-1] == 0.0 and np.all(np.diff(mp) <= 1e-18)   # monotone depletion
    # the area law lengthens the burn compared with the constant-area default
    assert s['pyrogen_duration'] > 1.5 * _runs['default']['summary']['pyrogen_duration']
