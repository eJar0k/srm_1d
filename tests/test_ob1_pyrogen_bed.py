"""
tests/test_ob1_pyrogen_bed.py -- OB1 opt-in pyrogen burn-model options
======================================================================

* ``Pyrogen.P_plateau`` : the burn RATE is evaluated at min(P, P_plateau)
  (Zhu et al. 2008 B/KNO3 closed-bomb law: r = r0*P up to 7 MPa, flat above).
* ``PyrogenChamber.basket_area_law`` : 'shrinking' (DEFAULT) = per-cell mass
  bookkeeping with the exact right-cylinder form function (all faces burning,
  D = D0-2s, H = H0-2s; A0 (m/m0)^(2/3) for L/D = 1); 'constant' = the pre-OB1
  pooled constant-area law (opt-in).
* ``PyrogenChamber.bed_flame_speed_m_s`` : diagnostic sequential lighting.
* the forward plenum's '0d' law uses the same exact cylinder form function.

``P_plateau`` and the bed speed are OFF by default. Hand-checked formulas,
shape-function checks, mass conservation, the library pellet geometry and the
superseded-snapshot upgrade.
"""
import dataclasses
import math

import numpy as np
import pytest

from srm_1d.igniter_plenum import (
    PyrogenChamber, _burn_area, _cyl_area_fraction, _cyl_area_fraction_general,
    _plenum_rates, chamber_params, pyrogen_params,
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
    assert ch.basket_area_law == 'shrinking'          # OB1 default (exact cylinder form function)
    assert ch.bed_flame_speed_m_s == 0.0
    pp, cp = pyrogen_params(ch.pyrogen), chamber_params(ch)
    assert pp.shape == (8,) and pp[7] == 0.0
    assert cp.shape == (8,) and cp[5] == 1.0 and cp[6] == 0.0 and cp[7] == ch.pyrogen.particle_LD_ratio
    ch2 = PyrogenChamber(pyrogen=_pyro(), m_pyrogen_initial=1e-3, A_burn_initial=1e-3,
                         A_throat=1e-5, V_plenum=1e-5, basket_area_law='constant')
    assert chamber_params(ch2)[5] == 0.0


def test_validation():
    base = dict(pyrogen=_pyro(), m_pyrogen_initial=1e-3, A_burn_initial=1e-3, A_throat=1e-5, V_plenum=1e-5)
    with pytest.raises(ValueError):
        PyrogenChamber(**base, basket_area_law='cubic')
    with pytest.raises(ValueError):
        PyrogenChamber(**base, bed_flame_speed_m_s=-1.0)
    with pytest.raises(ValueError):                      # sequential lighting needs the per-cell path
        PyrogenChamber(**base, basket_area_law='constant', bed_flame_speed_m_s=5.0)
    with pytest.raises(ValueError):
        PyrogenChamber(**{**base, 'pyrogen': _pyro(particle_LD_ratio=0.0)})
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
    out['constant'] = run_from_ric(ric, basket_area_law='constant', **kw)[0]
    out['shrinking'] = run_from_ric(ric, basket_area_law='shrinking', **kw)[0]
    return out


def test_run_default_is_the_shrinking_law(_runs):
    a, b = _runs['default'], _runs['shrinking']
    assert np.array_equal(np.asarray(a['P_head']), np.asarray(b['P_head']))
    assert a['summary']['steps'] == b['summary']['steps']


def test_run_shrinking_conserves_pyrogen_mass_and_changes_the_burn(_runs):
    s = _runs['shrinking']['summary']
    assert s['pyrogen_mass_burned'] == pytest.approx(4.5e-3, rel=1e-9)
    mp = np.asarray(_runs['shrinking']['m_pyrogen'])
    assert mp[0] > 0.99 * 4.5e-3 and mp[-1] == 0.0 and np.all(np.diff(mp) <= 1e-18)   # monotone depletion
    # the area law lengthens the burn compared with the opt-in constant-area law
    assert s['pyrogen_duration'] > 1.5 * _runs['constant']['summary']['pyrogen_duration']
    c = _runs['constant']['summary']
    assert c['pyrogen_mass_burned'] == pytest.approx(4.5e-3, rel=1e-9)


# ---------------------------------------------------------------- data round trips
def test_p_plateau_round_trips_through_yaml_and_igniter_block(tmp_path):
    import yaml
    from srm_1d.openmotor_adapter import (
        _block_to_pyrogen, _pyrogen_to_block, default_igniter_block, load_igniter,
    )
    p = _pyro(P_plateau=7.0e6)
    assert _block_to_pyrogen(_pyrogen_to_block(p)).P_plateau == 7.0e6
    f = tmp_path / 'zhu.yaml'
    f.write_text(yaml.safe_dump({'name': 'Z', 'a': 6.07375e-9, 'n': 1.0, 'rho': 1700.0, 'T_flame': 2719.0,
                                 'M': 0.0482, 'gamma': 1.17, 'P_plateau': 7.0e6}))
    q = load_pyrogen(str(f))
    assert q.P_plateau == 7.0e6 and q.n == 1.0
    # the motor-file igniter block carries the chamber option through to the sizing kwargs
    block = default_igniter_block('bpnv')
    assert block['basket_area_law'] == 'shrinking'
    block['basket_area_law'] = 'constant'
    _, sizing = load_igniter({'igniter': block})
    assert sizing['basket_area_law'] == 'constant'


# ---------------------------------------------------------------- exact cylinder form function
def _cyl_geometry(s_depth, lam):
    """Hand model (D0 = 1, H0 = lam): returns (V/V0, A/A0) at regression depth s."""
    D, H = 1.0 - 2.0 * s_depth, lam - 2.0 * s_depth
    V = math.pi * D * D * H / 4.0
    A = math.pi * D * H + math.pi * D * D / 2.0
    V0 = math.pi * lam / 4.0
    A0 = math.pi * lam + math.pi / 2.0
    return V / V0, A / A0


@pytest.mark.parametrize('lam', [0.4, 1.0, 1.504, 3.0])
def test_cylinder_form_function_matches_the_hand_geometry(lam):
    smax = 0.5 * min(1.0, lam)
    for frac_s in (0.0, 0.1, 0.3, 0.5, 0.8, 0.95):
        sd = frac_s * smax
        f, g = _cyl_geometry(sd, lam)
        assert _cyl_area_fraction_general(f, lam) == pytest.approx(g, rel=1e-9, abs=1e-12)


def test_cylinder_form_function_endpoints_and_burnout():
    for lam in (0.4, 1.0, 1.504, 3.0):
        assert _cyl_area_fraction(1.0, lam) == 1.0
        assert _cyl_area_fraction(0.0, lam) == 0.0
        # A decreases monotonically with the remaining mass (L/D >= 1; for L/D < 1 it tends to the disc area)
        fs = np.linspace(0.02, 1.0, 50)
        gs = np.array([_cyl_area_fraction(f, lam) for f in fs])
        assert np.all(np.diff(gs) >= -1e-12) or lam < 1.0
    # L/D = 1.504: burn-out at s = D0/2 (D -> 0 first); the area vanishes with the volume
    assert _cyl_area_fraction(1e-9, 1.504) < 1e-3
    # L/D < 1: H -> 0 first; the last remnant is a zero-thickness disc of area pi D^2/2
    lam = 0.4
    D = 1.0 - 2.0 * (lam / 2.0)
    disc = (math.pi * D * D / 2.0) / (math.pi * lam + math.pi / 2.0)
    assert _cyl_area_fraction(1e-12, lam) == pytest.approx(disc, rel=1e-3)


def test_cylinder_form_function_reduces_to_the_two_thirds_power_at_LD_1():
    for f in (0.999, 0.9, 0.5, 0.125, 0.01, 1e-4):
        assert _cyl_area_fraction_general(f, 1.0) == pytest.approx(f ** (2.0 / 3.0), rel=1e-9)
        assert _cyl_area_fraction(f, 1.0) == f ** (2.0 / 3.0)


def test_plenum_burn_area_uses_the_exact_shape_function():
    A0 = 3.0e-3
    # L/D = 1: identical to the old A0 (m/m0)^(2/3)
    for f in (0.9, 0.4, 0.05):
        assert _burn_area(f, 1.0, A0, 0, 1.0) == pytest.approx(A0 * f ** (2.0 / 3.0), rel=1e-12)
    # L/D = 1.504: the hand geometry
    lam = 1.504
    sd = 0.2
    f, g = _cyl_geometry(sd, lam)
    assert _burn_area(f, 1.0, A0, 0, lam) == pytest.approx(A0 * g, rel=1e-9)
    assert _burn_area(f, 1.0, A0, 1, lam) == A0          # end_burning unchanged
    assert _burn_area(0.0, 1.0, A0, 0, lam) == 0.0


def test_basket_cells_follow_the_shape_function_at_LD_1504():
    N, m0 = 2, 1.0e-3
    P = np.full(N, 3.0e6)
    lam = 1.504
    sd = 0.25
    f, g = _cyl_geometry(sd, lam)
    m_cell = np.full(N, f * m0)
    mdot = np.zeros(N)
    _compute_uncontained_pyrogen_mdot_cells(P, A_SR, N_SR, RHO, 1.0e-3, m_cell, m0, 0.0, 1e-9, 0, 1, N, mdot,
                                            1.0, 0.0, np.zeros(N), lam)
    assert mdot[0] == pytest.approx(RHO * A_SR * math.sqrt(3.0e6) * 1.0e-3 * g, rel=1e-9)


def test_single_cylinder_burn_time_matches_the_analytic_web_time():
    """L/D = 1.504 cylinder at constant P: D reaches 0 at s = D0/2, so t = (D0/2)/r (mass -> 0 there)."""
    d0, lam = 3.175e-3, 1.504
    h0 = lam * d0
    V0 = math.pi * d0 * d0 * h0 / 4.0
    m0 = RHO * V0
    A0 = math.pi * d0 * h0 + math.pi * d0 * d0 / 2.0
    P0 = 5.0e6
    r = A_SR * math.sqrt(P0)
    m_cell = np.array([m0])
    mdot = np.zeros(1)
    P = np.array([P0])
    t, dt = 0.0, 2.0e-6
    while m_cell[0] > 0.0:
        _compute_uncontained_pyrogen_mdot_cells(P, A_SR, N_SR, RHO, A0, m_cell, m0, t, dt, 0, 0, 1, mdot, 1.0, 0.0,
                                                np.zeros(1), lam)
        t += dt
        assert t < 1.0
    assert t == pytest.approx((d0 / 2.0) / r, rel=5e-3)


# ---------------------------------------------------------------- initial area and the library geometry
def test_initial_burn_area_is_the_cylinder_formula_for_any_LD():
    from srm_1d.openmotor_adapter import build_pyrogen_chamber
    from srm_1d.grain_geometry import build_snapped_geometry   # noqa: F401  (import check only)
    import types
    geo = types.SimpleNamespace(D_outer=0.05, L_motor=0.3, total_propellant_volume=lambda: 1e-4)
    nozzle = None
    for lam in (0.5, 1.0, 1.504, 3.0):
        pyro = _pyro(particle_diameter_m=3.0e-3, particle_LD_ratio=lam)
        ch = build_pyrogen_chamber(pyro, geo, nozzle, pyrogen_mass=4.5e-3, pyrogen_throat_area=1e-5,
                                   pyrogen_volume=1e-5)
        d, h = 3.0e-3, lam * 3.0e-3
        A_p = math.pi * d * h + math.pi * d * d / 2.0
        V_p = math.pi * d * d * h / 4.0
        assert ch.A_burn_initial == pytest.approx(4.5e-3 / (RHO * V_p) * A_p, rel=1e-12)
    # the former sphere branch (6 m/(rho d)) was wrong for L/D < 1; at L/D = 1 both agree
    ch1 = build_pyrogen_chamber(_pyro(particle_diameter_m=3e-3, particle_LD_ratio=1.0), geo, nozzle,
                                pyrogen_mass=4.5e-3, pyrogen_throat_area=1e-5, pyrogen_volume=1e-5)
    assert ch1.A_burn_initial == pytest.approx(6.0 * 4.5e-3 / (RHO * 3e-3), rel=1e-12)


def test_library_pellet_geometry_is_the_mil_p_46994b_type_I_A_pellet():
    for name in ('bpnv', 'mtv'):
        p = load_pyrogen(name)
        assert p.particle_diameter_m == 3.175e-3            # 0.125 in
        assert p.particle_LD_ratio == 1.504                 # 0.188 / 0.125
        assert p.particle_LD_ratio * p.particle_diameter_m == pytest.approx(0.188 * 0.0254, rel=1e-3)


# ---------------------------------------------------------------- superseded library snapshots
def _pre_ob1_block(name, with_p_plateau):
    blk = {
        'BPNV': {'name': 'BPNV', 'a': 2e-05, 'n': 0.5, 'rho': 1700.0, 'T_flame': 2719.0, 'M': 0.0482, 'gamma': 1.17,
                 'impetus_W': 5000.0, 'heat_flux_cal_cm2_s': 69.4, 'Cp_gas': None, 'gas_mass_fraction': 0.781,
                 'kappa_jet': 8.0, 'form': 'pellets', 'particle_diameter_m': 0.0032, 'particle_LD_ratio': 1.0,
                 'heat_delivery_mode': 'demar', 'pellet_emissivity': 0.7, 'radiation_absorption_length_m': 1.0},
        'MTV': {'name': 'MTV', 'a': 4.4e-05, 'n': 0.35, 'rho': 1800.0, 'T_flame': 3000.0, 'M': 0.032, 'gamma': 1.22,
                'impetus_W': 4700.0, 'heat_flux_cal_cm2_s': 110.0, 'Cp_gas': None, 'gas_mass_fraction': 1.0,
                'kappa_jet': 8.0, 'form': 'pellets', 'particle_diameter_m': 0.0032, 'particle_LD_ratio': 1.0,
                'heat_delivery_mode': 'demar', 'pellet_emissivity': 0.7, 'radiation_absorption_length_m': 0.5},
    }[name]
    if with_p_plateau:
        blk['P_plateau'] = 0.0
    return blk


@pytest.mark.parametrize('name', ['BPNV', 'MTV'])
@pytest.mark.parametrize('with_p_plateau', [False, True])
def test_pre_ob1_library_blocks_upgrade_to_the_current_entry(name, with_p_plateau):
    from srm_1d.openmotor_adapter import library_materials, upgrade_igniter_block
    block = {'pyrogen': _pre_ob1_block(name, with_p_plateau), 'mass': 0.0045}
    new, changed = upgrade_igniter_block(block)
    assert changed is True
    cur = next(m for m in library_materials() if m['name'] == name)
    assert new['pyrogen'] == cur and cur['particle_LD_ratio'] == 1.504
    assert new['mass'] == 0.0045                          # sizing untouched
    # a customised pellet (e.g. a user's own 4 mm) is NOT a superseded snapshot
    block['pyrogen']['particle_diameter_m'] = 4.0e-3
    assert upgrade_igniter_block(block)[1] is False


def test_committed_motor_igniter_block_upgrades_on_load():
    import warnings
    from pathlib import Path
    from srm_1d.openmotor_adapter import (
        IgniterMaterialUpgradedWarning, StaleIgniterMaterialWarning, load_igniter, load_ric,
    )
    ric = Path(__file__).resolve().parents[1] / 'motors' / 'machbusterNew.ric'
    d = load_ric(str(ric))
    ig = d['data']['igniter'] if 'data' in d and 'igniter' in d['data'] else d['igniter']
    assert ig['pyrogen']['particle_LD_ratio'] == 1.0         # the committed file embeds the pre-OB1 entry
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        pyro, _ = load_igniter({'igniter': ig})
    assert pyro.particle_LD_ratio == 1.504
    assert any(issubclass(w.category, IgniterMaterialUpgradedWarning) for w in caught)
    assert not any(issubclass(w.category, StaleIgniterMaterialWarning) for w in caught)
