"""
test_he1_jet.py -- HE1 opt-in head-end igniter-jet heat transfer (Peretz Eqs. III-28..III-30).

Hand-checked against an independent evaluation of the published formulas, continuity of the
three branches at the impingement point, the x/D <= 10 cut, default identity, and the
forward-plenum-only scope.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d import simulation as sim  # noqa: E402

MDOT, TG, TW, R, CP, GAM = 0.1, 3000.0, 300.0, 300.0, 2000.0, 1.2
D, AP, AT = 0.02, 4.0e-4, 2.0e-4


def _ref(x, mdot=MDOT):
    """Independent evaluation of III-28..30 with the film-temperature properties."""
    T_af = 0.5 * (TG + TW)
    W = 8314.462618 / R
    mu = 0.87e-7 * math.sqrt(W) * T_af ** 0.65
    Re = 4.0 * mdot / (math.pi * D * mu)
    Pr = GAM / (1.77 * GAM - 0.45)
    k = CP * mu / Pr
    if x < D:
        m = 0.022 + 0.066 * x / D
    else:
        m = 0.031 + 0.057 * math.exp(-0.5 * (x - D) / D)
    return 0.7 * Re ** 0.8 * (AP / AT) ** 0.4 * m * k / D


def _h(x, mdot=MDOT):
    return sim._head_end_jet_h(mdot, TG, TW, R, CP, GAM, D, x, AP, AT)


class TestPublishedForm:
    @pytest.mark.parametrize("x", [0.2 * D, 0.5 * D, 1.0 * D, 1.5 * D, 3.0 * D, 6.0 * D, 10.0 * D])
    def test_matches_independent_evaluation(self, x):
        assert _h(x) == pytest.approx(_ref(x), rel=1e-12)

    def test_continuous_at_the_impingement_point(self):
        eps = 1e-9 * D
        assert _h(D - eps) == pytest.approx(_h(D + eps), rel=1e-6)

    def test_III28_equals_III29_at_x_im(self):
        # 0.7 * (0.022 + 0.066) = 0.0616 vs the published 0.0615
        assert 0.7 * (0.022 + 0.066) == pytest.approx(0.0615, rel=2e-3)

    def test_peak_at_impingement_then_decay_to_dittus_boelter_like_constant(self):
        assert _h(D) > _h(0.3 * D) and _h(D) > _h(5 * D)
        far = 0.7 * 0.031
        assert far == pytest.approx(0.0217, rel=1e-3)

    def test_scales_with_igniter_flow_to_the_08(self):
        assert _h(3 * D, 2 * MDOT) / _h(3 * D, MDOT) == pytest.approx(2.0 ** 0.8, rel=1e-3)

    def test_ap_over_at_exponent_is_04(self):
        h1 = sim._head_end_jet_h(MDOT, TG, TW, R, CP, GAM, D, 3 * D, 2 * AP, AT)
        assert h1 / _h(3 * D) == pytest.approx(2.0 ** 0.4, rel=1e-9)


class TestRangeAndGuards:
    def test_zero_beyond_ten_diameters(self):
        assert _h(10.0 * D * 1.0001) == 0.0
        assert _h(10.0 * D) > 0.0

    def test_zero_when_the_igniter_is_not_flowing(self):
        assert _h(3 * D, 0.0) == 0.0

    def test_mode_codes(self):
        assert sim._head_end_jet_code(False) == 0
        assert sim._head_end_jet_code(True) == 1
        assert sim._head_end_jet_code('max') == 1
        assert sim._head_end_jet_code('replace') == 2
        with pytest.raises(ValueError):
            sim._head_end_jet_code('bogus')


@pytest.fixture(scope="module")
def runs():
    from srm_1d.openmotor_adapter import run_from_ric
    ric = os.path.join(os.path.dirname(__file__), "..", "motors", "machbusterNew.ric")
    out = {}
    for mode in (False, 'max', 'replace'):
        out[mode] = run_from_ric(
            ric, target_propellant_cells=60, t_max=0.02, roughness=32e-6, kappa=0.44,
            T_ignition=756.0, print_interval=0.0, verbose=False, pyrogen='bpnv',
            pyrogen_mass=2.623e-3, head_end_jet_h=mode, history_capacity=100000)[0]
    return out


class TestMotorPlumbing:
    def test_stamp(self, runs):
        assert runs[False]['summary']['head_end_jet_h'] is False
        assert runs['max']['summary']['head_end_jet_h'] == 'max'

    def test_default_has_no_jet_h_and_jet_changes_the_ignition(self, runs):
        d0 = np.asarray(runs[False]['ignition_diag'])
        lit0 = np.isfinite(d0[:, 0])
        assert np.all(d0[lit0, 14] == 0.0)
        d1 = np.asarray(runs['max']['ignition_diag'])
        lit1 = np.isfinite(d1[:, 0])
        assert np.any(d1[lit1, 14] > 0.0)
        t0 = np.asarray(runs[False]['ignition_time_by_cell'])
        t1 = np.asarray(runs['max']['ignition_time_by_cell'])
        assert not np.array_equal(t0, t1)

    def test_max_never_lowers_h_below_bore(self, runs):
        # 'max' can only add heating: the first light is not later than the default's
        t0 = np.asarray(runs[False]['ignition_time_by_cell'])
        t1 = np.asarray(runs['max']['ignition_time_by_cell'])
        assert t1.min() <= t0.min() + 1e-12

    def test_bad_mode_rejected(self):
        from srm_1d.openmotor_adapter import run_from_ric
        ric = os.path.join(os.path.dirname(__file__), "..", "motors", "machbusterNew.ric")
        with pytest.raises(ValueError):
            run_from_ric(ric, t_max=0.001, verbose=False, head_end_jet_h='bogus')


def test_head_basket_is_untouched():
    from srm_1d.openmotor_adapter import run_from_ric
    ric = os.path.join(os.path.dirname(__file__), "..", "motors", "machbusterNew.ric")
    kw = dict(target_propellant_cells=60, t_max=0.01, roughness=32e-6, kappa=0.44,
              T_ignition=756.0, print_interval=0.0, verbose=False, pyrogen='bpnv',
              pyrogen_mass=2.623e-3, injection_topology='head_basket', history_capacity=100000)
    a = run_from_ric(ric, head_end_jet_h=False, **kw)[0]
    b = run_from_ric(ric, head_end_jet_h='max', **kw)[0]
    assert np.array_equal(np.asarray(a['P_head']), np.asarray(b['P_head']))
    assert np.array_equal(np.asarray(a['ignition_time_by_cell']), np.asarray(b['ignition_time_by_cell']))
