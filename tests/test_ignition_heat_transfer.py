"""
test_ignition_heat_transfer.py -- IG2/IG3 ignition heat-transfer closures.

``run_simulation(ignition_heat_transfer=...)`` selects the h that heats the UNIGNITED
propellant surface in the Goodman kernel: 'peretz' (DEFAULT since IG3; Peretz et al. 1973,
Eq. III-22..III-26), 'gnielinski' (the pre-IG3 default, Ma/Gnielinski chain) or 'keller'
(Keller-Baer-Ryan 1966, Eq. 11). ``peretz_ap_reference`` ('local' | 'head' | 'aft') picks the
port area in Peretz's (A_p/A_t)^0.4 factor. The Ma erosive path is not affected.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d import simulation as sim  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def _peretz(P=2.0e6, u=300.0, Tg=2000.0, Tw=300.0, R=290.0, Cp=2000.0, gam=1.2,
            D=0.0127, x=0.05, Ap=1.27e-4, At=1.27e-4):
    return sim._peretz_heat_transfer_coeff.py_func(P, u, Tg, Tw, R, Cp, gam, D, x, Ap, At)


class TestPeretz:
    def test_matches_hand_evaluation(self):
        P, u, Tg, Tw, R, Cp, gam, D, x = 2.0e6, 300.0, 2000.0, 300.0, 290.0, 2000.0, 1.2, 0.0127, 0.05
        Ap_At = 1.5
        Taf = 0.5 * (Tg + Tw)
        W = 8314.462618 / R
        mu = 0.87e-7 * math.sqrt(W) * Taf ** 0.65
        rho = P / (R * Taf)
        Re = rho * u * D / mu
        Pr = gam / (1.77 * gam - 0.45)
        Nu = 2.56e-2 * (x / D) ** -0.1 * Re ** 0.8 * Pr ** 0.4 * Ap_At ** 0.4
        want = Nu * (Cp * mu / Pr) / D
        got = _peretz(P=P, u=u, Tg=Tg, Tw=Tw, R=R, Cp=Cp, gam=gam, D=D, x=x, Ap=1.5e-4, At=1.0e-4)
        assert got == pytest.approx(want, rel=1e-12)

    def test_coefficient_consistency_with_eq_iii_31(self):
        # a_h = (0.87e-6)^0.2 * a6 = 1.56e-3 (Peretz Eq. III-31 vs III-22)
        assert (0.87e-6) ** 0.2 * 2.56e-2 == pytest.approx(1.56e-3, rel=0.02)

    def test_1972_paper_constant_is_a_different_calibration_of_the_same_law(self):
        # AIAA Paper 72-1119 Eq. 14 prints 1.78e-3 for the SAME variables/exponents
        # (h = a Pr^-0.6 cp (pu/R)^0.8 W^0.1 Taf^-0.67 (xr dh)^-0.1 (Ap/At)^0.4, CGS).
        # The thesis (Eq. III-31) derives 1.56e-3 from a6 = 2.56e-2; the code carries
        # the thesis value. The paper value would correspond to a6 = 2.92e-2.
        assert 1.78e-3 / 1.56e-3 == pytest.approx(1.141, abs=1e-3)
        assert 1.78e-3 / (0.87e-6) ** 0.2 == pytest.approx(2.92e-2, rel=0.01)

    def test_mu_air_sanity(self):
        # Eq. III-24 for air (W=28.96) at 1500 K ~ 5.4e-5 Pa s
        W = 28.96
        assert 0.87e-7 * math.sqrt(W) * 1500.0 ** 0.65 == pytest.approx(5.4e-5, rel=0.03)

    def test_scalings(self):
        base = _peretz()
        # h ~ u^0.8 (above the laminar floor)
        assert _peretz(u=600.0) / base == pytest.approx(2.0 ** 0.8, rel=1e-9)
        # h ~ (Ap/At)^0.4
        assert _peretz(Ap=2.54e-4) / base == pytest.approx(2.0 ** 0.4, rel=1e-9)

    def test_laminar_floor_and_degenerate_inputs(self):
        assert _peretz(u=1e-9) > 0.0                      # Nu floor 3.66, no NaN
        assert _peretz(D=0.0) == 0.0 and _peretz(At=0.0) == 0.0


class TestKeller:
    def test_matches_eq_11_at_keller_diameter(self):
        rho, u, Tg = 5.0, 200.0, 2000.0
        G = rho * u * 0.1  # g/cm2 s
        want = 1.435e-4 * Tg ** 0.3 * G ** 0.905 * 41840.0
        got = sim._keller_heat_transfer_coeff.py_func(rho, u, Tg, sim._KELLER_D_H)
        assert got == pytest.approx(want, rel=1e-12)

    def test_diameter_scaling_exponent(self):
        a = sim._keller_heat_transfer_coeff.py_func(5.0, 200.0, 2000.0, 0.0127)
        b = sim._keller_heat_transfer_coeff.py_func(5.0, 200.0, 2000.0, 0.0254)
        assert b / a == pytest.approx(2.0 ** -0.095, rel=1e-12)


class TestModeResolution:
    def test_codes(self):
        assert sim._ignition_h_code('gnielinski') == 0
        assert sim._ignition_h_code('peretz') == 1
        assert sim._ignition_h_code('keller') == 2

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="ignition_heat_transfer"):
            sim._ignition_h_code('dittus_boelter')


def _short_run(**kw):
    from srm_1d.openmotor_adapter import run_from_ric
    res = run_from_ric(os.path.join(ROOT, "motors", "machbusterNew.ric"), pyrogen="bpnv",
                       pyrogen_mass=4.5e-3, injection_topology="head_basket", roughness=32e-6,
                       kappa=0.44, T_ignition=756.0, t_max=0.006, print_interval=0.0,
                       verbose=False, **kw)[0]
    return res


class TestDefaultIdentityAndEffect:
    def test_default_is_peretz_local(self):
        import inspect
        params = inspect.signature(sim.run_simulation).parameters
        assert params["ignition_heat_transfer"].default == "peretz"
        assert params["peretz_ap_reference"].default == "local"

    def test_default_equals_explicit_peretz_local(self):
        a = _short_run()
        b = _short_run(ignition_heat_transfer="peretz", peretz_ap_reference="local")
        assert np.array_equal(np.asarray(a["P_head"]), np.asarray(b["P_head"]))
        assert np.array_equal(np.asarray(a["ignition_time_by_cell"]),
                              np.asarray(b["ignition_time_by_cell"]))
        assert a["summary"]["ignition_heat_transfer"] == "peretz"
        assert a["summary"]["peretz_ap_reference"] == "local"

    def test_peretz_delays_first_ignition_vs_gnielinski(self):
        a = _short_run(ignition_heat_transfer="gnielinski")
        b = _short_run()
        assert a["summary"]["ignition_heat_transfer"] == "gnielinski"
        ta, tb = a["summary"]["first_ignition_time_s"], b["summary"]["first_ignition_time_s"]
        assert ta is not None and tb is not None and tb > ta

    def test_unknown_mode_raises_from_run(self):
        with pytest.raises(ValueError, match="ignition_heat_transfer"):
            _short_run(ignition_heat_transfer="nope")


class TestPeretzApReference:
    """IG3: port area of the (A_p/A_t)^0.4 factor: local | head | aft."""

    def test_codes_and_unknown(self):
        assert sim._peretz_ap_code('local') == 0
        assert sim._peretz_ap_code('head') == 1
        assert sim._peretz_ap_code('aft') == 2
        with pytest.raises(ValueError, match="peretz_ap_reference"):
            sim._peretz_ap_code('mid')

    def test_unknown_raises_from_run(self):
        with pytest.raises(ValueError, match="peretz_ap_reference"):
            _short_run(peretz_ap_reference="mid")

    @staticmethod
    def _kernel_Tsurf(ap_code, A_port, is_grain=None):
        """One call of the Goodman kernel in peretz mode; returns T_surf after the step."""
        N = len(A_port)
        dx = 0.01
        P = np.full(N, 2.0e6)
        T = np.full(N, 2000.0)
        T_surf = np.full(N, 300.0)
        delta = np.full(N, 1.0e-6)
        has_ignited = np.zeros(N, dtype=np.bool_)
        is_burning = np.zeros(N, dtype=np.bool_)
        if is_grain is None:
            is_grain = np.ones(N, dtype=np.bool_)
        args = (
            P, T, T_surf, delta, has_ignited, is_burning, is_grain,
            np.full(N, 1.0e10), np.zeros(N), np.zeros(N),
            np.zeros(N), np.zeros(N),
            np.ones(N), np.zeros(N), np.zeros(N),
            np.zeros(N), np.zeros(N), np.zeros(N, dtype=np.bool_),
            (np.arange(N) + 0.5) * dx, np.full(N, 5.0e4), np.full(N, 0.0127), np.full(N, 0.02),
            0.0, 1.0e-4, 1700.0, 3041.0, 293.0,
            0.5, 0.3685, 50e-6, 0.45, 1.0e-7, 0.3,
            10000.0, N, dx, 0.0, 300.0,
            2060.0, 1385.0,
            0.0, 0.45, False, 0.0,
            np.zeros((N, 3)),
            np.ones(N), False,
            np.concatenate([np.zeros((N, 2)), np.ones((N, 1))], axis=1), np.zeros(N),
            False, np.ones(N, dtype=np.bool_), 0,
            np.full(N, 5.0), np.full(N, 1500.0),
            1, np.full(N, 300.0), np.asarray(A_port, dtype=np.float64), 1.0e-4,
            np.full(N, 290.0), np.full(N, 1.2), ap_code,
            0, np.array([0.0, 1.0e-6]), np.full((1, 2), 293.0), np.zeros(N),   # SM1: goodman, dummy conduction state
            np.full(N, 3, dtype=np.int64), np.zeros((2, 2)), 0.0125,
            np.zeros((N, 15)), False, np.zeros((1, 1)), np.zeros((1, 1)), np.zeros(1), np.zeros(1),
            np.zeros(1, dtype=np.int64),   # SM2: diag, recorder off
            0, 0,                          # RT1: static driver; HE1: jet off
        )
        sim._goodman_ignition_sources_and_mass(*args)
        return T_surf

    def test_head_and_aft_select_the_end_cell_area(self):
        # port widening head -> aft (area ratio 1 -> 4); identical gas state in every cell
        A = np.array([1.0, 1.5, 2.0, 3.0, 4.0]) * 1.0e-4
        loc, head, aft = (self._kernel_Tsurf(c, A) for c in (0, 1, 2))
        # h ~ (Ap/At)^0.4 and the surface heats monotonically with h
        assert np.all(head[1:] < loc[1:]) and np.all(aft[:-1] > loc[:-1])
        # the reference cell itself: head == local at cell 0, aft == local at the last cell
        assert head[0] == loc[0] and aft[-1] == loc[-1]
        # per cell, the heating ratio follows (A_ref/A_local)^0.4 in the early-time limit
        # (cells differ only through the x-entrance factor, common to all three readings)
        for i in range(len(A)):
            assert (head[i] - 300.0) / (loc[i] - 300.0) == pytest.approx((A[0] / A[i]) ** 0.4, rel=0.03)
            assert (aft[i] - 300.0) / (loc[i] - 300.0) == pytest.approx((A[-1] / A[i]) ** 0.4, rel=0.03)

    def test_head_aft_use_grain_cells_only(self):
        # a non-grain (casing) cell fore and aft of the grain must not be the reference
        A = np.array([9.0, 1.0, 2.0, 3.0, 9.0]) * 1.0e-4
        g = np.array([False, True, True, True, False])
        head = self._kernel_Tsurf(1, A, g)
        aft = self._kernel_Tsurf(2, A, g)
        # grain cells (1..3): head reference = cell 1 (1e-4), aft = cell 3 (3e-4)
        assert aft[2] - 300.0 == pytest.approx((head[2] - 300.0) * 3.0 ** 0.4, rel=0.02)

    def test_run_head_equals_local_on_uniform_head_section_aft_lights_sooner(self):
        a = _short_run(peretz_ap_reference="local")
        b = _short_run(peretz_ap_reference="head")
        c = _short_run(peretz_ap_reference="aft")
        # Chunc's port is uniform over the head-end cells where it first lights, so head ~ local
        # there (later cells differ through the taper); aft is wider (A_p/A_t 1.80 vs 0.557)
        # -> larger h -> earlier first light
        ta, tb, tc = (r["summary"]["first_ignition_time_s"] for r in (a, b, c))
        assert ta is not None and tb is not None and tc is not None
        assert tb == pytest.approx(ta, rel=1e-3)
        assert tc < ta
