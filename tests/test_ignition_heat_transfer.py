"""
test_ignition_heat_transfer.py -- IG2 opt-in ignition heat-transfer closures.

``run_simulation(ignition_heat_transfer=...)`` selects the h that heats the UNIGNITED
propellant surface in the Goodman kernel: 'gnielinski' (default, unchanged), 'peretz'
(Peretz et al. 1973, Eq. III-22..III-26) or 'keller' (Keller-Baer-Ryan 1966, Eq. 11).
The Ma erosive path is not affected.
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
    def test_default_equals_explicit_gnielinski(self):
        a = _short_run()
        b = _short_run(ignition_heat_transfer="gnielinski")
        assert np.array_equal(np.asarray(a["P_head"]), np.asarray(b["P_head"]))
        assert np.array_equal(np.asarray(a["ignition_time_by_cell"]),
                              np.asarray(b["ignition_time_by_cell"]))
        assert a["summary"]["ignition_heat_transfer"] == "gnielinski"

    def test_peretz_delays_first_ignition(self):
        a = _short_run()
        b = _short_run(ignition_heat_transfer="peretz")
        assert b["summary"]["ignition_heat_transfer"] == "peretz"
        ta, tb = a["summary"]["first_ignition_time_s"], b["summary"]["first_ignition_time_s"]
        assert ta is not None and tb is not None and tb > ta

    def test_unknown_mode_raises_from_run(self):
        with pytest.raises(ValueError, match="ignition_heat_transfer"):
            _short_run(ignition_heat_transfer="nope")
