"""
test_solid_heating_models.py -- SM1: the three solid-heating models of the ignition gate
(``run_simulation(solid_heating_model='goodman'|'peretz'|'conduction')``).

Exact-solution checks of the kernels (semi-infinite solid, Robin / constant flux / step
changes in h, the latter against an independent Volterra solve) and the default-identity
gate. See docs/pyrogen_realism/SM1_REPORT.md.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d import simulation as sim  # noqa: E402
from srm_1d import solid_thermal as st  # noqa: E402
from srm_1d.tools import solid_verify as sv  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")

K = 0.271
AL = K / (1800.0 * 1500.0)
TI, TG = 293.0, 2800.0
Y = st.build_solid_grid(AL, 2.0, st.CONDUCTION_DY0, st.CONDUCTION_RATIO)
AF = st.CONDUCTION_AGE_FRACTION


def _robin(h, dt, tend=0.05):
    n = int(tend / dt)
    t = (np.arange(n) + 1) * dt
    return t, np.full(n, h), np.full(n, TG), sv.exact_robin(t, h, TG, TI, AL, K)


class TestKernelsAgainstExactSolutions:
    @pytest.mark.parametrize("h", [1000.0, 5000.0, 20000.0])
    @pytest.mark.parametrize("dt", [2e-6, 2e-5])
    def test_conduction_matches_robin(self, h, dt):
        t, hh, TT, ex = _robin(h, dt)
        r = sv.run_conduction(hh, TT, dt, TI, AL, K, Y, AF)
        win = ex <= 900.0
        assert np.abs(r - ex)[win].max() < 1.0                       # K
        t756 = sv.time_to(ex, dt, 756.0, TI)
        assert abs(sv.time_to(r, dt, 756.0, TI) - t756) / t756 < 0.005

    @pytest.mark.parametrize("h", [1000.0, 5000.0])
    def test_integral_models_within_goodman_accuracy(self, h):
        # the cubic-profile integral method is ~2 % high in T_s - T_i (Goodman 1958)
        t, hh, TT, ex = _robin(h, 2e-6)
        win = ex <= 900.0
        for r in (sv.run_goodman(hh, TT, 2e-6, TI, AL, K),
                  sv.run_peretz(hh, TT, 2e-6, TI, AL, K)):
            assert np.abs(r - ex)[win].max() < 0.03 * 607.0

    def test_peretz_equals_goodman_at_constant_h(self):
        t, hh, TT, ex = _robin(5000.0, 2e-6)
        g = sv.run_goodman(hh, TT, 2e-6, TI, AL, K)
        p = sv.run_peretz(hh, TT, 2e-6, TI, AL, K)
        # they differ only by the start-up seed (delta0 = 1e-6 m vs eps = 0.1 K), which
        # washes out; at late times the two ODEs are the same solution
        win = (ex <= 900.0) & (t > 5e-4)
        assert np.abs(g - p)[win].max() < 0.5

    def test_constant_flux_sqrt_t(self):
        q, h, dt = 1.0e7, 1.0, 2e-6
        n = int(2e-3 / dt)
        t = (np.arange(n) + 1) * dt
        ex = sv.exact_flux(t, q, TI, AL, K)
        r = sv.run_conduction(np.full(n, h), np.full(n, TI + q / h), dt, TI, AL, K, Y, AF)
        assert np.abs(r - ex).max() < 1.0

    @pytest.mark.parametrize("h1,fac,t1", [(5400.0, 1 / 14.0, 0.8e-3), (500.0, 10.0, 1.0e-3)])
    def test_step_changes_against_volterra(self, h1, fac, t1):
        dt, te = 2e-6, 4e-3
        hf = lambda t: np.where(t <= t1, h1, h1 * fac)  # noqa: E731
        tv, v = sv.volterra_richardson(hf, lambda t: np.full_like(t, TG), te, dt, TI, AL, K)
        n = int(te / dt)
        tt = (np.arange(n) + 1) * dt
        hh, TT = hf(tt), np.full(n, TG)
        c = sv.run_conduction(hh, TT, dt, TI, AL, K, Y, AF)
        g = sv.run_goodman(hh, TT, dt, TI, AL, K)
        p = sv.run_peretz(hh, TT, dt, TI, AL, K)
        after = tt > t1
        assert np.abs(c - v)[after].max() < 3.0
        # Goodman's algebraic T_s(h) and Peretz's held T_s are both wrong by 100s of K
        # right after the step, by construction (this is the SM1 finding)
        assert np.abs(g - v)[after].max() > 100.0
        assert np.abs(p - v)[after].max() > 20.0

    def test_peretz_holds_T_s_when_h_drops_goodman_does_not(self):
        n = 400
        dt = 2e-6
        h = np.where(np.arange(n) < 200, 5400.0, 5400.0 / 14.0)
        TT = np.full(n, TG)
        g = sv.run_goodman(h, TT, dt, TI, AL, K)
        p = sv.run_peretz(h, TT, dt, TI, AL, K)
        assert g[200] < g[199] - 50.0       # algebraic drop on the h change
        assert p[200] >= p[199]             # state carried, no drop

    def test_conduction_cools_slowly_after_h_drop(self):
        n = 400
        dt = 2e-6
        h = np.where(np.arange(n) < 200, 5400.0, 0.0)
        TT = np.full(n, TG)
        c = sv.run_conduction(h, TT, dt, TI, AL, K, Y, AF)
        assert c[200] < c[199]              # insulated surface relaxes into the solid
        assert c[-1] > TI + 0.4 * (c[199] - TI)

    def test_degenerate_inputs_hold(self):
        assert st._step_peretz_ode(400.0, 0.0, 2800.0, TI, AL, K, 1e-5) == 400.0
        assert st._step_peretz_ode(400.0, 5000.0, 300.0, TI, AL, K, 1e-5) == 400.0  # gas colder
        T = np.full(len(Y), TI)
        cp, dp = np.zeros(len(Y)), np.zeros(len(Y))
        Ts, m = st._step_conduction(T, Y, 3, 0.0, 0.0, 2800.0, TI, AL, K, 1e-5, AF, cp, dp)
        assert Ts == pytest.approx(TI)


class TestSimIntegration:
    @staticmethod
    def _run(**kw):
        from srm_1d.openmotor_adapter import run_from_ric
        return run_from_ric(os.path.join(ROOT, "motors", "machbusterNew.ric"), pyrogen="bpnv",
                            pyrogen_mass=4.5e-3, injection_topology="head_basket",
                            roughness=32e-6, kappa=0.44, T_ignition=756.0, t_max=0.006,
                            print_interval=0.0, verbose=False, **kw)[0]

    def test_default_is_goodman(self):
        import inspect
        assert inspect.signature(sim.run_simulation).parameters[
            "solid_heating_model"].default == "goodman"

    def test_default_equals_explicit_goodman_bitwise(self):
        a = self._run()
        b = self._run(solid_heating_model="goodman")
        assert np.array_equal(np.asarray(a["P_head"]), np.asarray(b["P_head"]))
        assert np.array_equal(np.asarray(a["ignition_time_by_cell"]),
                              np.asarray(b["ignition_time_by_cell"]))
        assert a["summary"]["solid_heating_model"] == "goodman"

    def test_codes_and_unknown(self):
        assert sim._solid_model_code("goodman") == sim.SOLID_GOODMAN
        assert sim._solid_model_code("peretz") == sim.SOLID_PERETZ
        assert sim._solid_model_code("conduction") == sim.SOLID_CONDUCTION
        with pytest.raises(ValueError, match="solid_heating_model"):
            sim._solid_model_code("exact")
        with pytest.raises(ValueError, match="solid_heating_model"):
            self._run(solid_heating_model="nope")

    @pytest.mark.parametrize("model", ["peretz", "conduction"])
    def test_models_run_and_ignite(self, model):
        a = self._run()
        b = self._run(solid_heating_model=model)
        assert b["summary"]["solid_heating_model"] == model
        ta = a["summary"]["first_ignition_time_s"]
        tb = b["summary"]["first_ignition_time_s"]
        assert ta is not None and tb is not None
        assert 0.3 * ta < tb < 3.0 * ta     # same physics, a different closure: same decade
