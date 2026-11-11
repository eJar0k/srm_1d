"""
test_as1_prompt_start.py -- AS1: the prompt-start success criterion of the auto igniter charge.

Success = at least ``ignition_prompt_fraction`` (0.95) of the grain cells ignited within ``ignition_prompt_window_s`` (0.5 s) of
igniter firing. An attempt that has not met it at the window end is failed there (termination_code 7) and retried at x1.25;
explicit charges are never resized (a warning reports the lit fraction at the window).
"""
import functools
import inspect
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

pytest.importorskip("numba")
from srm_1d import openmotor_adapter as oma  # noqa: E402
from srm_1d import simulation as sim  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
CHUNC = os.path.join(ROOT, "motors", "machbusterNew.ric")


def _run(**kw):
    kw.setdefault("t_max", 0.4)
    kw.setdefault("P_cutoff", 0.01e6)      # B0-row cutoff: a slow start is not ended by the 0.5 MPa blowdown cutoff
    kw.setdefault("print_interval", 0.0)
    kw.setdefault("verbose", False)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        res = oma.run_from_ric(CHUNC, **kw)[0]
    return res, [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]


class TestDefaults:
    def test_signature_defaults_are_the_users_numbers(self):
        p = inspect.signature(sim.run_simulation).parameters
        assert p["ignition_prompt_fraction"].default == 0.95
        assert p["ignition_prompt_window_s"].default == 0.5

    def test_invalid_values_raise(self):
        with pytest.raises(ValueError, match="ignition_prompt_fraction"):
            _run(ignition_prompt_fraction=0.0)
        with pytest.raises(ValueError, match="ignition_prompt_window_s"):
            _run(ignition_prompt_window_s=0.0)


class TestSummaryKeys:
    @pytest.fixture(scope="class")
    def run(self):
        return _run(pyrogen="bpnv", pyrogen_mass=4.5e-3, t_max=0.1)

    def test_prompt_keys(self, run):
        res, _ = run
        s = res["summary"]
        assert s["prompt_fraction_required"] == 0.95 and s["prompt_window_s"] == 0.5
        # the run ended at 0.1 s < the window: judged at its end
        assert s["prompt_window_effective_s"] == pytest.approx(float(res["time"][-1]))
        ign = np.asarray(res["ignition_time_by_cell"])
        isg = np.asarray(res["snapshots"][0]["is_grain"])
        want = float(np.mean(ign[isg] <= s["prompt_window_effective_s"]))
        assert s["prompt_lit_fraction"] == pytest.approx(want)
        assert s["prompt_start"] == (want >= 0.95)
        assert s["ignition_t90_s"] > s["ignition_t50_s"] > 0.0

    def test_demar_equivalent_pressure(self, run):
        res, _ = run
        s = res["summary"]
        V, W = s["pyrogen_free_volume_in3"], 5000.0       # BPNV impetus
        assert V > 0.0
        assert s["pyrogen_demar_equiv_psi"] == pytest.approx(4.5 * W / V, rel=1e-9)


class TestExplicitCharge:
    def test_failing_prompt_warns_with_the_lit_fraction_and_is_not_resized(self):
        # DeMar-size charge (0.8 g) lights the Chunc grain only after ~0.23 s: fails a 0.01 s window
        res, msgs = _run(pyrogen="bpnv", pyrogen_mass=0.815e-3, t_max=0.4, ignition_prompt_window_s=0.01)
        s = res["summary"]
        assert s["prompt_start"] is False and s["prompt_lit_fraction"] < 0.95
        assert s["ignition_attempts"] == 1 and not s["pyrogen_resized"] and not s["pyrogen_mass_auto"]
        assert s["termination_code"] != 7                  # an explicit run is never cut at the window
        assert s["ignited"] is True                        # it did light eventually (slow start)
        hit = [m for m in msgs if "prompt start not met" in m]
        assert hit and "explicitly set" in hit[0] and "never resized" in hit[0]
        assert f"{100 * s['prompt_lit_fraction']:.0f} %" in hit[0] and "slow start" in hit[0]


class TestAutoCharge:
    def test_attempts_failing_the_window_are_cut_there_and_retried(self):
        # a 10 ms window: the 0.8 g DeMar charge (99 % lit at 0.235 s) fails it; bigger charges meet it
        res, msgs = _run(ignition_prompt_window_s=0.01, t_max=0.3)
        s = res["summary"]
        assert s["pyrogen_mass_auto"] and s["ignition_attempts"] > 1 and s["pyrogen_resized"]
        terms = s["ignition_attempt_terminations"]
        assert 7 in terms[:-1] and terms[-1] != 7          # failed attempts stopped at the window; the last ran on
        fr = s["ignition_attempt_prompt_fractions"]
        assert all(f < 0.95 for f in fr[:-1]) and fr[-1] >= 0.95 and s["prompt_start"]
        m = s["ignition_attempt_masses_kg"]
        assert all(abs(m[k + 1] / m[k] - 1.25) < 1e-9 for k in range(len(m) - 1))
        assert any("auto-sized" in x and "DeMar-equivalent" in x and "ignites within 0.01 s" in x for x in msgs)
        # the cut attempts were cheap: ended at the window, not at t_max
        assert res["time"][-1] > 0.01

    def test_cap_reached_but_the_grain_lit_returns_a_full_run_with_a_warning(self):
        res, msgs = _run(ignition_prompt_window_s=1.0e-4, ignition_retry_max_attempts=2, t_max=0.4)
        s = res["summary"]
        assert s["ignition_attempts"] == 2 and s["prompt_start"] is False and s["ignited"] is True
        assert s["termination_code"] != 7 and res["time"][-1] > 0.3     # the last attempt was completed
        assert any("prompt start not met" in x and "auto-sizing cap" in x for x in msgs)

    def test_nothing_lit_after_the_cap_raises(self):
        with pytest.raises(sim.PyrogenAutoSizeError):
            _run(ignition_retry_max_attempts=1, pyrogen="bpnv", t_max=0.05, T_ignition=3500.0)


class TestPluginWarning:
    def test_slow_ignition_is_a_warning_not_an_error(self):
        pytest.importorskip("skfmm")
        from srm_1d.fmm_grain import _setup_openmotor_path
        _setup_openmotor_path()
        import srm_1d.srm1d_plugin as pl
        from motorlib import solvers, motor as om_motor  # type: ignore
        from motorlib.simResult import SimAlertLevel  # type: ignore
        d = oma.load_ric(CHUNC)
        d["igniter"] = dict(d["igniter"])
        d["igniter"]["mass"] = 0.815e-3                    # explicit: lights the grain after ~0.23 s
        cfg = pl._solver_config_defaults()
        cfg["t_max"] = 0.4
        cfg["P_cutoff"] = 1.0e4
        orig = sim.run_simulation
        wrapped = functools.wraps(orig)(lambda *a, **k: orig(*a, ignition_prompt_window_s=0.01, **k))
        sim.run_simulation = wrapped
        try:
            sr = solvers.get_solver("srm_1d-transient").simulate(om_motor.Motor(d), config=cfg)
        finally:
            sim.run_simulation = orig
        assert sr.success is True
        levels = [a.level for a in sr.alerts]
        assert SimAlertLevel.WARNING in levels and SimAlertLevel.ERROR not in levels
        texts = [str(getattr(a, "description", a)) for a in sr.alerts]
        assert any("slow ignition" in t and "within 0.01 s" in t for t in texts), texts
