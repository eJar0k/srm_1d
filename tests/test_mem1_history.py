"""
test_mem1_history.py -- MEM1: bounded history memory.

* no fixed 5,000,000-row floor (history is sized from the CFL estimate only);
* an auto-sized history that fills triggers ONE re-run with double the capacity (never a
  silent truncation); an explicit history_capacity is honoured as given;
* the in-place trim keeps the values (history size never changes the physics).
"""
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

pytest.importorskip("numba")
from srm_1d import simulation as sim  # noqa: E402
from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
CHUNC = os.path.join(ROOT, "motors", "machbusterNew.ric")
KW = dict(pyrogen="bpnv", pyrogen_mass=4.5e-3, t_max=0.01, print_interval=0.0, verbose=False)


def _run(**extra):
    kw = dict(KW)
    kw.update(extra)
    return run_from_ric(CHUNC, **kw)[0]


@pytest.fixture(scope="module")
def base():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _run()


class TestNoFloor:
    def test_capacity_is_the_estimate(self, base):
        s = base["summary"]
        assert s["history_capacity"] == s["history_estimated_steps"]
        assert s["history_capacity"] < 5_000_000
        assert s["steps"] <= s["history_capacity"]
        assert len(base["time"]) == s["steps"]

    def test_arrays_have_exactly_n_steps(self, base):
        n = base["summary"]["steps"]
        for key in ("time", "P_head", "mdot_ig", "mdot_grain", "n_burning", "dt"):
            assert len(base[key]) == n, key


class TestShrinkInPlace:
    def test_shrink_keeps_values_and_identity(self):
        a = np.arange(10.0)
        b = sim._shrink_hist(a, 4)
        assert b is a and b.shape == (4,)
        assert np.array_equal(b, [0.0, 1.0, 2.0, 3.0])

    def test_history_size_does_not_change_the_physics(self, base):
        big = _run(history_capacity=400_000)
        assert big["summary"]["history_capacity"] == 400_000
        assert big["summary"]["steps"] == base["summary"]["steps"]
        for key in ("time", "P_head", "mdot_grain", "T_ig"):
            assert np.array_equal(np.asarray(big[key]), np.asarray(base[key])), key


class TestFullHistory:
    def test_auto_history_that_fills_is_rerun_once_with_double(self, base, monkeypatch):
        steps = base["summary"]["steps"]
        est = int(0.7 * steps)                     # too small; 2x (1.4 steps) is enough
        monkeypatch.setattr(sim, "_estimate_history_steps", lambda *a, **k: est)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = _run()
        s = res["summary"]
        assert s["termination_code"] != 3
        assert s["history_capacity"] == 2 * est and s["steps"] == steps
        msgs = [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]
        assert any("re-running once" in m for m in msgs)
        assert not any("history arrays filled" in m for m in msgs)    # superseded
        assert np.array_equal(np.asarray(res["P_head"]), np.asarray(base["P_head"]))

    def test_second_overflow_still_truncates_with_warning(self, base, monkeypatch):
        steps = base["summary"]["steps"]
        est = int(0.3 * steps)                     # 2x = 0.6 steps: still too small
        monkeypatch.setattr(sim, "_estimate_history_steps", lambda *a, **k: est)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = _run()
        assert res["summary"]["termination_code"] == 3
        assert res["summary"]["history_capacity"] == 2 * est
        msgs = [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]
        assert any("history arrays filled" in m for m in msgs)

    def test_explicit_capacity_is_honoured_and_not_rerun(self):
        with pytest.warns(RuntimeWarning, match="history_capacity="):
            res = _run(history_capacity=300)
        assert res["summary"]["termination"] == "history array full"
        assert res["summary"]["history_capacity"] == 300 and len(res["time"]) == 300
