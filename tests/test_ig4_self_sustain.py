"""
test_ig4_self_sustain.py -- IG4: the ignition fallback's success test.

Success = after the igniter is spent the GRAIN's own gas-generation rate exceeds, at
some step, the igniter's peak gas-generation rate of the attempt. One lit cell that
cannot out-produce the igniter (pathfin54 at the DeMar charge, 0.585 g) is a hang-fire.
"""
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

pytest.importorskip("numba")
from srm_1d import openmotor_adapter as oma  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
PATHFIN = os.path.join(ROOT, "motors", "2025.12.25 54-2800 ST2.0.ric")
CHUNC = os.path.join(ROOT, "motors", "machbusterNew.ric")


def _pathfin_explicit(t_max=1.0):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        res = oma.run_from_ric(PATHFIN, pyrogen="bpnv", pyrogen_mass=0.585e-3,
                               t_max=t_max, print_interval=0.0, verbose=False)[0]
    return res, w


@pytest.fixture(scope="module")
def pathfin_dud():
    return _pathfin_explicit()


class TestHangFire:
    def test_explicit_dud_fails_the_test(self, pathfin_dud):
        res, w = pathfin_dud
        s = res["summary"]
        assert s["ignited"] is False
        assert s["grain_cells_ignited"] < s["grain_cells"]          # a cell or few lit
        assert s["ignition_attempts"] == 1 and not s["pyrogen_resized"]
        assert s["grain_mdot_after_igniter_kg_s"] <= s["igniter_peak_mdot_kg_s"]
        assert s["igniter_peak_mdot_kg_s"] > 0.0
        msgs = [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]
        hit = [m for m in msgs if "prompt start not met" in m and "hang-fire" in m]   # AS1: message reworded
        assert hit and "0.585 g" in hit[0] and "never resized" in hit[0]

    def test_pressure_cutoff_hang_fire_is_named(self, pathfin_dud):
        res, w = pathfin_dud
        if res["summary"]["termination_code"] == 2:                 # P_cutoff ended it
            assert any("hang-fire" in str(x.message) for x in w)

    def test_rates_in_summary_match_histories(self, pathfin_dud):
        res, _ = pathfin_dud
        s = res["summary"]
        assert s["igniter_peak_mdot_kg_s"] == pytest.approx(float(np.max(res["mdot_ig"])))
        t = np.asarray(res["time"])
        after = t > s["pyrogen_duration"]
        assert s["grain_mdot_after_igniter_kg_s"] == pytest.approx(
            float(np.max(np.asarray(res["mdot_grain"])[after])))

    def test_auto_charge_is_upsized_until_the_grain_drives_the_chamber(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = oma.run_from_ric(PATHFIN, t_max=0.4, print_interval=0.0, verbose=False)[0]
        s = res["summary"]
        assert s["ignited"] and s["pyrogen_resized"] and s["ignition_attempts"] > 1
        assert s["grain_mdot_after_igniter_kg_s"] > s["igniter_peak_mdot_kg_s"]
        assert s["pyrogen_mass_final_kg"] > s["ignition_attempt_masses_kg"][0]


class TestCanonicalStillSucceed:
    def test_chunc_auto(self):
        res = oma.run_from_ric(CHUNC, t_max=0.1, print_interval=0.0, verbose=False)[0]
        s = res["summary"]
        assert s["ignited"]
        assert s["grain_mdot_after_igniter_kg_s"] > s["igniter_peak_mdot_kg_s"]

    def test_all_lit_grain_with_big_charge(self):
        res = oma.run_from_ric(CHUNC, pyrogen="bpnv", pyrogen_mass=4.5e-3, t_max=0.1,
                               print_interval=0.0, verbose=False)[0]
        assert res["summary"]["ignited"] and res["summary"]["ignition_attempts"] == 1


class TestPluginAlert:
    def test_dud_raises_the_error_alert(self):
        from srm_1d.fmm_grain import _setup_openmotor_path
        pytest.importorskip("skfmm")
        _setup_openmotor_path()
        import srm_1d.srm1d_plugin as pl
        from motorlib import solvers, motor as om_motor  # type: ignore
        d = oma.load_ric(PATHFIN)
        d["igniter"] = dict(d["igniter"])
        d["igniter"]["mass"] = 0.585e-3                             # explicit: never resized
        cfg = pl._solver_config_defaults()
        cfg["t_max"] = 1.0
        sr = solvers.get_solver("srm_1d-transient").simulate(om_motor.Motor(d), config=cfg)
        assert sr.success is False
        texts = [str(getattr(a, "description", a)) for a in sr.alerts]
        assert any("self-sustaining" in t for t in texts), texts


class TestAutoSizeFailureInThePlugin:
    """The plugin returns a failed result + ERROR alert instead of letting PyrogenAutoSizeError reach the GUI."""

    def test_auto_size_error_becomes_a_failed_result(self, monkeypatch):
        from srm_1d.fmm_grain import _setup_openmotor_path
        pytest.importorskip("skfmm")
        _setup_openmotor_path()
        import srm_1d.simulation as sim
        import srm_1d.srm1d_plugin as pl
        from motorlib import solvers, motor as om_motor  # type: ignore

        def boom(*a, **k):
            raise sim.PyrogenAutoSizeError("did not light", masses_kg=[0.184e-3, 0.23e-3, 1.37e-3])

        cfg = pl._solver_config_defaults()      # (introspects the real run_simulation signature)
        monkeypatch.setattr(sim, "run_simulation", boom)
        motor = om_motor.Motor(oma.load_ric(PATHFIN))
        sr = solvers.get_solver("srm_1d-transient").simulate(motor, config=cfg)
        assert sr.success is False
        texts = [str(getattr(a, "description", a)) for a in sr.alerts]
        assert any("1.37 g" in t and "Set the pyrogen charge explicitly" in t and "0.184 g" in t for t in texts), texts
        assert any("3 attempts" in t for t in texts)

    def test_exception_carries_the_charges(self):
        import srm_1d.simulation as sim
        e = sim.PyrogenAutoSizeError("x", masses_kg=[1e-3, 2e-3])
        assert e.masses_kg == (1e-3, 2e-3)
