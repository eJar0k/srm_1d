"""
test_ig3_stage2.py -- IG3 stage 2: DeMar flux off by default, DeMar charge sizing,
auto-upsize ignition fallback with a physical early stop.
"""
import dataclasses
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d import openmotor_adapter as oma  # noqa: E402
from srm_1d import simulation as sim  # noqa: E402
from srm_1d.igniter_plenum import demar_pyrogen_mass, sutton_pyrogen_mass  # noqa: E402
from srm_1d.propellant import Pyrogen  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
CHUNC = os.path.join(ROOT, "motors", "machbusterNew.ric")


class TestDemarSizing:
    # DeMar, "Starting Your Up-Goers" slides: M[g] = V[in^3] P[psi] / W[psi in^3/g]
    @pytest.mark.parametrize("V, P, g", [(2.36, 300.0, 0.14), (17.7, 200.0, 0.71),
                                         (103.7, 400.0, 8.3)])
    def test_slide_examples(self, V, P, g):
        # W = 5000 psi in^3/g (BPNV, the slides' pyrogen)
        assert demar_pyrogen_mass(V, 5000.0, P) * 1e3 == pytest.approx(g, rel=0.02)

    def test_default_pressure_is_300_psi_and_returns_kg(self):
        assert demar_pyrogen_mass(10.0, 5000.0) == pytest.approx(10.0 * 300.0 / 5000.0 * 1e-3)

    @pytest.mark.parametrize("W", [0.0, -1.0, None])
    def test_missing_impetus_raises(self, W):
        with pytest.raises(ValueError, match="impetus_W"):
            demar_pyrogen_mass(10.0, W)

    def test_sutton_kept(self):
        assert sutton_pyrogen_mass(10.0) == pytest.approx(0.12 * 10.0 ** 0.7 * 1e-3)

    def test_auto_charge_uses_demar_and_marks_chamber(self):
        motor = oma.load_ric(CHUNC)
        args = oma.ric_to_sim_args(motor, gas_props=None)
        geo, nozzle = args["geo"], args["nozzle"]
        pyro = oma.load_pyrogen("bpnv")
        ch = oma.build_pyrogen_chamber(pyro, geo, nozzle)
        case_v = np.pi / 4.0 * geo.D_outer ** 2 * geo.L_motor
        v_in3 = max(case_v - geo.total_propellant_volume(), 0.0) / 0.0254 ** 3
        assert ch.m_pyrogen_initial == pytest.approx(v_in3 * 300.0 / pyro.impetus_W * 1e-3)
        assert ch.mass_auto and ch.resize is not None
        big = ch.resize(1.25 * ch.m_pyrogen_initial)
        assert big.mass_auto and big.m_pyrogen_initial == pytest.approx(1.25 * ch.m_pyrogen_initial)
        # derived sizes follow the new mass
        assert big.V_plenum == pytest.approx(1.25 * ch.V_plenum)
        assert big.A_burn_initial == pytest.approx(1.25 * ch.A_burn_initial)
        # an explicit charge is not auto
        ex = oma.build_pyrogen_chamber(pyro, geo, nozzle, pyrogen_mass=1.0e-3)
        assert not ex.mass_auto and ex.resize is None


class TestDefaultsOff:
    def test_dataclass_default_is_none(self):
        f = {f.name: f.default for f in dataclasses.fields(Pyrogen)}
        assert f["heat_delivery_mode"] == "none"

    def test_mapping_default_is_none(self):
        d = dict(name="x", a=1e-5, n=0.5, rho=1700.0, T_flame=2700.0, M=0.04, gamma=1.2)
        assert oma._pyrogen_from_mapping(d).heat_delivery_mode == "none"

    def test_every_library_yaml_is_none(self):
        for fname, pyro in oma._library_pyrogen_files():
            assert pyro.heat_delivery_mode == "none", fname

    def test_every_library_pyrogen_has_an_impetus(self):
        for fname, pyro in oma._library_pyrogen_files():
            assert pyro.impetus_W > 0.0, fname

    def test_pre_ig3_embedded_blocks_upgrade_to_current(self):
        for name in ("bpnv", "mtv", "mtv_fast", "thermite_cual"):
            cur = oma._pyrogen_to_block(oma.load_pyrogen(name))
            old = dict(cur)
            old["heat_delivery_mode"] = ("radiation" if name == "thermite_cual" else "demar")
            block = {"pyrogen": old, "mass": -1.0}
            new, changed = oma.upgrade_igniter_block(block)
            assert changed, name
            assert new["pyrogen"]["heat_delivery_mode"] == "none"
            assert new["mass"] == -1.0
            # idempotent on the current entry
            assert oma.upgrade_igniter_block(new)[1] is False


def _chunc_run(**kw):
    kw.setdefault("t_max", 0.05)
    kw.setdefault("print_interval", 0.0)
    return oma.run_from_ric(CHUNC, verbose=False, **kw)[0]


class TestFallback:
    def test_auto_charge_is_upsized_until_it_lights(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = _chunc_run()
        s = res["summary"]
        assert s["ignited"] and s["pyrogen_mass_auto"] and s["pyrogen_resized"]
        assert s["ignition_attempts"] > 1
        m = s["ignition_attempt_masses_kg"]
        assert m[1] == pytest.approx(1.25 * m[0], rel=1e-12)
        assert s["pyrogen_mass_final_kg"] == pytest.approx(m[-1])
        assert s["pyrogen_mass_initial"] == pytest.approx(m[-1])   # final attempt reported
        msgs = [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]
        assert any("auto-sized" in t and "g" in t for t in msgs)

    def test_explicit_charge_never_resized(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = _chunc_run(pyrogen="bpnv", pyrogen_mass=0.2e-3, t_max=0.5)
        s = res["summary"]
        assert not s["ignited"] and not s["pyrogen_mass_auto"]
        assert s["ignition_attempts"] == 1 and not s["pyrogen_resized"]
        assert s["pyrogen_mass_initial"] == pytest.approx(0.2e-3)
        assert s["termination_code"] == 6                    # stopped early
        assert res["time"][-1] < 0.5
        msgs = [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]
        assert any("explicitly set" in t for t in msgs)

    def test_attempt_cap_raises(self):
        with pytest.raises(sim.PyrogenAutoSizeError, match="did not light"):
            _chunc_run(ignition_retry_max_attempts=2)

    def test_no_retry_arg_validation(self):
        with pytest.raises(ValueError, match="ignition_retry_factor"):
            _chunc_run(ignition_retry_factor=1.0)
        with pytest.raises(ValueError, match="ignition_retry_max_attempts"):
            _chunc_run(ignition_retry_max_attempts=0)

    def test_signature_exposes_defaults(self):
        import inspect
        p = inspect.signature(sim.run_simulation).parameters
        assert p["ignition_retry_factor"].default == 1.25
        assert p["ignition_retry_max_attempts"].default == 10
        assert p["stop_if_unignitable"].default is True
        assert p["ignition_heat_transfer"].default == "peretz"


class TestEarlyStop:
    def test_lit_run_identical_with_and_without_the_stop(self):
        kw = dict(pyrogen="bpnv", pyrogen_mass=3.0e-3, t_max=0.02)
        a = _chunc_run(stop_if_unignitable=True, **kw)
        b = _chunc_run(stop_if_unignitable=False, **kw)
        assert a["summary"]["grain_cells_ignited"] > 0 and b["summary"]["grain_cells_ignited"] > 0
        assert np.array_equal(np.asarray(a["P_head"]), np.asarray(b["P_head"]))

    def test_stop_off_runs_to_t_max(self):
        res = _chunc_run(pyrogen="bpnv", pyrogen_mass=0.2e-3, t_max=0.1,
                         stop_if_unignitable=False)
        assert res["summary"]["termination_code"] == 0
        assert not res["summary"]["ignited"]

    def test_stop_condition_is_physical(self):
        # at the stop: no cell lit, pyrogen consumed, gas colder than T_ignition everywhere
        res = _chunc_run(pyrogen="bpnv", pyrogen_mass=0.2e-3, t_max=0.5, snapshot_interval=1e-4)
        s = res["summary"]
        assert s["termination_code"] == 6 and s["pyrogen_mass_remaining"] <= 1e-12
        snap = res["snapshots"][-1]
        # grain-cell gas colder than T_ignition (the last snapshot is <= 0.1 ms before the stop;
        # grain gas cools ~0.1 K per 0.1 ms here)
        assert snap["T"][snap["is_grain"]].max() < 756.0 * 1.01
        assert not np.any(np.asarray(res["ignition_time_by_cell"]) < 1e8)


class TestPartialIgnition:
    BALLS = os.path.join(ROOT, "motors", "BALLSstick.ric")

    def test_partly_lit_grain_is_a_failed_ignition(self):
        # BALLSstick at 4.83 g (explicit): 10 head-end cells (10 %) stay unlit
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = oma.run_from_ric(self.BALLS, pyrogen="bpnv", pyrogen_mass=4.829e-3,
                                   t_max=0.05, print_interval=0.0, verbose=False)[0]
        s = res["summary"]
        assert s["grain_cells_ignited"] < 0.99 * s["grain_cells"]
        assert not s["ignited"]
        assert any("did not fully ignite" in str(x.message) for x in w)

    def test_min_lit_fraction_validated(self):
        with pytest.raises(ValueError, match="ignition_min_lit_fraction"):
            _chunc_run(ignition_min_lit_fraction=0.0)
