"""
test_refresh_cadence.py -- split burn-rate / geometry refresh cadence.

Before solver-verification T5 a single ``burn_update_interval`` (default
``max(10, N // 5)``) gated BOTH ``update_cell_geometry`` and
``compute_burn_rates`` in ``simulation._run_time_loop``.  T5 splits the two:
``run_simulation(geometry_update_interval=None)`` resolves to the burn
interval (so defaults are bit-identical to before -- gate B0 enforces that on
the three canonical motors), and the summary now reports the intervals plus
the largest simulated-time gap between consecutive refreshes of each kind.

Everything here is a short Chunc run (``t_max = 0.05`` s, default 100-cell
grid -> 106 actual cells -> default refresh interval 21 steps).  The
``time`` history is the step-START time, and refreshes happen on steps
``s`` with ``s % interval == 0`` (``hist_idx == step``), so the expected gap
between refresh steps ``s1 < s2`` is exactly ``time[s2] - time[s1]`` -- the
tests recompute the gaps that way, independently of the loop's bookkeeping.
"""
import numpy as np
import pytest

pytest.importorskip("numba")

CANON = dict(
    roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
    T_ignition=756.0, t_max=0.05, print_interval=0.0, verbose=False,
)

SUMMARY_KEYS = (
    "burn_refresh_interval_steps",
    "geometry_refresh_interval_steps",
    "max_burn_refresh_gap_s",
    "max_geometry_refresh_gap_s",
)


def _chunc(**overrides):
    from srm_1d.openmotor_adapter import run_from_ric

    args = dict(CANON)
    args.update(overrides)
    result, _perf, _nz, _geo, _prop = run_from_ric(
        "motors/machbusterNew.ric", **args)
    return result


def _expected_max_gap(time, interval):
    """Max of time[s2]-time[s1] over consecutive refresh steps (s % k == 0)."""
    idx = np.arange(0, len(time), interval)
    if len(idx) < 2:
        return 0.0
    return float(np.max(np.diff(time[idx])))


@pytest.fixture(scope="module")
def default_run():
    return _chunc()


@pytest.fixture(scope="module")
def default_interval(default_run):
    return default_run["summary"]["burn_refresh_interval_steps"]


@pytest.fixture(scope="module")
def every_step_run():
    return _chunc(burn_update_interval=1, geometry_update_interval=1)


@pytest.fixture(scope="module")
def geometry_only_run(default_interval):
    """Burn rates on the default cadence, geometry EVERY step."""
    return _chunc(burn_update_interval=default_interval,
                  geometry_update_interval=1)


@pytest.fixture(scope="module")
def burn_only_run(default_interval):
    """Burn rates EVERY step, geometry on the default cadence."""
    return _chunc(burn_update_interval=1,
                  geometry_update_interval=default_interval)


class TestSummaryKeysAndDefaults:
    def test_keys_present(self, default_run):
        summary = default_run["summary"]
        for key in SUMMARY_KEYS:
            assert key in summary, key

    def test_default_intervals_unchanged_and_equal(self, default_run):
        """geometry_update_interval=None resolves to the burn interval, which
        is still the pre-T5 auto value max(10, N // 5)."""
        s = default_run["summary"]
        assert s["burn_refresh_interval_steps"] == max(10, s["cells"] // 5)
        assert s["geometry_refresh_interval_steps"] == s["burn_refresh_interval_steps"]
        assert isinstance(s["burn_refresh_interval_steps"], int)
        assert isinstance(s["geometry_refresh_interval_steps"], int)

    def test_default_gaps_match_time_history(self, default_run):
        s = default_run["summary"]
        k = s["burn_refresh_interval_steps"]
        expected = _expected_max_gap(np.asarray(default_run["time"]), k)
        assert expected > 0.0
        assert s["max_burn_refresh_gap_s"] == pytest.approx(expected, rel=1e-12)
        assert s["max_geometry_refresh_gap_s"] == pytest.approx(expected, rel=1e-12)

    def test_gap_bounded_by_interval_times_dt_extremes(self, default_run):
        """A gap is the sum of `interval` consecutive step sizes."""
        s = default_run["summary"]
        k = s["burn_refresh_interval_steps"]
        dt = np.asarray(default_run["dt"])
        gap = s["max_burn_refresh_gap_s"]
        assert k * float(np.min(dt)) <= gap * (1 + 1e-12)
        assert gap <= k * float(np.max(dt)) * (1 + 1e-12)


class TestNoneResolvesToBurnInterval:
    def test_explicit_equal_interval_is_bit_identical(self, default_run,
                                                      default_interval):
        """geometry_update_interval=None  ==  geometry_update_interval=<the
        resolved burn interval>, bit for bit."""
        explicit = _chunc(geometry_update_interval=default_interval)
        a, b = default_run, explicit
        assert a["summary"]["steps"] == b["summary"]["steps"]
        assert np.array_equal(a["time"], b["time"])
        assert np.array_equal(a["P_head"], b["P_head"])
        assert np.array_equal(a["min_pressure"], b["min_pressure"])
        for key in SUMMARY_KEYS:
            assert a["summary"][key] == b["summary"][key], key

    def test_none_follows_an_explicit_burn_interval(self):
        """With only burn_update_interval given, geometry follows it (the
        pre-T5 behaviour of a single shared cadence)."""
        s = _chunc(burn_update_interval=7, t_max=0.01)["summary"]
        assert s["burn_refresh_interval_steps"] == 7
        assert s["geometry_refresh_interval_steps"] == 7


class TestSplitIsWired:
    def test_geometry_only_every_step(self, geometry_only_run, default_interval,
                                      default_run):
        r = geometry_only_run
        s = r["summary"]
        time = np.asarray(r["time"])
        assert s["burn_refresh_interval_steps"] == default_interval
        assert s["geometry_refresh_interval_steps"] == 1
        # Every-step geometry gap == the largest single step the run took
        # (between refresh steps s and s+1 the gap is exactly dt_s).
        assert s["max_geometry_refresh_gap_s"] == pytest.approx(
            float(np.max(np.diff(time))), rel=1e-12)
        assert s["max_burn_refresh_gap_s"] == pytest.approx(
            _expected_max_gap(time, default_interval), rel=1e-12)
        # ...and the kwarg is really threaded: the trajectory moved.
        m = min(len(r["P_head"]), len(default_run["P_head"]))
        assert not np.array_equal(np.asarray(r["P_head"])[:m],
                                  np.asarray(default_run["P_head"])[:m])

    def test_burn_only_every_step(self, burn_only_run, default_interval,
                                  default_run):
        r = burn_only_run
        s = r["summary"]
        time = np.asarray(r["time"])
        assert s["burn_refresh_interval_steps"] == 1
        assert s["geometry_refresh_interval_steps"] == default_interval
        assert s["max_burn_refresh_gap_s"] == pytest.approx(
            float(np.max(np.diff(time))), rel=1e-12)
        assert s["max_geometry_refresh_gap_s"] == pytest.approx(
            _expected_max_gap(time, default_interval), rel=1e-12)
        m = min(len(r["P_head"]), len(default_run["P_head"]))
        assert not np.array_equal(np.asarray(r["P_head"])[:m],
                                  np.asarray(default_run["P_head"])[:m])

    def test_every_step_gap_is_one_step_and_smaller_than_default(
            self, every_step_run, default_run):
        s = every_step_run["summary"]
        time = np.asarray(every_step_run["time"])
        assert s["burn_refresh_interval_steps"] == 1
        assert s["geometry_refresh_interval_steps"] == 1
        one_step = float(np.max(np.diff(time)))
        assert s["max_burn_refresh_gap_s"] == pytest.approx(one_step, rel=1e-12)
        assert s["max_geometry_refresh_gap_s"] == pytest.approx(one_step, rel=1e-12)
        # a 21-step cadence stales the burn rate ~21x longer than every-step
        assert default_run["summary"]["max_burn_refresh_gap_s"] > 10.0 * one_step


class TestEveryStepVsDefault:
    def test_peak_within_half_percent(self, default_run, every_step_run):
        """Refreshing burn rate AND geometry every step moves the 50-ms Chunc
        peak head pressure by < 0.5 % (measured: see worklog/T5.md)."""
        p_def = float(np.max(default_run["P_head"]))
        p_all = float(np.max(every_step_run["P_head"]))
        min_def = float(np.min(default_run["min_pressure"]))
        min_all = float(np.min(every_step_run["min_pressure"]))
        print(f"\n[T5] peak P_head default={p_def:.6e} every-step={p_all:.6e} "
              f"rel={(p_all - p_def) / p_def:+.3e}; global min P "
              f"default={min_def:.6e} every-step={min_all:.6e}")
        assert abs(p_all - p_def) / p_def < 0.005


class TestValidation:
    @pytest.mark.parametrize("bad", [0, -1, -21, 2.5])
    def test_bad_geometry_interval_rejected(self, bad):
        with pytest.raises(ValueError, match="geometry_update_interval"):
            _chunc(geometry_update_interval=bad, t_max=0.01)
