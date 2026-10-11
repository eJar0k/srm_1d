"""
test_numerical_limits.py -- solver-verification PB (ports of Codex
d08bcc5 / d773b55 / e46562b / e749f71 / 86c587f, adapted to our tree).

(P1 hardening tests live in tests/test_hardening.py.)

P2 diagnostics: per-cell T-floor / T-ceiling / pressure-floor activation
counts, durations, first/last times, |clip energy| (T limits), per-cell max
floor deficit and optional port-Mach-cap face tallies. All BOOKKEEPING ONLY:
the tallies must never perturb the state (asserted bit-for-bit below).
"""
import numpy as np
import pytest

pytest.importorskip("numba")

from srm_1d import run_simulation  # noqa: E402
from srm_1d.solver import (  # noqa: E402
    LIMIT_PRESSURE_FLOOR,
    LIMIT_TEMPERATURE_CEILING,
    LIMIT_TEMPERATURE_FLOOR,
    N_LIMIT_DIAGNOSTICS,
    _piso_step_with_energy_diagnostics,
)
from srm_1d.tools.ignition_diagnostics import numerical_limits_summary  # noqa: E402
from tests.test_simulation_phase3 import _small_motor, _test_chamber  # noqa: E402

R_GAS = 300.0
CP_GAS = 2000.0
GAMMA = 1.2


def _limit_arrays(N):
    return dict(
        counts=np.zeros((N_LIMIT_DIAGNOSTICS, N), dtype=np.int64),
        dur=np.zeros((N_LIMIT_DIAGNOSTICS, N)),
        en=np.zeros((N_LIMIT_DIAGNOSTICS, N)),
        first=np.full((N_LIMIT_DIAGNOSTICS, N), -1.0),
        last=np.full((N_LIMIT_DIAGNOSTICS, N), -1.0),
        deficit=np.zeros(N),
    )


def _source_step(mass_rate, source_temperature, dt=1.0e-4, arrays=None,
                 step_start=0.0):
    """One closed cell with a mass/enthalpy source (T-clip driver)."""
    N = 1
    P0 = 1.0 * R_GAS * 300.0
    args = (
        np.array([1.0]), np.zeros(2), np.array([P0]), np.array([300.0]),
        np.array([1.0e-3]), np.array([0.035]),
        np.array([mass_rate]),
        np.array([mass_rate * source_temperature * CP_GAS]),
        np.zeros(2), np.zeros(1),
        0.01, dt, np.full(N, GAMMA), np.full(N, R_GAS), np.full(N, CP_GAS),
        np.full(N, 3000.0 * 1.01),
        0.0, P0, 300.0, N,
    )
    if arrays is None:
        return _piso_step_with_energy_diagnostics(*args)
    a = arrays
    return _piso_step_with_energy_diagnostics(
        *args, 0.0, a['counts'], a['dur'], a['en'], a['first'], a['last'],
        a['deficit'], None, None, None, None, step_start)


def _floor_step(pressure, floor_pa=None, arrays=None, step_start=0.0):
    """One venting cell whose raw post-correction pressure undershoots."""
    T = 300.0
    args = (
        np.array([pressure / (R_GAS * T)]), np.zeros(2), np.array([pressure]),
        np.array([T]), np.array([1.0e-3]), np.array([0.035]),
        np.zeros(1), np.zeros(1), np.zeros(2), np.zeros(1),
        0.01, 1.0e-5, np.array([GAMMA]), np.array([R_GAS]),
        np.array([CP_GAS]), np.array([3030.0]),
        0.0, 101325.0, T, 1,
    )
    kw = {}
    if floor_pa is not None:
        kw['pressure_floor_pa'] = floor_pa
    if arrays is None:
        return _piso_step_with_energy_diagnostics(*args, **kw)
    a = arrays
    return _piso_step_with_energy_diagnostics(
        *args, 0.0, a['counts'], a['dur'], a['en'], a['first'], a['last'],
        a['deficit'], None, None, None, None, step_start, **kw)


# ---------------------------------------------------------------- P2 (PISO)

class TestPisoLimitTallies:
    def test_temperature_ceiling_clip_is_localized(self):
        a = _limit_arrays(1)
        dt = 1.0e-4
        _source_step(100.0, 10_000.0, dt, a, step_start=0.5)
        assert a['counts'][LIMIT_TEMPERATURE_CEILING, 0] == 1
        assert a['dur'][LIMIT_TEMPERATURE_CEILING, 0] == pytest.approx(dt)
        assert a['en'][LIMIT_TEMPERATURE_CEILING, 0] > 0.0
        assert a['first'][LIMIT_TEMPERATURE_CEILING, 0] == 0.5
        assert a['last'][LIMIT_TEMPERATURE_CEILING, 0] == pytest.approx(0.5 + dt)
        assert a['counts'][LIMIT_TEMPERATURE_FLOOR, 0] == 0
        assert a['counts'][LIMIT_PRESSURE_FLOOR, 0] == 0

    def test_ceiling_energy_equals_abs_clip_power_times_dt(self):
        a = _limit_arrays(1)
        dt = 1.0e-4
        out = _source_step(100.0, 10_000.0, dt, a)
        clip_power = out[10]
        assert a['en'][LIMIT_TEMPERATURE_CEILING, 0] == pytest.approx(
            abs(clip_power) * dt, rel=1e-12)

    def test_temperature_floor_clip_is_localized(self):
        a = _limit_arrays(1)
        dt = 1.0e-4
        # A strong enthalpy SINK drives T_raw below T_ambient = 300 K.
        N = 1
        P0 = R_GAS * 300.0
        _piso_step_with_energy_diagnostics(
            np.array([1.0]), np.zeros(2), np.array([P0]), np.array([300.0]),
            np.array([1.0e-3]), np.array([0.035]),
            np.zeros(1), np.array([-5.0e6]), np.zeros(2), np.zeros(1),
            0.01, dt, np.full(N, GAMMA), np.full(N, R_GAS),
            np.full(N, CP_GAS), np.full(N, 3030.0),
            0.0, P0, 300.0, N, 0.0,
            a['counts'], a['dur'], a['en'], a['first'], a['last'],
            a['deficit'], None, None, None, None, 0.0)
        assert a['counts'][LIMIT_TEMPERATURE_FLOOR, 0] == 1
        assert a['en'][LIMIT_TEMPERATURE_FLOOR, 0] > 0.0
        assert a['counts'][LIMIT_TEMPERATURE_CEILING, 0] == 0

    def test_pressure_floor_override_records_pre_clamp_deficit(self):
        a = _limit_arrays(1)
        out = _floor_step(500.0, floor_pa=750.0, arrays=a)
        assert out[2][0] == pytest.approx(750.0)
        assert a['counts'][LIMIT_PRESSURE_FLOOR, 0] == 1
        assert a['deficit'][0] == pytest.approx(250.0)
        assert a['first'][LIMIT_PRESSURE_FLOOR, 0] == 0.0
        assert a['last'][LIMIT_PRESSURE_FLOOR, 0] == pytest.approx(1.0e-5)

    @pytest.mark.parametrize("which", ["ceiling", "floor"])
    def test_tallies_do_not_perturb_the_state(self, which):
        a = _limit_arrays(1)
        if which == "ceiling":
            bare = _source_step(100.0, 10_000.0)
            tallied = _source_step(100.0, 10_000.0, arrays=a)
        else:
            bare = _floor_step(500.0)
            tallied = _floor_step(500.0, arrays=a)
        assert len(bare) == len(tallied)
        for x, y in zip(bare, tallied):
            assert np.array_equal(np.asarray(x), np.asarray(y))
        assert a['counts'].sum() > 0

    def test_mach_cap_face_tallies(self):
        N = 3
        P = np.full(N, 2.0e5)
        T = np.full(N, 300.0)
        rho = P / (R_GAS * T)
        u = np.zeros(N + 1)
        A = np.full(N, 1.0e-3)
        cnt = np.zeros(N - 1, dtype=np.int64)
        dur = np.zeros(N - 1)
        first = np.full(N - 1, -1.0)
        last = np.full(N - 1, -1.0)
        a = _limit_arrays(N)
        dt = 1.0e-5
        _piso_step_with_energy_diagnostics(
            rho, u, P, T, A, np.full(N, 0.035),
            np.full(N, 0.5), np.full(N, 0.5 * 1500.0 * CP_GAS),
            np.zeros(N + 1), np.zeros(N),
            0.01, dt, np.full(N, GAMMA), np.full(N, R_GAS),
            np.full(N, CP_GAS), np.full(N, 3030.0),
            1.0e-3, 101325.0, 300.0, N,
            1.0e-6,  # absurdly low Mach cap -> every moving face clamps
            a['counts'], a['dur'], a['en'], a['first'], a['last'],
            a['deficit'], cnt, dur, first, last, 0.25)
        active = cnt > 0
        assert np.any(active)
        assert np.all(first[active] == 0.25)
        assert np.all(last[active] == pytest.approx(0.25 + dt))
        assert np.all(dur[active] == pytest.approx(dt))
        assert np.all(first[~active] == -1.0)


# ---------------------------------------------------------------- P2 (run)

def _short_run(**kw):
    geo, prop, nozzle = _small_motor()
    args = dict(T_ignition=294.0, t_max=0.001, P_cutoff=1.0,
                snapshot_interval=0.001, verbose=False)
    args.update(kw)
    return geo, run_simulation(geo, prop, nozzle, _test_chamber(), **args)


class TestRunLevel:
    def test_result_exposes_per_cell_limit_maps(self):
        geo, result = _short_run()
        lim = result['numerical_limits']
        N = geo.N_cells
        assert lim['x_m'].shape == (N,)
        for name in ('temperature_floor', 'temperature_ceiling',
                     'pressure_floor'):
            blk = lim[name]
            for key in ('activation_count_by_cell', 'duration_s_by_cell',
                        'first_activation_time_s_by_cell',
                        'last_activation_time_s_by_cell'):
                assert blk[key].shape == (N,)
            active = blk['activation_count_by_cell'] > 0
            assert np.all(np.isfinite(
                blk['first_activation_time_s_by_cell'][active]))
            assert np.all(np.isnan(
                blk['first_activation_time_s_by_cell'][~active]))
            assert np.all(blk['last_activation_time_s_by_cell'][active]
                          >= blk['first_activation_time_s_by_cell'][active])
        assert lim['temperature_floor']['threshold_k'] == pytest.approx(
            max(1.0, 294.0), rel=0.2)
        assert lim['pressure_floor']['threshold_pa'] == 1.0e3
        assert lim['pressure_floor']['maximum_raw_deficit_pa_by_cell'].shape == (N,)
        assert np.all(
            lim['pressure_floor']['maximum_raw_deficit_pa_by_cell'] >= 0.0)
        assert result['summary']['pressure_floor_pa'] == 1.0e3
        m = lim['port_mach_cap']
        assert m['threshold_mach'] == 0.0
        assert m['x_m'].shape == (N - 1,)
        assert np.count_nonzero(m['activation_count_by_face']) == 0
        assert np.all(np.isnan(m['first_activation_time_s_by_face']))

    def test_enabled_port_mach_cap_reports_face_activations(self):
        _geo, result = _short_run(port_mach_cap=1.0e-6)
        m = result['numerical_limits']['port_mach_cap']
        active = m['activation_count_by_face'] > 0
        assert np.any(active)
        assert np.all(np.isfinite(m['first_activation_time_s_by_face'][active]))
        assert np.all(m['last_activation_time_s_by_face'][active]
                      >= m['first_activation_time_s_by_face'][active])

    def test_pressure_floor_override_is_reported(self):
        _geo, result = _short_run(pressure_floor_pa=5.0e3)
        assert result['summary']['pressure_floor_pa'] == 5.0e3
        assert result['numerical_limits']['pressure_floor']['threshold_pa'] == 5.0e3

    def test_summary_helper_is_consistent_with_maps(self):
        _geo, result = _short_run()
        s = numerical_limits_summary(result)
        lim = result['numerical_limits']
        for name in ('temperature_floor', 'temperature_ceiling',
                     'pressure_floor'):
            assert s[name]['total_activations'] == int(
                lim[name]['activation_count_by_cell'].sum())
            assert s[name]['active_cells'] == int(np.count_nonzero(
                lim[name]['activation_count_by_cell']))
        assert s['port_mach_cap']['active_faces'] == 0


class TestChuncLimits:
    """Canonical Chunc, 50 ms, 406 cells: the pressure floor is known to
    fire (T3: 26 cell-steps) so the per-cell maps can be tied to the T3
    run-level counters and the T6 clip-power history."""

    @pytest.fixture(scope="class")
    def run(self):
        from srm_1d.openmotor_adapter import run_from_ric
        res = run_from_ric(
            "motors/machbusterNew.ric", roughness=32e-6, kappa=0.44,
            pyrogen="bpnv", pyrogen_mass=None, T_ignition=756.0,
            t_max=0.05, target_propellant_cells=400, print_interval=0.0,
            verbose=False,
            # T11: the new PISO defaults (K=1, S=True) remove the floor
            # events in this regime (P5/T11); this test needs them to tie the
            # per-cell map to the T3 counters, so it runs the pre-T11 step.
            piso_energy_coupling_iters=0, piso_corr2_storage=False)[0]
        return res

    def test_pressure_floor_map_ties_to_t3_counters(self, run):
        lim = run['numerical_limits']['pressure_floor']
        summ = run['summary']
        assert summ['pressure_floor_activations'] > 0
        assert int(lim['activation_count_by_cell'].sum()) == (
            summ['pressure_floor_activations'])
        first = np.nanmin(lim['first_activation_time_s_by_cell'])
        assert first == pytest.approx(summ['pressure_floor_first_time_s'])
        assert lim['maximum_raw_deficit_pa_by_cell'].max() == pytest.approx(
            1.0e3 - summ['pressure_floor_min_raw_P'], rel=1e-12)

    def test_t_clip_energy_bounds_the_net_clip_power_history(self, run):
        lim = run['numerical_limits']
        total_abs = (
            lim['temperature_floor']['absolute_correction_energy_j_by_cell'].sum()
            + lim['temperature_ceiling'][
                'absolute_correction_energy_j_by_cell'].sum())
        net = abs(np.sum(run['clipping_correction_power'] * run['dt']))
        assert total_abs > 0.0
        assert total_abs >= net * (1.0 - 1e-9)
