"""
test_pressure_floor.py -- accounting for the 1 kPa PISO pressure floor.

``solver._piso_step_with_energy_diagnostics`` STEP 4 clamps ``P_new`` to a
1 kPa floor and then rebuilds ``rho_new = P/(R*T)``.  That clamp injects mass
and (sensible) energy outside every ledger.  Solver-verification task T3 adds
ACCOUNTING ONLY -- counters appended to the PISO return tuple and threaded
into per-step histories + summary keys.  The floor value and behaviour are
unchanged (gate B0 enforces bit-identical default trajectories).

Part (b) forces a clamp in a tiny-N direct call; part (a) checks the keys /
history lengths / internal consistency on a short real Chunc run.
"""
import math

import numpy as np
import pytest

pytest.importorskip("numba")

from srm_1d.solver import _piso_step_with_energy_diagnostics  # noqa: E402

PRESSURE_FLOOR = 1.0e3  # Pa -- pinned: T3 must not change the floor value

R_GAS = 300.0
CP_GAS = 2000.0
GAMMA = 1.2


def _vent_step(N, P_initial, A_port_value=1.0e-3, dx=0.01, dt=1.0e-4,
               P_ambient=1.0, **piso_kw):
    """One PISO step of an N-cell duct venting through a throat the size of
    the port into a near-vacuum back pressure (no sources, no friction).

    A deliberately over-large step (dt = 1e-4 s against dx/a ~ 3e-5 s) lets
    the choked nozzle drain the 1.5 kPa duct within one step, so the raw
    post-correction pressure falls below the 1 kPa floor. (Measured, N=1,
    P0=1500 Pa: raw P_new = 172.9 Pa -> clamped to 1000 Pa.)
    """
    P = np.asarray(P_initial, dtype=float).copy()
    T = np.full(N, 300.0)
    rho = P / (R_GAS * T)
    u = np.zeros(N + 1)
    A_port = np.full(N, A_port_value)
    D_hyd = np.full(N, 0.035)
    zeros_N = np.zeros(N)
    return _piso_step_with_energy_diagnostics(
        rho, u, P, T, A_port, D_hyd,
        zeros_N, zeros_N.copy(), np.zeros(N + 1), zeros_N.copy(),
        dx, dt,
        np.full(N, GAMMA), np.full(N, R_GAS), np.full(N, CP_GAS),
        np.full(N, 3000.0 * 1.01),
        A_port_value, P_ambient, 300.0, N, **piso_kw,
    )


class TestPressureFloorCounters:
    """(b) Direct PISO call that forces the clamp."""

    def test_return_tuple_has_four_trailing_floor_values(self):
        out = _vent_step(1, [1500.0])
        # 12 pre-T3 values (rho, u, P, T, 3 energy, 3 power, clip, residual)
        # + (n_floor_cells, floor_mass_added, floor_energy_added,
        #    floor_min_raw_P)
        assert len(out) == 24

    def test_single_cell_clamp_is_counted_exactly(self):
        N = 1
        A_port_value, dx = 1.0e-3, 0.01
        out = _vent_step(N, [1500.0], A_port_value=A_port_value, dx=dx)
        rho_new, _u, P_new, T_new = out[0], out[1], out[2], out[3]
        n_cells, m_add, e_add, p_raw = out[12], out[13], out[14], out[15]

        # The clamp fired and left the carried state exactly at the floor,
        # with density rebuilt from the clamped pressure (unchanged behaviour).
        assert P_new[0] == PRESSURE_FLOOR
        assert rho_new[0] == pytest.approx(
            PRESSURE_FLOOR / (R_GAS * T_new[0]), rel=1e-14)

        assert n_cells == 1
        assert p_raw < PRESSURE_FLOOR
        expected_mass = (PRESSURE_FLOOR - p_raw) / (R_GAS * T_new[0]) \
            * A_port_value * dx
        assert m_add > 0.0
        assert m_add == pytest.approx(expected_mass, rel=1e-12)
        assert e_add == pytest.approx(m_add * CP_GAS * T_new[0], rel=1e-12)

    @pytest.mark.parametrize("P_initial, expected_cells", [
        ([1500.0, 1500.0, 1500.0], 3),   # every cell clamps
        ([3000.0, 1500.0, 1500.0], 2),   # mixed: the nozzle-side cell does not
    ])
    def test_multi_cell_count_matches_clamped_cells(self, P_initial,
                                                    expected_cells):
        N = 3
        # T11: the new PISO defaults (K=1, S=True) no longer drive this
        # over-large-dt duct below the floor in every cell (1 of 3 clamps);
        # this test pins the counter-vs-clamp identity, so it keeps exercising
        # the pre-T11 step (explicit K=0, S=False), where all 3 cells clamp.
        out = _vent_step(N, P_initial, energy_coupling_iters=0,
                         corr2_storage=False)
        P_new = out[2]
        n_cells, m_add, e_add, p_raw = out[12], out[13], out[14], out[15]

        # A cell ends exactly at the floor iff it was clamped.
        n_at_floor = int(np.sum(P_new == PRESSURE_FLOOR))
        assert n_at_floor >= 1
        assert n_cells == n_at_floor == expected_cells  # measured counts
        assert m_add > 0.0 and e_add > 0.0
        assert p_raw < PRESSURE_FLOOR
        assert np.all(P_new >= PRESSURE_FLOOR)

    def test_no_clamp_leaves_counters_zero(self):
        out = _vent_step(1, [1.0e5], dt=1.0e-6)
        P_new = out[2]
        n_cells, m_add, e_add, p_raw = out[12], out[13], out[14], out[15]

        assert P_new[0] > PRESSURE_FLOOR
        assert n_cells == 0
        assert m_add == 0.0
        assert e_add == 0.0
        # "no clamp" is a finite sentinel (fastmath assumes no infs), far above
        # any physical pressure so a running min over clamped cells ignores it.
        assert math.isfinite(p_raw)
        assert p_raw > 1.0e100


# ---------------------------------------------------------------------------
# (a) Short real run
# ---------------------------------------------------------------------------

SUMMARY_KEYS = (
    "pressure_floor_activations",
    "pressure_floor_steps",
    "pressure_floor_mass_added_kg",
    "pressure_floor_energy_added_J",
    "pressure_floor_first_time_s",
    "pressure_floor_min_raw_P",
)


@pytest.fixture(scope="module")
def chunc_short():
    from srm_1d.openmotor_adapter import run_from_ric

    result, _perf, _nz, _geo, _prop = run_from_ric(
        "motors/machbusterNew.ric",
        roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
        T_ignition=756.0, t_max=0.05, print_interval=0.0, verbose=False,
    )
    return result


class TestPressureFloorRun:
    def test_summary_keys_and_history_lengths(self, chunc_short):
        summary = chunc_short["summary"]
        for key in SUMMARY_KEYS:
            assert key in summary, key
        n = len(chunc_short["time"])
        assert n > 0
        assert len(chunc_short["pressure_floor_cells"]) == n
        assert len(chunc_short["pressure_floor_mass"]) == n

    def test_summary_matches_histories(self, chunc_short):
        summary = chunc_short["summary"]
        cells = np.asarray(chunc_short["pressure_floor_cells"])
        mass = np.asarray(chunc_short["pressure_floor_mass"])
        time = np.asarray(chunc_short["time"])

        assert np.all(cells >= 0.0)
        assert np.all(mass >= 0.0)
        assert summary["pressure_floor_activations"] == int(round(cells.sum()))
        assert summary["pressure_floor_steps"] == int(np.count_nonzero(cells))
        assert summary["pressure_floor_mass_added_kg"] == pytest.approx(
            float(mass.sum()), rel=1e-12, abs=0.0)

    def test_zero_or_consistent(self, chunc_short):
        """Zero clamps -> every counter is its 'never' value; otherwise the
        counters are mutually consistent (activations >= steps >= 1, mass and
        energy >= 0, raw P below the floor, first time inside the run)."""
        s = chunc_short["summary"]
        time = np.asarray(chunc_short["time"])
        if s["pressure_floor_activations"] == 0:
            assert s["pressure_floor_steps"] == 0
            assert s["pressure_floor_mass_added_kg"] == 0.0
            assert s["pressure_floor_energy_added_J"] == 0.0
            assert s["pressure_floor_first_time_s"] == -1.0
            assert math.isnan(s["pressure_floor_min_raw_P"])
        else:
            assert s["pressure_floor_activations"] >= s["pressure_floor_steps"] >= 1
            assert s["pressure_floor_mass_added_kg"] >= 0.0
            assert s["pressure_floor_energy_added_J"] >= 0.0
            assert s["pressure_floor_min_raw_P"] < PRESSURE_FLOOR
            assert time[0] <= s["pressure_floor_first_time_s"] <= time[-1]

    def test_carried_pressure_never_below_floor(self, chunc_short):
        """The clamp is applied before the state is carried, so the recorded
        minimum cell pressure can never dip under the floor (pins the floor
        value/behaviour this task must not change)."""
        assert float(np.min(chunc_short["min_pressure"])) >= PRESSURE_FLOOR
