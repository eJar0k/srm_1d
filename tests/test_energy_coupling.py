"""
test_energy_coupling.py -- P5 opt-in energy-coupled PISO iteration.

``piso_energy_coupling_iters`` (K) re-runs the pressure correctors + energy
update with the T-updated compressibility; ``piso_corr2_storage`` (S) adds the
storage term the original second corrector omits.  K = 0, S = False must be the
original step bit-for-bit (the B0 gate also enforces this on full runs).
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d.solver import _piso_step_with_energy_diagnostics as piso  # noqa: E402

N, DX, A = 20, 0.02, 1.0e-3
R_GAS, CP, GAMMA = 300.0, 1200.0, 1.25


def _step(K=0, S=False, T_src=3000.0, dt=1.0e-6, tol=1.0e-9, extra=True):
    P = np.full(N, 2.0e5)
    T = np.full(N, 600.0)
    rho = P / (R_GAS * T)
    u = np.zeros(N + 1)
    src = np.zeros(N)
    src[0:3] = 5.0
    args = (rho, u, P, T, np.full(N, A), np.full(N, 0.035), src, src * CP * T_src,
            np.zeros(N + 1), np.full(N, 0.02), DX, dt, np.full(N, GAMMA),
            np.full(N, R_GAS), np.full(N, CP), np.full(N, 3500.0), 1.0e-4,
            101325.0, 300.0, N, 0.0)
    out = piso(*args, K, tol, S) if extra else piso(*args)
    mass_before = out[16]
    carried_rate = (np.sum(out[0]) * A * DX - mass_before) / dt
    expected_rate = np.sum(src) * DX - out[20]
    return out, carried_rate, expected_rate


class TestK0Identity:
    def test_defaults_equal_explicit_zero(self):
        a, *_ = _step(extra=False)
        b, *_ = _step(K=0, S=False)
        for x, y in zip(a[:21], b[:21]):
            assert np.array_equal(np.asarray(x), np.asarray(y))
        assert b[21] == 0  # passes

    def test_tuple_has_p5_outputs(self):
        out, *_ = _step(K=2)
        assert len(out) == 24


class TestEosJumpShrinks:
    def test_original_step_double_counts_storage(self):
        # no T change at all: carried mass rate is ~2x the true rate (P'2
        # re-applies the storage the first corrector took)
        _, carried, expected = _step(T_src=600.0)
        assert carried == pytest.approx(2.0 * expected, rel=1e-2)

    def test_literal_iteration_removes_T_lag_only(self):
        out0, *_ = _step(K=0, S=False)
        out3, *_ = _step(K=3, S=False)
        lag0 = abs(out0[18] - out0[23])   # eos jump minus its P'2 part
        lag3 = abs(out3[18] - out3[23])
        assert lag3 < 1e-3 * lag0

    def test_full_iteration_closes_mass_and_eos_jump(self):
        out, carried, expected = _step(K=3, S=True)
        assert carried == pytest.approx(expected, rel=1e-3)
        assert abs(out[18]) < 1e-6 * abs(out[17])
        assert out[22] < 1e-6  # converged relative T change

    def test_K0_with_storage_term_only_fixes_double_count(self):
        _, carried, expected = _step(K=0, S=True, T_src=600.0)
        assert carried == pytest.approx(expected, rel=1e-3)


@pytest.fixture(scope="module")
def short_runs():
    from srm_1d.openmotor_adapter import run_from_ric
    ric = os.path.join(os.path.dirname(__file__), "..", "motors", "machbusterNew.ric")
    out = {}
    for K, S in ((0, False), (3, True)):
        out[(K, S)] = run_from_ric(
            ric, target_propellant_cells=100, t_max=0.01, roughness=32e-6,
            kappa=0.44, T_ignition=756.0, print_interval=0.0, verbose=False,
            pyrogen='bpnv', pyrogen_mass=None,
            piso_energy_coupling_iters=K, piso_corr2_storage=S)[0]
    return out


class TestShortRun:
    def test_eos_jump_shrinks_on_short_run(self, short_runs):
        j0 = np.abs(np.asarray(short_runs[(0, False)]['storage_jump_eos_J'])[:-1]).sum()
        j3 = np.abs(np.asarray(short_runs[(3, True)]['storage_jump_eos_J'])[:-1]).sum()
        assert j3 < 0.05 * j0

    def test_pass_counters_recorded(self, short_runs):
        r0, r3 = short_runs[(0, False)], short_runs[(3, True)]
        assert np.all(np.asarray(r0['piso_coupling_passes']) == 0)
        p3 = np.asarray(r3['piso_coupling_passes'])
        assert p3.max() <= 3 and p3.mean() > 0.5
        assert np.all(np.isfinite(r3['piso_coupling_rel_dT']))
