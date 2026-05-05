"""Tests for the PISO solver numerical building blocks."""
import numpy as np
import pytest
from srm_1d.solver import thomas_solve, compute_dt_cfl


class TestThomasSolve:
    def test_identity_system(self):
        """b=1, a=c=0: solution equals RHS."""
        N = 5
        x = thomas_solve(np.zeros(N), np.ones(N), np.zeros(N),
                         np.array([1., 2., 3., 4., 5.]), N)
        np.testing.assert_allclose(x, [1., 2., 3., 4., 5.])

    def test_known_tridiagonal(self):
        """Solve a known 4x4 system and check against numpy.linalg.solve."""
        N = 4
        a = np.array([0.0, -1.0, -1.0, -1.0])
        b = np.array([4.0, 4.0, 4.0, 4.0])
        c = np.array([-1.0, -1.0, -1.0, 0.0])
        d = np.array([1.0, 2.0, 3.0, 4.0])

        # Build full matrix for reference
        A = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
        x_ref = np.linalg.solve(A, d)
        x_tdma = thomas_solve(a, b, c, d, N)
        np.testing.assert_allclose(x_tdma, x_ref, atol=1e-12)

    def test_large_diagonally_dominant(self):
        """100-element diagonally dominant system."""
        N = 100
        a = np.full(N, -1.0); a[0] = 0.0
        b = np.full(N, 3.0)
        c = np.full(N, -1.0); c[-1] = 0.0
        d = np.ones(N)

        A = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
        x_ref = np.linalg.solve(A, d)
        x_tdma = thomas_solve(a, b, c, d, N)
        np.testing.assert_allclose(x_tdma, x_ref, atol=1e-10)


class TestCFL:
    def test_basic_cfl(self):
        """CFL condition should scale dt inversely with wave speed."""
        u = np.array([10.0, 20.0, 30.0])
        dt = compute_dt_cfl(u, 340.0, 0.01, 3, 0.5, 1.0)
        # max wave speed = 30 + 340 = 370
        expected = 0.5 * 0.01 / 370.0
        assert dt == pytest.approx(expected, rel=1e-6)

    def test_dt_max_cap(self):
        """dt should not exceed dt_max."""
        u = np.array([0.01])
        dt = compute_dt_cfl(u, 340.0, 0.01, 1, 0.5, 1e-6)
        assert dt == pytest.approx(1e-6)

    def test_zero_velocity(self):
        """With zero flow, dt is set by sound speed alone."""
        u = np.zeros(5)
        dt = compute_dt_cfl(u, 340.0, 0.01, 5, 0.5, 1.0)
        expected = 0.5 * 0.01 / 340.0
        assert dt == pytest.approx(expected, rel=1e-6)
