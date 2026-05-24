"""
tests/test_uniform_band_weights.py — v0.7.3 Phase A.1 kernel tests
====================================================================

Direct-kernel tests for `_compute_uniform_band_weights`, the top-hat
axial weighting that unifies submerged-igniter pyrogen distribution
(head_cartridge, aft_cartridge_zero_axial, aft_cartridge_fore_firing).

See srm_1d/docs/v0_7_2/candidates_post_phaseA.md (v0.7.3 design).
"""
import numpy as np
import pytest

from srm_1d.simulation import _compute_uniform_band_weights


def _uniform_grid(N, dx_uniform=0.01):
    dx = np.full(N, dx_uniform)
    return dx


# ---------------------------------------------------------------------
# Basic correctness
# ---------------------------------------------------------------------

def test_head_basket_first_3_of_10_cells():
    """i_start=0, i_end=2 → uniform weights across cells 0-2, zero
    elsewhere. Each gets 1/3 with uniform dx."""
    N = 10
    dx = _uniform_grid(N)
    w = np.zeros(N)
    _compute_uniform_band_weights(dx, 0, 2, N, w)
    expected = np.zeros(N)
    expected[0:3] = 1.0 / 3.0
    np.testing.assert_allclose(w, expected, atol=1e-12)
    assert np.isclose(w.sum(), 1.0, atol=1e-12)


def test_aft_cartridge_last_2_cells():
    """i_start=N-2, i_end=N-1 → uniform weights across cells N-2 and
    N-1, zero elsewhere. Each gets 0.5."""
    N = 8
    dx = _uniform_grid(N)
    w = np.zeros(N)
    _compute_uniform_band_weights(dx, N - 2, N - 1, N, w)
    expected = np.zeros(N)
    expected[N - 2:] = 0.5
    np.testing.assert_allclose(w, expected, atol=1e-12)
    assert np.isclose(w.sum(), 1.0, atol=1e-12)


def test_single_cell_range_collapses_to_cell_0_equivalent():
    """i_start=i_end=0 puts all weight at cell 0 (equivalent to the
    forward_plenum cell-0-only injection)."""
    N = 5
    dx = _uniform_grid(N)
    w = np.zeros(N)
    _compute_uniform_band_weights(dx, 0, 0, N, w)
    assert w[0] == 1.0
    assert np.all(w[1:] == 0.0)


def test_mid_bore_range():
    """Cells [3, 5] in a 10-cell motor get weight; others zero."""
    N = 10
    dx = _uniform_grid(N)
    w = np.full(N, 99.9)  # sentinel — kernel must overwrite
    _compute_uniform_band_weights(dx, 3, 5, N, w)
    expected = np.zeros(N)
    expected[3:6] = 1.0 / 3.0
    np.testing.assert_allclose(w, expected, atol=1e-12)


# ---------------------------------------------------------------------
# Conservation across parametrized geometries
# ---------------------------------------------------------------------

@pytest.mark.parametrize('lo,hi,N', [
    (0, 0, 10),
    (0, 4, 10),
    (3, 7, 10),
    (8, 9, 10),
    (0, 49, 50),
    (25, 49, 50),
])
def test_sum_to_one_uniform_grid(lo, hi, N):
    """Conservation guarantee across a range of injection windows."""
    dx = _uniform_grid(N)
    w = np.zeros(N)
    _compute_uniform_band_weights(dx, lo, hi, N, w)
    assert np.isclose(w.sum(), 1.0, atol=1e-10), (
        f"sum(w) != 1 for range [{lo}, {hi}] in N={N}"
    )


def test_sum_to_one_non_uniform_grid():
    """Non-uniform dx (cells linearly widen) — weights account for
    cell-width contribution: w[i] = dx[i] / sum_{i in range}(dx[i])."""
    N = 6
    dx = np.array([0.005, 0.010, 0.015, 0.020, 0.025, 0.030])
    w = np.zeros(N)
    _compute_uniform_band_weights(dx, 1, 4, N, w)
    total = dx[1] + dx[2] + dx[3] + dx[4]
    expected = np.zeros(N)
    for i in range(1, 5):
        expected[i] = dx[i] / total
    np.testing.assert_allclose(w, expected, atol=1e-12)
    assert np.isclose(w.sum(), 1.0, atol=1e-10)


# ---------------------------------------------------------------------
# Defensive edge cases
# ---------------------------------------------------------------------

def test_invalid_range_routes_to_cell_zero():
    """i_start > i_end is an invalid range; kernel falls back to
    cell-0-only injection (defensive)."""
    N = 10
    dx = _uniform_grid(N)
    w = np.zeros(N)
    _compute_uniform_band_weights(dx, 5, 3, N, w)
    assert w[0] == 1.0
    assert np.all(w[1:] == 0.0)


def test_out_of_bounds_indices_clamped():
    """i_start < 0 and i_end >= N get clamped to [0, N-1]."""
    N = 10
    dx = _uniform_grid(N)
    w = np.zeros(N)
    # Request [-3, 15] — should clamp to [0, 9] = all cells
    _compute_uniform_band_weights(dx, -3, 15, N, w)
    expected = np.full(N, 1.0 / N)
    np.testing.assert_allclose(w, expected, atol=1e-12)


def test_full_bore_range():
    """[0, N-1] gives uniform 1/N weight to every cell."""
    N = 20
    dx = _uniform_grid(N)
    w = np.zeros(N)
    _compute_uniform_band_weights(dx, 0, N - 1, N, w)
    expected = np.full(N, 1.0 / N)
    np.testing.assert_allclose(w, expected, atol=1e-12)
