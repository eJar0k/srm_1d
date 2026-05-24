"""
tests/test_pyrogen_axial_weights.py — v0.7.2 Phase A.1 kernel tests
====================================================================

Direct-kernel tests for `_compute_pyrogen_axial_weights`, the
exponential-decay axial weighting that splits pyrogen injection across
multiple bore cells. These tests gate the Phase A regression criterion
(L_jet=0 recovers v0.7.1 behavior byte-for-byte), the conservation
guarantee (sum == 1), and the functional form (matches a hand
calculation for a known geometry).

See srm_1d/docs/v0_7_2/candidates/03_pyrogen_spatial_distribution.md.
"""

import numpy as np
import pytest

from srm_1d.simulation import _compute_pyrogen_axial_weights


def _uniform_grid(N, dx_uniform=0.01):
    """Build a uniform grid: N cells of width dx_uniform.

    x_centers[i] = (i + 0.5) * dx_uniform.
    """
    dx = np.full(N, dx_uniform)
    x_centers = (np.arange(N) + 0.5) * dx_uniform
    return x_centers, dx


# ---------------------------------------------------------------------
# Regression gate: L_jet = 0 recovers v0.7.1 cell-0-only behavior
# ---------------------------------------------------------------------

def test_L_jet_zero_recovers_cell_zero_only():
    """The Phase A regression gate.

    L_jet <= 0 must put all weight in cell 0, byte-for-byte recovering
    v0.7.1's cell-0-only pyrogen injection.
    """
    N = 50
    x_centers, dx = _uniform_grid(N)

    w = _compute_pyrogen_axial_weights(x_centers, dx, 0.0, N)

    assert w[0] == 1.0
    assert np.all(w[1:] == 0.0)
    assert np.isclose(w.sum(), 1.0)


def test_L_jet_negative_recovers_cell_zero_only():
    """Defensive: negative L_jet also routes to cell 0."""
    N = 20
    x_centers, dx = _uniform_grid(N)

    w = _compute_pyrogen_axial_weights(x_centers, dx, -1.0, N)

    assert w[0] == 1.0
    assert np.all(w[1:] == 0.0)


# ---------------------------------------------------------------------
# Conservation: weights sum to exactly 1
# ---------------------------------------------------------------------

@pytest.mark.parametrize('L_jet', [0.005, 0.01, 0.05, 0.1, 1.0])
def test_weights_sum_to_one_uniform_grid(L_jet):
    """Conservation guarantee across a range of L_jet on a uniform grid."""
    N = 100
    x_centers, dx = _uniform_grid(N, dx_uniform=0.005)

    w = _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)

    # Float-precision tolerance is more generous than 1e-12 because
    # we do N sums of small floats; 1e-14 per term scales to ~1e-12
    # for N=100. We test against 1e-10 to be safe.
    assert np.isclose(w.sum(), 1.0, atol=1e-10), (
        f"sum(w)={w.sum()} for L_jet={L_jet}, expected 1.0"
    )


@pytest.mark.parametrize('L_jet', [0.01, 0.1])
def test_weights_sum_to_one_non_uniform_grid(L_jet):
    """Conservation also holds for non-uniform dx."""
    N = 30
    # Non-uniform: cells get linearly wider toward the aft end
    dx = np.linspace(0.005, 0.02, N)
    x_centers = np.cumsum(dx) - 0.5 * dx  # cell centers from the head

    w = _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)

    assert np.isclose(w.sum(), 1.0, atol=1e-10)


# ---------------------------------------------------------------------
# Functional form: matches hand calculation for known geometry
# ---------------------------------------------------------------------

def test_weights_match_hand_calculation_5_cell_uniform():
    """For a 5-cell uniform geometry with L_jet = 2 * dx, verify the
    weight pattern matches the textbook exponential-decay formula
    computed by hand.

    Geometry: 5 cells of 0.01 m each. x_centers = [0.005, 0.015, 0.025,
    0.035, 0.045]. L_jet = 0.02 (= 2 * dx).

    Raw weights: exp(-x_i / L_jet) * dx_i, all dx_i identical.
    Then normalized to sum to 1.
    """
    N = 5
    x_centers, dx = _uniform_grid(N, dx_uniform=0.01)
    L_jet = 0.02

    w = _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)

    # Hand-calc raw exponentials
    raw = np.exp(-x_centers / L_jet) * dx
    expected = raw / raw.sum()

    np.testing.assert_allclose(w, expected, atol=1e-12)
    # Sanity check on the expected ordering: cell 0 has highest weight
    assert expected[0] > expected[1] > expected[2] > expected[3] > expected[4]


def test_weights_monotonically_decreasing_for_uniform_grid():
    """For any L_jet > 0 on a uniform grid, weights must monotonically
    decrease with cell index (because exp(-x/L) is monotonically
    decreasing and dx is constant).
    """
    N = 25
    x_centers, dx = _uniform_grid(N)
    L_jet = 0.05

    w = _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)

    diffs = np.diff(w)
    assert np.all(diffs <= 0.0), (
        f"weights not monotonically decreasing: max diff = {diffs.max()}"
    )


def test_large_L_jet_gives_nearly_uniform_weights():
    """In the limit L_jet >> bore length, exp(-x/L_jet) → 1 so weights
    should approach dx[i] / sum(dx) (uniform on a uniform grid → 1/N).
    """
    N = 20
    x_centers, dx = _uniform_grid(N, dx_uniform=0.01)
    L_jet = 1000.0  # >> total bore length 0.20 m

    w = _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)

    expected_uniform = 1.0 / N
    np.testing.assert_allclose(w, expected_uniform, atol=1e-3)


def test_small_L_jet_concentrates_at_head_end():
    """L_jet << dx[0] should put most weight in cell 0 (approaching
    the L_jet=0 limit) without being exactly the cell-0-only case.
    """
    N = 10
    x_centers, dx = _uniform_grid(N, dx_uniform=0.01)
    L_jet = 1e-4  # 100x smaller than dx

    w = _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)

    assert w[0] > 0.99   # cell 0 holds essentially all the weight
    assert np.isclose(w.sum(), 1.0, atol=1e-10)
