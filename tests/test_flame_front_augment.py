"""
tests/test_flame_front_augment.py — v0.7.2 Phase B-v2 kernel tests
====================================================================

Direct-kernel tests for `_compute_flame_front_augment`, the
reformulated Phase B kernel that gates h_c augmentation on "cell i is
the immediate downstream neighbor of a cell that ignited within the
last tau_window seconds" rather than the original commit-065d193
cumulative-G magnitude (which double-counted with PISO's local-Re
tracking).

See srm_1d/docs/v0_7_2/candidates/02_spatial_ignition_front_coupling.md.
"""

import numpy as np

from srm_1d.simulation import _compute_flame_front_augment


def _zeros_bool(N):
    return np.zeros(N, dtype=np.bool_)


# ---------------------------------------------------------------------
# Default state: no cells burning → no augmentation anywhere
# ---------------------------------------------------------------------

def test_no_cells_burning_returns_all_unity():
    """With nothing burning, every cell stays at the no-op default 1.0."""
    N = 10
    is_burning = _zeros_bool(N)
    has_ignited = _zeros_bool(N)
    ignition_time = np.full(N, 1e10)  # never ignited
    augment = np.full(N, 99.9)  # sentinel — kernel must overwrite

    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 0.0,
        1e-3, 3.0, N, augment,
    )

    assert np.all(augment == 1.0), f"expected all 1.0, got {augment}"


# ---------------------------------------------------------------------
# Strict-sequential: cell 0 never gets boost (no upstream)
# ---------------------------------------------------------------------

def test_cell_0_never_augmented():
    """Even when cell 0 itself has 'just ignited' (impossible but
    defensive), cell 0's augment stays 1.0 because there's no
    upstream cell j with j+1=0.
    """
    N = 5
    is_burning = np.array([True, False, False, False, False])
    has_ignited = np.array([True, False, False, False, False])
    ignition_time = np.array([0.5e-3, 1e10, 1e10, 1e10, 1e10])
    augment = np.zeros(N)

    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 1e-3,
        1.5e-3, 3.0, N, augment,
    )

    # Cell 0 stays at 1.0; cell 1 gets the boost; cells 2-4 unchanged.
    assert augment[0] == 1.0
    assert augment[1] == 3.0
    assert augment[2] == 1.0
    assert augment[3] == 1.0
    assert augment[4] == 1.0


# ---------------------------------------------------------------------
# Front propagation: only the immediate-downstream-of-burning cell
# ---------------------------------------------------------------------

def test_only_immediate_downstream_neighbor_boosted():
    """Cells 0-2 burning, cells 3-9 unignited. Only cell 3 should get
    the boost (the leading-edge downstream neighbor).
    """
    N = 10
    is_burning = np.array([True, True, True, False, False, False, False, False, False, False])
    has_ignited = is_burning.copy()
    # All three cells ignited recently.
    ignition_time = np.array([0.0, 0.2e-3, 0.4e-3, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10])
    augment = np.zeros(N)

    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 0.5e-3,
        1e-3, 3.0, N, augment,
    )

    # Cells 0, 1 sit upstream of burning cells → augment stays 1.0
    # (the kernel boosts cell j+1 only if cell j is burning AND j+1
    # is unignited; cells 1 and 2 are ignited so they don't get
    # boosted).
    assert augment[0] == 1.0
    assert augment[1] == 1.0
    assert augment[2] == 1.0
    # Cell 3 IS downstream of burning cell 2 AND unignited → boost.
    assert augment[3] == 3.0
    # Cells 4-9 are downstream of unignited cells → no boost
    assert np.all(augment[4:] == 1.0)


# ---------------------------------------------------------------------
# Time window: boost expires after tau_window
# ---------------------------------------------------------------------

def test_boost_expires_after_tau_window():
    """Cell 0 ignited at t=0. At t = tau_window the boost on cell 1
    should expire (transition from boost to 1.0).
    """
    N = 3
    is_burning = np.array([True, False, False])
    has_ignited = np.array([True, False, False])
    ignition_time = np.array([0.0, 1e10, 1e10])
    augment = np.zeros(N)
    tau_window = 1e-3
    boost = 3.0

    # At t = tau_window/2: cell 1 should get boost.
    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 0.5e-3,
        tau_window, boost, N, augment,
    )
    assert augment[1] == boost, f"at t=tau/2 expected boost; got {augment[1]}"

    # At t = tau_window (boundary): boost should expire (strictly less than).
    augment.fill(0.0)
    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, tau_window,
        tau_window, boost, N, augment,
    )
    assert augment[1] == 1.0, f"at t=tau (boundary) expected 1.0; got {augment[1]}"

    # At t = 2*tau_window: boost is well past expiry.
    augment.fill(0.0)
    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 2.0 * tau_window,
        tau_window, boost, N, augment,
    )
    assert augment[1] == 1.0


# ---------------------------------------------------------------------
# Multi-cell scenario: independent boosts on multiple fronts
# ---------------------------------------------------------------------

def test_multiple_burning_fronts_each_boost_their_neighbor():
    """Suppose cells 0 and 4 both recently ignited (e.g., one from
    pyrogen, one from radiation), and cells 1, 2, 3, 5 are
    unignited. Cells 1 and 5 should both get boost; cells 2, 3 not.
    """
    N = 6
    is_burning = np.array([True, False, False, False, True, False])
    has_ignited = is_burning.copy()
    ignition_time = np.array([0.3e-3, 1e10, 1e10, 1e10, 0.4e-3, 1e10])
    augment = np.zeros(N)

    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 0.5e-3,
        1e-3, 3.0, N, augment,
    )

    assert augment[0] == 1.0
    assert augment[1] == 3.0  # downstream of cell 0
    assert augment[2] == 1.0
    assert augment[3] == 1.0
    assert augment[4] == 1.0
    assert augment[5] == 3.0  # downstream of cell 4


# ---------------------------------------------------------------------
# Final-cell edge case: cell N-1 burning doesn't boost anything
# (no cell N to boost)
# ---------------------------------------------------------------------

def test_last_cell_burning_no_phantom_boost():
    """If the last cell N-1 has ignited, there's no cell N to boost.
    The kernel must not write past the array end.
    """
    N = 4
    is_burning = np.array([False, False, False, True])
    has_ignited = is_burning.copy()
    ignition_time = np.array([1e10, 1e10, 1e10, 0.4e-3])
    augment = np.zeros(N)

    # Should not crash; should leave all augment at 1.0.
    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 0.5e-3,
        1e-3, 3.0, N, augment,
    )

    assert np.all(augment == 1.0)


# ---------------------------------------------------------------------
# Boost = 1.0 (no-op): behavior equivalent to flame_spread_enabled=False
# ---------------------------------------------------------------------

def test_boost_of_unity_is_noop():
    """Boost value of 1.0 means no augmentation. Equivalent to disabling
    flame spread entirely.
    """
    N = 5
    is_burning = np.array([True, False, False, False, False])
    has_ignited = np.array([True, False, False, False, False])
    ignition_time = np.array([0.3e-3, 1e10, 1e10, 1e10, 1e10])
    augment = np.zeros(N)

    _compute_flame_front_augment(
        is_burning, has_ignited, ignition_time, 0.5e-3,
        1e-3, 1.0, N, augment,
    )

    assert np.all(augment == 1.0)
