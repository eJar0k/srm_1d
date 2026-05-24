"""
tests/test_uncontained_pyrogen.py — v0.7.3 Phase A.2 kernel tests
====================================================================

Direct-kernel tests for `_compute_uncontained_pyrogen_mdot`, the
per-cell Saint-Robert burn rate computation for the uncontained
submerged-igniter topologies (head_basket, aft_basket). These tests
verify the math in isolation; integration tests against
_run_time_loop come later in Phase A.2.

See srm_1d/igniter_plenum.py PyrogenChamber docstring for the
architecture (uncontained-bore-P vs plenum-with-orifice split) and
srm_1d/docs/v0_7_2/candidates_post_phaseA.md for the v0.7.3 design.
"""
import numpy as np
import pytest

from srm_1d.simulation import _compute_uncontained_pyrogen_mdot


# Realistic BPNV-class scalars used across multiple tests
BPNV_A = 2.0e-5    # m/s / Pa^n
BPNV_N = 0.5
BPNV_RHO = 1700.0  # kg/m^3


# ---------------------------------------------------------------------
# Basic Saint-Robert math
# ---------------------------------------------------------------------

def test_constant_pressure_uniform_mdot_per_cell():
    """All cartridge cells at the same bore P → all get the same mdot.
    Uses 5 cartridge cells at P=1 MPa with A_burn_per_cell=1e-4 m^2.
    """
    N = 10
    P_bore = np.full(N, 1.0e6)  # 1 MPa everywhere
    A_burn_per_cell = 1.0e-4
    m_remain = 1.0  # 1 kg pyrogen (plenty)
    dt = 1.0e-4
    mdot = np.zeros(N)

    new_remain = _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, A_burn_per_cell,
        m_remain, dt, 0, 4, N, mdot,
    )

    # r_b at 1 MPa: 2e-5 * (1e6)^0.5 = 2e-5 * 1000 = 0.02 m/s
    # mdot per cell: 1700 * 0.02 * 1e-4 = 0.0034 kg/s
    expected_per_cell = 0.0034
    np.testing.assert_allclose(mdot[0:5], expected_per_cell, rtol=1e-9)
    assert np.all(mdot[5:] == 0.0), "cells outside [0,4] must be zero"
    # Total consumption over dt: 5 * 0.0034 * 1e-4 = 1.7e-6 kg
    assert new_remain == pytest.approx(m_remain - 5 * 0.0034 * 1.0e-4)


def test_varying_pressure_higher_mdot_at_higher_P():
    """Linearly varying P_bore → mdot[i] follows a * P[i]^n."""
    N = 5
    # P_bore varies from 1 MPa at cell 0 to 5 MPa at cell 4
    P_bore = np.linspace(1.0e6, 5.0e6, N)
    A_burn_per_cell = 1.0e-4
    mdot = np.zeros(N)

    _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, A_burn_per_cell,
        1.0, 1.0e-4, 0, 4, N, mdot,
    )

    # mdot[i] should follow sqrt(P[i]) ratio
    for i in range(N):
        expected = BPNV_RHO * BPNV_A * (P_bore[i] ** BPNV_N) * A_burn_per_cell
        assert mdot[i] == pytest.approx(expected, rel=1e-9)
    # Monotonically increasing (since P is monotonic and n > 0)
    assert np.all(np.diff(mdot) > 0)


def test_negative_pressure_clamped_to_zero():
    """Defensive: numerical artifact P < 0 produces mdot=0, not NaN."""
    N = 3
    P_bore = np.array([-1.0e3, 1.0e6, 5.0e5])
    mdot = np.zeros(N)

    _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 1e-4,
        1.0, 1e-4, 0, 2, N, mdot,
    )

    assert mdot[0] == 0.0
    assert mdot[1] > 0.0  # 1 MPa
    assert mdot[2] > 0.0  # 0.5 MPa
    assert np.all(np.isfinite(mdot))


# ---------------------------------------------------------------------
# Mass conservation / depletion
# ---------------------------------------------------------------------

def test_depletion_scales_mdot_uniformly_in_last_step():
    """When sum(mdot) * dt would exceed m_pyrogen_remaining, all
    per-cell mdots are scaled down so the last step consumes exactly
    the remaining mass.
    """
    N = 5
    P_bore = np.full(N, 1.0e6)
    A_burn_per_cell = 1.0e-4
    # Pyrogen almost gone — only 1e-7 kg left
    m_remain = 1.0e-7
    dt = 1.0e-4
    mdot = np.zeros(N)

    # Without conservation: would burn 5 * 0.0034 * 1e-4 = 1.7e-6 kg
    # remaining 1e-7 kg → scale = 1e-7 / 1.7e-6 ≈ 0.0588
    new_remain = _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, A_burn_per_cell,
        m_remain, dt, 0, 4, N, mdot,
    )

    # All scaled equally (uniform P → uniform mdot ratio)
    expected_total = m_remain / dt
    assert sum(mdot) == pytest.approx(expected_total, rel=1e-9)
    for i in range(5):
        assert mdot[i] == pytest.approx(expected_total / 5.0, rel=1e-9)
    assert new_remain == 0.0


def test_zero_remaining_returns_all_zero():
    """After pyrogen depleted, subsequent steps return all-zero mdot."""
    N = 4
    P_bore = np.full(N, 2.0e6)
    mdot = np.zeros(N)

    new_remain = _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 1e-4,
        0.0, 1e-4, 0, 3, N, mdot,
    )

    assert np.all(mdot == 0.0)
    assert new_remain == 0.0


def test_normal_step_decrements_remaining_correctly():
    """Standard non-depleting step: new_remain = old - sum(mdot)*dt."""
    N = 5
    P_bore = np.full(N, 1.0e6)
    m_remain = 0.05  # 50g — plenty for many steps at this rate
    dt = 1.0e-4
    mdot = np.zeros(N)

    new_remain = _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 1e-4,
        m_remain, dt, 0, 4, N, mdot,
    )

    consumed = sum(mdot) * dt
    assert new_remain == pytest.approx(m_remain - consumed, rel=1e-12)
    assert new_remain > 0.0
    assert consumed > 0.0


# ---------------------------------------------------------------------
# Range / topology variants
# ---------------------------------------------------------------------

def test_head_basket_range_at_start():
    """i_start=0, i_end=2 fills only cells 0-2 (head topology)."""
    N = 10
    P_bore = np.full(N, 1.0e6)
    mdot = np.full(N, 99.9)  # sentinel — kernel must overwrite to 0

    _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 1e-4,
        1.0, 1e-4, 0, 2, N, mdot,
    )

    assert np.all(mdot[0:3] > 0.0)
    assert np.all(mdot[3:] == 0.0)


def test_aft_basket_range_at_end():
    """i_start=N-3, i_end=N-1 fills only the last 3 cells (aft topology)."""
    N = 10
    P_bore = np.full(N, 1.0e6)
    mdot = np.zeros(N)

    _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 1e-4,
        1.0, 1e-4, N - 3, N - 1, N, mdot,
    )

    assert np.all(mdot[:N - 3] == 0.0)
    assert np.all(mdot[N - 3:] > 0.0)


# ---------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------

def test_invalid_range_returns_unchanged():
    """i_start > i_end → no consumption, remaining unchanged."""
    N = 5
    P_bore = np.full(N, 1.0e6)
    m_remain = 0.1
    mdot = np.zeros(N)

    new_remain = _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 1e-4,
        m_remain, 1e-4, 4, 2, N, mdot,
    )

    assert np.all(mdot == 0.0)
    assert new_remain == m_remain


def test_zero_A_burn_returns_unchanged():
    """Zero burning area → no consumption."""
    N = 5
    P_bore = np.full(N, 1.0e6)
    m_remain = 0.1
    mdot = np.zeros(N)

    new_remain = _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 0.0,
        m_remain, 1e-4, 0, 4, N, mdot,
    )

    assert np.all(mdot == 0.0)
    assert new_remain == m_remain


def test_out_of_bounds_indices_clamped():
    """i_start=-3 and i_end=N+5 clamp to [0, N-1]."""
    N = 5
    P_bore = np.full(N, 1.0e6)
    mdot = np.zeros(N)

    _compute_uncontained_pyrogen_mdot(
        P_bore, BPNV_A, BPNV_N, BPNV_RHO, 1e-4,
        1.0, 1e-4, -3, N + 5, N, mdot,
    )

    # All N cells should get mdot
    assert np.all(mdot > 0.0)
