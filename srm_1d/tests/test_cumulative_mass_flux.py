"""
tests/test_cumulative_mass_flux.py — v0.7.2 Phase B.1 kernel tests
====================================================================

Direct-kernel tests for `_compute_cumulative_mass_flux` and
`_blowing_augmentation`, the two Numba kernels that build the
spatial-ignition-coupling mechanism. Kernel-level tests verify the
math in isolation; Phase B.3 will add integration tests after the
kernels are wired into the Goodman call site.

See srm_1d/docs/v0_7_2/candidates/02_spatial_ignition_front_coupling.md.
"""

import numpy as np

from srm_1d.simulation import (
    _compute_cumulative_mass_flux,
    _blowing_augmentation,
)


# ---------------------------------------------------------------------
# _blowing_augmentation: math identities
# ---------------------------------------------------------------------

def test_blowing_augmentation_zero_returns_unity():
    """G_cum=0 → factor=1.0 (no augmentation, recovers Phase A)."""
    assert _blowing_augmentation(0.0, 1.0) == 1.0


def test_blowing_augmentation_negative_returns_unity():
    """Defensive: negative G_cum routes to unity (kernel safety)."""
    assert _blowing_augmentation(-1.0, 1.0) == 1.0


def test_blowing_augmentation_zero_G_ref_returns_unity():
    """Defensive: G_ref=0 is a divide-by-zero hazard; route to unity."""
    assert _blowing_augmentation(1.0, 0.0) == 1.0


def test_blowing_augmentation_unity_ratio_matches_formula():
    """G_cum/G_ref = 1 → factor = 2^0.8 ≈ 1.7411."""
    factor = _blowing_augmentation(1.0, 1.0)
    expected = 2.0 ** 0.8
    assert abs(factor - expected) < 1e-12


def test_blowing_augmentation_monotonic_in_G_cum():
    """For fixed G_ref, factor must strictly increase with G_cum."""
    G_ref = 1.0
    G_values = np.linspace(0.0, 10.0, 25)
    factors = np.array([_blowing_augmentation(g, G_ref) for g in G_values])
    diffs = np.diff(factors)
    assert np.all(diffs >= 0.0), (
        f"factor not monotonically non-decreasing; min diff = {diffs.min()}"
    )


def test_blowing_augmentation_dittus_boelter_scaling():
    """Verify (1 + G_cum/G_ref)^0.8 across a sweep — full math check."""
    G_ref = 2.5
    for G_cum in [0.1, 0.5, 1.0, 2.5, 5.0, 10.0]:
        factor = _blowing_augmentation(G_cum, G_ref)
        expected = (1.0 + G_cum / G_ref) ** 0.8
        assert abs(factor - expected) < 1e-12, (
            f"G_cum={G_cum}, G_ref={G_ref}: got {factor}, expected {expected}"
        )


# ---------------------------------------------------------------------
# _compute_cumulative_mass_flux: structural tests
# ---------------------------------------------------------------------

def _setup_uniform_bore(N=10, A_port=1.0e-4, dx=0.01):
    """Build a uniform-bore test geometry: constant A_port, constant dx."""
    A_port_arr = np.full(N, A_port)
    return A_port_arr, dx


def test_G_cum_uniform_no_burning_equals_G_igniter():
    """With no cells burning, G_cum[i] = G_igniter / A_port[i] for all i
    (the pyrogen mass flux just propagates downstream unchanged).
    """
    N = 8
    A_port, dx = _setup_uniform_bore(N)
    G_igniter = 0.01  # kg/s
    rho_prop = 1700.0
    r_b = np.zeros(N)
    P_burn = np.full(N, 0.1)
    is_burning = np.zeros(N, dtype=np.bool_)
    G_cum = np.zeros(N)

    _compute_cumulative_mass_flux(
        G_igniter, rho_prop, r_b, P_burn, A_port,
        is_burning, dx, N, G_cum,
    )

    expected = G_igniter / A_port[0]
    assert np.allclose(G_cum, expected), (
        f"G_cum should equal G_igniter/A_port[i] uniformly; got {G_cum}, "
        f"expected {expected}"
    )


def test_G_cum_monotonic_when_upstream_cells_burning():
    """G_cum must increase monotonically with i when upstream cells
    are burning (mass addition only goes up).
    """
    N = 10
    A_port, dx = _setup_uniform_bore(N)
    G_igniter = 0.001
    rho_prop = 1700.0
    r_b = np.full(N, 0.005)        # 5 mm/s burning rate
    P_burn = np.full(N, 0.05)      # 5 cm perimeter
    is_burning = np.ones(N, dtype=np.bool_)  # all burning
    G_cum = np.zeros(N)

    _compute_cumulative_mass_flux(
        G_igniter, rho_prop, r_b, P_burn, A_port,
        is_burning, dx, N, G_cum,
    )

    diffs = np.diff(G_cum)
    assert np.all(diffs > 0.0), (
        f"G_cum should strictly increase when all upstream cells "
        f"burning; got {G_cum}, diffs {diffs}"
    )


def test_G_cum_hand_calculation_3_cell():
    """Hand calc for a 3-cell motor with cell 0 burning.

    G_igniter = 0.001 kg/s, A_port = 1e-4 m^2, dx = 0.01 m,
    rho_p = 1700 kg/m^3, r_b = 0.005 m/s, P_burn[0] = 0.05 m.
    Only cell 0 is burning.

    Expected:
        G_cum[0] = G_igniter / A_port = 0.001 / 1e-4 = 10.0 kg/(m^2*s)
        After cell 0:
            running_mdot = G_igniter + rho_p * r_b * P_burn * dx
                         = 0.001 + 1700 * 0.005 * 0.05 * 0.01
                         = 0.001 + 0.00425 = 0.00525 kg/s
        G_cum[1] = 0.00525 / 1e-4 = 52.5 kg/(m^2*s)
        G_cum[2] = 52.5 (cell 1 not burning so no further accumulation)
    """
    N = 3
    A_port = np.full(N, 1e-4)
    dx = 0.01
    G_igniter = 0.001
    rho_prop = 1700.0
    r_b = np.full(N, 0.005)
    P_burn = np.full(N, 0.05)
    is_burning = np.array([True, False, False])
    G_cum = np.zeros(N)

    _compute_cumulative_mass_flux(
        G_igniter, rho_prop, r_b, P_burn, A_port,
        is_burning, dx, N, G_cum,
    )

    expected = np.array([10.0, 52.5, 52.5])
    np.testing.assert_allclose(G_cum, expected, atol=1e-10)


def test_G_cum_zero_A_port_clamps_to_zero():
    """Defensive: cells with A_port=0 (closed end?) should clamp G_cum
    to 0 rather than divide-by-zero or NaN.
    """
    N = 4
    A_port = np.array([1e-4, 0.0, 1e-4, 1e-4])
    dx = 0.01
    G_igniter = 0.001
    rho_prop = 1700.0
    r_b = np.zeros(N)
    P_burn = np.zeros(N)
    is_burning = np.zeros(N, dtype=np.bool_)
    G_cum = np.zeros(N)

    _compute_cumulative_mass_flux(
        G_igniter, rho_prop, r_b, P_burn, A_port,
        is_burning, dx, N, G_cum,
    )

    assert G_cum[1] == 0.0, f"A_port=0 cell should clamp to 0; got {G_cum[1]}"
    assert np.all(np.isfinite(G_cum)), f"All G_cum values must be finite; got {G_cum}"


def test_G_cum_zero_G_igniter_propagates_to_burning_contributions_only():
    """With no pyrogen (G_igniter=0) but cells already burning, G_cum
    should still build up purely from the burning-cell contributions.
    """
    N = 5
    A_port = np.full(N, 1e-4)
    dx = 0.01
    G_igniter = 0.0
    rho_prop = 1700.0
    r_b = np.full(N, 0.005)
    P_burn = np.full(N, 0.05)
    is_burning = np.array([True, True, True, False, False])
    G_cum = np.zeros(N)

    _compute_cumulative_mass_flux(
        G_igniter, rho_prop, r_b, P_burn, A_port,
        is_burning, dx, N, G_cum,
    )

    # G_cum[0] should be 0 (nothing upstream of cell 0 and G_igniter=0)
    assert G_cum[0] == 0.0
    # G_cum[1] should be from cell 0's burning contribution
    expected_g1 = rho_prop * r_b[0] * P_burn[0] * dx / A_port[1]
    assert abs(G_cum[1] - expected_g1) < 1e-10
    # G_cum[4] should equal G_cum[3] (cell 3 not burning, no addition)
    assert G_cum[4] == G_cum[3]
