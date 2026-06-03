"""
test_zn_burn_rate.py — v0.7.4 Phase Z Zeldovich-Novozhilov relaxation
=====================================================================

Covers the `_advance_zn_burn_rate` kernel (analytic relaxation, the
tau = kappa_zn*alpha_s/r^2 timescale, first-touch seeding, non-burning
reset, plateau preservation / no overshoot) plus a full-sim check that
enabling Z-N does not collapse the run.

GOTCHA #1: delete srm_1d/__pycache__/ (+ .nbi/.nbc) after editing the
@njit kernels before running.
"""

import math
import numpy as np
import pytest

from srm_1d.simulation import _advance_zn_burn_rate, ZN_R_FLOOR


ALPHA = 1.0e-7      # solid thermal diffusivity [m^2/s]
KAPPA_ZN = 1.0
DT = 1.0e-5


def test_single_step_matches_analytic_relaxation():
    N = 1
    r_dyn = np.array([0.005])
    r_qs = np.array([0.010])
    is_burning = np.array([True])

    r0 = r_dyn[0]
    tau = KAPPA_ZN * ALPHA / (r0 * r0)        # r0 > ZN_R_FLOOR
    expected = r_qs[0] + (r0 - r_qs[0]) * math.exp(-DT / tau)

    _advance_zn_burn_rate(r_dyn, r_qs, is_burning, ALPHA, KAPPA_ZN,
                          ZN_R_FLOOR, DT, N)
    assert r_dyn[0] == pytest.approx(expected, rel=1e-12)
    # tau ≈ 4 ms at 5 mm/s (Lengelle 1975) — one µs step barely moves it.
    assert tau == pytest.approx(4.0e-3, rel=1e-6)


def test_first_touch_seeds_to_quasi_steady():
    """A just-ignited cell (r_dyn==0) jumps straight to r_qs (no ramp-up)."""
    r_dyn = np.array([0.0])
    r_qs = np.array([0.008])
    is_burning = np.array([True])
    _advance_zn_burn_rate(r_dyn, r_qs, is_burning, ALPHA, KAPPA_ZN,
                          ZN_R_FLOOR, DT, 1)
    assert r_dyn[0] == pytest.approx(0.008)


def test_non_burning_cell_resets():
    r_dyn = np.array([0.006])
    r_qs = np.array([0.010])
    is_burning = np.array([False])
    _advance_zn_burn_rate(r_dyn, r_qs, is_burning, ALPHA, KAPPA_ZN,
                          ZN_R_FLOOR, DT, 1)
    assert r_dyn[0] == 0.0


def test_plateau_preserved_after_many_steps():
    """Constant r_qs for >> tau ⇒ r_dyn converges to r_qs (no plateau droop)."""
    r_dyn = np.array([0.005])
    r_qs = np.array([0.010])
    is_burning = np.array([True])
    for _ in range(5000):   # 5000 * 1e-5 s = 50 ms >> tau ≈ 4 ms
        _advance_zn_burn_rate(r_dyn, r_qs, is_burning, ALPHA, KAPPA_ZN,
                              ZN_R_FLOOR, DT, 1)
    assert r_dyn[0] == pytest.approx(0.010, rel=1e-3)


def test_relaxation_never_overshoots():
    """r_dyn stays monotonically between its start and r_qs (analytic form)."""
    r_dyn = np.array([0.004])
    r_qs = np.array([0.012])
    is_burning = np.array([True])
    prev = r_dyn[0]
    for _ in range(2000):
        _advance_zn_burn_rate(r_dyn, r_qs, is_burning, ALPHA, KAPPA_ZN,
                              ZN_R_FLOOR, DT, 1)
        assert 0.004 - 1e-12 <= r_dyn[0] <= 0.012 + 1e-12  # bounded
        assert r_dyn[0] >= prev - 1e-15                     # monotone up
        prev = r_dyn[0]


def test_zero_target_holds_zero():
    """Burning cell whose QS target is still 0 (pre-burn-update lag) stays 0."""
    r_dyn = np.array([0.0])
    r_qs = np.array([0.0])
    is_burning = np.array([True])
    _advance_zn_burn_rate(r_dyn, r_qs, is_burning, ALPHA, KAPPA_ZN,
                          ZN_R_FLOOR, DT, 1)
    assert r_dyn[0] == 0.0


# ---------------------------------------------------------------------------
# Integration: enabling Z-N must not collapse the run nor depress the plateau.
# ---------------------------------------------------------------------------

def test_zn_full_sim_does_not_collapse():
    from srm_1d.openmotor_adapter import run_from_ric

    result, _perf, _nz, _geo, _prop = run_from_ric(
        'motors/hasegawa_a.ric',
        roughness=37.1e-6, kappa=0.45, pyrogen='bpnv', pyrogen_mass=None,
        T_ignition=850.0, P_cutoff=0.05e6,
        snapshot_interval=2.0, print_interval=20.0, verbose=False,
        zn_enabled=True, kappa_zn=1.0,
    )
    s = result['summary']
    assert s['termination_code'] != 4, "run collapsed under Z-N"
    assert s['t_burn'] > 0.5, "burn ended implausibly early"
    # Plateau preserved: a healthy Hasegawa A peaks in a sane band.
    assert 2.0e6 < s['P_peak'] < 30.0e6
