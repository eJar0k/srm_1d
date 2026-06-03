"""Tests for the Ma et al. erosive burning model."""
import numpy as np
import pytest
from srm_1d.burn_rate import (
    haaland_friction, gnielinski_nusselt, transpiration_correction,
    burn_rate_cell, compute_burn_rates,
)


# ---- Haaland friction factor ----

class TestHaalandFriction:
    def test_laminar_exact(self):
        """Laminar regime: f = 64/Re (Hagen-Poiseuille, exact)."""
        assert haaland_friction(1000.0, 20e-6, 0.04) == pytest.approx(0.064, rel=1e-6)
        assert haaland_friction(500.0, 20e-6, 0.04) == pytest.approx(0.128, rel=1e-6)

    def test_turbulent_range(self):
        """Turbulent regime: Darcy factor should be 0.01-0.05 for typical Re."""
        f = haaland_friction(50000.0, 20e-6, 0.04)
        assert 0.01 < f < 0.05

    def test_smooth_pipe_limit(self):
        """Zero roughness should give the smooth-pipe Haaland result."""
        f = haaland_friction(100000.0, 0.0, 0.04)
        assert 0.01 < f < 0.03

    def test_transition_continuity(self):
        """Transition region (2300-4000): blend should be between endpoints."""
        f_lam = haaland_friction(2299.0, 20e-6, 0.04)
        f_mid = haaland_friction(3150.0, 20e-6, 0.04)
        f_turb = haaland_friction(4001.0, 20e-6, 0.04)
        f_lo = min(f_lam, f_turb)
        f_hi = max(f_lam, f_turb)
        assert f_lo <= f_mid <= f_hi

    def test_minimum_Re(self):
        """Very low Re should not crash or return negative."""
        f = haaland_friction(0.1, 20e-6, 0.04)
        assert f > 0


# ---- Gnielinski Nusselt ----

class TestGnielinskiNusselt:
    def test_laminar_floor(self):
        """Nu should never drop below 3.66 (fully developed laminar limit)."""
        Nu = gnielinski_nusselt(100.0, 0.5, 0.04, 0.5, 0.03, 3000.0, 1000.0, 0.45)
        assert Nu >= 3.66

    def test_turbulent_increases_with_Re(self):
        """Higher Re should give higher Nu in the turbulent regime."""
        f1 = haaland_friction(10000.0, 20e-6, 0.04)
        f2 = haaland_friction(50000.0, 20e-6, 0.04)
        Nu1 = gnielinski_nusselt(10000.0, 0.5, 0.04, 0.5, f1, 3000.0, 1000.0, 0.45)
        Nu2 = gnielinski_nusselt(50000.0, 0.5, 0.04, 0.5, f2, 3000.0, 1000.0, 0.45)
        assert Nu2 > Nu1

    def test_kappa_zero_no_correction(self):
        """kappa=0 should disable the temperature ratio correction (K=1)."""
        f = haaland_friction(50000.0, 20e-6, 0.04)
        Nu_k0 = gnielinski_nusselt(50000.0, 0.5, 0.04, 0.5, f, 3000.0, 1000.0, 0.0)
        Nu_k45 = gnielinski_nusselt(50000.0, 0.5, 0.04, 0.5, f, 3000.0, 1000.0, 0.45)
        # With T_gas > T_surface, kappa > 0 increases Nu
        assert Nu_k45 > Nu_k0


# ---- Transpiration correction ----

class TestTranspirationCorrection:
    def test_zero_blowing(self):
        """beta=0 should give ratio=1 (no correction)."""
        assert transpiration_correction(0.0) == pytest.approx(1.0, abs=1e-6)

    def test_small_blowing_taylor(self):
        """Small beta should match Taylor expansion: 1 - beta/2."""
        assert transpiration_correction(1e-8) == pytest.approx(1.0 - 1e-8/2, rel=1e-4)

    def test_large_blowing_suppression(self):
        """Large beta should suppress heat transfer (ratio → 0)."""
        assert transpiration_correction(100.0) < 0.01

    def test_bounded_zero_one(self):
        """Result should always be in [0, 1]."""
        for beta in [0.0, 0.1, 1.0, 10.0, 100.0, 500.0]:
            r = transpiration_correction(beta)
            assert 0.0 <= r <= 1.0


# ---- Single-cell burn rate ----

class TestBurnRateCell:
    """Test burn_rate_cell with Hasegawa Propellant 1 parameters."""

    # Hasegawa Propellant 1 parameters
    P_REF = 4.9e6
    R_REF = 4.9e-3
    N_SR = 0.3
    A_SR = R_REF / (P_REF ** N_SR)
    RHO_P = 1700.0
    CPS = 1500.0
    T_SURFACE = 1000.0
    T_INITIAL = 293.0
    T_FLAME = 3041.0
    PR = 0.4943
    K_THERMAL = 0.3685
    CP_GAS = 2060.0
    KAPPA = 0.45

    # Single-tab spanning the operating range
    TAB_MIN_P = np.array([0.0])
    TAB_MAX_P = np.array([2e7])
    TAB_A = np.array([A_SR])
    TAB_N = np.array([N_SR])
    N_TABS = 1

    def _call(self, P, Re, D=0.04, x=0.5):
        return burn_rate_cell(
            P, Re, D, x, 20e-6,
            self.PR, self.K_THERMAL, self.CP_GAS,
            self.T_FLAME, self.T_SURFACE,
            self.RHO_P, self.CPS, self.T_INITIAL,
            self.TAB_MIN_P, self.TAB_MAX_P, self.TAB_A, self.TAB_N, self.N_TABS,
            self.KAPPA,
        )

    def test_low_Re_no_erosion(self):
        """At low Re, erosive component should be zero."""
        r_total, r_erosive = self._call(5e6, 50.0)
        assert r_erosive == 0.0
        assert r_total == pytest.approx(self.A_SR * 5e6**self.N_SR, rel=1e-6)

    def test_saint_robert_at_reference(self):
        """At the reference point, r0 should match the measured rate."""
        r_total, _ = self._call(self.P_REF, 50.0)
        assert r_total == pytest.approx(self.R_REF, rel=0.01)

    def test_erosion_positive_at_high_Re(self):
        """At high Re, there should be positive erosive burning."""
        r_total, r_erosive = self._call(5e6, 500000.0)
        assert r_erosive > 0.0
        assert r_total > self.A_SR * 5e6**self.N_SR

    def test_erosion_ratio_reasonable(self):
        """Erosion ratio should not exceed ~5x for typical conditions."""
        r_total, _ = self._call(5e6, 500000.0)
        r0 = self.A_SR * 5e6**self.N_SR
        assert r_total / r0 < 5.0


# ---- Vectorized wrapper ----

class TestComputeBurnRates:
    def test_non_burning_cells_zero(self):
        """Cells where is_burning=False should have zero burn rate."""
        N = 10
        P = np.full(N, 5e6)
        Re = np.full(N, 50000.0)
        D_hyd = np.full(N, 0.04)
        x = np.linspace(0.01, 0.5, N)
        is_burning = np.zeros(N, dtype=np.bool_)
        is_burning[3:7] = True

        a = 4.9e-3 / (4.9e6 ** 0.3)
        r_total, r_erosive = compute_burn_rates(
            P, Re, D_hyd, x, is_burning, 20e-6,
            0.49, 0.37, 2060.0, 3041.0, 1000.0,
            1700.0, 1500.0, 293.0,
            np.array([0.0]), np.array([2e7]),
            np.array([a]), np.array([0.3]), 1,
            0.45, N,
        )
        assert np.all(r_total[:3] == 0.0)
        assert np.all(r_total[7:] == 0.0)
        assert np.all(r_total[3:7] > 0.0)
