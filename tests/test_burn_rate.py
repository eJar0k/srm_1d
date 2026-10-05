"""Tests for the Ma et al. erosive burning model."""
import numpy as np
import pytest
from srm_1d.burn_rate import (
    haaland_friction, gnielinski_nusselt, transpiration_correction,
    burn_rate_cell, compute_burn_rates, mukunda_paul_eta,
    saint_robert_from_tabs,
    ErosiveModel, EROSIVE_MODELS, resolve_erosive_model, register_erosive_model,
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

class _HasegawaPropellant1Fixture:
    """Shared Hasegawa Propellant 1 parameters + burn_rate_cell caller.

    NOT a pytest test class (no ``Test`` prefix, not collected) --
    subclassed by TestBurnRateCell and the erosive-model-registry test
    classes below so they all exercise the same reference propellant
    without pytest re-collecting/re-running TestBurnRateCell's own test
    methods on every subclass.
    """

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

    MU_GAS = 9.0e-5

    def _call(self, P, Re, D=0.04, x=0.5, erosive_mp_cap=0.0, erosive_model_code=0):
        return burn_rate_cell(
            P, Re, D, x, 20e-6,
            self.PR, self.K_THERMAL, self.CP_GAS,
            self.T_FLAME, self.T_SURFACE,
            self.RHO_P, self.CPS, self.T_INITIAL,
            self.TAB_MIN_P, self.TAB_MAX_P, self.TAB_A, self.TAB_N, self.N_TABS,
            self.KAPPA, self.MU_GAS, erosive_mp_cap, erosive_model_code,
        )


class TestBurnRateCell(_HasegawaPropellant1Fixture):
    """Test burn_rate_cell with Hasegawa Propellant 1 parameters."""

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

    def test_mp_cap_off_by_default_identical(self):
        """erosive_mp_cap=0.0 must reproduce the uncapped Ma result."""
        for Re in (50.0, 50000.0, 500000.0):
            assert self._call(5e6, Re) == self._call(5e6, Re,
                                                     erosive_mp_cap=0.0)

    def test_mp_cap_bounds_eta_at_mukunda_paul(self):
        """With the cap on, eta = r_total/r0 never exceeds scale × the MP
        Eq. 12 universal law; at high Re (high g) the cap actively bites,
        and the 1.1 (MP ±10% band) variant sits between 1.0 and uncapped."""
        from srm_1d.burn_rate import mukunda_paul_eta
        P, Re, D = 5e6, 500000.0, 0.04
        r0 = self.A_SR * P**self.N_SR
        G = Re * self.MU_GAS / D
        eta_mp = mukunda_paul_eta(G, r0, D, self.MU_GAS, self.RHO_P)
        r_unc, _ = self._call(P, Re, D=D)
        r_cap, r_ero_cap = self._call(P, Re, D=D, erosive_mp_cap=1.0)
        r_cap11, _ = self._call(P, Re, D=D, erosive_mp_cap=1.1)
        assert r_cap / r0 <= eta_mp * (1.0 + 1e-12)
        assert r_cap11 / r0 <= 1.1 * eta_mp * (1.0 + 1e-12)
        assert r_ero_cap == pytest.approx(r_cap - r0)
        if r_unc / r0 > eta_mp:      # Ma above MP here -> cap must bite
            assert r_cap < r_unc
            assert r_cap == pytest.approx(r0 * eta_mp)
            assert r_cap <= r_cap11 <= r_unc


# ---- Ma bisection: guaranteed upper bracket (solver-verification T2) ----

class TestMaBisectionBracket(_HasegawaPropellant1Fixture):
    """The code-0 (Ma 2020) bisection used a fixed 20*r0 upper bracket and
    silently returned ~20*r0 whenever the true root lay above it (small
    port / high Re). burn_rate_cell now checks F(20*r0) and, if negative,
    widens the bracket to the analytic r0 + ebf*h0 (transpiration factor
    beta/(exp(beta)-1) <= 1  =>  h(r) <= h0  =>  F(r0 + ebf*h0) >= 0).

    The residual is evaluated here with the SAME module-level correlation
    helpers burn_rate_cell calls (saint_robert_from_tabs -> haaland_friction
    -> gnielinski_nusselt), in the same order -- burn_rate.py code-0 path,
    "Normal burn rate" ... "Precompute alpha" -- plus the three lines of
    inline algebra it uses for h0 / ebf / alpha (copied below).
    """

    ROUGHNESS = 20e-6   # the value _call hard-codes
    X_FROM_HEAD = 0.5   # _call's default x
    D_HYD_MM = (2, 4, 8, 16, 32)
    P_MPA = (2, 7, 15)

    def _ma_terms(self, P, Re, D):
        """(r0, h0, ebf, alpha) exactly as burn_rate_cell builds them."""
        r0 = saint_robert_from_tabs(
            P, self.TAB_MIN_P, self.TAB_MAX_P, self.TAB_A, self.TAB_N,
            self.N_TABS)
        f = haaland_friction(Re, self.ROUGHNESS, D)
        Nu = gnielinski_nusselt(Re, self.PR, D, self.X_FROM_HEAD, f,
                                self.T_FLAME, self.T_SURFACE, self.KAPPA)
        h0 = Nu * self.K_THERMAL / D
        ebf = (self.T_FLAME - self.T_SURFACE) / (
            self.RHO_P * self.CPS * (self.T_SURFACE - self.T_INITIAL))
        alpha = self.RHO_P * self.CP_GAS / h0
        return r0, h0, ebf, alpha

    @staticmethod
    def _F(r, r0, h0, ebf, alpha):
        """Ma residual F(r) = r - r0 - ebf*h(r), with burn_rate_cell's
        h(beta) branches (< 1e-6 Taylor / > 500 zero / exp)."""
        beta = alpha * r
        if beta < 1e-6:
            h = h0 * (1.0 - beta / 2.0)
        elif beta > 500.0:
            h = 0.0
        else:
            h = h0 * beta / (np.exp(beta) - 1.0)
        return r - r0 - ebf * h

    def _envelope(self):
        """Yield (D_hyd_mm, Re, P_MPa, D, P) over the 150-state envelope."""
        for D_mm in self.D_HYD_MM:
            for Re in np.logspace(4.0, 7.0, 10):
                for P_MPa in self.P_MPA:
                    yield D_mm, Re, P_MPa, D_mm * 1e-3, P_MPa * 1e6

    def test_envelope_finite_above_r0_and_residual_closed(self):
        """150 states (D_hyd 2-32 mm x Re 1e4-1e7 x P 2/7/15 MPa): finite,
        r_total >= r0, and the Ma residual |F(r_total)| <= 1e-7 m/s."""
        n = 0
        bad = []
        for D_mm, Re, P_MPa, D, P in self._envelope():
            n += 1
            r_total, r_erosive = self._call(P, Re, D=D)
            r0, h0, ebf, alpha = self._ma_terms(P, Re, D)
            F = self._F(r_total, r0, h0, ebf, alpha)
            ok = (np.isfinite(r_total) and np.isfinite(r_erosive)
                  and r_total >= r0 and r_erosive >= 0.0 and abs(F) <= 1e-7)
            if not ok:
                bad.append((D_mm, Re, P_MPa, r_total / r0, F))
        assert n == 150
        assert not bad, f"{len(bad)} envelope states failed: {bad[:5]}"

    def test_clamp_regression_root_above_20_r0(self):
        """At least one envelope state has its root above 20*r0 -- the old
        fixed bracket would have silently clamped it to ~20*r0."""
        above = [(D_mm, Re, P_MPa)
                 for D_mm, Re, P_MPa, D, P in self._envelope()
                 if self._call(P, Re, D=D)[0]
                 > 20.0 * self._ma_terms(P, Re, D)[0] * (1.0 + 1e-6)]
        assert len(above) >= 1

    def test_widened_bracket_is_a_guaranteed_bracket(self):
        """The analytic reason for the fix: F(r0 + ebf*h0) >= 0 in every
        envelope state (h(r) <= h0 because beta/(exp(beta)-1) <= 1)."""
        for D_mm, Re, P_MPa, D, P in self._envelope():
            r0, h0, ebf, alpha = self._ma_terms(P, Re, D)
            assert self._F(r0 + ebf * h0, r0, h0, ebf, alpha) >= 0.0


# ---- Erosive model registry (docs/erosive_model_registry) ----

class TestErosiveModelRegistry:
    def test_resolve_known_names(self):
        assert resolve_erosive_model('ma2020') == 0
        assert resolve_erosive_model('mukunda_paul') == 1
        assert resolve_erosive_model('saint_robert') == 2

    def test_registry_descriptors(self):
        assert EROSIVE_MODELS['ma2020'].code == 0
        assert EROSIVE_MODELS['mukunda_paul'].code == 1
        assert EROSIVE_MODELS['saint_robert'].code == 2
        assert EROSIVE_MODELS['ma2020'].name == 'ma2020'

    def test_resolve_unknown_raises(self):
        with pytest.raises(ValueError, match="unknown erosive_model"):
            resolve_erosive_model('nope')

    def test_resolve_unknown_lists_registered_names(self):
        with pytest.raises(ValueError, match="ma2020"):
            resolve_erosive_model('nope')

    def test_register_duplicate_name_rejected(self):
        dup = ErosiveModel('ma2020', 99, 'dup', 'dup name', ())
        with pytest.raises(ValueError, match="name already registered"):
            register_erosive_model(dup)

    def test_register_duplicate_code_rejected(self):
        dup = ErosiveModel('brand_new_model', 0, 'dup code', 'dup code', ())
        with pytest.raises(ValueError, match="code already registered"):
            register_erosive_model(dup)
        # must not have been inserted despite the name being novel
        assert 'brand_new_model' not in EROSIVE_MODELS


# ---- G2: mukunda_paul kernel branch (erosive_model_code=1) ----

class TestMukundaPaulModel(_HasegawaPropellant1Fixture):
    """Reuses TestBurnRateCell's Hasegawa Propellant 1 fixture."""

    def test_high_Re_matches_closed_form(self):
        """At high Re (g>35), mukunda_paul must return EXACTLY
        r0 * mukunda_paul_eta(...) -- no bisection, no iteration."""
        P, Re, D = 5e6, 500000.0, 0.04
        r0 = self.A_SR * P ** self.N_SR
        G = Re * self.MU_GAS / D
        eta_mp = mukunda_paul_eta(G, r0, D, self.MU_GAS, self.RHO_P)
        r_total, r_erosive = self._call(P, Re, D=D, erosive_model_code=1)
        assert r_total == pytest.approx(r0 * eta_mp, rel=1e-12)
        assert r_erosive == pytest.approx(r0 * eta_mp - r0, rel=1e-12)

    def test_g_below_threshold_gives_r0(self):
        """Below MP's g=35 threshold, eta=1 -> mukunda_paul returns r0
        exactly (same as the Ma path at this low-erosion condition)."""
        P, Re, D = 5e6, 50000.0, 0.04
        r0 = self.A_SR * P ** self.N_SR
        r_total, r_erosive = self._call(P, Re, D=D, erosive_model_code=1)
        assert r_total == pytest.approx(r0, rel=1e-12)
        assert r_erosive == 0.0

    def test_Re_below_100_gives_r0(self):
        """The shared Re_local < 100 early-return applies to mukunda_paul
        too (falls through before the model branch is ever reached)."""
        P = 5e6
        r0 = self.A_SR * P ** self.N_SR
        r_total, r_erosive = self._call(P, 50.0, erosive_model_code=1)
        assert r_total == pytest.approx(r0, rel=1e-12)
        assert r_erosive == 0.0


# ---- G2b: saint_robert kernel branch (erosive_model_code=2) ----

class TestSaintRobertModel(_HasegawaPropellant1Fixture):
    """Reuses TestBurnRateCell's Hasegawa Propellant 1 fixture."""

    def test_low_Re_returns_r0(self):
        P = 5e6
        r0 = self.A_SR * P ** self.N_SR
        r_total, r_erosive = self._call(P, 50.0, erosive_model_code=2)
        assert r_total == pytest.approx(r0, rel=1e-12)
        assert r_erosive == 0.0

    def test_very_high_Re_still_returns_r0(self):
        """No flow dependence at all -- even at an Re that would drive
        strong erosion under ma2020, saint_robert returns bare r0."""
        P = 5e6
        r0 = self.A_SR * P ** self.N_SR
        r_total, r_erosive = self._call(P, 5e6, erosive_model_code=2)
        assert r_total == pytest.approx(r0, rel=1e-12)
        assert r_erosive == 0.0
        # Sanity: ma2020 at the same condition is NOT bare r0 (confirms
        # the test is actually exercising the erosive-off branch).
        r_total_ma, _ = self._call(P, 5e6, erosive_model_code=0)
        assert r_total_ma > r0


# ---- G3: mukunda_paul == ma2020 with erosive_mp_cap=1.0, where Ma exceeds MP ----

class TestMPModelMatchesCappedMa(_HasegawaPropellant1Fixture):
    """Reuses TestBurnRateCell's Hasegawa Propellant 1 fixture."""

    def test_mukunda_paul_equals_ma_cap_1_0_where_ma_exceeds_mp(self):
        P, Re, D = 5e6, 5e5, 0.04  # eta_Ma (1.713) > eta_MP (1.621) here
        r0 = self.A_SR * P ** self.N_SR
        G = Re * self.MU_GAS / D
        eta_mp = mukunda_paul_eta(G, r0, D, self.MU_GAS, self.RHO_P)
        r_ma_uncapped, _ = self._call(P, Re, D=D, erosive_model_code=0)
        assert r_ma_uncapped / r0 > eta_mp  # confirm the premise holds
        r_mp, r_ero_mp = self._call(P, Re, D=D, erosive_model_code=1)
        r_ma_capped, r_ero_ma_capped = self._call(
            P, Re, D=D, erosive_mp_cap=1.0, erosive_model_code=0)
        assert r_mp == r_ma_capped
        assert r_ero_mp == r_ero_ma_capped


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
            0.45, 9.0e-5, 0.0, 0, N,
        )
        assert np.all(r_total[:3] == 0.0)
        assert np.all(r_total[7:] == 0.0)
        assert np.all(r_total[3:7] > 0.0)


# ---- run_simulation kwarg validation (T1 task 6 integration guard) ----
#
# Both raise before geo/propellant/nozzle/pyrogen_chamber are ever
# touched (docs/erosive_model_registry/DESIGN.md §3), so plain sentinel
# objects stand in -- no real motor construction, no time loop.

class TestRunSimulationErosiveModelValidation:
    def test_unknown_model_raises_before_sim_setup(self):
        from srm_1d.simulation import run_simulation
        with pytest.raises(ValueError, match="unknown erosive_model"):
            run_simulation(object(), object(), object(), object(),
                            erosive_model='nope')

    def test_cap_with_non_ma2020_model_raises(self):
        from srm_1d.simulation import run_simulation
        with pytest.raises(ValueError, match="ma2020-only"):
            run_simulation(object(), object(), object(), object(),
                            erosive_model='mukunda_paul', erosive_mp_cap=1.0)

    def test_cap_with_ma2020_model_does_not_raise_here(self):
        """erosive_mp_cap > 0 with the default ma2020 model must pass this
        particular validation (it will proceed on to touch the sentinel
        pyrogen_chamber next, which is the expected AttributeError --
        confirms our ValueError did NOT fire)."""
        from srm_1d.simulation import run_simulation
        with pytest.raises(AttributeError):
            run_simulation(object(), object(), object(), object(),
                            erosive_model='ma2020', erosive_mp_cap=1.0)
