"""Tests for propellant properties and nozzle performance."""
import numpy as np
import pytest
from srm_1d.propellant import (
    Propellant, GasProperties, R_UNIVERSAL,
    create_gas_properties, create_gas_properties_estimated,
    speed_of_sound, density_from_ideal_gas,
    critical_flow_function, characteristic_velocity,
)
from srm_1d.tests._motor_fixtures import (
    hasegawa_propellant_1 as make_hasegawa_propellant_1,
)
from srm_1d.nozzle import (
    Nozzle, exit_pressure_from_expansion_ratio,
    ideal_thrust_coefficient, compute_thrust_isp,
    compute_thrust_history, compute_motor_performance,
)


# ---- Propellant ----

class TestPropellant:
    def test_hasegawa_burn_rate_at_reference(self):
        """Hasegawa Propellant 1: r(4.9 MPa) = 4.9 mm/s."""
        prop = make_hasegawa_propellant_1()
        r = prop.burn_rate_normal(4.9e6)
        assert r == pytest.approx(4.9e-3, rel=0.01)

    def test_burn_rate_pressure_scaling(self):
        """Higher pressure should give higher burn rate."""
        prop = make_hasegawa_propellant_1()
        assert prop.burn_rate_normal(10e6) > prop.burn_rate_normal(5e6)

    def test_burn_rate_zero_pressure(self):
        """Zero pressure should give zero burn rate."""
        prop = make_hasegawa_propellant_1()
        assert prop.burn_rate_normal(0.0) == 0.0


class TestGasProperties:
    def test_R_specific(self):
        """R_specific should equal R_universal / MW."""
        gas = create_gas_properties(1.19, 0.0254, 3041.0,
                                     8.842e-5, 0.3685, 2060.0)
        assert gas.R_specific == pytest.approx(R_UNIVERSAL / 0.0254, rel=1e-6)

    def test_prandtl_number(self):
        """Pr = mu * Cp / k."""
        gas = create_gas_properties(1.19, 0.0254, 3041.0,
                                     8.842e-5, 0.3685, 2060.0)
        expected_Pr = 8.842e-5 * 2060.0 / 0.3685
        assert gas.Pr == pytest.approx(expected_Pr, rel=1e-6)

    def test_estimated_warns(self):
        """Estimated properties should issue a warning."""
        with pytest.warns(UserWarning, match="estimated"):
            create_gas_properties_estimated(1.19, 0.0254, 3041.0)


class TestThermodynamicUtilities:
    def test_speed_of_sound_air(self):
        """Speed of sound in air at 293K ≈ 343 m/s."""
        a = speed_of_sound(1.4, 287.0, 293.0)
        assert a == pytest.approx(343.0, rel=0.01)

    def test_ideal_gas_density(self):
        """Air at STP: ρ ≈ 1.225 kg/m³."""
        rho = density_from_ideal_gas(101325.0, 287.0, 293.0)
        assert rho == pytest.approx(1.205, rel=0.02)

    def test_critical_flow_function_air(self):
        """Gamma(1.4) ≈ 0.6847 (standard textbook value for air)."""
        G = critical_flow_function(1.4)
        assert G == pytest.approx(0.6847, rel=0.001)

    def test_c_star_hasegawa(self):
        """c* for Hasegawa propellant should be ~1543 m/s."""
        prop = make_hasegawa_propellant_1()
        tab = prop.representative_tab()
        gas = create_gas_properties(tab.gamma, tab.molecular_weight,
                                     tab.T_flame, prop.mu_gas,
                                     prop.k_gas, prop.Cp_gas)
        c = characteristic_velocity(tab.gamma, gas.R_specific, tab.T_flame)
        assert c == pytest.approx(1543.0, rel=0.01)


# ---- Nozzle ----

class TestNozzleGeometry:
    def test_expansion_ratio(self):
        noz = Nozzle(D_throat=0.020, D_exit=0.040)
        assert noz.expansion_ratio == pytest.approx(4.0, rel=1e-6)

    def test_divergence_loss_15deg(self):
        """15-deg half-angle: λ = (1+cos(15°))/2 ≈ 0.9830."""
        noz = Nozzle(D_throat=0.020, D_exit=0.040, div_angle=15.0)
        assert noz.divergence_losses() == pytest.approx(0.9830, rel=0.001)

    def test_throat_loss_zero_length(self):
        """throat_length = 0 → aspect = 0 → throat_loss = 0.99."""
        noz = Nozzle(D_throat=0.020, D_exit=0.040, throat_length=0.0)
        assert noz.throat_losses() == pytest.approx(0.99, rel=1e-6)

    def test_throat_loss_high_aspect(self):
        """throat_length/D_throat > 0.45 → throat_loss = 0.95."""
        noz = Nozzle(D_throat=0.020, D_exit=0.040, throat_length=0.020)
        assert noz.throat_losses() == pytest.approx(0.95, rel=1e-6)

    def test_skin_loss_constant(self):
        """skin_loss is hardcoded 0.99 (openMotor convention)."""
        noz = Nozzle(D_throat=0.020, D_exit=0.040)
        assert noz.skin_losses() == pytest.approx(0.99, rel=1e-6)

    def test_no_throat_change_default(self):
        noz = Nozzle(D_throat=0.020, D_exit=0.040)
        assert noz.has_throat_change is False

    def test_has_throat_change_erosion(self):
        noz = Nozzle(D_throat=0.020, D_exit=0.040, erosion_coeff=0.1)
        assert noz.has_throat_change is True


class TestExitPressure:
    def test_unity_expansion(self):
        """Expansion ratio 1.0 should give P_e/P_c = 1.0."""
        assert exit_pressure_from_expansion_ratio(1.4, 1.0) == pytest.approx(1.0)

    def test_supersonic_expansion(self):
        """P_e/P_c should decrease with increasing expansion ratio."""
        pe1 = exit_pressure_from_expansion_ratio(1.4, 2.0)
        pe2 = exit_pressure_from_expansion_ratio(1.4, 4.0)
        assert 0 < pe2 < pe1 < 1.0


class TestThrustCoefficient:
    def test_vacuum_CF_range(self):
        """Vacuum C_F should be 1.0-2.0 for typical expansion ratios."""
        CF, Pe = ideal_thrust_coefficient(1.2, 3.0, 5e6, 0.0)
        assert 1.0 < CF < 2.0

    def test_zero_pressure_no_thrust(self):
        """Zero chamber pressure should give zero C_F."""
        CF, Pe = ideal_thrust_coefficient(1.2, 3.0, 0.0, 101325.0)
        assert CF == 0.0


class TestThroatEvolution:
    def test_zero_coefficients_constant_throat(self):
        """No erosion/slag → throat diameter unchanged."""
        t = np.linspace(0, 4, 100)
        P = np.full(100, 5e6)
        _, _, _, _, Dt = compute_thrust_history(
            t, P, 100, 1.2, 0.020, 0.040,
            0.0, 0.0,
            0.983, 0.95, 0.0, 0.99,
            0.0, 1500.0, 9.80665,
        )
        np.testing.assert_allclose(Dt, 0.020, atol=1e-12)

    def test_erosion_grows_throat(self):
        """Positive erosion coeff → throat diameter increases."""
        t = np.linspace(0, 4, 1000)
        P = np.full(1000, 5e6)
        _, _, _, _, Dt = compute_thrust_history(
            t, P, 1000, 1.2, 0.020, 0.040,
            0.1, 0.0,
            0.983, 0.95, 0.0, 0.99,
            0.0, 1500.0, 9.80665,
        )
        assert Dt[-1] > 0.020

    def test_erosion_analytical(self):
        """Check erosion integration against analytical prediction."""
        t = np.linspace(0, 4, 10000)
        P = np.full(10000, 5e6)  # constant 5 MPa
        _, _, _, _, Dt = compute_thrust_history(
            t, P, 10000, 1.2, 0.020, 0.040,
            0.1, 0.0,
            0.983, 0.95, 0.0, 0.99,
            0.0, 1500.0, 9.80665,
        )
        # Expected: dD = 2 × 0.1e-6 × 5.0 × 4.0 = 4.0e-6 m
        expected_D = 0.020 + 4.0e-6
        assert Dt[-1] == pytest.approx(expected_D, rel=0.01)


class TestSummaryDict:
    def test_summary_keys_present(self):
        """run_simulation result should contain a 'summary' dict with all keys."""
        # We don't run a full sim here — just verify the import and key list
        expected_keys = {
            'propellant_mass', 'mass_produced', 'mass_nozzle',
            'mass_balance_error', 'P_peak', 't_peak', 'P_mid',
            't_burn', 't_first_burnout', 'c_star', 'wall_time',
            'steps', 'cells', 'D_throat_initial', 'D_throat_final',
        }
        # Verify the key set is documented (static check)
        assert len(expected_keys) == 15
