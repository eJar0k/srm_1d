"""
propellant.py — Propellant and Combustion Gas Properties
=========================================================

PURPOSE:
    Defines propellant physical properties, baseline burn rate model,
    and combustion gas thermodynamic/transport properties. This is the
    foundation that every other module depends on.

SAINT-ROBERT'S LAW (multi-tab):
    The baseline (non-erosive) burn rate is:

        r₀ = a(P) · P^n(P)

    where a, n, gamma, T_flame, and molecular_weight are stored
    per-pressure-range as `PropellantTab` entries (matches openMotor's
    motorlib.propellant.PropellantTab schema). Hard-switchover lookup:
    the tab whose (min_pressure, max_pressure) range contains P is used;
    pressures outside any range fall back to the tab with the closest
    boundary (matches openMotor's getCombustionProperties behavior).

    For srm_1d v0.3.x, only `a` and `n` vary in the hot loop. Gas
    thermo (gamma, T_flame, molecular_weight) is evaluated once at
    simulation start from a representative tab — see TODO note for
    per-step gas-thermo lookup.

GAS TRANSPORT PROPERTIES:
    The combustion gas properties (viscosity, thermal conductivity,
    Prandtl number) MUST come from a chemical equilibrium solver —
    NASA CEA, RPA, or similar. Use the EFFECTIVE (not frozen) values
    at the nozzle inlet / chamber conditions.

    The frozen-vs-effective distinction matters enormously:
        k_frozen  ≈ 0.37 W/(m·K)    k_effective  ≈ 0.65 W/(m·K)
        Pr_frozen ≈ 0.49             Pr_effective ≈ 0.38
    Using frozen values would under-predict heat transfer by ~40%.

UNITS (SI throughout):
    Temperature: K          Pressure: Pa
    Density: kg/m³          Velocity: m/s
    Length: m               Mass flux: kg/(m²·s)
    Viscosity: Pa·s         Thermal conductivity: W/(m·K)
    Specific heat: J/(kg·K) Molecular weight: kg/mol
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional


# Universal gas constant [J/(mol·K)]
R_UNIVERSAL = 8.314462


# ================================================================
# Propellant Tab (per-pressure-range combustion properties)
# ================================================================

@dataclass
class PropellantTab:
    """
    Per-pressure-range combustion properties (a, n, gamma, T_flame, MW).

    Mirrors openMotor's `motorlib.propellant.PropellantTab` schema
    (snake_case'd; units kept in our SI conventions). Multiple tabs
    can be stacked on a `Propellant` to model pressure-dependent
    combustion behavior with hard switchover.

    Attributes
    ----------
    min_pressure, max_pressure : float
        Inclusive operating-pressure range for this tab [Pa].
        (openMotor: minPressure, maxPressure.)
    a : float
        Saint-Robert burn rate coefficient [m/s / Pa^n]. Units depend
        on n. (openMotor: a.)
    n : float
        Saint-Robert pressure exponent [-]. (openMotor: n.)
    gamma : float
        Ratio of specific heats [-]. (openMotor: k.)
    T_flame : float
        Adiabatic flame temperature [K]. (openMotor: t.)
    molecular_weight : float
        Mean molecular weight of combustion gas [kg/mol]. CEA reports
        g/mol — divide by 1000. (openMotor: m, in g/mol.)
    """
    min_pressure: float
    max_pressure: float
    a: float
    n: float
    gamma: float
    T_flame: float
    molecular_weight: float


# ================================================================
# Pyrogen (single-tab igniter propellant)
# ================================================================

@dataclass
class Pyrogen:
    """
    Pyrogen propellant properties for the v0.7.0 igniter chamber.

    This mirrors the scalar combustion fields on PropellantTab, but is
    intentionally single-tab: small pyrogen charges burn over a narrow
    pressure range and do not need openMotor-style pressure bands.
    """
    name: str
    a: float             # Saint-Robert coefficient [m/s / Pa^n]
    n: float             # Saint-Robert pressure exponent [-]
    rho: float           # Solid pyrogen density [kg/m^3]
    T_flame: float       # Adiabatic flame temperature [K]
    M: float             # Product molecular weight [kg/mol]
    gamma: float         # Product gas specific heat ratio [-]
    impetus_W: float = 0.0  # Optional DeMar impetus [psi*in^3/g]


# ================================================================
# Propellant
# ================================================================

@dataclass
class Propellant:
    """
    Complete propellant specification for internal ballistics simulation.

    Combustion properties (a, n, gamma, T_flame, molecular_weight) live
    in `tabs` per openMotor's schema. Solid properties and gas transport
    are stored at the propellant level (openMotor doesn't model heat
    transfer, so transport properties are srm_1d-specific additions).

    Attributes
    ----------
    name : str
        Human-readable propellant identifier.
    tabs : list[PropellantTab]
        One or more pressure-dependent combustion tabs. Use a single
        tab spanning [0, large] for pressure-independent burn rate.
    rho_propellant : float
        Solid propellant density [kg/m³].
    Cps : float
        Specific heat of the solid propellant [J/(kg·K)].
    T_surface : float
        Propellant burning surface temperature [K].
    T_initial : float
        Initial (ambient) propellant temperature [K].
    mu_gas, k_gas, Cp_gas : float
        Combustion gas viscosity / thermal conductivity / Cp from CEA/RPA
        (EFFECTIVE values, not frozen). Stored at propellant level
        because srm_1d's heat-transfer model assumes single-tab gas
        transport — gas thermo (gamma, T, MW) varies per-tab but
        transport stays fixed.
    """
    name: str
    tabs: List[PropellantTab]

    # Solid propellant physical properties
    rho_propellant: float
    Cps: float
    T_surface: float
    T_initial: float

    # Combustion gas transport properties (from CEA/RPA, EFFECTIVE values)
    mu_gas: float
    k_gas: float
    Cp_gas: float
    k_solid: float = 0.3  # Solid conductivity for Goodman ignition [W/(m*K)]

    def select_tab(self, P) -> PropellantTab:
        """
        Hard-switchover tab lookup. Mirrors openMotor's
        getCombustionProperties: if any tab strictly contains P
        (min < P < max), use it; otherwise fall back to the tab with
        the closest boundary.
        """
        for tab in self.tabs:
            if tab.min_pressure < P < tab.max_pressure:
                return tab
        # Fallback: closest boundary
        best = self.tabs[0]
        best_dist = min(abs(P - best.min_pressure), abs(P - best.max_pressure))
        for tab in self.tabs[1:]:
            d = min(abs(P - tab.min_pressure), abs(P - tab.max_pressure))
            if d < best_dist:
                best = tab
                best_dist = d
        return best

    def representative_tab(self, P_expected: Optional[float] = None) -> PropellantTab:
        """
        Pick one tab to serve as the source for sim-start gas thermo
        (gamma, T_flame, molecular_weight). If `P_expected` is given,
        select the tab covering it; otherwise pick the tab with the
        widest pressure range. (Single-tab propellants always return
        their only tab.)
        """
        if len(self.tabs) == 1:
            return self.tabs[0]
        if P_expected is not None:
            return self.select_tab(P_expected)
        return max(self.tabs,
                   key=lambda t: t.max_pressure - t.min_pressure)

    def burn_rate_normal(self, P):
        """
        Normal (non-erosive) burn rate r₀ = a(P)·P^n(P) using tab lookup.
        Scalar P only (vectorize via np.vectorize if needed).
        """
        tab = self.select_tab(P)
        return tab.a * max(P, 0.0) ** tab.n

    def tab_arrays(self):
        """
        Return parallel numpy arrays of (min_p, max_p, a, n) suitable
        for passing into Numba-compiled burn-rate code. Length = len(tabs).
        """
        n_tabs = len(self.tabs)
        min_p = np.empty(n_tabs)
        max_p = np.empty(n_tabs)
        a_arr = np.empty(n_tabs)
        n_arr = np.empty(n_tabs)
        for i, t in enumerate(self.tabs):
            min_p[i] = t.min_pressure
            max_p[i] = t.max_pressure
            a_arr[i] = t.a
            n_arr[i] = t.n
        return min_p, max_p, a_arr, n_arr


# ================================================================
# Gas Properties
# ================================================================

@dataclass
class GasProperties:
    """
    Derived combustion gas properties for the flow solver.

    All transport properties (mu, k_thermal, Cp, Pr) should come from
    CEA/RPA at chamber conditions, using EFFECTIVE (not frozen) values.

    This object is created once at simulation start from the Propellant
    data and passed to the solver as a convenient container.
    """
    gamma: float            # Ratio of specific heats [-]
    molecular_weight: float # Mean molecular weight [kg/mol]
    T_flame: float          # Gas stagnation temperature [K]
    Cp: float               # Specific heat at constant pressure [J/(kg·K)]
    R_specific: float       # Specific gas constant [J/(kg·K)]
    mu: float               # Dynamic viscosity [Pa·s]
    k_thermal: float        # Thermal conductivity [W/(m·K)]
    Pr: float               # Prandtl number [-]


def create_gas_properties(gamma, molecular_weight, T_flame,
                          mu, k_thermal, Cp):
    """
    Create a GasProperties object from CEA/RPA output.

    This is the preferred way to create gas properties. All transport
    properties come directly from the equilibrium solver, with no
    kinetic-theory approximations.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats from CEA/RPA chamber conditions.
    molecular_weight : float
        Mean molecular weight [kg/mol].
    T_flame : float
        Gas stagnation temperature [K].
    mu : float
        Dynamic viscosity [Pa·s].
    k_thermal : float
        EFFECTIVE thermal conductivity [W/(m·K)].
    Cp : float
        EFFECTIVE specific heat [J/(kg·K)].

    Returns
    -------
    GasProperties
    """
    R_specific = R_UNIVERSAL / molecular_weight
    Pr = mu * Cp / k_thermal

    return GasProperties(
        gamma=gamma,
        molecular_weight=molecular_weight,
        T_flame=T_flame,
        Cp=Cp,
        R_specific=R_specific,
        mu=mu,
        k_thermal=k_thermal,
        Pr=Pr,
    )


def create_gas_properties_estimated(gamma, molecular_weight, T_flame):
    """
    FALLBACK: estimate gas properties when CEA/RPA data is unavailable.

    Uses Sutherland's law for viscosity and the modified Eucken
    correlation for thermal conductivity. Expect 30-50% errors in
    thermal conductivity and Prandtl number compared to equilibrium
    values.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats.
    molecular_weight : float
        Mean molecular weight [kg/mol].
    T_flame : float
        Gas stagnation temperature [K].

    Returns
    -------
    GasProperties
    """
    import warnings
    warnings.warn(
        "Using estimated gas transport properties (Sutherland/Eucken). "
        "For accurate erosive burning predictions, provide CEA/RPA-derived "
        "values using create_gas_properties() instead.",
        UserWarning,
        stacklevel=2,
    )

    R_specific = R_UNIVERSAL / molecular_weight
    Cp = gamma * R_specific / (gamma - 1.0)

    # Sutherland's law for viscosity
    mu_ref = 1.73e-5     # Reference viscosity [Pa·s] at T_ref
    T_ref = 300.0         # Reference temperature [K]
    S_viscosity = 240.0   # Sutherland constant for combustion gas [K]
    mu = mu_ref * (T_flame / T_ref) ** 1.5 * (T_ref + S_viscosity) / (T_flame + S_viscosity)

    # Modified Eucken for thermal conductivity (underestimates by ~2-3x)
    k_thermal = mu * (Cp + 1.25 * R_specific)

    Pr = mu * Cp / k_thermal

    return GasProperties(
        gamma=gamma,
        molecular_weight=molecular_weight,
        T_flame=T_flame,
        Cp=Cp,
        R_specific=R_specific,
        mu=mu,
        k_thermal=k_thermal,
        Pr=Pr,
    )


# ================================================================
# Thermodynamic Utilities
# ================================================================

def speed_of_sound(gamma, R_specific, T):
    """
    Local speed of sound: a = √(γ · R · T).

    Exact for an ideal gas, not an approximation.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats [-].
    R_specific : float
        Specific gas constant [J/(kg·K)].
    T : float
        Local gas temperature [K].

    Returns
    -------
    float
        Speed of sound [m/s].
    """
    return np.sqrt(gamma * R_specific * T)


def density_from_ideal_gas(P, R_specific, T):
    """
    Gas density from the ideal gas equation of state: ρ = P / (R · T).

    For the pressures and temperatures in amateur solid rocket motors
    (1-10 MPa, 1500-3500 K), the ideal gas assumption introduces
    errors of less than 1%.

    Parameters
    ----------
    P : float or ndarray
        Pressure [Pa].
    R_specific : float
        Specific gas constant [J/(kg·K)].
    T : float or ndarray
        Temperature [K].

    Returns
    -------
    float or ndarray
        Density [kg/m³].
    """
    return P / (R_specific * T)


def critical_flow_function(gamma):
    """
    Choked-flow function for isentropic nozzle flow.

        Γ = √γ · (2/(γ+1))^((γ+1) / (2·(γ-1)))

    This relates chamber pressure to mass flow through the throat:
        ṁ = P · A_throat · Γ / √(R · T)

    Used as the nozzle boundary condition in the PISO solver and
    for computing characteristic velocity c*.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats [-].

    Returns
    -------
    float
        Γ, the critical flow function [-].
    """
    return np.sqrt(gamma) * (
        2.0 / (gamma + 1.0)
    ) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))


def characteristic_velocity(gamma, R_specific, T_flame):
    """
    Characteristic velocity c* [m/s].

        c* = √(R · T_flame) / Γ(γ)

    c* is a measure of the combustion efficiency — it depends only on
    the thermochemistry (flame temperature and gas composition), not on
    the nozzle geometry.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats [-].
    R_specific : float
        Specific gas constant [J/(kg·K)].
    T_flame : float
        Flame temperature [K].

    Returns
    -------
    float
        Characteristic velocity [m/s].
    """
    Gamma = critical_flow_function(gamma)
    return np.sqrt(R_specific * T_flame) / Gamma


# ================================================================
# Propellant Definitions
# ================================================================
#
# Named propellants live in ``srm_1d/motors/<motor>.ric`` (combustion
# data) and ``srm_1d/motors/<motor>.transport.yaml`` (transport
# properties). Use ``run_from_ric(...)`` from openmotor_adapter or
# ``convert_propellant`` directly.
