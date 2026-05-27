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
    heat_flux_cal_cm2_s: Optional[float] = None  # DeMar heat flux [cal/(cm^2*s)]
    Cp_gas: Optional[float] = None  # Optional explicit product-gas Cp [J/(kg*K)]
    # v0.7.2 Phase A: pyrogen-jet axial distribution. L_jet = kappa_jet * d_throat
    # sets the characteristic exponential-decay length over which pyrogen
    # mass / enthalpy / momentum are deposited into the bore (rather than the
    # v0.7.1 cell-0-only model). Defensible defaults from coaxial-jet
    # theory + sonic-injection literature: kappa_jet ~ 6-10 for choked
    # sonic axial-vent pyrogens (Witze 1974 potential-core length;
    # Hersch/Rieser 1971 NTRS 19710018794); ~ 2-4 for predominantly
    # radial-vent pyrogens (JICF at q < 5). kappa_jet = 0 recovers the
    # v0.7.1 cell-0-only behavior byte-for-byte. See
    # srm_1d/docs/v0_7_2/candidates/03_pyrogen_spatial_distribution.md.
    kappa_jet: float = 8.0

    # v0.7.3 Phase B.3: physical form of the pyrogen charge.
    # v0.7.4 Phase C.1: form is now informational only — actual A_burn
    # is computed in build_pyrogen_chamber from the explicit physical
    # particle dimensions below (`particle_diameter_m`,
    # `particle_LD_ratio`). The form field is retained as a YAML
    # convention marker (e.g., 'pellets' for BKNO3/MTV) and may be
    # consumed by future preset / lookup logic, but the legacy
    # ×1/×5/×20 multiplier dispatch is removed.
    # See `[[pyrogen-form-archetypes]]` memory.
    form: str = 'pellets'

    # v0.7.4 Phase C.1: explicit physical particle dimensions for the
    # pyrogen charge. Replaces the v0.7.3 form-archetype A_burn
    # multipliers with first-principles geometry:
    #   sphere    (LD_ratio <= 1):  A = 6 · m / (ρ · d)
    #   cylinder  (LD_ratio  > 1):  A = m · (4·λ + 2) / (ρ · λ · d)
    # where λ = particle_LD_ratio. The user supplies a characteristic
    # diameter d for each pyrogen, plus L/D = 1.0 (sphere) for
    # powders and L/D >= 2.0 for pellets / slivers. Defensible defaults
    # (per amateur HPR practice):
    #   - BKNO3 (BPNV) pellets: d ≈ 5 mm, L/D ≈ 3 (cylinder)
    #   - MTV pellets:          d ≈ 5 mm, L/D ≈ 3 (cylinder)
    #   - Cu/Al thermite:       d ≈ 100-200 µm sphere (mesh-dependent;
    #                                                  agglomeration uncertain)
    # Override via explicit `pyrogen_burn_area=` kwarg always wins.
    particle_diameter_m: float = 5.0e-3
    particle_LD_ratio: float = 3.0

    # v0.7.3 Phase B.4: pyrogen-to-propellant heat delivery mode for
    # uncontained topologies (head_basket, aft_basket). Mutually
    # exclusive modes (no double-counting of the radiative component):
    #   'demar'      — apply Pyrogen.heat_flux_cal_cm2_s as a
    #                  time-averaged total flux distributed across
    #                  cartridge cells. Empirically grounded for
    #                  pellet pyrogens (BKNO3 BPNV, MTV) where DeMar
    #                  2021 measured. Doesn't model distance falloff.
    #   'radiation'  — Stefan-Boltzmann pellet emission:
    #                  q = sigma * pellet_emissivity * T_flame^4 *
    #                      F_view * exp(-d / radiation_absorption_length_m)
    #                  where F_view is the geometric view factor
    #                  A_port_j / (4*pi*(x_j - x_i)^2 + A_port_j).
    #                  Physically modeled; extensible to powder/chunks.
    #   'none'       — neither pathway. Recovers v0.7.3-phaseA
    #                  byte-for-byte; backward-compat regression target.
    # Default 'demar' for pellet pyrogens (best empirical grounding);
    # powder/chunks should override to 'radiation' in their YAMLs.
    heat_delivery_mode: str = 'demar'

    # v0.7.3 Phase B.4 (radiation mode only): pellet emissivity at
    # T_flame. Lit range 0.5-0.9 for ceramic/metal-oxide combustion
    # products; 0.7 is a defensible mid-range default. Highly tunable.
    pellet_emissivity: float = 0.7

    # v0.7.3 Phase B.4 (radiation mode only): aggregate radiation
    # absorption length capturing combined gas + particle attenuation.
    # Clean (non-aluminized) pyrogen exhaust: ~1 m. Aluminized exhaust
    # (Al2O3 droplet fog): ~0.1 m. Default 1.0 m matches clean
    # pellet pyrogens (BKNO3, MTV without metal).
    radiation_absorption_length_m: float = 1.0

    @property
    def species(self) -> "GasSpecies":
        """
        Return this pyrogen's combustion-product gas as a GasSpecies.

        If Cp_gas is not given explicitly, derive it from the ideal-gas
        identity Cp = gamma * R_specific / (gamma - 1).
        """
        if self.Cp_gas is not None:
            cp = self.Cp_gas
        else:
            R_specific = R_UNIVERSAL / self.M
            cp = self.gamma * R_specific / (self.gamma - 1.0)
        return GasSpecies(
            name=f"{self.name}_gas",
            gamma=self.gamma,
            Cp=cp,
            molecular_weight=self.M,
            T_flame=self.T_flame,
        )


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
    radiation_emissivity : float
        Effective particle/flame emissivity for adjacent-cell ignition
        radiation. This is a material property, not a heat-transfer
        multiplier; use 0.0 to disable the radiation path.
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
    radiation_emissivity: float = 0.0
    # v0.7.2 Phase B (DEFAULT DISABLED, opt-in only): flame-front h_c
    # augmentation. When True, an unignited cell's Bartz h_c is
    # multiplied by flame_spread_boost IF its immediate upstream
    # neighbor ignited within the last flame_spread_tau seconds.
    #
    # **Empirically disabled by default** after Phase B-v1 (cumulative-G)
    # and Phase B-v2 (flame-front gating) both showed the augmentation
    # AMPLIFIES the ignition spike rather than smoothing it in this
    # codebase — because PISO's local-Re tracking already captures
    # upstream-mass-flux contributions to h_c at unignited cells, the
    # extra augmentation is double-counting. Phase A (pyrogen axial
    # distribution) is the load-bearing v0.7.2 ship; Phase B remains
    # available as diagnostic infrastructure for experimentation /
    # future iteration. See docs/v0_7_2/candidates/02_*.md and the
    # negative findings at commits 065d193 (v1) and (v2 commit pending).
    flame_spread_enabled: bool = False
    flame_spread_tau: float = 1.0e-3   # [s] window after upstream ignition
    flame_spread_boost: float = 3.0    # h_c multiplier when boost active

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

    def species(self, P_expected: Optional[float] = None) -> "GasSpecies":
        """
        Return this propellant's combustion-product gas as a GasSpecies,
        built from the representative tab. Multi-tab Cp(p) variation is
        not yet modeled (gas thermo is taken at the representative tab;
        burn-rate APN coefficients still switch per-tab).
        """
        tab = self.representative_tab(P_expected)
        return GasSpecies(
            name=f"{self.name}_gas",
            gamma=tab.gamma,
            Cp=self.Cp_gas,
            molecular_weight=tab.molecular_weight,
            T_flame=tab.T_flame,
        )


# ================================================================
# Gas Species (per-cell mixture component)
# ================================================================

@dataclass
class GasSpecies:
    """
    Lightweight bulk-flow thermo for one gas species in an N-species
    mixture.

    Introduced in v0.7.1 to support per-cell variable gamma/Cp/R via the
    SPINBALL-style "infinite-gases mixture" formulation (see
    ``docs/v0_7_1/DESIGN.md``). Each cell carries mass fractions
    ``Y[i, s]`` for s = 0..S-1; per-cell (gamma, Cp, R, M) derive from
    these via standard ideal-gas mixing rules.

    A GasSpecies captures only the **bulk-flow** properties (gamma, Cp,
    M, T_flame). Burn-rate APN coefficients stay on Pyrogen /
    PropellantTab. Transport properties (k_thermal, mu_gas, Pr) remain
    scalar / single-source in v0.7.1 and are NOT carried per species;
    they live at the Propellant level.

    Attributes
    ----------
    name : str
        Human-readable identifier (e.g. "hasegawa_prop1_gas",
        "bpnv_pyrogen_gas"). Used in diagnostics only.
    gamma : float
        Ratio of specific heats [-].
    Cp : float
        Specific heat at constant pressure [J/(kg*K)].
    molecular_weight : float
        Mean molecular weight [kg/mol].
    T_flame : float
        Adiabatic flame temperature [K]. Used only as the source-term
        injection temperature when this species is added to a cell by
        combustion; it is NOT a bulk-mixture state variable.
    """
    name: str
    gamma: float
    Cp: float
    molecular_weight: float
    T_flame: float

    @property
    def R_specific(self) -> float:
        return R_UNIVERSAL / self.molecular_weight


def ambient_air_species(T_ambient: float = 298.15) -> "GasSpecies":
    """
    Default ambient-air GasSpecies for bore pre-fill.

    Used as the v0.7.1 default for species index 2 (pre-fill / ambient).
    The species has no continuing source after t=0: it is an initial
    condition only, purged through the nozzle during chamber fill.

    Defaults are dry-air values at standard conditions:
        gamma = 1.40
        Cp    = 1005 J/(kg*K)
        M     = 0.02897 kg/mol
        T_flame = T_ambient (used only as a "source temperature" if this
                  species were ever re-injected; nominally unused)

    Parameters
    ----------
    T_ambient : float
        Ambient temperature for T_flame assignment [K]. Default 298.15 K.

    Returns
    -------
    GasSpecies
    """
    return GasSpecies(
        name="ambient_air",
        gamma=1.40,
        Cp=1005.0,
        molecular_weight=0.02897,
        T_flame=float(T_ambient),
    )


def species_array(species_list):
    """
    Pack a list of GasSpecies into a 2D numpy array suitable for passing
    into Numba kernels. Column layout:
        column 0: gamma
        column 1: Cp
        column 2: molecular_weight
        column 3: T_flame

    Parameters
    ----------
    species_list : list[GasSpecies]
        Ordered species. Index 0 conventionally reserved for the igniter
        species; index 1 for the main grain combustion species; higher
        indices for future sources (head-end motor, ablation, ...).

    Returns
    -------
    np.ndarray[S, 4]
        Float64 array. Row s carries (gamma, Cp, M, T_flame) for species s.
    """
    n = len(species_list)
    if n == 0:
        raise ValueError("species_array: at least one species required")
    arr = np.empty((n, 4), dtype=np.float64)
    for s, sp in enumerate(species_list):
        if sp.gamma <= 1.0:
            raise ValueError(
                f"species[{s}] '{sp.name}': gamma must be > 1 (got {sp.gamma})")
        if sp.Cp <= 0.0:
            raise ValueError(
                f"species[{s}] '{sp.name}': Cp must be > 0 (got {sp.Cp})")
        if sp.molecular_weight <= 0.0:
            raise ValueError(
                f"species[{s}] '{sp.name}': molecular_weight must be > 0 "
                f"(got {sp.molecular_weight})")
        arr[s, 0] = sp.gamma
        arr[s, 1] = sp.Cp
        arr[s, 2] = sp.molecular_weight
        arr[s, 3] = sp.T_flame
    return arr


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
