"""
igniter_plenum.py -- Standalone v0.7.0 pyrogen plenum model
============================================================

Forward 0D pyrogen chamber with a Saint-Robert burn law and
choked/subsonic vent to a downstream chamber. v0.7.0 Phase 3 wires this
module into ``simulation.py`` as the hot-gas igniter source.
"""

import numpy as np
from dataclasses import dataclass

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper

from .propellant import Pyrogen, R_UNIVERSAL


BURN_LAW_0D = 0
BURN_LAW_END_BURNING = 1


_VALID_TOPOLOGIES = (
    'forward_plenum',     # v0.7.0+ default: plenum upstream of bore, choked vent
    'head_basket',        # uncontained pyrogen in head-end bore cells [0, i_end]
    'aft_basket',         # uncontained pyrogen in aft bore cells [i_start, N-1]
)

# Numba-compatible topology codes (passed into the time loop as an int)
TOPOLOGY_FORWARD_PLENUM = 0
TOPOLOGY_HEAD_BASKET    = 1
TOPOLOGY_AFT_BASKET     = 2


def _topology_code(topology):
    if topology == 'forward_plenum':
        return TOPOLOGY_FORWARD_PLENUM
    if topology == 'head_basket':
        return TOPOLOGY_HEAD_BASKET
    if topology == 'aft_basket':
        return TOPOLOGY_AFT_BASKET
    raise ValueError(f"unknown injection_topology: {topology!r}")


@dataclass
class PyrogenChamber:
    """
    Pyrogen chamber / igniter configuration. v0.7.3 Phase A supports
    two physical containment models selected by ``injection_topology``:

    **Plenum-containment model** (``'forward_plenum'``, v0.7.0+ default):
    pyrogen burns inside an enclosed volume with its own internal
    pressure ``P_ig`` (decoupled from bore pressure by the orifice),
    Saint-Robert burn rate evaluated at ``P_ig``, choked/subsonic
    venting through ``A_throat`` into bore cell 0 (distributed via
    Phase A exponential decay if ``pyrogen.kappa_jet > 0``). Axial
    momentum injected at face 1 with downstream sign.

    **Uncontained model** (``'head_basket'`` | ``'aft_basket'``):
    pyrogen pellets sit physically inside the bore — no plenum
    chamber wall, no internal pressure separation, no defined orifice
    or burst threshold. Each pellet burns at the LOCAL bore pressure
    of its host cell: ``r_b[i] = a · P_bore[i]^n``. Per-cell mass
    enters the bore directly via ``mdot[i] = rho_p · r_b[i] ·
    A_burn_per_cell`` where ``A_burn_per_cell = A_burn_initial /
    n_cartridge_cells`` (initial burning surface distributed
    uniformly across the cartridge's axial extent). No axial momentum
    injection — PISO handles flow naturally via pressure gradient
    between high-P pyrogen cells and surrounding bore. Plenum-state
    fields (``A_throat``, ``V_plenum``, ``A_burn_initial``) are
    repurposed: ``A_burn_initial`` becomes total burning surface
    distributed across cells; ``V_plenum`` and ``A_throat`` are
    ignored for the uncontained-burn computation but still validated
    at the Python boundary so existing motor configs don't break.

    This split represents general physical practice for amateur
    high-power-rocketry pyrogen configurations. The ``head_basket``
    topology models a head-end BKNO3 / MTV pellet pack glued or taped
    to the forward bulkhead — a plausible amateur build. A future
    ``'aft_fore_firing'`` topology (deferred from v0.7.3 Phase A)
    would add opt-in upstream momentum injection at face i_start to
    model nozzle-inserted igniter cartridges (the classic Super Loki
    factory configuration per NASA CR-61238 cross-section, the
    AeroTech FirstFire / Firestar line, and amateur Loki Research-
    style HPR pyrogens).

    **Provenance correction (2026-05-25)**: previous versions of
    this docstring cited "NASA CR-61238, MIT Super Loki Report,
    Smithsonian/NASM" as support for modeling ISP Super Loki as
    head_basket. The lit dive in
    ``docs/v0_7_3/references/super_loki_igniter_lit_dive.md`` could
    not corroborate that citation chain — the public Super Loki
    sources describe the factory igniter as "separable and installed
    at the launch site" (i.e. nozzle-inserted, not head-end basket).
    The RCS Rocket Motor Components Super Loki recreation that
    ``examples/ISP_Super_Loki.py`` models does NOT ship with an
    igniter; the user assembles their own. head_basket is therefore
    a defensible amateur-build choice, but not "the" factory
    topology.

    Topology options (``injection_topology``):

    - ``'forward_plenum'`` (default, v0.7.0-v0.7.2 behavior): plenum
      upstream of bore.

    - ``'head_basket'``: uncontained pyrogen in cells ``[0, i_end]``
      where ``i_end`` is set by ``cartridge_length_m``. Cells burn at
      local bore P; mass + enthalpy enter their host cells. No
      momentum injection. Models amateur head-end pellet packs.

    - ``'aft_basket'``: uncontained pyrogen in cells
      ``[i_start, N-1]``. Same physics as head_basket, just at the
      aft end of the bore. Diagnostic value: tests whether the
      simultaneous-ignition artifact is driven by mass-injection
      position (head vs aft). Empirically inadequate as a startup
      mechanism on its own — see v0.7.3 Phase B validation findings.

    Cartridge length: if ``cartridge_length_m < 0`` (default), derived
    from pyrogen mass and bore geometry as
    ``L_cart = m_pyrogen_initial / (rho_pyrogen * A_port_avg)`` where
    A_port_avg is the bore-volume-weighted port area. User override
    via explicit positive ``cartridge_length_m`` always wins.

    See srm_1d/docs/v0_7_2/candidates_post_phaseA.md (v0.7.3 design)
    and the Super Loki igniter lit-check notes in the v0.7.3
    DESIGN doc.
    """
    pyrogen: Pyrogen
    m_pyrogen_initial: float
    A_burn_initial: float
    A_throat: float
    V_plenum: float
    burn_law: str = "0d"
    # v0.7.3 Phase A — submerged-igniter topology fields
    injection_topology: str = 'forward_plenum'
    cartridge_length_m: float = -1.0  # < 0 => derive from pyrogen mass
    # v0.7.4 — realistic basket cartridge geometry (replaces the solid-puck
    # L_cart). Used by head_basket / aft_basket only when cartridge_length_m
    # < 0. A real floating igniter / basket fills only a FRACTION of the bore
    # cross-section, and loose pellets pack at a fraction of solid density:
    #   L_cart = m / ( (pellet_packing_fraction·ρ) · (basket_fill_fraction·A_port) )
    basket_fill_fraction: float = 0.5      # basket cross-section ÷ bore cross-section
    pellet_packing_fraction: float = 0.60  # bulk(tap) density ÷ solid density
    #   Packing is shape-dependent and ~size-INDEPENDENT for monodisperse
    #   convex particles (random loose packing ≈ 0.60 for spheres / L≈D
    #   cylinders; size enters only weakly via wall effects). The strong
    #   size-dependence lives in the burn AREA (A/m = 6/ρd), set separately
    #   by particle_diameter_m. So particle diameter is the single unifying
    #   knob: it sets specific surface area, while fill/packing are fixed
    #   structural assumptions.

    def __post_init__(self):
        if self.m_pyrogen_initial <= 0.0:
            raise ValueError("m_pyrogen_initial must be positive")
        if self.A_burn_initial < 0.0:
            raise ValueError("A_burn_initial must be nonnegative")
        if self.A_throat < 0.0:
            raise ValueError("A_throat must be nonnegative")
        if self.V_plenum <= 0.0:
            raise ValueError("V_plenum must be positive")
        _burn_law_code(self.burn_law)
        if self.injection_topology not in _VALID_TOPOLOGIES:
            raise ValueError(
                f"injection_topology must be one of {_VALID_TOPOLOGIES}; "
                f"got '{self.injection_topology}'"
            )
        if not (0.0 < self.basket_fill_fraction <= 1.0):
            raise ValueError(
                "basket_fill_fraction must be in (0, 1]; "
                f"got {self.basket_fill_fraction}"
            )
        if not (0.0 < self.pellet_packing_fraction <= 1.0):
            raise ValueError(
                "pellet_packing_fraction must be in (0, 1]; "
                f"got {self.pellet_packing_fraction}"
            )

    def resolve_cartridge_length(self, A_port_avg):
        """Return cartridge length in meters, deriving from pyrogen mass
        if not user-specified.

        v0.7.4 — realistic packed-bed cartridge (replaces the old solid-puck
        ``L_cart = m / (ρ · A_port_avg)`` which assumed the pyrogen filled the
        entire bore cross-section at full material density):

            L_cart = m / ( (φ_pack·ρ) · (f_area·A_port_avg) )

        where ``φ_pack = pellet_packing_fraction`` (loose-pellet bulk density
        ÷ solid density) and ``f_area = basket_fill_fraction`` (basket
        cross-section ÷ bore cross-section). Both < 1 → the cartridge is
        longer (and spans more cells) than the old puck assumption.

        A_port_avg is bore-volume-weighted port area (caller supplies).
        """
        if self.cartridge_length_m > 0.0:
            return float(self.cartridge_length_m)
        if A_port_avg <= 0.0:
            raise ValueError(
                "Cannot derive cartridge length from pyrogen mass: "
                "A_port_avg must be positive"
            )
        rho_bulk = self.pellet_packing_fraction * self.pyrogen.rho
        A_basket = self.basket_fill_fraction * A_port_avg
        return float(self.m_pyrogen_initial / (rho_bulk * A_basket))

    def resolve_injection_cells(self, x_centers, dx, N, A_port):
        """Compute (i_start, i_end) cell-index range for pyrogen
        injection given bore geometry. Inclusive of both endpoints.

        Returns (i_start, i_end) for submerged topologies, or (0, 0)
        for forward_plenum (legacy default; the actual Phase A axial
        decay is applied elsewhere).
        """
        if self.injection_topology == 'forward_plenum':
            return (0, 0)
        # Bore-volume-weighted A_port for cartridge length derivation
        total_volume = 0.0
        for i in range(N):
            total_volume += A_port[i] * dx
        bore_length = N * dx
        A_port_avg = total_volume / bore_length if bore_length > 0 else 0.0
        L_cart = self.resolve_cartridge_length(A_port_avg)
        # Snap to whole cells: find smallest n such that n*dx >= L_cart
        n_cart = max(1, int(np.ceil(L_cart / dx)))
        n_cart = min(n_cart, N)
        if self.injection_topology == 'head_basket':
            return (0, n_cart - 1)
        # aft_basket
        return (N - n_cart, N - 1)


def _burn_law_code(burn_law):
    if burn_law == "0d":
        return BURN_LAW_0D
    if burn_law == "end_burning":
        return BURN_LAW_END_BURNING
    raise ValueError("burn_law must be '0d' or 'end_burning'")


def pyrogen_params(pyrogen):
    """Return pyrogen scalars as a Numba-friendly array.

    Index 6 (gas_mass_fraction) is the condensed-phase split: only this
    fraction of the burned SOLID mass becomes pressure-generating gas.
    """
    return np.array([
        pyrogen.a,
        pyrogen.n,
        pyrogen.rho,
        pyrogen.T_flame,
        pyrogen.M,
        pyrogen.gamma,
        pyrogen.gas_mass_fraction,
    ], dtype=np.float64)


def chamber_params(chamber):
    """Return chamber scalars as a Numba-friendly array."""
    return np.array([
        chamber.m_pyrogen_initial,
        chamber.A_burn_initial,
        chamber.A_throat,
        chamber.V_plenum,
        float(_burn_law_code(chamber.burn_law)),
    ], dtype=np.float64)


def initial_plenum_state(chamber, P_initial=101325.0, T_initial=300.0):
    """
    Build the initial state vector ``[m_pyrogen, m_gas, T_gas]``.

    The initial gas mass is the ullage gas needed to produce
    ``P_initial`` in the fixed plenum volume.
    """
    if P_initial < 0.0:
        raise ValueError("P_initial must be nonnegative")
    if T_initial <= 0.0:
        raise ValueError("T_initial must be positive")
    m_gas = P_initial * chamber.pyrogen.M * chamber.V_plenum / (
        R_UNIVERSAL * T_initial
    )
    return np.array([
        chamber.m_pyrogen_initial,
        m_gas,
        T_initial,
    ], dtype=np.float64)


def step_plenum(chamber, state, dt, P_main):
    """
    Python wrapper around the Numba RK4 plenum step.

    Returns ``(new_state, mdot_out, mdot_generated, P_ig)``.
    """
    return _step_plenum_ode(
        state,
        pyrogen_params(chamber.pyrogen),
        chamber_params(chamber),
        dt,
        P_main,
    )


def sutton_pyrogen_mass(V_free_in3):
    """
    Sutton Eq. 15-4 default igniter charge mass.

    ``V_free_in3`` is the motor free volume in cubic inches. Returns kg.
    """
    if V_free_in3 < 0.0:
        raise ValueError("V_free_in3 must be nonnegative")
    grams = 0.12 * V_free_in3 ** 0.7
    return grams * 1e-3


@njit(cache=True)
def _plenum_pressure(m_gas, T_gas, M, V_plenum):
    if m_gas <= 0.0 or T_gas <= 0.0 or M <= 0.0 or V_plenum <= 0.0:
        return 0.0
    return m_gas * R_UNIVERSAL * T_gas / (M * V_plenum)


@njit(cache=True)
def _critical_pressure_ratio(gamma):
    if gamma <= 1.0:
        return 0.0
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


@njit(cache=True)
def _burn_area(m_pyrogen, m_pyrogen_initial, A_burn_initial, burn_law_code):
    if m_pyrogen <= 0.0 or m_pyrogen_initial <= 0.0:
        return 0.0
    if burn_law_code == BURN_LAW_0D:
        frac = m_pyrogen / m_pyrogen_initial
        if frac < 0.0:
            frac = 0.0
        return A_burn_initial * frac ** (2.0 / 3.0)
    if burn_law_code == BURN_LAW_END_BURNING:
        return A_burn_initial
    return 0.0


@njit(cache=True)
def _choked_orifice_mdot(P_ig, T_ig, A_t, gamma, R, M, P_main):
    """
    Choked or subsonic ideal-gas orifice mass flow.

    Returns zero for invalid states, zero area, or non-positive pressure
    drop. The gas constant ``R`` is the universal gas constant.
    """
    if A_t <= 0.0 or P_ig <= 0.0 or T_ig <= 0.0:
        return 0.0
    if gamma <= 1.0 or R <= 0.0 or M <= 0.0:
        return 0.0
    if P_main >= P_ig:
        return 0.0

    pressure_ratio = P_main / P_ig
    if pressure_ratio < 0.0:
        pressure_ratio = 0.0

    crit = _critical_pressure_ratio(gamma)
    if pressure_ratio < crit:
        Gamma = (gamma * (2.0 / (gamma + 1.0)) **
                 ((gamma + 1.0) / (gamma - 1.0))) ** 0.5
        return P_ig * A_t * Gamma / ((R * T_ig / M) ** 0.5)

    term = pressure_ratio ** (2.0 / gamma) - pressure_ratio ** (
        (gamma + 1.0) / gamma
    )
    if term <= 0.0:
        return 0.0
    coeff = (2.0 * gamma / (gamma - 1.0)) * M / (R * T_ig)
    return P_ig * A_t * (coeff ** 0.5) * (term ** 0.5)


@njit(cache=True)
def _plenum_rates(
    m_pyrogen, m_gas, T_gas,
    a, n, rho_pyrogen, T_flame, M, gamma,
    m_pyrogen_initial, A_burn_initial, A_throat, V_plenum, burn_law_code,
    P_main,
):
    P_ig = _plenum_pressure(m_gas, T_gas, M, V_plenum)
    A_burn = _burn_area(
        m_pyrogen, m_pyrogen_initial, A_burn_initial, int(burn_law_code)
    )
    if A_burn <= 0.0 or a <= 0.0 or rho_pyrogen <= 0.0:
        mdot_generated = 0.0
    else:
        mdot_generated = rho_pyrogen * A_burn * a * max(P_ig, 0.0) ** n
    mdot_out = _choked_orifice_mdot(
        P_ig, T_gas, A_throat, gamma, R_UNIVERSAL, M, P_main
    )
    return mdot_generated, mdot_out, P_ig


@njit(cache=True)
def _plenum_rhs(
    m_pyrogen, m_gas, T_gas,
    a, n, rho_pyrogen, T_flame, M, gamma,
    m_pyrogen_initial, A_burn_initial, A_throat, V_plenum, burn_law_code,
    dt_limit, P_main, gas_mass_fraction,
):
    if m_pyrogen < 0.0:
        m_pyrogen = 0.0
    if m_gas < 0.0:
        m_gas = 0.0
    if T_gas < 1.0:
        T_gas = 1.0

    mdot_generated, mdot_out, _P_ig = _plenum_rates(
        m_pyrogen, m_gas, T_gas,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, P_main,
    )
    if dt_limit > 0.0 and mdot_generated * dt_limit > m_pyrogen:
        mdot_generated = m_pyrogen / dt_limit

    # Condensed-phase split: the SOLID depletes at the full burn rate
    # (mdot_generated), but only gas_mass_fraction of it becomes gas that
    # accumulates in / vents from the plenum. (1 - gas_mass_fraction) is
    # condensed product that carries no pressure. mdot_out is the orifice
    # vent of the (gas-only) plenum inventory and needs no extra factor.
    mdot_gas_in = gas_mass_fraction * mdot_generated

    dm_pyrogen_dt = -mdot_generated
    dm_gas_dt = mdot_gas_in - mdot_out

    m_eff = m_gas
    if m_eff < 1e-12:
        m_eff = 1e-12
    dT_dt = (
        gamma * (mdot_gas_in * T_flame - mdot_out * T_gas)
        - T_gas * dm_gas_dt
    ) / m_eff

    return dm_pyrogen_dt, dm_gas_dt, dT_dt


@njit(cache=True)
def _step_plenum_ode(state, pyrogen_params_arr, chamber_params_arr, dt, P_main):
    """
    RK4 step for the standalone pyrogen-plenum ODE.

    State is ``[m_pyrogen_remaining, m_gas, T_gas]``. Returns
    ``(new_state, mdot_out, mdot_generated, P_ig)`` evaluated at the
    updated state.
    """
    a = pyrogen_params_arr[0]
    n = pyrogen_params_arr[1]
    rho_pyrogen = pyrogen_params_arr[2]
    T_flame = pyrogen_params_arr[3]
    M = pyrogen_params_arr[4]
    gamma = pyrogen_params_arr[5]
    gas_mass_fraction = pyrogen_params_arr[6]

    m_pyrogen_initial = chamber_params_arr[0]
    A_burn_initial = chamber_params_arr[1]
    A_throat = chamber_params_arr[2]
    V_plenum = chamber_params_arr[3]
    burn_law_code = chamber_params_arr[4]

    m0 = state[0]
    g0 = state[1]
    T0 = state[2]

    if dt <= 0.0:
        new_state = np.empty(3, dtype=np.float64)
        new_state[0] = max(m0, 0.0)
        new_state[1] = max(g0, 0.0)
        new_state[2] = max(T0, 1.0)
        mdot_generated, mdot_out, P_ig = _plenum_rates(
            new_state[0], new_state[1], new_state[2],
            a, n, rho_pyrogen, T_flame, M, gamma,
            m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
            burn_law_code, P_main,
        )
        return new_state, mdot_out, mdot_generated, P_ig

    k1_m, k1_g, k1_T = _plenum_rhs(
        m0, g0, T0,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, dt, P_main, gas_mass_fraction,
    )
    k2_m, k2_g, k2_T = _plenum_rhs(
        m0 + 0.5 * dt * k1_m,
        g0 + 0.5 * dt * k1_g,
        T0 + 0.5 * dt * k1_T,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, dt, P_main, gas_mass_fraction,
    )
    k3_m, k3_g, k3_T = _plenum_rhs(
        m0 + 0.5 * dt * k2_m,
        g0 + 0.5 * dt * k2_g,
        T0 + 0.5 * dt * k2_T,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, dt, P_main, gas_mass_fraction,
    )
    k4_m, k4_g, k4_T = _plenum_rhs(
        m0 + dt * k3_m,
        g0 + dt * k3_g,
        T0 + dt * k3_T,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, dt, P_main, gas_mass_fraction,
    )

    new_m = m0 + (dt / 6.0) * (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m)
    new_g = g0 + (dt / 6.0) * (k1_g + 2.0 * k2_g + 2.0 * k3_g + k4_g)
    new_T = T0 + (dt / 6.0) * (k1_T + 2.0 * k2_T + 2.0 * k3_T + k4_T)

    if new_m < 0.0:
        new_g += new_m
        new_m = 0.0
    if new_g < 0.0:
        new_g = 0.0
    if new_T < 1.0:
        new_T = 1.0

    new_state = np.empty(3, dtype=np.float64)
    new_state[0] = new_m
    new_state[1] = new_g
    new_state[2] = new_T

    mdot_generated, mdot_out, P_ig = _plenum_rates(
        new_m, new_g, new_T,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, P_main,
    )
    return new_state, mdot_out, mdot_generated, P_ig
