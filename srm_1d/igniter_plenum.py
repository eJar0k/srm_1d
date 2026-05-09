"""
igniter_plenum.py -- Standalone v0.7.0 pyrogen plenum model
============================================================

Phase 1 implementation only: a forward 0D pyrogen chamber with a
Saint-Robert burn law and choked/subsonic vent to a downstream chamber.
This module is intentionally not wired into simulation.py yet.
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


@dataclass
class PyrogenChamber:
    """
    Standalone pyrogen chamber configuration.

    Parameters are SI. ``burn_law`` is validated at the Python boundary
    and passed into Numba kernels as an integer code.
    """
    pyrogen: Pyrogen
    m_pyrogen_initial: float
    A_burn_initial: float
    A_throat: float
    V_plenum: float
    burn_law: str = "0d"

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


def _burn_law_code(burn_law):
    if burn_law == "0d":
        return BURN_LAW_0D
    if burn_law == "end_burning":
        return BURN_LAW_END_BURNING
    raise ValueError("burn_law must be '0d' or 'end_burning'")


def pyrogen_params(pyrogen):
    """Return pyrogen scalars as a Numba-friendly array."""
    return np.array([
        pyrogen.a,
        pyrogen.n,
        pyrogen.rho,
        pyrogen.T_flame,
        pyrogen.M,
        pyrogen.gamma,
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
    dt_limit, P_main,
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

    dm_pyrogen_dt = -mdot_generated
    dm_gas_dt = mdot_generated - mdot_out

    m_eff = m_gas
    if m_eff < 1e-12:
        m_eff = 1e-12
    dT_dt = (
        gamma * (mdot_generated * T_flame - mdot_out * T_gas)
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
        burn_law_code, dt, P_main,
    )
    k2_m, k2_g, k2_T = _plenum_rhs(
        m0 + 0.5 * dt * k1_m,
        g0 + 0.5 * dt * k1_g,
        T0 + 0.5 * dt * k1_T,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, dt, P_main,
    )
    k3_m, k3_g, k3_T = _plenum_rhs(
        m0 + 0.5 * dt * k2_m,
        g0 + 0.5 * dt * k2_g,
        T0 + 0.5 * dt * k2_T,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, dt, P_main,
    )
    k4_m, k4_g, k4_T = _plenum_rhs(
        m0 + dt * k3_m,
        g0 + dt * k3_g,
        T0 + dt * k3_T,
        a, n, rho_pyrogen, T_flame, M, gamma,
        m_pyrogen_initial, A_burn_initial, A_throat, V_plenum,
        burn_law_code, dt, P_main,
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
