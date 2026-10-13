"""EN1 -- verification of ``solver.py`` against exact solutions, both energy forms.

Drives ``_piso_step_with_energy_diagnostics`` directly (the solver module is pure
numerics: no propellant, no grain).  Four problems, each run with
``total_energy=False`` (the shipped sensible-enthalpy step) and ``total_energy=True``
(the EN1 conservative total-energy step):

  1. Sod shock tube (ideal gas, gamma = 1.4).   Exact Riemann solution (Toro, ch. 4).
  2. Steady isentropic quasi-1-D converging-diverging flow, supersonic exit.
  3. Steady normal shock in the divergent section (back pressure chosen so that the
     shock stands at a prescribed station); exact quasi-1-D shock relations.
  4. Closed-vessel filling: exact energy balance (adiabatic, uniform).

Boundary conditions and how they map onto the motor
---------------------------------------------------
* Head end = the motor's head closure: wall, ``u[0] = 0`` (hard-wired in the kernel).
  The stagnation reservoir of problems 2-3 is realised exactly like propellant burning:
  a mass + enthalpy SOURCE in the first cell (``mass_source``, ``thermal_source =
  mdot*Cp*T0``; injected gas is at rest, so its enthalpy is the stagnation enthalpy).
* Nozzle end = the motor's nozzle boundary ``_nozzle_boundary_flow`` (unchanged code):
    - problems 1 and 4: ``A_throat = 0`` -> the boundary passes no flow (closed end);
    - problems 2 and 3: ``A_throat = 1000*A(x_last)`` -> the same subsonic-orifice
      branch of the boundary acts as a pressure outlet (cell pressure = ``P_ambient`` to
      ~1e-6 relative; the pressure drop that carries the flow is mdot^2/(2 rho A_big^2)).
      This is the standard way to hold a back pressure; the choked branch of the motor
      boundary (a lumped throat, not a resolved duct) is not exercised by a
      resolved-duct problem and is verified in the motor tables instead.
* The momentum equation, pressure correctors, P5 coupling (K = 1, S = True) are the
  shipped defaults for BOTH forms; only the energy step differs.
"""
import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from srm_1d.solver import (  # noqa: E402
    _piso_step_with_energy_diagnostics as piso, compute_dt_cfl, _nozzle_boundary_flow)

GAMMA = 1.4
R_GAS = 287.0
CP = GAMMA * R_GAS / (GAMMA - 1.0)
CV = CP - R_GAS


# ----------------------------------------------------------------------------
# generic stepping helper
# ----------------------------------------------------------------------------
class Tube:
    """Minimal 1-D duct state + time stepping through the shipped PISO kernel."""

    def __init__(self, A_cell, dx, rho, u, P, T, gamma=GAMMA, R=R_GAS,
                 mass_source=None, thermal_source=None, A_throat=0.0,
                 P_ambient=1.0e5, T_ambient=1.0, total_energy=False,
                 cfl=0.3, K=1, S=True, dt_fixed=None, cons_mom=False):
        N = len(A_cell)
        self.N, self.dx = N, dx
        self.A = np.asarray(A_cell, float)
        self.gamma, self.R = gamma, R
        self.Cp = gamma * R / (gamma - 1.0)
        self.rho, self.u = rho.copy(), u.copy()
        self.P, self.T = P.copy(), T.copy()
        z = np.zeros(N)
        self.msrc = z.copy() if mass_source is None else mass_source
        self.tsrc = z.copy() if thermal_source is None else thermal_source
        self.A_throat, self.P_amb, self.T_amb = A_throat, P_ambient, T_ambient
        self.total, self.cfl, self.K, self.S = total_energy, cfl, K, S
        self.cons_mom = cons_mom
        self.dt_fixed = dt_fixed
        self.t = 0.0
        self.n = 0
        self.diverged = False
        self.max_wall_s = 60.0           # per advance_to call; flags runaway runs
        self.p_amb_fn = None            # optional back-pressure schedule P_amb(t)
        self._gam = np.full(N, gamma)
        self._R = np.full(N, R)
        self._Cp = np.full(N, self.Cp)
        self._Tc = np.full(N, 1.0e5)
        self._D = np.full(N, 0.1)
        self._f = np.zeros(N)
        self._mom = np.zeros(N + 1)
        self.last = None

    def dt(self):
        if self.dt_fixed is not None:
            return self.dt_fixed
        a = math.sqrt(self.gamma * self.R * float(self.T.max()))
        return compute_dt_cfl(self.u, a, self.dx, self.N + 1, self.cfl, 1.0)

    def step(self, dt=None):
        dt = self.dt() if dt is None else dt
        if self.p_amb_fn is not None:
            self.P_amb = self.p_amb_fn(self.t + dt)
        out = piso(self.rho, self.u, self.P, self.T, self.A, self._D,
                   self.msrc, self.tsrc, self._mom, self._f, self.dx, dt,
                   self._gam, self._R, self._Cp, self._Tc, self.A_throat,
                   self.P_amb, self.T_amb, self.N, 0.0,
                   pressure_floor_pa=1.0, energy_coupling_iters=self.K,
                   energy_coupling_tol=1.0e-6, corr2_storage=self.S,
                   total_energy=self.total,
                   conservative_momentum=self.cons_mom)
        self.rho, self.u, self.P, self.T = out[0], out[1], out[2], out[3]
        self.t += dt
        self.n += 1
        self.last = out
        return out

    def advance_to(self, t_end):
        t_wall = time.time()
        while self.t < t_end - 1e-12:
            if time.time() - t_wall > self.max_wall_s:   # runaway guard
                self.diverged = True
                return
            dt_cfl = self.dt()
            if not (dt_cfl > 1.0e-10):   # NaN or collapsed time step: the run diverged
                self.diverged = True
                return
            dt = min(dt_cfl, t_end - self.t)
            self.step(dt)
            if not np.isfinite(self.T).all() or float(self.T.max()) > 2.0e4:
                self.diverged = True
                return

    def u_cell(self):
        return 0.5 * (self.u[:-1] + self.u[1:])

    def mach(self):
        return self.u_cell() / np.sqrt(self.gamma * self.R * self.T)


# ----------------------------------------------------------------------------
# 1. Sod -- exact Riemann solver (Toro, "Riemann Solvers and Numerical Methods for
#    Fluid Dynamics", ch. 4; ideal gas)
# ----------------------------------------------------------------------------
def _f_and_df(p, rho_k, p_k, gamma):
    a_k = math.sqrt(gamma * p_k / rho_k)
    if p > p_k:  # shock
        A = 2.0 / ((gamma + 1.0) * rho_k)
        B = (gamma - 1.0) / (gamma + 1.0) * p_k
        f = (p - p_k) * math.sqrt(A / (p + B))
        df = math.sqrt(A / (p + B)) * (1.0 - 0.5 * (p - p_k) / (p + B))
    else:  # rarefaction
        f = 2.0 * a_k / (gamma - 1.0) * ((p / p_k) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)
        df = (1.0 / (rho_k * a_k)) * (p / p_k) ** (-(gamma + 1.0) / (2.0 * gamma))
    return f, df


def riemann_exact(left, right, gamma, xi):
    """Exact solution (rho, u, p) at similarity coordinate xi = (x - x0)/t.

    left/right = (rho, u, p).  Returns arrays plus a dict of the star-state
    quantities and wave speeds."""
    rL, uL, pL = left
    rR, uR, pR = right
    aL, aR = math.sqrt(gamma * pL / rL), math.sqrt(gamma * pR / rR)
    p = 0.5 * (pL + pR)
    for _ in range(100):
        fL, dfL = _f_and_df(p, rL, pL, gamma)
        fR, dfR = _f_and_df(p, rR, pR, gamma)
        g = fL + fR + (uR - uL)
        pn = max(p - g / (dfL + dfR), 1e-8 * p)
        if abs(pn - p) < 1e-14 * p:
            p = pn
            break
        p = pn
    fL, _ = _f_and_df(p, rL, pL, gamma)
    fR, _ = _f_and_df(p, rR, pR, gamma)
    ps = p
    us = 0.5 * (uL + uR) + 0.5 * (fR - fL)
    gm = gamma
    # left is a rarefaction, right a shock for Sod-like data; handle generally
    info = {'p_star': ps, 'u_star': us}
    xi = np.asarray(xi, float)
    rho = np.empty_like(xi)
    u = np.empty_like(xi)
    pr = np.empty_like(xi)

    def side(sign, rk, uk, pk, ak):
        # sign=-1: left wave, +1: right wave. returns (star density, wave descriptor)
        if ps > pk:  # shock
            S = uk + sign * ak * math.sqrt((gm + 1.0) / (2.0 * gm) * ps / pk
                                          + (gm - 1.0) / (2.0 * gm))
            rs = rk * ((ps / pk + (gm - 1.0) / (gm + 1.0))
                       / ((gm - 1.0) / (gm + 1.0) * ps / pk + 1.0))
            return rs, ('shock', S)
        rs = rk * (ps / pk) ** (1.0 / gm)
        a_s = ak * (ps / pk) ** ((gm - 1.0) / (2.0 * gm))
        return rs, ('fan', uk + sign * ak, us + sign * a_s)

    rsL, wL = side(-1, rL, uL, pL, aL)
    rsR, wR = side(+1, rR, uR, pR, aR)
    info.update(rho_star_L=rsL, rho_star_R=rsR, wave_L=wL, wave_R=wR)
    for k, x in enumerate(xi):
        if x <= us:  # left of contact
            if wL[0] == 'shock':
                if x <= wL[1]:
                    rho[k], u[k], pr[k] = rL, uL, pL
                else:
                    rho[k], u[k], pr[k] = rsL, us, ps
            else:
                head, tail = wL[1], wL[2]
                if x <= head:
                    rho[k], u[k], pr[k] = rL, uL, pL
                elif x >= tail:
                    rho[k], u[k], pr[k] = rsL, us, ps
                else:
                    uf = 2.0 / (gm + 1.0) * (aL + 0.5 * (gm - 1.0) * uL + x)
                    af = 2.0 / (gm + 1.0) * (aL + 0.5 * (gm - 1.0) * (uL - x))
                    rho[k] = rL * (af / aL) ** (2.0 / (gm - 1.0))
                    u[k], pr[k] = uf, pL * (af / aL) ** (2.0 * gm / (gm - 1.0))
        else:
            if wR[0] == 'shock':
                if x >= wR[1]:
                    rho[k], u[k], pr[k] = rR, uR, pR
                else:
                    rho[k], u[k], pr[k] = rsR, us, ps
            else:
                head, tail = wR[1], wR[2]
                if x >= head:
                    rho[k], u[k], pr[k] = rR, uR, pR
                elif x <= tail:
                    rho[k], u[k], pr[k] = rsR, us, ps
                else:
                    uf = 2.0 / (gm + 1.0) * (-aR + 0.5 * (gm - 1.0) * uR + x)
                    af = 2.0 / (gm + 1.0) * (aR - 0.5 * (gm - 1.0) * (uR - x))
                    rho[k] = rR * (af / aR) ** (2.0 / (gm - 1.0))
                    u[k], pr[k] = uf, pR * (af / aR) ** (2.0 * gm / (gm - 1.0))
    return rho, u, pr, info


SOD_L = (1.0, 0.0, 1.0e5)       # rho [kg/m3], u, P [Pa]  (Sod scaled so that sqrt(P/rho) = 316 m/s)
SOD_R = (0.125, 0.0, 1.0e4)


def sod_case(N, total, cfl=0.3, K=1, S=True, t_nondim=0.2, L=1.0, x0=0.5, snapshots=(0.6,),
             cons_mom=False):
    """Run Sod to t = t_nondim * L/sqrt(P_L/rho_L); returns dict of results."""
    dx = L / N
    x = (np.arange(N) + 0.5) * dx
    rho = np.where(x < x0, SOD_L[0], SOD_R[0])
    P = np.where(x < x0, SOD_L[2], SOD_R[2])
    T = P / (rho * R_GAS)
    u = np.zeros(N + 1)
    tube = Tube(np.full(N, 1.0e-3), dx, rho, u, P, T, total_energy=total, cfl=cfl,
                K=K, S=S, A_throat=0.0, T_ambient=1.0, cons_mom=cons_mom)
    t_end = t_nondim * L / math.sqrt(SOD_L[2] / SOD_L[0])
    snaps = {}
    for fr in snapshots:
        tube.advance_to(fr * t_end)
        snaps[fr] = (tube.t, tube.P.copy())
    tube.advance_to(t_end)
    uc = tube.u_cell()
    xi = (x - x0) / tube.t
    rho_e, u_e, p_e, info = riemann_exact(SOD_L, SOD_R, GAMMA, xi)
    T_e = p_e / (rho_e * R_GAS)
    out = dict(N=N, total=total, steps=tube.n, t=tube.t, x=x, rho=tube.rho,
               u=uc, P=tube.P, T=tube.T, rho_e=rho_e, u_e=u_e, P_e=p_e, T_e=T_e,
               info=info, snaps=snaps, x0=x0, tube=tube)
    # energy bookkeeping
    return out


def _cross_from_right(x, f, level):
    """x where f (decreasing towards the right) crosses ``level`` scanning from the right."""
    for k in range(len(x) - 1, 0, -1):
        if f[k - 1] >= level > f[k]:
            return x[k - 1] + (f[k - 1] - level) / (f[k - 1] - f[k]) * (x[k] - x[k - 1])
    return float('nan')


def sod_metrics(res):
    x, info, t = res['x'], res['info'], res['t']
    dx = x[1] - x[0]
    wR = info['wave_R']
    S_shock = wR[1]
    xs_exact = res['x0'] + S_shock * t
    P_e = res['P_e']
    level = 0.5 * (info['p_star'] + SOD_R[2])
    xs_num = _cross_from_right(x, res['P'], level)
    # shock speed from the 60 % snapshot
    t1, P1 = res['snaps'][0.6]
    xs1 = _cross_from_right(x, P1, level)
    S_num = (xs_num - xs1) / (t - t1)
    # plateaus: mid-points of (fan tail, contact) and (contact, shock)
    xc = res['x0'] + info['u_star'] * t
    xtail = res['x0'] + info['wave_L'][2] * t
    plats = {}
    for name, a, b in (('star-L', xtail, xc), ('star-R', xc, xs_exact)):
        mid = 0.5 * (a + b)
        w = 0.15 * (b - a)
        m = (x > mid - w) & (x < mid + w)
        plats[name] = {k: (float(res[k][m].mean()), float(res[k + '_e'][m].mean()))
                       for k in ('rho', 'u', 'P', 'T')}
    # L1 errors, normalised by the exact-solution scale
    scales = {'rho': SOD_L[0], 'u': max(abs(info['u_star']), 1.0), 'P': SOD_L[2],
              'T': res['T_e'].max()}
    l1 = {k: float(np.mean(np.abs(res[k] - res[k + '_e'])) / scales[k])
          for k in ('rho', 'u', 'P', 'T')}
    # total-energy conservation of the discrete state (closed tube): sum m(e + KE)
    return dict(xs_exact=xs_exact, xs_num=xs_num, S_exact=S_shock, S_num=S_num,
                plateaus=plats, l1=l1, dx=dx)


def sod_total_energy_drift(res):
    """Relative drift of the discrete total energy  sum m (Cv T + u_cell^2/2)  of the
    closed tube vs its initial value (should be 0 for a conservative scheme)."""
    N = res['N']
    dx = 1.0 / N
    A = 1.0e-3
    rho0 = np.where((np.arange(N) + 0.5) * dx < res['x0'], SOD_L[0], SOD_R[0])
    P0 = np.where((np.arange(N) + 0.5) * dx < res['x0'], SOD_L[2], SOD_R[2])
    E0 = np.sum(rho0 * A * dx * (P0 / (rho0 * (GAMMA - 1.0))))
    u_f = res['tube'].u
    ke = 0.25 * (u_f[:-1] ** 2 + u_f[1:] ** 2)
    E1 = np.sum(res['tube'].rho * A * dx * (CV * res['T'] + ke))
    M0 = np.sum(rho0) * A * dx
    M1 = np.sum(res['tube'].rho) * A * dx
    return dict(dE_rel=(E1 - E0) / E0, dM_rel=(M1 - M0) / M0)


# ----------------------------------------------------------------------------
# 2./3. quasi-1-D converging-diverging duct
# ----------------------------------------------------------------------------
def area_ratio_of_M(M, gamma=GAMMA):
    return (1.0 / M) * ((2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * M * M)) ** (
        (gamma + 1.0) / (2.0 * (gamma - 1.0)))


def mach_of_area_ratio(ar, supersonic, gamma=GAMMA):
    lo, hi = (1.0, 20.0) if supersonic else (1.0e-6, 1.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = area_ratio_of_M(mid, gamma) - ar
        if (f > 0) == supersonic:   # supersonic branch: A/A* increases with M
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def total_to_static_T(M, gamma=GAMMA):
    return 1.0 + 0.5 * (gamma - 1.0) * M * M


def total_to_static_P(M, gamma=GAMMA):
    return total_to_static_T(M, gamma) ** (gamma / (gamma - 1.0))


def flow_function(gamma=GAMMA):
    return math.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))


class Duct:
    """Converging-diverging profile A(x)/A* with a smooth throat at x_t."""

    def __init__(self, L=1.0, x_t=0.4, A_star=1.0e-3, ar_in=8.0, ar_exit=2.0):
        self.L, self.x_t, self.As = L, x_t, A_star
        self.ar_in, self.ar_exit = ar_in, ar_exit

    def ar(self, x):
        x = np.asarray(x, float)
        conv = 1.0 + (self.ar_in - 1.0) * (1.0 - x / self.x_t) ** 2
        div = 1.0 + (self.ar_exit - 1.0) * ((x - self.x_t) / (self.L - self.x_t)) ** 2
        return np.where(x <= self.x_t, conv, div)

    def A(self, x):
        return self.As * self.ar(x)

    def x_of_ar_div(self, ar):
        return self.x_t + (self.L - self.x_t) * math.sqrt((ar - 1.0) / (self.ar_exit - 1.0))


def isentropic_exact(duct, x, P0, T0, A_star_eff=None, supersonic_from=None,
                     gamma=GAMMA):
    """Isentropic quasi-1-D state at stations x.  Subsonic for x < supersonic_from
    (default: throat), supersonic after.  Returns M, T/T0, P/P0 arrays."""
    ar = duct.ar(x) * duct.As / (duct.As if A_star_eff is None else A_star_eff)
    sup_from = duct.x_t if supersonic_from is None else supersonic_from
    M = np.array([mach_of_area_ratio(a, xi > sup_from, gamma) for a, xi in zip(ar, x)])
    return M, 1.0 / total_to_static_T(M, gamma), 1.0 / total_to_static_P(M, gamma)


def normal_shock_p0_ratio(M1, gamma=GAMMA):
    g = gamma
    a = ((g + 1.0) * M1 ** 2 / ((g - 1.0) * M1 ** 2 + 2.0)) ** (g / (g - 1.0))
    b = ((g + 1.0) / (2.0 * g * M1 ** 2 - (g - 1.0))) ** (1.0 / (g - 1.0))
    return a * b


def post_shock_M(M1, gamma=GAMMA):
    g = gamma
    return math.sqrt(((g - 1.0) * M1 ** 2 + 2.0) / (2.0 * g * M1 ** 2 - (g - 1.0)))


def duct_exact(duct, x, P0, T0, shock_ar=None, gamma=GAMMA, throat_M=None):
    """Exact steady state (supersonic exit, or a normal shock where A/A* = shock_ar in
    the divergent section).  Returns M, T, P arrays at x and a dict with the
    shock station, the post-shock stagnation ratio and the exit back pressure of the
    LAST station."""
    x = np.asarray(x, float)
    ar = duct.ar(x)
    info = {}
    if throat_M is not None:
        # fully subsonic (venturi): the duct's minimum area is A_t = A*_eff * ar(M_t)
        ar_eff = ar * area_ratio_of_M(throat_M, gamma)
        M = np.array([mach_of_area_ratio(a, False, gamma) for a in ar_eff])
        info['A_star_eff_over_A_t'] = 1.0 / area_ratio_of_M(throat_M, gamma)
        T = T0 / total_to_static_T(M, gamma)
        return M, T, P0 / total_to_static_P(M, gamma), info
    M = np.array([mach_of_area_ratio(a, xi > duct.x_t, gamma) for a, xi in zip(ar, x)])
    if shock_ar is not None:
        x_s = duct.x_of_ar_div(shock_ar)
        M1 = mach_of_area_ratio(shock_ar, True, gamma)
        p02_p01 = normal_shock_p0_ratio(M1, gamma)
        ar2 = ar * p02_p01       # A / A2*  with A2* = A* / (p02/p01)
        # downstream stations: subsonic on the A2* branch
        for k, xi in enumerate(x):
            if xi > x_s:
                M[k] = mach_of_area_ratio(ar2[k], False, gamma)
        info.update(x_shock=x_s, M1=M1, M2=post_shock_M(M1, gamma), p02_p01=p02_p01)
        P_over_P0 = np.where(x > x_s, p02_p01, 1.0) / total_to_static_P(M, gamma)
    else:
        P_over_P0 = 1.0 / total_to_static_P(M, gamma)
    T = T0 / total_to_static_T(M, gamma)
    P = P0 * P_over_P0
    return M, T, P, info


def _back_pressure(P_c, T_c, A_bc, mdot, T_amb):
    """Back pressure ``Pa`` for which the motor's nozzle-boundary function passes exactly
    ``mdot`` from a last-cell state (P_c, T_c) through the area ``A_bc`` (bisection on
    the subsonic-orifice branch).  This is how a prescribed exit pressure is realised
    with the unchanged boundary code (see the module docstring)."""
    lo, hi = 0.5 * P_c, P_c
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        m = _nozzle_boundary_flow(P_c, T_c, A_bc, GAMMA, R_GAS, mid, T_amb)[0]
        if m > mdot:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def duct_case(N, total, shock_ar=None, throat_M=None, cfl=0.3, K=1, S=True, P0=5.0e5, T0=1000.0,
              t_end=0.05, duct=None, big=3.0, check_every=0.005, verbose=False,
              cons_mom=False, init_shift_ar=-0.25, ramp_s=3.0e-3):
    """Run the duct to (near) steady state and return numerical + exact profiles.

    Initial condition: the exact isentropic fields (for the shock case: the exact shocked
    fields with the shock initially at the area ratio ``shock_ar + init_shift_ar``, so it
    has to travel to its exact station).  The sensible-enthalpy form then drifts away
    from this state to its own steady state; a form that is consistent with the
    conservation laws must stay (or return) there.
    """
    duct = duct or Duct()
    dx = duct.L / N
    x = (np.arange(N) + 0.5) * dx
    A = duct.A(x)
    M_ex, T_ex, P_ex, info = duct_exact(duct, x, P0, T0, shock_ar, throat_M=throat_M)
    a_eff = duct.As * info.get('A_star_eff_over_A_t', 1.0)
    mdot = a_eff * P0 * flow_function() / math.sqrt(R_GAS * T0)
    msrc = np.zeros(N)
    tsrc = np.zeros(N)
    msrc[0] = mdot / dx
    tsrc[0] = msrc[0] * CP * T0
    A_bc = big * A[-1]
    T_amb = 1.0      # also the temperature floor of the kernel (max(1, T_ambient)); inflow never occurs at steady state
    P_back = _back_pressure(P_ex[-1], T_ex[-1], A_bc, mdot, T_amb)
    # initial fields: the exact solution.  Isentropic case: the target solution.  Shock
    # case: the exact solution for a shock initially at the area ratio
    # ``shock_ar + init_shift_ar`` (weaker, further upstream, so the back pressure only has to FALL: no reverse flow), held at the matching
    # back pressure; the back pressure is then ramped (``ramp_s``) to the target value
    # so the shock has to travel to its exact station.
    if shock_ar is None:        # isentropic cases start from the exact solution
        M_i, T_i, P_i, _ = M_ex, T_ex, P_ex, None
        P_back0 = P_back
    else:
        M_i, T_i, P_i, _ = duct_exact(duct, x, P0, T0, shock_ar + init_shift_ar)
        P_back0 = _back_pressure(P_i[-1], T_i[-1], A_bc, mdot, T_amb)
    T = T_i.copy()
    P = P_i.copy()
    rho = P / (R_GAS * T)
    uc = M_i * np.sqrt(GAMMA * R_GAS * T_i)           # exact cell velocity
    u = np.zeros(N + 1)
    u[1:N] = 0.5 * (uc[:-1] + uc[1:])
    u[N] = uc[-1]
    tube = Tube(A, dx, rho, u, P, T, mass_source=msrc, thermal_source=tsrc,
                A_throat=A_bc, P_ambient=P_back0, T_ambient=T_amb,
                total_energy=total, cfl=cfl, K=K, S=S, cons_mom=cons_mom)
    last_P = tube.P.copy()
    resid = []
    if shock_ar is not None:
        tube.p_amb_fn = lambda t: P_back0 + (P_back - P_back0) * min(1.0, t / ramp_s)
    while tube.t < t_end:
        tube.advance_to(min(tube.t + check_every, t_end))
        r = float(np.max(np.abs(tube.P - last_P) / tube.P))
        resid.append((tube.t, r))
        last_P = tube.P.copy()
        if verbose:
            print(f"   t={tube.t:.4f} resid={r:.3e} Tmax={tube.T.max():.0f}", flush=True)
        if tube.diverged:
            resid.append((tube.t, float('nan')))
            break
    return dict(N=N, total=total, x=x, A=A, duct=duct, tube=tube, P0=P0, T0=T0,
                M=tube.mach(), T=tube.T, P=tube.P, M_e=M_ex, T_e=T_ex, P_e=P_ex,
                info=info, mdot=mdot, resid=resid, steps=tube.n, P_back=P_back,
                diverged=tube.diverged)


def duct_metrics(res, exclude_shock_cells=6):
    x = res['x']
    N = res['N']
    info = res['info']
    out = {}
    mask = np.ones(N, bool)
    if 'x_shock' in info:
        dx = x[1] - x[0]
        mask = np.abs(x - info['x_shock']) > exclude_shock_cells * dx
    for name in ('M', 'T', 'P'):
        num, ex = res[name], res[name + '_e']
        scale = {'M': np.max(ex), 'T': res['T0'], 'P': res['P0']}[name]
        out['L1_' + name] = float(np.mean(np.abs(num - ex)) / scale)
        out['L1ex_' + name] = float(np.mean(np.abs(num - ex)[mask]) / scale)
        out['Linf_' + name] = float(np.max(np.abs(num - ex)[mask]) / scale)
        if 'x_shock' in info:
            up = x < info['x_shock'] - 8.0 * (x[1] - x[0])      # isentropic, supersonic part
            out['L1up_' + name] = float(np.mean(np.abs(num - ex)[up]) / scale)
    # throat values
    kt = int(np.argmin(np.abs(x - res['duct'].x_t)))
    out['M_throat'] = (float(res['M'][kt]), float(res['M_e'][kt]))
    out['T_over_T0_throat'] = (float(res['T'][kt] / res['T0']), float(res['T_e'][kt] / res['T0']))
    out['M_exit'] = (float(res['M'][-1]), float(res['M_e'][-1]))
    out['T_over_T0_exit'] = (float(res['T'][-1] / res['T0']), float(res['T_e'][-1] / res['T0']))
    if 'x_shock' in info:
        P_pre = res['P_e'][np.where(x < info['x_shock'])[0][-1]]
        P_post = res['P_e'][np.where(x > info['x_shock'])[0][0]]
        level = 0.5 * (P_pre + P_post)
        # scan for the downward crossing of the pressure jump -> UPWARD jump
        xs = float('nan')
        P = res['P']
        for k in range(1, N):
            if P[k - 1] < level <= P[k] and abs(x[k] - info['x_shock']) < 0.2:
                xs = x[k - 1] + (level - P[k - 1]) / (P[k] - P[k - 1]) * (x[k] - x[k - 1])
                break
        out['x_shock'] = (xs, info['x_shock'])
        # post-shock plateau: 8 cells past the smeared shock vs the exact value there
        kk = int(np.searchsorted(x, info['x_shock'])) + 8
        kk = min(kk, N - 1)
        out['post_shock_P_ratio'] = (float(res['P'][kk] / res['P0']), float(res['P_e'][kk] / res['P0']))
        out['post_shock_T_ratio'] = (float(res['T'][kk] / res['T0']), float(res['T_e'][kk] / res['T0']))
        out['post_shock_M'] = (float(res['M'][kk]), float(res['M_e'][kk]))
        jj = int(np.searchsorted(x, info['x_shock'])) - 8
        out['pre_shock_M'] = (float(res['M'][jj]), float(res['M_e'][jj]))
        # total-pressure loss across the shock, from the numerical state
        def p0(k):
            return res['P'][k] * total_to_static_P(res['M'][k])
        out['p0_ratio_across_shock'] = (float(p0(kk) / p0(jj)), float(info['p02_p01']))
        # stagnation temperature from the numerical state: T (1 + (g-1)/2 M^2)/T0
        def t0(k):
            return res['T'][k] * total_to_static_T(res['M'][k]) / res['T0']
        out['T0_pre_post'] = (float(t0(jj)), float(t0(kk)))
    return out


# ----------------------------------------------------------------------------
# 4. closed-vessel filling
# ----------------------------------------------------------------------------
def vessel_exact(t, m0, T0, mdot, T_in, V):
    m = m0 + mdot * t
    E = m0 * CV * T0 + mdot * t * CP * T_in
    T = E / (m * CV)
    return m, T, m * R_GAS * T / V


def vessel_case(total, dt, t_end=0.01, N=20, P0=1.0e5, T0=300.0, T_in=1200.0,
                mass_ratio=2.0, K=1, S=True):
    L, A = 1.0, 1.0e-3
    dx = L / N
    V = A * L
    rho0 = P0 / (R_GAS * T0)
    m0 = rho0 * V
    mdot = mass_ratio * m0 / t_end
    msrc = np.full(N, mdot / L)
    tsrc = msrc * CP * T_in
    tube = Tube(np.full(N, A), dx, np.full(N, rho0), np.zeros(N + 1), np.full(N, P0),
                np.full(N, T0), mass_source=msrc, thermal_source=tsrc, A_throat=0.0,
                total_energy=total, K=K, S=S, dt_fixed=dt)
    n = int(round(t_end / dt))
    for _ in range(n):
        tube.step()
    m_e, T_e, P_e = vessel_exact(tube.t, m0, T0, mdot, T_in, V)
    m_n = float(np.sum(tube.rho) * A * dx)
    return dict(total=total, dt=dt, t=tube.t, T=(float(tube.T.mean()), T_e),
                P=(float(tube.P.mean()), P_e), m=(m_n, m_e),
                umax=float(np.abs(tube.u).max()),
                P_spread=float(tube.P.max() - tube.P.min()),
                T_sensible_limit=(m0 * T0 + mdot * t_end * T_in) / (m0 + mdot * t_end),
                P_sensible_limit=(m0 + mdot * t_end) * R_GAS
                * ((m0 * T0 + mdot * t_end * T_in) / (m0 + mdot * t_end)) / V)


# ----------------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------------
def order(e_coarse, e_fine, ratio=2.0):
    if e_fine <= 0 or e_coarse <= 0:
        return float('nan')
    return math.log(e_coarse / e_fine) / math.log(ratio)


def run_all(args):
    results = {}
    ncells = [100, 200, 400, 800, 1600] if not args.quick else [100, 200, 400]
    # ---- Sod
    sod = {}
    for total in (False, True):
        tag = 'total' if total else 'sensible'
        sod[tag] = {}
        for N in ncells:
            t0 = time.time()
            r = sod_case(N, total)
            m = sod_metrics(r)
            m['drift'] = sod_total_energy_drift(r)
            m['steps'] = r['steps']
            m['wall'] = time.time() - t0
            sod[tag][N] = m
            print(f"[sod {tag} N={N}] xs {m['xs_num']:.5f} (exact {m['xs_exact']:.5f}) "
                  f"S {m['S_num']:.1f} (exact {m['S_exact']:.1f}) L1 {m['l1']} "
                  f"drift {m['drift']} wall {m['wall']:.1f}s", flush=True)
    results['sod'] = sod
    # ---- vessel
    ves = {}
    for total in (False, True):
        tag = 'total' if total else 'sensible'
        ves[tag] = [vessel_case(total, dt) for dt in (4e-5, 2e-5, 1e-5, 5e-6, 2.5e-6)]
        for v in ves[tag]:
            print(f"[vessel {tag} dt={v['dt']:.1e}] T {v['T']} P {v['P']}", flush=True)
    results['vessel'] = ves
    # ---- duct
    nd = [100, 200, 400] if args.quick else [100, 200, 400, 800]
    for name, shock_ar in (('nozzle', None), ('shock', 1.5)):
        results[name] = {}
        for total in (False, True):
            tag = 'total' if total else 'sensible'
            results[name][tag] = {}
            for N in nd:
                t0 = time.time()
                try:
                    r = duct_case(N, total, shock_ar=shock_ar)
                    m = duct_metrics(r)
                    m['resid_last'] = r['resid'][-1][1]
                    m['steps'] = r['steps']
                except Exception as exc:   # report, do not hide
                    m = dict(error=repr(exc))
                m['wall'] = time.time() - t0
                results[name][tag][N] = m
                print(f"[{name} {tag} N={N}] {m}", flush=True)
    return results


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    print(f"[provenance] srm_1d from {Path(sys.modules['srm_1d'].__file__).parent}")
    res = run_all(a)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(_jsonable(res), f, indent=1)
