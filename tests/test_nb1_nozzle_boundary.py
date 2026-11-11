"""NB1: the nozzle boundary uses the true (isentropic) stagnation state.

Exact-answer tests of ``solver._nozzle_boundary_flow`` and of the PISO step with the
boundary in both modes (``nozzle_isentropic`` True = isentropic, False = the pre-NB1
static-as-stagnation boundary).  Spec: docs/solver_verification/NB1_TASK.md;
worklog: docs/solver_verification/worklog/NB1.md.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from srm_1d.solver import (  # noqa: E402
    _nozzle_boundary_flow,
    _nozzle_choke_mach,
    _subsonic_mach_from_area_ratio,
    _critical_flow_function,
    _critical_pressure_ratio,
    NOZZLE_STATE_CHOKED_OUT,
    NOZZLE_STATE_SUBSONIC_OUT,
    NOZZLE_STATE_SUBSONIC_IN,
    NOZZLE_STATE_CHOKED_IN,
    compute_dt_cfl,
    piso_step,
)

GAMMAS = (1.1, 1.2, 1.4)
R = 300.0
T = 2800.0
PA = 101325.0
TA = 293.0


def _g(M, gamma):
    return (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))


def _area_ratio(M, gamma):
    x = 0.5 * (gamma - 1.0)
    e = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return ((1.0 + x * M * M) / (1.0 + x)) ** e / M


def _m_sub_bisect(r, gamma):
    """Independent subsonic area-Mach root (plain bisection on A/A*, decreasing in M)."""
    lo, hi = 1.0e-9, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _area_ratio(mid, gamma) > r:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _flow(P, A_t, A_e, gamma, isentropic=True, Pa=PA):
    return _nozzle_boundary_flow(P, T, A_t, gamma, R, Pa, TA, A_e,
                                 _nozzle_choke_mach(A_e, A_t, gamma, isentropic))


class TestAreaMachRoot:
    @pytest.mark.parametrize("gamma", GAMMAS)
    @pytest.mark.parametrize("r", [1.0000001, 1.01, 1.2, 2.0, 5.0, 40.0, 1000.0])
    def test_matches_bisection(self, gamma, r):
        M = _subsonic_mach_from_area_ratio(r, gamma)
        assert 0.0 < M < 1.0
        assert _area_ratio(M, gamma) == pytest.approx(r, rel=1e-12)
        assert M == pytest.approx(_m_sub_bisect(r, gamma), rel=1e-6)

    def test_no_subsonic_root_below_one(self):
        assert _subsonic_mach_from_area_ratio(1.0, 1.2) == 1.0
        assert _subsonic_mach_from_area_ratio(0.5, 1.2) == 1.0


class TestChoked:
    @pytest.mark.parametrize("gamma", GAMMAS)
    @pytest.mark.parametrize("r", [1.05, 1.5, 3.0, 12.0])
    def test_p0_over_p_and_textbook_mdot(self, gamma, r):
        A_t = 1.0e-4
        A_e = r * A_t
        P = 4.0e6
        mdot, dmdp, T_up, state = _flow(P, A_t, A_e, gamma)
        assert state == NOZZLE_STATE_CHOKED_OUT
        M_e = _m_sub_bisect(r, gamma)
        P0 = P * _g(M_e, gamma)
        T0 = T * (1.0 + 0.5 * (gamma - 1.0) * M_e * M_e)
        gamma_fn = _critical_flow_function(gamma)
        # textbook choked flow from the true stagnation state
        assert mdot == pytest.approx(gamma_fn * A_t * P0 / math.sqrt(R * T0), rel=1e-10)
        # and equal to the continuity flow at the port exit, rho u A_e
        assert mdot == pytest.approx(A_e * P * M_e * math.sqrt(gamma / (R * T)), rel=1e-9)
        # the implied P0/P (mdot / the static-as-stagnation flow, with T0) is g(M_e)
        P0_implied = mdot * math.sqrt(R * T0) / (gamma_fn * A_t)
        assert P0_implied / P == pytest.approx(_g(M_e, gamma), rel=1e-10)
        assert dmdp == pytest.approx(mdot / P, rel=1e-12)
        # the advected temperature is the STATIC one (no double-counted kinetic energy)
        assert T_up == pytest.approx(T)

    def test_static_mode_unchanged(self):
        gamma = 1.2
        A_t = 1.0e-4
        mdot, dmdp, T_up, state = _flow(2.0e6, A_t, 5.0 * A_t, gamma, isentropic=False)
        gf = _critical_flow_function(gamma)
        assert mdot == pytest.approx(2.0e6 * A_t * gf / math.sqrt(R * T), rel=1e-14)
        assert dmdp == pytest.approx(A_t * gf / math.sqrt(R * T), rel=1e-14)
        assert state == NOZZLE_STATE_CHOKED_OUT

    def test_large_exit_area_recovers_static(self):
        gamma = 1.2
        A_t = 1.0e-4
        a = _flow(2.0e6, A_t, 1.0e6 * A_t, gamma)[0]
        b = _flow(2.0e6, A_t, 1.0e6 * A_t, gamma, isentropic=False)[0]
        assert a == pytest.approx(b, rel=1e-6)

    def test_chunc_like_mach_gives_about_seven_percent(self):
        # M_e ~ 0.35, gamma 1.2: P0/P - 1 ~ 7 %, and the choked flow is higher by
        # g(M_e)/sqrt(1 + (gamma-1)/2 M_e^2) = (1 + x M^2)^(k - 1/2).
        gamma = 1.2
        r = _area_ratio(0.35, gamma)
        A_t = 1.0e-4
        new = _flow(6.0e6, A_t, r * A_t, gamma)[0]
        old = _flow(6.0e6, A_t, r * A_t, gamma, isentropic=False)[0]
        assert _g(0.35, gamma) == pytest.approx(1.0757, abs=2e-3)
        assert new / old == pytest.approx((1.0 + 0.1 * 0.35 ** 2) ** (6.0 - 0.5), rel=1e-9)


class TestSubsonicAndSwitch:
    @pytest.mark.parametrize("gamma", GAMMAS)
    @pytest.mark.parametrize("r", [1.3, 2.0, 6.0])
    def test_dmdp_matches_finite_difference_both_branches(self, gamma, r):
        A_t = 1.0e-4
        A_e = r * A_t
        crit = _critical_pressure_ratio(gamma)
        M_c = _m_sub_bisect(r, gamma)
        P_switch = PA / (crit * _g(M_c, gamma))
        seen = set()
        for P in (PA * 1.0005, PA * 1.01, 0.5 * (PA + P_switch), 0.999 * P_switch,
                  1.001 * P_switch, 1.5 * P_switch, 40.0 * P_switch):
            mdot, dmdp, _, state = _flow(P, A_t, A_e, gamma)
            seen.add(state)
            h = 1.0e-6 * P
            fd = (_flow(P + h, A_t, A_e, gamma)[0] - _flow(P - h, A_t, A_e, gamma)[0]) / (2 * h)
            assert dmdp == pytest.approx(fd, rel=2e-5)
            assert dmdp >= 0.0 and mdot > 0.0
        assert seen == {NOZZLE_STATE_SUBSONIC_OUT, NOZZLE_STATE_CHOKED_OUT}

    @pytest.mark.parametrize("gamma", GAMMAS)
    @pytest.mark.parametrize("r", [1.3, 2.0, 6.0])
    def test_subsonic_consistency_and_switch_continuity(self, gamma, r):
        A_t = 1.0e-4
        A_e = r * A_t
        crit = _critical_pressure_ratio(gamma)
        M_c = _m_sub_bisect(r, gamma)
        P_switch = PA / (crit * _g(M_c, gamma))
        # self-consistency of a subsonic point: continuity at the port exit and at the
        # throat (static P_amb), with one isentropic stagnation state.
        P = 0.5 * (PA + P_switch)
        mdot, _, _, state = _flow(P, A_t, A_e, gamma)
        assert state == NOZZLE_STATE_SUBSONIC_OUT
        M_e = mdot / (A_e * P * math.sqrt(gamma / (R * T)))
        P0 = P * _g(M_e, gamma)
        T0 = T * (1.0 + 0.5 * (gamma - 1.0) * M_e ** 2)
        pr = PA / P0
        assert pr > crit
        M_t = math.sqrt(2.0 / (gamma - 1.0) * (pr ** (-(gamma - 1.0) / gamma) - 1.0))
        T_t = T0 / (1.0 + 0.5 * (gamma - 1.0) * M_t ** 2)
        mdot_t = A_t * PA * M_t * math.sqrt(gamma / (R * T_t))
        assert mdot_t == pytest.approx(mdot, rel=1e-9)
        # continuity of mdot and of dmdp across the switch
        lo = _flow(P_switch * (1.0 - 1e-9), A_t, A_e, gamma)
        hi = _flow(P_switch * (1.0 + 1e-9), A_t, A_e, gamma)
        assert lo[3] == NOZZLE_STATE_SUBSONIC_OUT and hi[3] == NOZZLE_STATE_CHOKED_OUT
        assert lo[0] == pytest.approx(hi[0], rel=1e-7)
        assert lo[1] == pytest.approx(hi[1], rel=1e-4)

    def test_switch_uses_stagnation_not_static_ratio(self):
        # P_amb / P is below the critical ratio but P_amb / P0 is not -> still subsonic.
        gamma, r = 1.2, 1.5
        A_t = 1.0e-4
        crit = _critical_pressure_ratio(gamma)
        M_c = _m_sub_bisect(r, gamma)
        P_switch = PA / (crit * _g(M_c, gamma))
        P_static_switch = PA / crit            # where the OLD boundary chokes
        assert P_switch < P_static_switch      # isentropic chokes earlier (lower static P)
        P = 0.5 * (P_switch + P_static_switch)
        assert _flow(P, A_t, r * A_t, gamma)[3] == NOZZLE_STATE_CHOKED_OUT
        assert _flow(P, A_t, r * A_t, gamma, isentropic=False)[3] == NOZZLE_STATE_SUBSONIC_OUT


class TestPortExitNotLargerThanThroat:
    @pytest.mark.parametrize("gamma", GAMMAS)
    @pytest.mark.parametrize("ratio", [1.0, 0.9, 0.3])
    def test_port_exit_is_the_sonic_section(self, gamma, ratio):
        A_t = 1.0e-4
        A_e = ratio * A_t
        P = 3.0e6
        mdot, dmdp, T_up, state = _flow(P, A_t, A_e, gamma)
        assert np.isfinite(mdot) and np.isfinite(dmdp) and mdot > 0.0
        assert state == NOZZLE_STATE_CHOKED_OUT
        # M_e = 1: the exit cell flows at its local sound speed, mdot = rho a A_e
        rho = P / (R * T)
        assert mdot == pytest.approx(rho * math.sqrt(gamma * R * T) * A_e, rel=1e-12)
        assert dmdp == pytest.approx(mdot / P, rel=1e-12)
        assert T_up == pytest.approx(T)

    def test_continuous_in_area_ratio_through_one(self):
        gamma = 1.2
        A_t = 1.0e-4
        below = _flow(3.0e6, A_t, A_t * (1.0 - 1e-9), gamma)[0]
        at = _flow(3.0e6, A_t, A_t, gamma)[0]
        above = _flow(3.0e6, A_t, A_t * (1.0 + 1e-9), gamma)[0]
        assert below == pytest.approx(at, rel=1e-7)
        assert above == pytest.approx(at, rel=1e-4)   # sqrt(2 eps) sensitivity at the sonic point


class TestReverseFlowUnchanged:
    @pytest.mark.parametrize("P", [60000.0, 95000.0, 100000.0])
    def test_inflow_identical_in_both_modes(self, P):
        a = _flow(P, 1.0e-4, 3.0e-4, 1.2, isentropic=True)
        b = _flow(P, 1.0e-4, 3.0e-4, 1.2, isentropic=False)
        assert a == b
        assert a[0] < 0.0
        assert a[3] in (NOZZLE_STATE_SUBSONIC_IN, NOZZLE_STATE_CHOKED_IN)
        assert a[2] == pytest.approx(TA)

    def test_balanced_is_zero(self):
        assert _flow(PA, 1.0e-4, 3.0e-4, 1.2)[0] == 0.0


# ---------------------------------------------------------------------------------------
# Steady duct with uniform wall mass addition feeding a choked throat (exact solution).
#
# Constant area A, closed head end, uniform mass source m' [kg/(m s)] injected at rest
# with stagnation temperature T0 (thermal_source = m' Cp T0), no friction.  Steady 1-D
# flow with mass addition at zero axial velocity:
#     momentum:  P + rho u^2 = P_head          (P (1 + gamma M^2) = const)
#     energy:    T0 = const                    (injection at rest, stagnation enthalpy)
# Exit (x = L): mdot = m' L, area ratio A/A_t gives the subsonic M_e of the choked
# throat, P0_e = mdot sqrt(R T0)/(Gamma A_t), P_e = P0_e / g(M_e), and
#     P_head = P_e (1 + gamma M_e^2).
# ---------------------------------------------------------------------------------------
DUCT_GAMMA = 1.2
DUCT_R = 300.0
DUCT_T0 = 3000.0
DUCT_A = 1.0e-3
DUCT_L = 0.5
DUCT_RATIO = 1.6            # A / A_t   (exit Mach ~ 0.4)
DUCT_MDOT = 0.8             # kg/s total


def _duct_exact():
    A_t = DUCT_A / DUCT_RATIO
    M_e = _m_sub_bisect(DUCT_RATIO, DUCT_GAMMA)
    gamma_fn = _critical_flow_function(DUCT_GAMMA)
    P0_e = DUCT_MDOT * math.sqrt(DUCT_R * DUCT_T0) / (gamma_fn * A_t)
    P_e = P0_e / _g(M_e, DUCT_GAMMA)
    return dict(A_t=A_t, M_e=M_e, P0_e=P0_e, P_e=P_e,
                P_head=P_e * (1.0 + DUCT_GAMMA * M_e ** 2))


def _run_duct(N, isentropic, t_end=0.03, max_steps=60000):
    ex = _duct_exact()
    dx = DUCT_L / N
    Cp = DUCT_GAMMA * DUCT_R / (DUCT_GAMMA - 1.0)
    A = np.full(N, DUCT_A)
    P = np.full(N, 0.8 * ex['P_head'])
    T = np.full(N, DUCT_T0)
    rho = P / (DUCT_R * T)
    u = np.zeros(N + 1)
    msrc = np.full(N, DUCT_MDOT / DUCT_L)
    tsrc = msrc * Cp * DUCT_T0
    gam = np.full(N, DUCT_GAMMA)
    Rr = np.full(N, DUCT_R)
    Cpa = np.full(N, Cp)
    Tc = np.full(N, 1.0e5)
    D = np.full(N, 0.05)
    f = np.zeros(N)
    mom = np.zeros(N + 1)
    t = 0.0
    n = 0
    a = math.sqrt(DUCT_GAMMA * DUCT_R * DUCT_T0)
    while t < t_end and n < max_steps:
        dt = compute_dt_cfl(u, a, dx, N + 1, 0.3, 1.0)
        rho, u, P, T = piso_step(
            rho, u, P, T, A, D, msrc, tsrc, mom, f, dx, dt, gam, Rr, Cpa, Tc,
            ex['A_t'], PA, 300.0, N,
            total_energy=True, conservative_momentum=True,
            nozzle_isentropic=isentropic)
        t += dt
        n += 1
    return ex, P, T, u


@pytest.mark.parametrize("N", [25, 50, 100])
def test_steady_duct_head_pressure(N):
    ex, P_new, _, _ = _run_duct(N, True)
    _, P_old, _, _ = _run_duct(N, False)
    err_new = P_new[0] / ex['P_head'] - 1.0
    err_old = P_old[0] / ex['P_head'] - 1.0
    # The old boundary understates P0 by about f^(k-1/2) = 1.09 at M_e ~ 0.4 and the
    # chamber pressure settles that much too HIGH; the new one matches to discretisation
    # error (the exit cell centre sits dx/2 short of the exit).
    assert err_old > 0.05
    assert abs(err_new) < 0.02 + 0.5 / N
    assert abs(err_new) < 0.3 * abs(err_old)
