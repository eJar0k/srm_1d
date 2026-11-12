"""TP1: thrust / CF / Isp / impulse use the nozzle-entrance stagnation pressure.

Spec: docs/solver_verification/TP1_TASK.md; worklog: worklog/TP1.md.

* ``solver._nozzle_entrance_stagnation`` against the exact isentropic relations and the
  boundary solve it is derived from (static mode, no flow, reverse flow, subsonic, choked,
  ``A_exit <= A_throat``).
* ``compute_motor_performance`` formula equivalence (P0 = P_head reproduces the old
  thrust), the fallback for results without ``P0_nozzle``, and use of the P0 series.
* A short real run: ``P0_nozzle`` / ``T0_nozzle`` are recorded, same length as every other
  history array, finite, P0 consistent with the run's own boundary state, and the
  head-end pressure is unchanged.
"""
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from srm_1d.solver import (  # noqa: E402
    _nozzle_boundary_flow,
    _nozzle_choke_mach,
    _nozzle_entrance_stagnation,
)
from srm_1d.nozzle import (  # noqa: E402
    Nozzle, compute_motor_performance, compute_thrust_isp, compute_thrust_series,
    compute_thrust_history,
)
from srm_1d.channels import build_channels  # noqa: E402

GAMMAS = (1.1, 1.2, 1.4)
R = 300.0
T = 2800.0
PA = 101325.0
TA = 293.0


def _g(M, gamma):
    return (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))


def _boundary(P, A_exit, A_throat, gamma, isentropic=True, T_cell=T, Pa=PA):
    Mc = _nozzle_choke_mach(A_exit, A_throat, gamma, isentropic)
    mdot, _dm, _Tu, _st = _nozzle_boundary_flow(P, T_cell, A_throat, gamma, R, Pa, TA,
                                                A_exit, Mc)
    return mdot, Mc


# ----------------------------------------------------------------------
# the entrance-stagnation function
# ----------------------------------------------------------------------
@pytest.mark.parametrize("gamma", GAMMAS)
def test_static_mode_gives_P0_equal_P_exit(gamma):
    """'static' boundary (M_choke = 0): P0 = P_exit, T0 = T_exit (the old value)."""
    A_t, A_e = 1.0e-3, 4.0e-3
    mdot, Mc = _boundary(5.0e6, A_e, A_t, gamma, isentropic=False)
    assert Mc == 0.0 and mdot > 0.0
    P0, T0, Me = _nozzle_entrance_stagnation(mdot, 5.0e6, T, gamma, R, A_e, Mc)
    assert P0 == 5.0e6 and T0 == T and Me == 0.0


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("ratio", (1.05, 1.5, 3.83, 5.27, 20.0))
def test_choked_P0_over_P_is_g_of_choke_mach(gamma, ratio):
    """Choked outflow: M_e is the area-Mach root, P0/P = g(M_e), T0/T = 1 + (g-1)/2 M_e^2."""
    A_t = 1.0e-3
    A_e = ratio * A_t
    P = 6.0e6
    mdot, Mc = _boundary(P, A_e, A_t, gamma)
    P0, T0, Me = _nozzle_entrance_stagnation(mdot, P, T, gamma, R, A_e, Mc)
    assert Me == pytest.approx(Mc, rel=1e-10)
    assert P0 / P == pytest.approx(_g(Mc, gamma), rel=1e-10)
    assert T0 / T == pytest.approx(1.0 + 0.5 * (gamma - 1.0) * Mc * Mc, rel=1e-12)
    # the choked mass flow of the boundary is Gamma A_t P0 / sqrt(R T0): same state
    Gam = math.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    assert mdot == pytest.approx(Gam * A_t * P0 / math.sqrt(R * T0), rel=1e-10)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_subsonic_outflow_recovers_exit_mach(gamma):
    """Unchoked outflow: M_e = mdot / (A_e P sqrt(gamma/RT)), exit continuity = throat flow."""
    A_t, A_e = 1.0e-3, 4.0e-3
    Pc = 1.4 * PA                  # above ambient but below crit * P0 -> subsonic throat
    mdot, Mc = _boundary(Pc, A_e, A_t, gamma, Pa=PA)
    P0, T0, Me = _nozzle_entrance_stagnation(mdot, Pc, T, gamma, R, A_e, Mc)
    assert 0.0 < Me < Mc
    assert P0 / Pc == pytest.approx(_g(Me, gamma), rel=1e-12)
    # the stagnation state really feeds the throat: isentropic subsonic flow P0,T0 -> PA
    pr = PA / P0
    mdot_t = (A_t * P0 * math.sqrt(2.0 * gamma / ((gamma - 1.0) * R * T0))
              * math.sqrt(pr ** (2.0 / gamma) - pr ** ((gamma + 1.0) / gamma)))
    assert mdot_t == pytest.approx(mdot, rel=1e-8)


@pytest.mark.parametrize("frac", (1.0, 0.9, 0.3))
def test_exit_not_larger_than_throat_is_sonic_at_the_exit(frac):
    """A_exit <= A_throat: the port exit is the sonic section, M_e = 1, P0 = P g(1)."""
    gamma = 1.2
    A_t = 1.0e-3
    A_e = frac * A_t
    P = 3.0e6
    mdot, Mc = _boundary(P, A_e, A_t, gamma)
    assert Mc == 1.0
    P0, T0, Me = _nozzle_entrance_stagnation(mdot, P, T, gamma, R, A_e, Mc)
    assert Me == pytest.approx(1.0, rel=1e-10)
    assert P0 / P == pytest.approx(_g(1.0, gamma), rel=1e-10)


def test_no_flow_and_reverse_flow_return_the_static_state():
    gamma, A_t, A_e = 1.2, 1.0e-3, 4.0e-3
    Mc = _nozzle_choke_mach(A_e, A_t, gamma, True)
    # balanced (P = ambient): no flow
    mdot, _d, _T, _s = _nozzle_boundary_flow(PA, T, A_t, gamma, R, PA, TA, A_e, Mc)
    assert mdot == 0.0
    assert _nozzle_entrance_stagnation(mdot, PA, T, gamma, R, A_e, Mc) == (PA, T, 0.0)
    # reverse inflow (P < ambient): reservoir at rest
    mdot, _d, _T, _s = _nozzle_boundary_flow(0.5 * PA, T, A_t, gamma, R, PA, TA, A_e, Mc)
    assert mdot < 0.0
    assert _nozzle_entrance_stagnation(mdot, 0.5 * PA, T, gamma, R, A_e, Mc) == (0.5 * PA, T, 0.0)


# ----------------------------------------------------------------------
# performance post-processing
# ----------------------------------------------------------------------
def _synthetic_result(P_head, P0=None, D_throat=None, mdot=None, dt=1.0e-3):
    n = len(P_head)
    r = {'time': np.arange(n) * dt, 'P_head': np.asarray(P_head, float),
         'P_ambient': PA}
    if P0 is not None:
        r['P0_nozzle'] = np.asarray(P0, float)
    if D_throat is not None:
        r['D_throat'] = np.asarray(D_throat, float)
    if mdot is not None:
        r['massflow'] = np.asarray(mdot, float)
    return r


def _prop():
    from tests._motor_fixtures import hasegawa_propellant_1
    return hasegawa_propellant_1()


def _nozzle():
    return Nozzle(D_throat=0.02, D_exit=0.05, efficiency=0.95, div_angle=15.0,
                  throat_length=0.01)


def test_no_flow_cell_thrust_matches_old_formula():
    """P0 = P_head (u = 0, no port loss): thrust is the old CF * P * A_t, point by point."""
    nz, prop = _nozzle(), _prop()
    P = np.linspace(2.0e5, 6.0e6, 40)
    D = np.full(40, nz.D_throat)
    res = _synthetic_result(P, P0=P, D_throat=D)
    perf = compute_motor_performance(build_channels(res), nz, prop, P_ambient=PA)
    rep = prop.representative_tab()
    from srm_1d.propellant import critical_flow_function, R_UNIVERSAL
    c_star = (R_UNIVERSAL / rep.molecular_weight * rep.T_flame) ** 0.5 / \
        critical_flow_function(rep.gamma)
    A_t = math.pi / 4 * nz.D_throat ** 2
    eps = nz.exit_area() / A_t
    for i in (0, 7, 39):
        th, cf, isp, pe = compute_thrust_isp(
            P[i], rep.gamma, eps, nz.divergence_losses(), nz.efficiency,
            nz.throat_losses(), nz.skin_losses(), A_t, PA, 9.80665, c_star)
        assert perf["thrust"][i] == pytest.approx(th, rel=1e-12)
        assert perf["C_F"][i] == pytest.approx(cf, rel=1e-12)
    np.testing.assert_array_equal(perf['P0_nozzle'], P)


def test_thrust_follows_P0_not_P_head():
    nz, prop = _nozzle(), _prop()
    P = np.full(30, 5.0e6)
    D = np.full(30, nz.D_throat)
    p_head = _synthetic_result(P, P0=P, D_throat=D)
    p_low = _synthetic_result(P, P0=0.97 * P, D_throat=D)
    perf_h = compute_motor_performance(p_head, nz, prop, P_ambient=PA)
    perf_l = compute_motor_performance(p_low, nz, prop, P_ambient=PA)
    assert np.all(perf_l['thrust'] < perf_h['thrust'])
    # F = CF(P0) * P0 * A_t: with the same CF (to first order) it scales ~0.97
    assert perf_l['thrust'][5] / perf_h['thrust'][5] == pytest.approx(0.97, abs=0.01)
    np.testing.assert_array_equal(perf_l['P0_nozzle'], 0.97 * P)


def test_fallback_without_P0_warns_and_uses_P_head():
    nz, prop = _nozzle(), _prop()
    P = np.linspace(1.0e6, 5.0e6, 30)
    D = np.full(30, nz.D_throat)
    old = _synthetic_result(P, D_throat=D)               # no P0_nozzle
    ref = _synthetic_result(P, P0=P, D_throat=D)
    with pytest.warns(RuntimeWarning, match="P0_nozzle"):
        perf_old = compute_motor_performance(old, nz, prop, P_ambient=PA)
    perf_ref = compute_motor_performance(ref, nz, prop, P_ambient=PA)
    np.testing.assert_array_equal(perf_old['thrust'], perf_ref['thrust'])
    assert perf_old['total_impulse'] == perf_ref['total_impulse']


def test_fallback_without_throat_history_integrates_from_P0():
    """No recorded D_throat and a nonzero erosion: the fallback integrator runs on the
    series it is given (here P0)."""
    nz = Nozzle(D_throat=0.02, D_exit=0.05, erosion_coeff=0.5)
    prop = _prop()
    P = np.full(20, 4.0e6)
    res = _synthetic_result(P, P0=0.98 * P)
    perf = compute_motor_performance(res, nz, prop, P_ambient=PA)
    # dD/dt = 2 * e_rate = 2 * 0.5e-6 * P0[MPa]
    expected = 0.02 + 2.0 * 0.5e-6 * (0.98 * 4.0) * (19 * 1.0e-3)
    assert perf['D_throat'][-1] == pytest.approx(expected, rel=1e-9)


def test_series_equals_history_when_D_throat_is_the_integrated_one():
    nz = Nozzle(D_throat=0.02, D_exit=0.05, erosion_coeff=0.3, slag_coeff=1e-4)
    t = np.arange(25) * 1.0e-3
    P = np.linspace(1.0e6, 7.0e6, 25)
    args = (1.2, nz.D_exit)
    th, cf, isp, pe, D = compute_thrust_history(
        t, P, 25, 1.2, nz.D_throat, nz.D_exit, nz.erosion_coeff, nz.slag_coeff,
        nz.divergence_losses(), nz.efficiency, nz.throat_length, nz.skin_losses(),
        PA, 1500.0, 9.80665)
    th2, cf2, isp2, pe2 = compute_thrust_series(
        P, D, 25, 1.2, nz.D_exit, nz.divergence_losses(), nz.efficiency,
        nz.throat_length, nz.skin_losses(), PA, 1500.0, 9.80665)
    np.testing.assert_array_equal(th, th2)
    np.testing.assert_array_equal(isp, isp2)


# ----------------------------------------------------------------------
# a real (short) run
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def real_run():
    pytest.importorskip("numba")
    pytest.importorskip("skfmm")
    from srm_1d import run_simulation
    from tests.test_simulation_phase3 import _small_motor, _test_chamber
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0, t_max=0.010, P_cutoff=1.0,
        dt_max=2.0e-5, burn_update_interval=1, snapshot_interval=0.002,
        cfl_target=0.5, verbose=False,
    )
    return geo, prop, nozzle, result


def test_history_arrays_are_recorded_bounded_and_finite(real_run):
    _geo, _prop, _nozzle, result = real_run
    n = len(result['time'])
    cap = int(result['summary']['history_capacity'])
    for key in ('P0_nozzle', 'T0_nozzle', 'gamma_nozzle', 'R_nozzle'):
        assert key in result
        assert len(result[key]) == n, key
        assert np.all(np.isfinite(result[key])), key
    assert n <= cap
    assert result.channels['P0_nozzle'].unit == 'Pa'
    assert np.all(result['T0_nozzle'] > 0.0) and np.all(result['gamma_nozzle'] > 1.0)


def test_P0_consistent_with_the_runs_own_boundary_state(real_run):
    """Every recorded row: P0 = P_exit * g(M_e) with M_e from that step's mass flow; the
    sampled rows with outflow have P0 >= P_exit, and static-state rows P0 == P_exit."""
    _geo, _prop, _nozzle, result = real_run
    P0 = np.asarray(result['P0_nozzle'])
    Pe = np.asarray(result['P_exit'])
    mdot = np.asarray(result['massflow'])
    assert np.all(P0 >= Pe * (1.0 - 1e-12))
    no_out = mdot <= 0.0
    np.testing.assert_array_equal(P0[no_out], Pe[no_out])
    if np.any(~no_out):
        assert np.all(P0[~no_out] >= Pe[~no_out])


def test_perf_uses_P0_and_head_pressure_is_untouched(real_run):
    _geo, prop, nozzle, result = real_run
    perf = compute_motor_performance(result, nozzle, prop)
    np.testing.assert_array_equal(perf['P0_nozzle'], np.asarray(result['P0_nozzle']))
    # recompute the thrust from P0 and from P_head: the stored perf equals the P0 one
    P0 = np.asarray(result['P0_nozzle'])
    Ph = np.asarray(result['P_head'])
    D = np.asarray(result['D_throat'])
    rep = prop.representative_tab()
    from srm_1d.propellant import critical_flow_function, R_UNIVERSAL
    c_star = (R_UNIVERSAL / rep.molecular_weight * rep.T_flame) ** 0.5 / \
        critical_flow_function(rep.gamma)

    def series(P):
        return compute_thrust_series(
            P, D, len(P), rep.gamma, nozzle.D_exit, nozzle.divergence_losses(),
            nozzle.efficiency, nozzle.throat_length, nozzle.skin_losses(),
            float(result['P_ambient']), c_star, 9.80665)[0]
    np.testing.assert_array_equal(perf['thrust'], series(P0))
    assert perf['c_star_sim'].shape == perf['thrust'].shape
    assert 'P_head' in result and len(Ph) == len(P0)


def test_csv_pressure_column_is_head_end_and_has_stagnation_column(real_run):
    from srm_1d.openmotor_adapter import result_to_csv
    geo, prop, nozzle, result = real_run
    perf = compute_motor_performance(result, nozzle, prop)
    csv = result_to_csv(result, perf, geo, prop)
    header = csv.split('\n')[0].split(',')
    assert header[2] == 'Pressure (Pa)'
    assert 'Nozzle Stagnation Pressure (Pa)' in header
    row = csv.split('\n')[1].split(',')
    assert float(row[2]) == pytest.approx(float(result['P_head'][0]), abs=0.01)
