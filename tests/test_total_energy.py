"""
test_total_energy.py -- EN1 opt-in conservative total-energy PISO step and the
opt-in conservative momentum predictor, checked against exact solutions.

The probe ``docs/solver_verification/probes/en1_verification.py`` drives the
shipped kernel on four problems with known solutions (Sod shock tube, steady
isentropic quasi-1-D nozzle, normal shock in the divergent section, closed-vessel
filling).  These tests pin the headline results on coarse grids, and pin that the
defaults (``total_energy=False``, ``conservative_momentum=False``) are the shipped
step.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "docs", "solver_verification", "probes"))

import en1_verification as ev  # noqa: E402
from srm_1d.solver import _piso_step_with_energy_diagnostics as piso  # noqa: E402


# ---------------------------------------------------------------- defaults
def _vessel_step_args():
    N, dx, A = 12, 0.05, 1.0e-3
    P = np.full(N, 2.0e5)
    T = np.full(N, 600.0)
    rho = P / (287.0 * T)
    u = np.linspace(0.0, 30.0, N + 1)
    src = np.zeros(N)
    src[0:2] = 4.0
    return (rho, u, P, T, np.full(N, A), np.full(N, 0.035), src,
            src * 1004.5 * 2000.0, np.zeros(N + 1), np.full(N, 0.02), dx, 1.0e-6,
            np.full(N, 1.4), np.full(N, 287.0), np.full(N, 1004.5),
            np.full(N, 5000.0), 1.0e-4, 101325.0, 300.0, N, 0.0)


class TestDefaultsAreTheShippedStep:
    def test_explicit_false_flags_equal_omitted(self):
        a = piso(*_vessel_step_args())
        b = piso(*_vessel_step_args(), total_energy=False,
                 conservative_momentum=False)
        # the flow state is bit-identical; the diagnostics may differ at round-off
        # (omitted vs explicit arguments compile to different fastmath specialisations)
        for x, y in zip(a[:4], b[:4]):
            assert np.array_equal(np.asarray(x), np.asarray(y))
        for x, y in zip(a[4:], b[4:]):
            assert np.allclose(np.asarray(x), np.asarray(y), rtol=1e-9,
                               atol=1e-9 * abs(a[5]))

    def test_total_energy_changes_the_step(self):
        a = piso(*_vessel_step_args())
        b = piso(*_vessel_step_args(), total_energy=True)
        assert not np.array_equal(a[3], b[3])

    def test_conservative_momentum_changes_the_step(self):
        a = piso(*_vessel_step_args())
        b = piso(*_vessel_step_args(), conservative_momentum=True)
        assert not np.array_equal(a[1], b[1])

    def test_total_energy_ledger_closes(self):
        out = piso(*_vessel_step_args(), total_energy=True)
        E_before, E_after = out[4], out[5]
        # residual = dE/dt - flux - source - clip: round-off only
        assert abs(out[11]) < 1e-9 * abs(E_before)
        assert E_after > 0.0


# ---------------------------------------------------------------- closed vessel
class TestClosedVessel:
    def test_total_form_is_exact(self):
        r = ev.vessel_case(True, 2.0e-5)
        assert r['T'][0] == pytest.approx(r['T'][1], rel=1e-9)
        assert r['P'][0] == pytest.approx(r['P'][1], rel=1e-9)
        assert r['m'][0] == pytest.approx(r['m'][1], rel=1e-9)

    def test_sensible_form_misses_the_pressure_work(self):
        r = ev.vessel_case(False, 2.0e-5)
        # the shipped form leaves out V dp/dt: T -> mass-weighted mean of the inlet
        # and initial temperatures instead of the energy-balance value.
        assert r['T'][0] == pytest.approx(r['T_sensible_limit'], rel=1e-9)
        assert r['T'][0] < 0.8 * r['T'][1]
        assert r['P'][0] < 0.8 * r['P'][1]


# ---------------------------------------------------------------- Sod
@pytest.fixture(scope="module")
def sod():
    out = {}
    for total, cm in ((False, False), (True, False), (False, True), (True, True)):
        r = ev.sod_case(200, total, cons_mom=cm)
        out[(total, cm)] = (r, ev.sod_metrics(r))
    return out


class TestSod:
    def test_total_plus_conservative_momentum_matches_exact(self, sod):
        _, m = sod[(True, True)]
        assert m['S_num'] == pytest.approx(m['S_exact'], rel=0.04)
        for k in ('rho', 'u', 'P', 'T'):
            pl = m['plateaus']['star-L'][k]
            assert pl[0] == pytest.approx(pl[1], rel=0.05)
            pr = m['plateaus']['star-R'][k]
            assert pr[0] == pytest.approx(pr[1], rel=0.05)
        assert all(v < 0.04 for v in m['l1'].values())

    def test_legacy_momentum_has_wrong_jump_conditions_even_with_total_energy(self, sod):
        _, m = sod[(True, False)]
        assert abs(m['S_num'] / m['S_exact'] - 1.0) > 0.20

    def test_sensible_energy_has_wrong_temperature_even_with_conservative_momentum(self, sod):
        _, m = sod[(False, True)]
        Tn, Te = m['plateaus']['star-L']['T']
        assert Tn > 1.25 * Te      # no expansion cooling: T stays at the initial value
        assert m['l1']['T'] > 0.08

    def test_conservative_scheme_error_shrinks_with_dt(self):
        errs = []
        for cfl in (0.3, 0.1):
            r = ev.sod_case(200, True, cfl=cfl, cons_mom=True)
            m = ev.sod_metrics(r)
            errs.append(abs(m['plateaus']['star-L']['P'][0] / m['plateaus']['star-L']['P'][1] - 1.0))
        assert errs[1] < 0.6 * errs[0]


# ---------------------------------------------------------------- quasi-1-D duct
@pytest.fixture(scope="module")
def venturi():
    return {k: ev.duct_case(100, k[0], throat_M=0.8, cons_mom=k[1], t_end=0.03)
            for k in ((False, False), (True, False), (True, True))}


class TestIsentropicVenturi:
    """Steady subsonic converging-diverging flow (throat Mach 0.8): T/T0 = 1/(1+(g-1)/2 M^2)."""

    def test_total_energy_gets_the_static_temperature(self, venturi):
        mt = ev.duct_metrics(venturi[(True, False)])
        ms = ev.duct_metrics(venturi[(False, False)])
        assert not venturi[(True, False)]['diverged']
        assert mt['L1_T'] < 0.005
        assert ms['L1_T'] > 0.03
        Tn, Te = mt['T_over_T0_exit']
        assert Tn == pytest.approx(Te, rel=0.005)
        Tn_s, Te_s = ms['T_over_T0_exit']
        assert Tn_s == pytest.approx(1.0, abs=1e-6)     # no cooling at all
        assert Te_s < 0.99

    def test_momentum_form_does_not_matter_in_steady_flow(self, venturi):
        a = ev.duct_metrics(venturi[(True, False)])
        b = ev.duct_metrics(venturi[(True, True)])
        for k in ('L1_M', 'L1_T', 'L1_P'):
            assert a[k] == pytest.approx(b[k], rel=0.15)


@pytest.fixture(scope="module")
def shock():
    return {k: ev.duct_case(100, k[0], shock_ar=1.5, cons_mom=k[1], t_end=0.035)
            for k in ((False, False), (True, False), (True, True))}


class TestNormalShockInDivergentSection:
    def test_total_energy_shock_station_and_post_shock_state(self, shock):
        r = shock[(True, False)]
        assert not r['diverged']
        m = ev.duct_metrics(r)
        dx = r['x'][1] - r['x'][0]
        xs_num, xs_ex = m['x_shock']
        assert abs(xs_num - xs_ex) < 2.0 * dx
        assert m['post_shock_P_ratio'][0] == pytest.approx(m['post_shock_P_ratio'][1], rel=0.01)
        assert m['post_shock_T_ratio'][0] == pytest.approx(m['post_shock_T_ratio'][1], rel=0.01)
        assert m['post_shock_M'][0] == pytest.approx(m['post_shock_M'][1], rel=0.02)
        # supersonic isentropic part upstream of the shock
        assert m['L1up_T'] < 0.01 and m['L1up_P'] < 0.02

    def test_sensible_form_has_the_wrong_shock(self, shock):
        m = ev.duct_metrics(shock[(False, False)])
        assert abs(m['x_shock'][0] - m['x_shock'][1]) > 0.05
        assert abs(m['post_shock_P_ratio'][0] / m['post_shock_P_ratio'][1] - 1.0) > 0.3
        assert m['post_shock_T_ratio'][0] == pytest.approx(1.0, abs=1e-4)   # T stays T0
        assert m['L1up_T'] > 0.05

    def test_steady_shock_is_independent_of_the_momentum_form(self, shock):
        a = ev.duct_metrics(shock[(True, False)])
        b = ev.duct_metrics(shock[(True, True)])
        assert a['x_shock'][0] == pytest.approx(b['x_shock'][0], abs=2e-3)


# ---------------------------------------------------------------- run_simulation plumbing
def test_run_simulation_rejects_unknown_formulations():
    from srm_1d.openmotor_adapter import run_from_ric
    ric = os.path.join(ROOT, "motors", "machbusterNew.ric")
    with pytest.raises(ValueError):
        run_from_ric(ric, t_max=0.001, verbose=False, energy_formulation='bogus')
    with pytest.raises(ValueError):
        run_from_ric(ric, t_max=0.001, verbose=False, momentum_formulation='bogus')


@pytest.fixture(scope="module")
def short_motor_runs():
    from srm_1d.openmotor_adapter import run_from_ric
    ric = os.path.join(ROOT, "motors", "machbusterNew.ric")
    out = {}
    for ef, mf in (('sensible', 'legacy'), ('total', 'legacy'), ('total', 'conservative')):
        out[(ef, mf)] = run_from_ric(
            ric, target_propellant_cells=60, t_max=0.012, roughness=32e-6, kappa=0.44,
            T_ignition=756.0, print_interval=0.0, verbose=False, pyrogen='bpnv',
            pyrogen_mass=2.623e-3, energy_formulation=ef,
            momentum_formulation=mf)[0]
    return out


class TestMotorPlumbing:
    def test_stamps_in_summary(self, short_motor_runs):
        for (ef, mf), r in short_motor_runs.items():
            assert r['summary']['energy_formulation'] == ef
            assert r['summary']['momentum_formulation'] == mf

    def test_total_energy_ledger_closes_in_a_motor_run(self, short_motor_runs):
        r = short_motor_runs[('total', 'legacy')]
        E = np.asarray(r['gas_sensible_energy'])
        res = np.asarray(r['energy_residual'])
        # the PISO residual is round-off relative to the stored energy
        assert np.max(np.abs(res)) < 1e-3 * max(np.max(np.abs(E)), 1.0)
        # source ledger is formulation independent
        assert np.max(np.abs(np.asarray(r['source_ledger_mismatch_W']))) < 1e-4

    def test_runs_stay_finite_and_bounded(self, short_motor_runs):
        for r in short_motor_runs.values():
            P = np.asarray(r['P_head'])
            assert np.all(np.isfinite(P)) and P.max() < 60e6
            assert r['summary']['termination'] == 't_max reached'
