"""
test_verification_limits.py -- solver-verification T4 (docs/solver_verification).

Pure verification tests (no source changes) pinning four things the Codex
review of srm_1d reported and that our tree should provably satisfy:

  (a) Tab-boundary parity with openMotor: srm_1d's hard-switchover tab lookup
      (``select_tab_idx`` / ``saint_robert_from_tabs``, plus the Python-level
      ``Propellant.select_tab``) picks the SAME tab as openMotor's
      ``Propellant.getCombustionProperties`` for boundary, shared-boundary,
      in-gap, below-all and above-all pressures (and 1-ulp neighbours of
      every boundary), in every tab ordering.
  (b) Zero crossflow (Re below the ``Re < 100`` guard): the burn rate is
      exactly Saint-Robert and the erosive increment is exactly zero, for
      every registry erosive-model code.
  (c) Normal-component invariance: the erosive models only ADD to the
      Saint-Robert rate (r_total - r_erosive == r0, r_erosive >= 0); the
      saint_robert code is erosive-free everywhere.
  (d) Mode invariants (integration): across the three registry erosive
      models a full-burnout motor delivers the same propellant mass and
      (to a few %) the same total impulse -- erosive burning redistributes
      the burn in time, it does not create or destroy propellant.

Part (a) needs the openMotor checkout (``motorlib``) and is skipped cleanly
if it is unavailable; parts (b)-(d) do not.
"""

import os

import numpy as np
import pytest

pytest.importorskip("numba")

from srm_1d.burn_rate import (  # noqa: E402
    burn_rate_cell, select_tab_idx, saint_robert_from_tabs,
    EROSIVE_MODELS, resolve_erosive_model,
)

MOTORS = os.path.join(os.path.dirname(__file__), '..', 'motors')

MPA = 1.0e6

# ================================================================
# Shared constants: Hasegawa Propellant 1 -- the same values the
# existing direct ``burn_rate_cell`` call sites in tests/test_burn_rate.py
# (``_HasegawaPropellant1Fixture``) use.
# ================================================================

_P_REF = 4.9e6
_R_REF = 4.9e-3
_N_SR = 0.3
_A_SR = _R_REF / (_P_REF ** _N_SR)
RHO_P = 1700.0
CPS = 1500.0
T_SURFACE = 1000.0
T_INITIAL = 293.0
T_FLAME = 3041.0
PR = 0.4943
K_THERMAL = 0.3685
CP_GAS = 2060.0
KAPPA = 0.45
MU_GAS = 9.0e-5
ROUGHNESS = 20e-6

# Tab sets as (min_p, max_p, a, n) arrays, the layout ``Propellant.tab_arrays``
# produces and ``burn_rate_cell`` consumes.
SINGLE_TAB = (
    np.array([0.0]), np.array([2e7]), np.array([_A_SR]), np.array([_N_SR]),
)

# Three tabs with a gap (the part-(a) openMotor propellant, a/n in the
# srm_1d / openMotor shared convention: m/s per Pa^n, Pa).
_TAB_SPEC = (
    # (minPressure [Pa], maxPressure [Pa], a [m/(s*Pa^n)], n)
    (0.1 * MPA, 2.0 * MPA, 3.0e-5, 0.40),
    (2.0 * MPA, 6.0 * MPA, 5.5e-5, 0.32),
    (8.0 * MPA, 15.0 * MPA, 1.2e-4, 0.20),
)
THREE_TAB = tuple(np.array(col) for col in zip(*_TAB_SPEC))


def _cell(P, Re, D, code, tabs=SINGLE_TAB, x=0.5, cap=0.0):
    """Direct ``burn_rate_cell`` call with the Hasegawa-1 constants."""
    mn, mx, a, n = tabs
    return burn_rate_cell(
        P, Re, D, x, ROUGHNESS,
        PR, K_THERMAL, CP_GAS, T_FLAME, T_SURFACE,
        RHO_P, CPS, T_INITIAL,
        mn, mx, a, n, len(mn),
        KAPPA, MU_GAS, cap, code,
    )


# ================================================================
# (a) Tab-boundary parity vs openMotor
# ================================================================

# Tab-list orderings: the strict-containment-then-closest-boundary tie-break
# (first tab wins an exact tie) is order dependent in BOTH implementations,
# so parity must hold in every ordering, not just the ascending one.
_ORDERINGS = {
    'ascending': (0, 1, 2),
    'descending': (2, 1, 0),
    'shuffled': (1, 2, 0),
}


def _named_pressures():
    """(label, P[Pa]) cases the task spec calls out."""
    cases = []
    for k, (lo, hi, _, _) in enumerate(_TAB_SPEC):
        cases.append((f'tab{k}_min', lo))
        cases.append((f'tab{k}_max', hi))
        cases.append((f'tab{k}_mid', 0.5 * (lo + hi)))
    cases += [
        ('shared_boundary_2MPa', 2.0 * MPA),
        ('gap_7MPa_equidistant', 7.0 * MPA),
        ('gap_6.9MPa', 6.9 * MPA),
        ('gap_7.1MPa', 7.1 * MPA),
        ('below_all_0.05MPa', 0.05 * MPA),
        ('above_all_20MPa', 20.0 * MPA),
    ]
    # 1-ulp neighbours of every distinct boundary: pins strict '<' (an
    # '<=' implementation would flip the selection exactly at the edge).
    for edge in sorted({v for lo, hi, _, _ in _TAB_SPEC for v in (lo, hi)}):
        cases.append((f'{edge / MPA:g}MPa_minus_ulp', float(np.nextafter(edge, 0.0))))
        cases.append((f'{edge / MPA:g}MPa_plus_ulp', float(np.nextafter(edge, np.inf))))
    return cases


_NAMED_PRESSURES = _named_pressures()


@pytest.fixture(scope="module")
def om_motorlib():
    """openMotor's ``motorlib`` (skip the part-(a) tests cleanly if the
    checkout is unavailable). Imported the way tests/test_srm1d_plugin.py
    does."""
    try:
        from srm_1d.fmm_grain import _setup_openmotor_path
        _setup_openmotor_path()
        from motorlib.propellant import Propellant as OMPropellant  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"openMotor checkout unavailable: {exc}")
    return OMPropellant


def _propellant_dict(order):
    """A propellant dict in openMotor's own schema, tabs in ``order``.

    Per-tab gas transport is provided so the srm_1d adapter reads it from the
    propellant itself (no ``gas_props`` override): one dict feeds BOTH
    openMotor and the srm_1d adapter.
    """
    tabs = []
    for k in order:
        lo, hi, a, n = _TAB_SPEC[k]
        tabs.append({
            'minPressure': lo, 'maxPressure': hi, 'a': a, 'n': n,
            'k': 1.22, 't': 3000.0, 'm': 25.0,
            'mu': 9.0e-5,
            'kThermalFrozen': 0.35, 'cpFrozen': 2000.0,
            'kThermalEffective': 0.65, 'cpEffective': 2700.0,
        })
    return {
        'name': 'T4 three-tab gap propellant',
        'density': 1700.0,
        'transportVariant': 'frozen',
        'tabs': tabs,
    }


@pytest.fixture(scope="module", params=sorted(_ORDERINGS))
def tab_parity_pair(request, om_motorlib):
    """(openMotor Propellant, srm_1d Propellant) built from the SAME data.

    The srm_1d propellant is produced by the adapter's propellant conversion
    (``convert_propellant``) applied to openMotor's own serialization
    (``Propellant.getProperties()`` -- exactly what a ``.ric`` file stores),
    so unit handling is the adapter's, not the test's.
    """
    from srm_1d.openmotor_adapter import convert_propellant
    order = _ORDERINGS[request.param]
    om_prop = om_motorlib(_propellant_dict(order))
    # openMotor's FloatProperty silently IGNORES out-of-range values; make
    # sure it actually holds the data we think it does (no vacuous test).
    held = om_prop.getProperties()['tabs']
    assert len(held) == 3
    for tab_held, k in zip(held, order):
        lo, hi, a, n = _TAB_SPEC[k]
        assert (tab_held['minPressure'], tab_held['maxPressure'],
                tab_held['a'], tab_held['n']) == (lo, hi, a, n)
    srm_prop = convert_propellant(om_prop.getProperties())
    return om_prop, srm_prop


def _assert_same_tab(om_prop, srm_prop, P, label=''):
    om_a, om_n, om_k, _, _ = om_prop.getCombustionProperties(P)
    mn, mx, a_arr, n_arr = srm_prop.tab_arrays()
    n_tabs = len(mn)

    # srm_1d numba lookup selects the same tab (a and n identify it).
    idx = select_tab_idx(P, mn, mx, n_tabs)
    assert (a_arr[idx], n_arr[idx]) == (om_a, om_n), (
        f"{label} P={P!r}: select_tab_idx -> tab {idx} "
        f"(a={a_arr[idx]}, n={n_arr[idx]}) but openMotor picked "
        f"(a={om_a}, n={om_n})")

    # Python-level Propellant.select_tab agrees too.
    tab = srm_prop.select_tab(P)
    assert (tab.a, tab.n) == (om_a, om_n), f"{label} P={P!r}: select_tab"
    assert tab.gamma == om_k

    # The burn rate the kernel uses is openMotor's burn rate at this P.
    r_srm = saint_robert_from_tabs(P, mn, mx, a_arr, n_arr, n_tabs)
    assert r_srm == pytest.approx(om_prop.getBurnRate(P), rel=1e-12), label


@pytest.mark.parametrize('label,P', _NAMED_PRESSURES,
                         ids=[c[0] for c in _NAMED_PRESSURES])
def test_tab_lookup_matches_openmotor_named_pressures(tab_parity_pair, label, P):
    om_prop, srm_prop = tab_parity_pair
    _assert_same_tab(om_prop, srm_prop, P, label)


def test_tab_lookup_matches_openmotor_dense_sweep(tab_parity_pair):
    """2501 pressures spanning below-all .. above-all (includes exact tab
    edges and the gap) -- no disagreement anywhere."""
    om_prop, srm_prop = tab_parity_pair
    for P in np.linspace(0.01 * MPA, 25.0 * MPA, 2501):
        _assert_same_tab(om_prop, srm_prop, float(P), 'sweep')


def test_tab_selection_semantics_pinned(tab_parity_pair):
    """Pin the semantics themselves (not just agreement): strict containment
    first, closest-boundary fallback, shared boundary and the in-gap
    equidistant tie both resolve to the FIRST listed tab that attains the
    minimum distance."""
    om_prop, srm_prop = tab_parity_pair
    mn, mx, a_arr, n_arr = srm_prop.tab_arrays()
    spec_by_a = {a: (lo, hi) for lo, hi, a, _ in _TAB_SPEC}

    def picked_range(P):
        return spec_by_a[a_arr[select_tab_idx(P, mn, mx, len(mn))]]

    # Strict containment: midpoints land in their own tab regardless of order.
    for lo, hi, _, _ in _TAB_SPEC:
        assert picked_range(0.5 * (lo + hi)) == (lo, hi)
    # Closest-boundary fallback outside every tab / in the gap.
    assert picked_range(0.05 * MPA) == (0.1 * MPA, 2.0 * MPA)
    assert picked_range(20.0 * MPA) == (8.0 * MPA, 15.0 * MPA)
    assert picked_range(6.9 * MPA) == (2.0 * MPA, 6.0 * MPA)
    assert picked_range(7.1 * MPA) == (8.0 * MPA, 15.0 * MPA)
    # Exact tab edge is NOT strictly contained: 6 MPa and 8 MPa resolve by
    # distance (0) to the tab that owns that edge.
    assert picked_range(6.0 * MPA) == (2.0 * MPA, 6.0 * MPA)
    assert picked_range(8.0 * MPA) == (8.0 * MPA, 15.0 * MPA)


# ================================================================
# (b) Zero crossflow
# ================================================================

_ALL_CODES = sorted(m.code for m in EROSIVE_MODELS.values())
_TABSETS = {'single_tab': SINGLE_TAB, 'three_tab_gap': THREE_TAB}
# Includes an in-gap pressure (7 MPa) and one above all three-tab ranges.
_PRESSURES_B = (0.5 * MPA, 1.0 * MPA, 5.0 * MPA, 7.0 * MPA, 12.0 * MPA, 20.0 * MPA)


def test_registry_codes_are_the_expected_three():
    """The parametrisation below covers every registered code; fail loudly
    if a new erosive model is registered without being added here."""
    assert _ALL_CODES == [0, 1, 2]
    assert [resolve_erosive_model(n) for n in ('ma2020', 'mukunda_paul', 'saint_robert')] == [0, 1, 2]


@pytest.mark.parametrize('tabset', sorted(_TABSETS))
@pytest.mark.parametrize('code', _ALL_CODES)
@pytest.mark.parametrize('Re', [0.0, 50.0, 99.0])
def test_zero_crossflow_is_exactly_saint_robert(code, Re, tabset):
    """Below the Re < 100 guard no erosive model may add anything: the
    result is bit-exactly the Saint-Robert tab-lookup rate and r_erosive is
    exactly 0 (the guard returns the very same r0 expression)."""
    tabs = _TABSETS[tabset]
    mn, mx, a, n = tabs
    for P in _PRESSURES_B:
        for D in (0.01, 0.04):
            r0 = saint_robert_from_tabs(P, mn, mx, a, n, len(mn))
            assert r0 > 0.0
            r_total, r_erosive = _cell(P, Re, D, code, tabs)
            assert r_total == r0, (code, Re, tabset, P, D, r_total, r0)
            assert r_erosive == 0.0, (code, Re, tabset, P, D, r_erosive)


# ================================================================
# (c) Normal-component invariance
# ================================================================

_PRESSURES_C = (1.0 * MPA, 5.0 * MPA, 12.0 * MPA)
_RE_C = (1e4, 1e5, 1e6)
_D_C = (0.010, 0.040)


@pytest.mark.parametrize('D', _D_C)
@pytest.mark.parametrize('Re', _RE_C)
@pytest.mark.parametrize('P', _PRESSURES_C)
@pytest.mark.parametrize('code', [0, 1])
def test_erosive_models_only_add_to_normal_rate(code, P, Re, D):
    """Ma 2020 (code 0) and Mukunda-Paul (code 1): the erosive component is
    non-negative and the normal part is untouched, r_total - r_erosive == r0."""
    mn, mx, a, n = SINGLE_TAB
    r0 = saint_robert_from_tabs(P, mn, mx, a, n, len(mn))
    r_total, r_erosive = _cell(P, Re, D, code)
    assert r_erosive >= 0.0
    assert r_total >= r0
    assert (r_total - r_erosive) == pytest.approx(r0, rel=1e-12)


@pytest.mark.parametrize('D', _D_C)
@pytest.mark.parametrize('Re', _RE_C)
@pytest.mark.parametrize('P', _PRESSURES_C)
def test_saint_robert_mode_has_no_erosive_component(P, Re, D):
    mn, mx, a, n = SINGLE_TAB
    r0 = saint_robert_from_tabs(P, mn, mx, a, n, len(mn))
    r_total, r_erosive = _cell(P, Re, D, 2)
    assert r_erosive == 0.0
    assert r_total == pytest.approx(r0, rel=1e-12)


@pytest.mark.parametrize('code', [0, 1])
def test_invariance_grid_is_not_vacuous(code):
    """The part-(c) grid must contain states where the model really erodes
    (otherwise ``r_erosive >= 0`` / ``r_total - r_erosive == r0`` could pass
    trivially on r_erosive == 0 everywhere)."""
    mn, mx, a, n = SINGLE_TAB
    n_eroding = 0
    for P in _PRESSURES_C:
        for Re in _RE_C:
            for D in _D_C:
                _, r_erosive = _cell(P, Re, D, code)
                n_eroding += r_erosive > 0.0
    assert n_eroding >= 3, f"code {code}: only {n_eroding} eroding states"


# ================================================================
# (d) Mode invariants (integration)
# ================================================================
#
# Motor choice: ``hlbltest.ric`` -- the cheapest committed motor that runs
# to the end of its burn (termination 'pressure cutoff', not 't_max') in
# well under 30 s per run warm; see docs/solver_verification/worklog/T4.md
# for the candidate timing table (example_bates 90 s, L1843 27 s, ...).

INTEGRATION_MOTOR = 'hlbltest.ric'
_MODES = ('ma2020', 'mukunda_paul', 'saint_robert')
MASS_TOL = 0.005       # propellant mass delivered, across modes
IMPULSE_TOL = 0.03     # total impulse vs ma2020
_FULL_BURN = ('complete burnout', 'pressure cutoff')


@pytest.fixture(scope="module")
def mode_runs():
    """One canonical-knob run of ``INTEGRATION_MOTOR`` per erosive model."""
    from srm_1d.openmotor_adapter import run_from_ric
    from srm_1d.nozzle import compute_motor_performance

    path = os.path.join(MOTORS, INTEGRATION_MOTOR)
    runs = {}
    for mode in _MODES:
        result, perf, nozzle, geo, prop = run_from_ric(
            path, pyrogen='bpnv', pyrogen_mass=None,
            roughness=32e-6, kappa=0.44, T_ignition=756.0,
            print_interval=0.0, verbose=False,
            erosive_model=mode,
        )
        s = result['summary']
        # run_from_ric already returns this same call's result as ``perf``;
        # recompute explicitly (the task's named source for total impulse)
        # and make sure the two agree.
        perf2 = compute_motor_performance(
            result, nozzle, prop, P_ambient=perf['P_ambient'])
        assert perf2['total_impulse'] == perf['total_impulse']
        runs[mode] = {
            'termination': s['termination'],
            'termination_code': s['termination_code'],
            'erosive_model': s['erosive_model'],
            'mass_produced': s['mass_produced'],
            'propellant_mass': s['propellant_mass'],
            'impulse': perf2['total_impulse'],
            'burn_time': perf2['burn_time'],      # thrust > 5% of peak
            't_end': s['t_burn'],                 # simulated time at termination
            'P_peak': s['P_peak'],
        }
    return runs


def _rel(x, ref):
    return (x - ref) / ref


def test_modes_burn_to_completion(mode_runs):
    """Each mode ran the motor out (burnout / pressure cutoff), not into the
    t_max cap, and the summary carries the model that produced it."""
    for mode, r in mode_runs.items():
        assert r['erosive_model'] == mode
        assert r['termination'] in _FULL_BURN, (mode, r['termination'])


def test_mode_invariant_propellant_mass_delivered(mode_runs):
    """Erosive burning redistributes WHEN the propellant burns, not how much:
    mass_produced agrees within 0.5% across all three modes."""
    m_ref = mode_runs['ma2020']['mass_produced']
    masses = [r['mass_produced'] for r in mode_runs.values()]
    spread = (max(masses) - min(masses)) / m_ref
    deltas = {m: _rel(r['mass_produced'], m_ref) for m, r in mode_runs.items()}
    assert spread <= MASS_TOL, f"mass spread {spread:.3%}; per-mode vs ma2020 {deltas}"


def test_mode_invariant_total_impulse(mode_runs):
    """Total impulse is within 3% of ma2020's for every mode."""
    i_ref = mode_runs['ma2020']['impulse']
    for mode, r in mode_runs.items():
        d = _rel(r['impulse'], i_ref)
        assert abs(d) <= IMPULSE_TOL, f"{mode}: impulse {r['impulse']:.1f} vs ma2020 {i_ref:.1f} ({d:+.2%})"


def test_modes_actually_differ_in_the_burn(mode_runs):
    """Non-vacuity guard: switching the erosive model must change the burn
    history (peak pressure or burn duration), otherwise the invariants above
    would hold trivially because the knob did nothing."""
    ma, sr = mode_runs['ma2020'], mode_runs['saint_robert']
    d_peak = abs(_rel(sr['P_peak'], ma['P_peak']))
    d_burn = abs(_rel(sr['t_end'], ma['t_end']))
    assert max(d_peak, d_burn) > 1e-3, (d_peak, d_burn)


def test_mode_deltas_recorded(mode_runs, capsys):
    """Not an assertion on physics: prints the per-mode deltas vs ma2020 (run
    with ``-s`` to see them; the numbers go in worklog/T4.md)."""
    ref = mode_runs['ma2020']
    lines = [f"[T4 mode invariants] motor={INTEGRATION_MOTOR}"]
    for mode in _MODES:
        r = mode_runs[mode]
        lines.append(
            f"  {mode:13s} term={r['termination']!r} "
            f"m={r['mass_produced']:.5f}kg ({_rel(r['mass_produced'], ref['mass_produced']):+.4%}) "
            f"I={r['impulse']:.2f}Ns ({_rel(r['impulse'], ref['impulse']):+.4%}) "
            f"t_burn(5%F)={r['burn_time']:.4f}s ({_rel(r['burn_time'], ref['burn_time']):+.4%}) "
            f"t_end={r['t_end']:.4f}s ({_rel(r['t_end'], ref['t_end']):+.4%}) "
            f"Ppk={r['P_peak'] / 1e6:.3f}MPa ({_rel(r['P_peak'], ref['P_peak']):+.4%})")
    with capsys.disabled():
        print("\n" + "\n".join(lines))
    assert len(lines) == 1 + len(_MODES)
