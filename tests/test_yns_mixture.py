"""
v0.7.1 Phase 2c — N-species mixing-rule tests.

Tests ``_compute_mixture_cell`` and ``_refresh_mixture_arrays``:
- pure-species limits (Y = one-hot)
- 50/50 binary mix matches hand calculation
- thermo consistency: gamma = Cp / (Cp - R)
- runs through a full Hasegawa A short sim; checks the returned
  ``gamma_mix_final / Cp_mix_final / R_mix_final / M_mix_final`` arrays
  for physical bounds.
"""
import numpy as np
import pytest

from srm_1d.simulation import (
    _compute_mixture_cell, _refresh_mixture_arrays,
    _compute_T_ceiling_arr,
)
from srm_1d.propellant import (
    GasSpecies, species_array, ambient_air_species, Pyrogen,
    R_UNIVERSAL,
)


def _bpnv_species():
    """Synthetic BPNV-like pyrogen species (matches load_pyrogen('bpnv'))."""
    p = Pyrogen(name='BPNV', a=0.0, n=0.0, rho=1800.0,
                T_flame=2800.0, M=0.030, gamma=1.25)
    return p.species  # Cp_gas derives to gamma*R/(gamma-1)*1/M ~ 1386


def _hasegawa_species():
    """Synthetic Hasegawa Prop 1 species (matches the real motor params)."""
    return GasSpecies(
        name='hasegawa_prop1_gas',
        gamma=1.19, Cp=2060.0, molecular_weight=0.0254, T_flame=3041.0,
    )


# ================================================================
# Pure-species limits
# ================================================================

def test_mixture_pure_pyrogen_limit():
    """Y = (1, 0, 0) returns species-0 thermo exactly."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    arr = species_array([sp0, sp1, sp2])

    Y = np.array([1.0, 0.0, 0.0])
    gamma, Cp, R, M = _compute_mixture_cell(Y, arr)

    assert Cp == pytest.approx(sp0.Cp)
    assert M == pytest.approx(sp0.molecular_weight)
    assert R == pytest.approx(R_UNIVERSAL / sp0.molecular_weight)
    # gamma derives: gamma = Cp / (Cp - R). For ideal gas built from gamma_0,
    # this round-trips back to the original gamma to machine precision.
    gamma_expected = sp0.Cp / (sp0.Cp - R_UNIVERSAL / sp0.molecular_weight)
    assert gamma == pytest.approx(gamma_expected, rel=1.0e-10)


def test_mixture_pure_propellant_limit():
    """Y = (0, 1, 0) returns species-1 thermo exactly."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    arr = species_array([sp0, sp1, sp2])

    Y = np.array([0.0, 1.0, 0.0])
    gamma, Cp, R, M = _compute_mixture_cell(Y, arr)

    assert Cp == pytest.approx(sp1.Cp)
    assert M == pytest.approx(sp1.molecular_weight)


def test_mixture_pure_ambient_limit():
    """Y = (0, 0, 1) returns species-2 (air) thermo exactly."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    arr = species_array([sp0, sp1, sp2])

    Y = np.array([0.0, 0.0, 1.0])
    gamma, Cp, R, M = _compute_mixture_cell(Y, arr)

    assert Cp == pytest.approx(1005.0)
    assert M == pytest.approx(0.02897)


# ================================================================
# Binary mixing
# ================================================================

def test_mixture_50_50_pyrogen_propellant():
    """Hand-calc check on a 50/50 BPNV/Hasegawa mix."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    arr = species_array([sp0, sp1, sp2])

    Y = np.array([0.5, 0.5, 0.0])
    gamma, Cp, R, M = _compute_mixture_cell(Y, arr)

    expected_Cp = 0.5 * sp0.Cp + 0.5 * sp1.Cp
    expected_inv_M = 0.5 / sp0.molecular_weight + 0.5 / sp1.molecular_weight
    expected_M = 1.0 / expected_inv_M
    expected_R = R_UNIVERSAL / expected_M
    expected_gamma = expected_Cp / (expected_Cp - expected_R)

    assert Cp == pytest.approx(expected_Cp)
    assert M == pytest.approx(expected_M)
    assert R == pytest.approx(expected_R)
    assert gamma == pytest.approx(expected_gamma)


def test_mixture_consistency_gamma_from_Cp_R():
    """For any non-degenerate mix, gamma = Cp / (Cp - R) must hold."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    arr = species_array([sp0, sp1, sp2])

    for Y in [
        np.array([0.3, 0.5, 0.2]),
        np.array([0.1, 0.1, 0.8]),
        np.array([0.7, 0.0, 0.3]),
        np.array([0.99, 0.005, 0.005]),
    ]:
        gamma, Cp, R, M = _compute_mixture_cell(Y, arr)
        assert gamma == pytest.approx(Cp / (Cp - R), rel=1.0e-12)


def test_mixture_degenerate_zero_Y_returns_species_0():
    """Y all zero (pathological) falls back to species 0 thermo."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    arr = species_array([sp0, sp1, sp2])

    Y = np.zeros(3)
    gamma, Cp, R, M = _compute_mixture_cell(Y, arr)

    assert Cp == pytest.approx(sp0.Cp)
    assert M == pytest.approx(sp0.molecular_weight)


# ================================================================
# Per-cell array refresh
# ================================================================

def test_refresh_mixture_arrays_matches_cell_calls():
    """``_refresh_mixture_arrays`` is just _compute_mixture_cell looped
    over cells; per-cell results must match a direct loop."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    arr = species_array([sp0, sp1, sp2])
    N = 5
    Y = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.5, 0.0],
        [0.3, 0.4, 0.3],
    ])

    gamma_arr = np.empty(N)
    Cp_arr = np.empty(N)
    R_arr = np.empty(N)
    M_arr = np.empty(N)
    _refresh_mixture_arrays(Y, arr, gamma_arr, Cp_arr, R_arr, M_arr, N)

    for i in range(N):
        g, Cp, R, M = _compute_mixture_cell(Y[i, :], arr)
        assert gamma_arr[i] == pytest.approx(g)
        assert Cp_arr[i] == pytest.approx(Cp)
        assert R_arr[i] == pytest.approx(R)
        assert M_arr[i] == pytest.approx(M)


# ================================================================
# Integration: Hasegawa A short run
# ================================================================

def test_hasegawa_a_mixture_arrays_in_result():
    """A short Hasegawa A run must expose mixture arrays in the result;
    they must have sensible physical bounds across all cells."""
    pytest.importorskip("numba")
    from srm_1d.openmotor_adapter import run_from_ric
    result, _perf, _nz, _geo, _prop = run_from_ric(
        'motors/hasegawa_a.ric',
        pyrogen='bpnv',
        pyrogen_mass=12.3e-3,
        pyrogen_throat_area=38.5e-6,
        pyrogen_volume=3.2e-6,
        T_ignition=927.0,
        roughness=37.5e-6,
        kappa=0.429,
        t_max=0.05,
        cfl_target=0.3,
        snapshot_interval=0.01,
        verbose=False,
    )

    gamma_arr = result['gamma_mix_final']
    Cp_arr = result['Cp_mix_final']
    R_arr = result['R_mix_final']
    M_arr = result['M_mix_final']

    # All arrays must be same length and finite
    assert gamma_arr.shape == Cp_arr.shape == R_arr.shape == M_arr.shape
    assert np.all(np.isfinite(gamma_arr))
    assert np.all(np.isfinite(Cp_arr))

    # Physical bounds: gamma in (1, 2), Cp > 800 J/(kg·K) for combustion
    # products + air, M in (0.020, 0.040) for our 3-species range.
    assert np.all(gamma_arr > 1.0)
    assert np.all(gamma_arr < 2.0)
    assert np.all(Cp_arr > 800.0)
    assert np.all(Cp_arr < 3000.0)
    assert np.all(M_arr > 0.015)
    assert np.all(M_arr < 0.050)

    # gamma = Cp / (Cp - R) must hold cell-by-cell (machine epsilon)
    for i in range(len(gamma_arr)):
        expected = Cp_arr[i] / (Cp_arr[i] - R_arr[i])
        assert gamma_arr[i] == pytest.approx(expected, rel=1.0e-10)


# ================================================================
# Strict T_ceiling kernel (DESIGN §5 with IC guard)
# ================================================================

def _three_species_arr():
    """Standard 3-species fixture matching the Hasegawa A registry."""
    sp0 = _bpnv_species()
    sp1 = _hasegawa_species()
    sp2 = ambient_air_species(298.15)
    return species_array([sp0, sp1, sp2]), sp0, sp1, sp2


def test_t_ceiling_strict_pure_pyrogen_cell_bounds_by_pyrogen_tflame():
    """A cell at Y_pyrogen = 1 caps at T_flame_pyrogen · 1.01 — the
    strict-form benefit vs the prior relaxed (max-over-all-species)
    ceiling that would have used T_flame_propellant · 1.01."""
    arr, sp0, sp1, sp2 = _three_species_arr()
    Y = np.array([[1.0, 0.0, 0.0]])
    T_ceiling = np.empty(1)
    _compute_T_ceiling_arr(Y, arr, T_ceiling, 1, T_initial_gas=300.0)
    assert T_ceiling[0] == pytest.approx(sp0.T_flame * 1.01)
    # And tighter than the propellant-dominated ceiling
    assert T_ceiling[0] < sp1.T_flame * 1.01


def test_t_ceiling_strict_pure_propellant_cell_bounds_by_propellant_tflame():
    """A cell at Y_propellant = 1 caps at T_flame_propellant · 1.01."""
    arr, _sp0, sp1, _sp2 = _three_species_arr()
    Y = np.array([[0.0, 1.0, 0.0]])
    T_ceiling = np.empty(1)
    _compute_T_ceiling_arr(Y, arr, T_ceiling, 1, T_initial_gas=300.0)
    assert T_ceiling[0] == pytest.approx(sp1.T_flame * 1.01)


def test_t_ceiling_strict_pyrogen_plus_propellant_takes_max():
    """A 50/50 mix of pyrogen + propellant uses max(T_flame_pyrogen,
    T_flame_propellant) · 1.01 — same as v0.7.0's scalar ceiling when
    propellant is the hottest species in the cell."""
    arr, _sp0, sp1, _sp2 = _three_species_arr()
    Y = np.array([[0.5, 0.5, 0.0]])
    T_ceiling = np.empty(1)
    _compute_T_ceiling_arr(Y, arr, T_ceiling, 1, T_initial_gas=300.0)
    assert T_ceiling[0] == pytest.approx(sp1.T_flame * 1.01)


def test_t_ceiling_strict_filters_below_Y_min():
    """A cell with Y_propellant = 0.04 (below Y_min=0.05) and Y_pyrogen
    = 0.96 should cap at T_flame_pyrogen · 1.01, ignoring the
    sub-threshold propellant fraction."""
    arr, sp0, _sp1, _sp2 = _three_species_arr()
    Y = np.array([[0.96, 0.04, 0.0]])
    T_ceiling = np.empty(1)
    _compute_T_ceiling_arr(Y, arr, T_ceiling, 1, T_initial_gas=300.0)
    assert T_ceiling[0] == pytest.approx(sp0.T_flame * 1.01)


def test_t_ceiling_ic_guard_activates_for_pure_ambient_with_hot_seed():
    """Pre-fill IC: Y = 100% ambient, T_initial_gas = T_flame_prop. The
    strict §5 alone would cap at T_ambient · 1.01 (~300 K) and clip
    the bore gas to ambient on step 0. The IC guard raises the ceiling
    to T_initial_gas · 1.01 — preserving the documented v0.7.0 IC."""
    arr, _sp0, sp1, _sp2 = _three_species_arr()
    Y = np.array([[0.0, 0.0, 1.0]])
    T_ceiling = np.empty(1)
    T_initial_gas = sp1.T_flame  # the v0.7.1 default IC
    _compute_T_ceiling_arr(Y, arr, T_ceiling, 1, T_initial_gas=T_initial_gas)
    # The strict per-species max would be T_ambient · 1.01; the IC guard
    # promotes it to T_initial_gas · 1.01 = T_flame_prop · 1.01.
    assert T_ceiling[0] == pytest.approx(T_initial_gas * 1.01)


def test_t_ceiling_ic_guard_inactive_once_combustion_dominates():
    """Once pyrogen displaces ambient (Y_pyrogen > Y_min), the per-
    species T_flame_pyrogen · 1.01 binds even if T_initial_gas is
    cooler. Verifies the IC guard doesn't permanently inflate the
    ceiling above the active-species bound."""
    arr, sp0, _sp1, _sp2 = _three_species_arr()
    Y = np.array([[1.0, 0.0, 0.0]])  # pure pyrogen
    T_ceiling = np.empty(1)
    # IC seed cooler than pyrogen T_flame — typical of v0.7.0 IC where
    # T_initial_gas is the propellant T_flame (3041 K) but pyrogen
    # T_flame (sp0) is 2800 K. Wait, pyrogen 2800 < propellant 3041,
    # so the IC guard at T_initial_gas would actually be HIGHER than
    # the per-species pyrogen ceiling. Use a low IC instead to test
    # the inactive-guard case.
    T_initial_gas = 500.0  # well below pyrogen T_flame
    _compute_T_ceiling_arr(Y, arr, T_ceiling, 1, T_initial_gas=T_initial_gas)
    assert T_ceiling[0] == pytest.approx(sp0.T_flame * 1.01)


def test_t_ceiling_strict_array_length():
    """_compute_T_ceiling_arr fills all N cells (no off-by-one)."""
    arr, sp0, sp1, sp2 = _three_species_arr()
    N = 4
    Y = np.array([
        [1.0, 0.0, 0.0],  # pyrogen
        [0.0, 1.0, 0.0],  # propellant
        [0.0, 0.0, 1.0],  # ambient
        [0.3, 0.4, 0.3],  # mixed (all above Y_min)
    ])
    T_ceiling = np.empty(N)
    _compute_T_ceiling_arr(Y, arr, T_ceiling, N, T_initial_gas=300.0)
    assert T_ceiling[0] == pytest.approx(sp0.T_flame * 1.01)
    assert T_ceiling[1] == pytest.approx(sp1.T_flame * 1.01)
    # ambient cell: per-species max is T_ambient · 1.01 ≈ 301; IC guard
    # at T_initial_gas=300 raises to 303. Take the larger of the two.
    expected_ambient = max(sp2.T_flame * 1.01, 300.0 * 1.01)
    assert T_ceiling[2] == pytest.approx(expected_ambient)
    # mixed cell: max over pyrogen, propellant, ambient → propellant.
    assert T_ceiling[3] == pytest.approx(sp1.T_flame * 1.01)
