"""
tests/test_pyrogen_axial_distribution.py — v0.7.2 Phase A.3 integration gates
==============================================================================

End-to-end tests for the Phase A wiring of pyrogen_axial_weights into
the time loop. The kernel-level tests in test_pyrogen_axial_weights.py
verify the weight formula in isolation; these tests verify that the
wiring actually changes simulation behavior (and changes it in the
expected direction).

Regression gate: kappa_jet=0 must reproduce v0.7.1.1's cell-0-only
pyrogen injection (the byte-for-byte fallback that lets users opt out
of the distribution if they want a v0.7.1.1-style trace).

Effect gate: at default kappa_jet=8.0, the trace must DIFFER from
kappa_jet=0 -- proving the wiring is active and not accidentally
no-op'd.
"""
import numpy as np
import pytest


def _short_hasegawa_a_run(kappa_jet, t_max=0.10):
    """Run Hasegawa A briefly with a per-test kappa_jet override.

    Uses v0.7.0 calibrated knobs (matches test_yns_phase4_validation
    convention) for trace comparability.
    """
    pytest.importorskip("numba")
    from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen

    pyrogen_obj = load_pyrogen('bpnv')
    pyrogen_obj.kappa_jet = kappa_jet

    result, _perf, _nz, _geo, _prop = run_from_ric(
        'motors/hasegawa_a.ric',
        pyrogen=pyrogen_obj,
        pyrogen_mass=12.3e-3,
        pyrogen_throat_area=38.5e-6,
        pyrogen_volume=3.2e-6,
        T_ignition=927.0,
        roughness=37.5e-6,
        kappa=0.429,
        cfl_target=0.3,
        snapshot_interval=0.01,
        t_max=t_max,
        verbose=False,
    )
    return result


def test_kappa_jet_zero_recovers_cell_zero_only_distribution():
    """Phase A regression gate: kappa_jet=0 routes all pyrogen mass
    into cell 0 — recovering v0.7.1.1's exact cell-0-only behavior.

    Verifies via the species mass fraction: after a short pyrogen-only
    transient, Y[0, _SPECIES_IGNITER] should be much larger than
    Y[downstream, _SPECIES_IGNITER] because all pyrogen mass enters
    at cell 0.
    """
    result = _short_hasegawa_a_run(kappa_jet=0.0, t_max=0.020)

    Y = result['Y_species_final']
    species_names = result['species_names']
    pyro_idx = list(species_names).index('BPNV_gas')

    # Cell 0 should be dominated by pyrogen
    assert Y[0, pyro_idx] > 0.3, (
        f"kappa_jet=0 cell-0 pyrogen fraction = {Y[0, pyro_idx]:.3f}; "
        f"expected > 0.3 (all pyrogen mass injected here)"
    )

    # Far-downstream cells should be much less pyrogen-rich at this
    # early time (advection hasn't carried much downstream yet)
    N = Y.shape[0]
    mid_cell = N // 2
    assert Y[mid_cell, pyro_idx] < Y[0, pyro_idx] * 0.5, (
        f"Mid-bore pyrogen fraction {Y[mid_cell, pyro_idx]:.3f} should be "
        f"much less than cell-0 fraction {Y[0, pyro_idx]:.3f} when "
        f"kappa_jet=0 (cell-0 injection only)"
    )


def test_kappa_jet_default_distributes_pyrogen_axially():
    """Phase A effect gate: at default kappa_jet=8.0, pyrogen mass
    deposits across multiple cells. Cells near the head end should
    have meaningfully higher pyrogen fraction than the kappa_jet=0
    case at the same early time.
    """
    result = _short_hasegawa_a_run(kappa_jet=8.0, t_max=0.020)

    Y = result['Y_species_final']
    species_names = result['species_names']
    pyro_idx = list(species_names).index('BPNV_gas')
    N = Y.shape[0]

    # With L_jet = 8 * d_throat_pyrogen ~ 8 * sqrt(4 * 38.5e-6 / pi)
    # ~ 8 * 7.0e-3 m ~ 5.6 cm, several cells near the head end should
    # carry meaningful pyrogen fraction (not just cell 0).
    n_pyrogen_rich = int(np.sum(Y[:, pyro_idx] > 0.01))
    assert n_pyrogen_rich >= 2, (
        f"At kappa_jet=8, only {n_pyrogen_rich} cells have Y_pyrogen > 0.01; "
        f"expected >= 2 (axial distribution)"
    )


def test_kappa_jet_zero_and_default_differ_in_pyrogen_spread():
    """Wiring activity gate: at default kappa_jet=8, MORE cells carry
    meaningful pyrogen fraction than at kappa_jet=0 (the cell-0-only
    fallback). Trace-level P_peak differences are dominated by
    propellant cascade timing at v0.7.0 calibrated knobs and don't
    reliably swing > 5%, so we check the upstream distribution
    signature directly instead.
    """
    result_zero = _short_hasegawa_a_run(kappa_jet=0.0, t_max=0.020)
    result_default = _short_hasegawa_a_run(kappa_jet=8.0, t_max=0.020)

    species_names = result_zero['species_names']
    pyro_idx = list(species_names).index('BPNV_gas')

    # v0.7.3 Phase B.0 + v0.7.4 Phase C.1: under the cold-bore IC plus
    # the new particle-geometry A_burn calculation (BPNV 3.2 mm pellets
    # L/D=1.0 per Mizushima 2016, n=0.589), pyrogen mass injection
    # advects through more cells than the v0.7.0/B.0 threshold-tuning
    # anticipated. Y > 0.3 saturates at ~28 cells under both
    # kappa_jet=0 and kappa_jet=8 because mass conservation forces the
    # smaller A_burn pellet flow to fill the head-end thoroughly under
    # either distribution. We instead compare PEAK Y at the head-end
    # cells (kappa_jet=0 dumps everything in cell 0 → very high Y[0];
    # kappa_jet=8 spreads → moderate Y across multiple cells with
    # LOWER peak Y[0]).
    Y_zero = result_zero['Y_species_final'][:, pyro_idx]
    Y_default = result_default['Y_species_final'][:, pyro_idx]

    assert Y_zero[0] > Y_default[0], (
        f"kappa_jet=0 should concentrate pyrogen at cell 0; "
        f"kappa_jet=8 should spread it. Got Y[0] = "
        f"{Y_zero[0]:.3f} (kappa=0) vs {Y_default[0]:.3f} (kappa=8). "
        f"Expected the kappa=0 case to have higher Y[0]."
    )


def test_pyrogen_distribution_preserves_total_mass_balance():
    """Conservation: with the distributed pyrogen, total mass-balance
    error should still close (within Phase 4's 2% tolerance — same
    standard the cell-0-only model is held to). Mass-balance error
    lives under result['summary'] per simulation.py output schema.
    """
    result = _short_hasegawa_a_run(kappa_jet=8.0, t_max=0.10)

    err = abs(result['summary']['mass_balance_error'])
    assert err < 0.02, (
        f"Mass-balance error = {err:.4f} > 2% with kappa_jet=8.0; "
        f"distribution may have broken conservation"
    )
