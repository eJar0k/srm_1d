"""
tests/test_spatial_ignition_coupling.py — v0.7.2 Phase B.3 integration gates
=============================================================================

End-to-end tests for the Phase B wiring of cumulative-G augmented h_c
into the Goodman ignition kernel. Kernel-level tests in
test_cumulative_mass_flux.py verify the math in isolation; these tests
verify that the wiring actually changes simulation behavior and that
flame_spread_enabled=False reproduces Phase A behavior.

See srm_1d/docs/v0_7_2/candidates/02_spatial_ignition_front_coupling.md.
"""
import numpy as np
import pytest


def _short_hasegawa_run(flame_spread_enabled, t_max=0.10):
    """Run Hasegawa A briefly with a flame_spread_enabled override.

    Uses v0.7.0 calibrated knobs (matches test_pyrogen_axial_distribution
    convention) so the Phase A / B test families are directly
    comparable. The propellant.flame_spread_enabled flag is set by
    intercepting the adapter's ric_to_sim_args output and mutating
    propellant in place before run_simulation runs.
    """
    pytest.importorskip("numba")
    from srm_1d.openmotor_adapter import (
        load_ric, ric_to_sim_args, build_pyrogen_chamber, load_pyrogen,
        load_transport,
    )
    from srm_1d.simulation import run_simulation

    motor = load_ric('srm_1d/motors/hasegawa_a.ric')
    gas_props = load_transport('srm_1d/motors/hasegawa_a.transport.yaml')
    args = ric_to_sim_args(motor, gas_props=gas_props)

    geo = args.pop('geo')
    propellant = args.pop('propellant')
    propellant.flame_spread_enabled = flame_spread_enabled
    nozzle = args.pop('nozzle')

    pyrogen = load_pyrogen('bpnv')
    pyrogen_chamber = build_pyrogen_chamber(
        pyrogen, geo, nozzle,
        pyrogen_mass=12.3e-3,
        pyrogen_throat_area=38.5e-6,
        pyrogen_volume=3.2e-6,
    )

    result = run_simulation(
        geo=geo, propellant=propellant, nozzle=nozzle,
        pyrogen_chamber=pyrogen_chamber,
        roughness=37.5e-6,
        kappa=0.429,
        T_ignition=927.0,
        cfl_target=0.3,
        snapshot_interval=0.01,
        t_max=t_max,
        verbose=False,
        **args,
    )
    return result


def test_flame_spread_enabled_differs_from_disabled():
    """Phase B-v2 effect gate: at flame_spread_enabled=True (the default),
    the flame-front augmentation boosts h_c at cells immediately
    downstream of a recently-ignited cell. The trace must differ from
    the disabled case (=False) by at least 0.3% at P_peak — direction
    is intentionally NOT asserted because Phase B-v2's physical effect
    in this codebase is empirically observed rather than predicted
    (the Phase B-v1 negative finding showed how the codebase's local-Re
    tracking interacts with augmentation in non-obvious ways).
    """
    result_off = _short_hasegawa_run(
        flame_spread_enabled=False, t_max=0.10
    )
    result_on = _short_hasegawa_run(
        flame_spread_enabled=True, t_max=0.10
    )

    p_off = float((result_off['P_head'] / 1e6).max())
    p_on = float((result_on['P_head'] / 1e6).max())

    rel_diff = abs(p_on - p_off) / max(p_off, 1.0)
    assert rel_diff > 0.003, (
        f"flame_spread_enabled=True/False should yield distinguishable "
        f"traces; got off={p_off:.3f} MPa, on={p_on:.3f} MPa, "
        f"rel diff={rel_diff:.3%}; expected > 0.3% (wiring may be no-op'd)"
    )


def test_flame_spread_disabled_preserves_phase_a_window():
    """Phase B regression gate: with flame_spread_enabled=False, the
    trace must land in the Phase A baseline P_peak window (~9.3 MPa
    per Phase A.4 validation runs with the same calibrated knobs).
    This proves the augmentation can be cleanly turned off for
    diagnostic A/B comparisons.
    """
    result = _short_hasegawa_run(
        flame_spread_enabled=False, t_max=0.10
    )
    p_peak = float((result['P_head'] / 1e6).max())

    assert 7.0 <= p_peak <= 11.0, (
        f"flame_spread_enabled=False P_peak={p_peak:.3f} MPa outside "
        f"Phase A baseline window [7.0, 11.0] MPa"
    )


def test_flame_spread_preserves_mass_balance():
    """Conservation gate: Phase B affects PRE-ignition h_c only — no
    impact on mass sources / sinks. Mass-balance error must still close
    within Phase 4's 2% tolerance with flame_spread_enabled=True.
    """
    result = _short_hasegawa_run(
        flame_spread_enabled=True, t_max=0.10
    )
    err = abs(result['summary']['mass_balance_error'])
    assert err < 0.02, (
        f"Mass-balance error = {err:.4f} > 2% with flame_spread_enabled=True; "
        f"Phase B should not affect conservation"
    )
