"""
test_canonical_examples.py — regression gates for canonical example
runs (v0.7.3.2+).

Purpose
-------
Catches the v0.7.3-phaseB / v0.7.3.1 silent-failure mode where the
canonical ``hasegawa_motor_a.py`` example was producing collapsed
traces (numerical-collapse termination_code = 4, t_burn ≈ 1.5 ms)
while the existing test suite passed (because the existing tests use
explicit calibrated pyrogen sizing that bypasses the Sutton-default
instability).

The new pytest tests below run the canonical CONFIG end-to-end and
assert termination_code in {0, 1, 2, 3} (anything except numerical
collapse) AND t_burn > 0.5 s (rules out the < 1 ms collapse-trip).

Background
----------
- Sutton-default `pyrogen_throat_area = 0.01·A_main` gave 9 mm² for
  Hasegawa A — too small under Phase B.0 cold-bore IC. Drove plenum
  P_ig past 600 MPa equilibrium → ignition runaway → collapse.
- v0.7.3.2 fix: Kn-based throat sizing
  (A_throat = pyrogen_burn_area / Kn_design with Kn=100) +
  cfl_target=0.3 (was 0.5) + source_cfl_factor=0.05 (was 0.10).
- See srm_1d/docs/v0_7_4/references/frozen_collapse_investigation.md
  for the regression timeline.
"""
import pytest


def test_canonical_hasegawa_motor_a_does_not_collapse():
    """Canonical hasegawa_motor_a.py config (Sutton-default pyrogen +
    v0.7.0-calibrated knobs) must run to completion without tripping
    numerical collapse. Was silently broken from v0.7.3-phaseB through
    v0.7.3.1 by Phase B.0's cold-bore IC; fixed in v0.7.3.2 via
    Kn-based throat sizing + tighter CFL defaults.
    """
    pytest.importorskip("numba")
    from srm_1d.openmotor_adapter import run_from_ric

    result, _perf, _nz, _geo, _prop = run_from_ric(
        'motors/hasegawa_a.ric',
        roughness=37.1e-6,
        kappa=0.45,
        pyrogen='bpnv',
        pyrogen_mass=None,         # Sutton sizing
        T_ignition=850.0,
        P_cutoff=0.05e6,
        snapshot_interval=2.0,
        print_interval=20.0,
        verbose=False,
    )

    summary = result['summary']
    term_code = summary['termination_code']
    t_burn = float(summary['t_burn'])
    p_peak = float(summary.get('P_peak', 0.0)) / 1e6

    # Termination code 4 = numerical collapse (the failure mode this
    # gates). Codes 0/1/2/3 are all clean.
    assert term_code != 4, (
        f"Canonical Hasegawa A collapsed (termination_code=4); "
        f"t_burn={t_burn:.4f}s, P_peak={p_peak:.2f} MPa. "
        f"v0.7.3.2 architecture fix likely regressed."
    )

    # t_burn > 0.5s confirms the sim actually ignited and burned for a
    # meaningful duration. The < 2 ms collapse case would fail here too.
    assert t_burn > 0.5, (
        f"Canonical Hasegawa A produced suspiciously short burn "
        f"(t_burn={t_burn:.4f}s); P_peak={p_peak:.2f} MPa. "
        f"Likely silent ignition failure even without collapse trip."
    )

    # P_peak should be on the order of the Hasegawa A experimental
    # plateau (~6 MPa). Sanity range: 2 < P_peak < 30 MPa.
    assert 2.0 < p_peak < 30.0, (
        f"Canonical Hasegawa A P_peak={p_peak:.2f} MPa outside "
        f"the (2, 30) MPa sanity band — likely degenerate."
    )


def test_canonical_chunc_does_not_collapse():
    """Canonical Chunc (machbusterNew, head_basket + mtv radiation, frozen
    transport) must run to completion without numerical collapse. Chunc is
    the v0.7.4 ignition-spike diagnostic; this gates that the default
    (flame-front + Z-N OFF) path stays healthy as those features land.
    """
    pytest.importorskip("numba")
    from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen

    pyro = load_pyrogen('mtv')
    pyro.heat_delivery_mode = 'radiation'
    result, _perf, _nz, _geo, _prop = run_from_ric(
        'motors/machbusterNew.ric',
        transport_path='motors/machbusterNew.frozen.transport.yaml',
        pyrogen=pyro,
        injection_topology='head_basket',
        roughness=37.1e-6, kappa=0.45, T_ignition=850.0,
        P_cutoff=0.05e6, cfl_target=0.3, t_max=3.0,
        snapshot_interval=0.5, print_interval=20.0, verbose=False,
    )

    summary = result['summary']
    term_code = summary['termination_code']
    t_burn = float(summary['t_burn'])
    p_peak = float(summary.get('P_peak', 0.0)) / 1e6

    assert term_code != 4, (
        f"Canonical Chunc collapsed (termination_code=4); "
        f"t_burn={t_burn:.4f}s, P_peak={p_peak:.2f} MPa."
    )
    assert t_burn > 0.5, (
        f"Canonical Chunc produced suspiciously short burn "
        f"(t_burn={t_burn:.4f}s); P_peak={p_peak:.2f} MPa."
    )
    # Chunc plateau is ~8.8 MPa; the (uncorrected) baseline spike reaches
    # ~16-17 MPa, so a wide sanity band catches only collapse/degeneracy.
    assert 4.0 < p_peak < 40.0, (
        f"Canonical Chunc P_peak={p_peak:.2f} MPa outside the (4, 40) MPa "
        f"sanity band — likely degenerate."
    )


def test_verify_run_health_passes_clean_runs():
    """verify_run_health should return True for a normally-terminated
    run with reasonable t_burn."""
    from srm_1d.run_artifacts import verify_run_health
    fake_result = {
        'summary': {
            'termination_code': 1,  # complete burnout
            't_burn': 4.16,
            'P_peak': 6.2e6,
        }
    }
    assert verify_run_health(fake_result, motor_name='test_fake') is True


def test_verify_run_health_flags_collapsed_runs():
    """verify_run_health should return False (and not raise by default)
    for a collapsed run."""
    from srm_1d.run_artifacts import verify_run_health
    fake_result = {
        'summary': {
            'termination_code': 4,  # numerical collapse
            't_burn': 0.0015,
            'P_peak': 1.44e6,
        }
    }
    assert verify_run_health(fake_result, motor_name='test_collapsed') is False


def test_verify_run_health_raises_when_asked():
    """verify_run_health with raise_on_fail=True should raise on a
    collapsed run — useful for CI / pytest contexts."""
    from srm_1d.run_artifacts import verify_run_health
    fake_result = {
        'summary': {
            'termination_code': 4,
            't_burn': 0.0015,
            'P_peak': 1.44e6,
        }
    }
    with pytest.raises(RuntimeError, match="health check failed"):
        verify_run_health(fake_result, raise_on_fail=True)


def test_verify_run_health_flags_too_short_burn():
    """A run that terminated cleanly but had too-short t_burn (e.g.,
    P_cutoff tripped immediately) is still flagged as unhealthy."""
    from srm_1d.run_artifacts import verify_run_health
    fake_result = {
        'summary': {
            'termination_code': 2,  # pressure cutoff (clean)
            't_burn': 0.05,         # but suspiciously short
            'P_peak': 0.1e6,
        }
    }
    assert verify_run_health(fake_result, min_t_burn_s=0.1) is False
