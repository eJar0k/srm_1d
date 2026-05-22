import numpy as np
import pytest

from srm_1d.tools.sensitivity import (
    pressure_trace_metrics,
    segmented_pressure_fitness,
    mse_fitness,
    _windowed_peak_time,
    _peak_alignment_offset,
)


def test_pressure_trace_metrics_reports_segments_and_pyrogen_summary():
    t_exp = np.array([0.04, 0.08, 0.20, 0.80, 3.00])
    p_exp = np.array([5.0, 6.0, 5.2, 5.3, 2.0])
    result = {
        'time': np.array([0.0, 0.04, 0.08, 0.20, 0.80, 3.00]),
        'P_head': np.array([0.1, 5.1, 6.2, 5.4, 5.1, 2.4]) * 1e6,
        'summary': {
            't_burn': 3.0,
            'pyrogen_duration': 0.150,
            'pyrogen_peak_P': 8.0e6,
            'pyrogen_mass_burned': 0.004,
        },
    }

    metrics = pressure_trace_metrics(t_exp, p_exp, t_min=0.01)(result)

    assert metrics['mse_all'] > 0.0
    assert metrics['mse_spike'] == pytest.approx(0.025)
    assert metrics['mse_plateau'] == pytest.approx(0.04)
    assert metrics['P_peak_sim_MPa'] == pytest.approx(6.2)
    assert metrics['peak_error_pct'] == pytest.approx((6.2 - 6.0) / 6.0 * 100.0)
    assert metrics['pyrogen_duration_ms'] == pytest.approx(150.0)
    assert metrics['pyrogen_peak_P_MPa'] == pytest.approx(8.0)
    assert metrics['pyrogen_mass_burned_g'] == pytest.approx(4.0)


# ================================================================
# Peak-alignment hook (v0.7.1 Phase 5)
# ================================================================

def test_windowed_peak_time_finds_argmax_inside_window():
    t = np.array([0.0, 0.02, 0.05, 0.10, 0.20, 0.50])
    p = np.array([0.1, 4.0, 6.0, 5.5, 5.2, 3.0])
    # Window contains all but the first sample; peak at t=0.05.
    assert _windowed_peak_time(t, p, (0.01, 0.50)) == pytest.approx(0.05)


def test_windowed_peak_time_returns_nan_on_empty_window():
    t = np.array([0.0, 0.5, 1.0])
    p = np.array([1.0, 2.0, 0.5])
    out = _windowed_peak_time(t, p, (0.6, 0.9))
    assert not np.isfinite(out)


def test_peak_alignment_offset_matches_sim_to_exp_peak():
    # Sim peaks at t=0.08; experimental peaks at t=0.05. Offset should
    # be 0.05 - 0.08 = -0.03 (shift sim earlier by 30 ms).
    t_sim = np.array([0.0, 0.04, 0.08, 0.20, 1.0])
    p_sim = np.array([0.1, 4.0, 6.0, 5.0, 2.0])
    t_exp = np.array([0.0, 0.03, 0.05, 0.20, 1.0])
    p_exp = np.array([0.1, 5.0, 6.0, 5.0, 2.0])
    offset = _peak_alignment_offset(
        t_sim, p_sim, t_exp, p_exp, peak_align_window=(0.0, 0.30),
    )
    assert offset == pytest.approx(-0.03)


def test_peak_alignment_offset_no_window_is_zero():
    t = np.array([0.0, 0.05, 0.10])
    p = np.array([1.0, 2.0, 1.0])
    assert _peak_alignment_offset(t, p, t, p, peak_align_window=None) == 0.0


def test_mse_fitness_alignment_removes_ignition_phasing_residual():
    """A sim trace identical in shape but shifted in time should score
    near zero with peak_align_window enabled and a large MSE without it.

    The synthetic profile is entirely a function of ``t - t_peak`` so a
    pure time-shift is exactly invertible: the aligned MSE should be at
    interpolation-noise floor (~1e-3 MPa^2 at this sampling density)."""
    base_t = np.linspace(0.0, 1.0, 1001)  # 1 ms grid resolution

    def make_profile(t_peak):
        # Sharp spike at t_peak then exponential decay relative to it.
        # Entire profile is f(t - t_peak), so shifting in time is
        # exactly equivalent to changing t_peak.
        tau = base_t - t_peak
        spike = 5.0 * np.exp(-(tau / 0.03) ** 2)
        decay = 3.0 * np.exp(-tau * 0.5) * (tau >= 0)
        return spike + decay

    # Experimental peak at t=0.05, sim peak at t=0.08 (30 ms late).
    p_exp = make_profile(0.05)
    p_sim = make_profile(0.08)
    result = {
        'time': base_t.copy(),
        'P_head': p_sim * 1e6,
        'summary': {},
    }

    unaligned = mse_fitness(base_t, p_exp, t_min=0.01)(result)
    aligned = mse_fitness(
        base_t, p_exp, t_min=0.01,
        peak_align_window=(0.0, 0.20),
    )(result)

    assert unaligned > 0.05, "Unaligned MSE should be substantial"
    # With perfect shift-invariance and 1 ms sampling, the aligned MSE
    # is bounded by interpolation noise (~1e-3 MPa^2). The crucial check
    # is the order-of-magnitude reduction vs unaligned.
    assert aligned < unaligned * 0.05, (
        f"Aligned MSE ({aligned:.4f}) should be << unaligned "
        f"({unaligned:.4f}) — alignment failing to shift sim trace"
    )


def test_pressure_trace_metrics_reports_t_offset_when_alignment_enabled():
    t_exp = np.array([0.0, 0.04, 0.08, 0.20, 1.0])
    p_exp = np.array([0.1, 5.0, 4.0, 3.0, 1.0])
    result = {
        'time': np.array([0.0, 0.06, 0.10, 0.20, 1.0]),
        'P_head': np.array([0.1, 5.0, 4.0, 3.0, 1.0]) * 1e6,
        'summary': {},
    }
    # Sim peak at t=0.06, exp peak at t=0.04 -> offset should be -0.02.
    metrics = pressure_trace_metrics(
        t_exp, p_exp, t_min=0.01,
        peak_align_window=(0.0, 0.20),
    )(result)
    assert 't_offset_applied_s' in metrics
    assert metrics['t_offset_applied_s'] == pytest.approx(-0.02)


def test_pressure_trace_metrics_t_offset_zero_when_alignment_disabled():
    t_exp = np.array([0.0, 0.04, 0.08, 0.20, 1.0])
    p_exp = np.array([0.1, 5.0, 4.0, 3.0, 1.0])
    result = {
        'time': np.array([0.0, 0.06, 0.10, 0.20, 1.0]),
        'P_head': np.array([0.1, 5.0, 4.0, 3.0, 1.0]) * 1e6,
        'summary': {},
    }
    metrics = pressure_trace_metrics(t_exp, p_exp, t_min=0.01)(result)
    assert metrics['t_offset_applied_s'] == 0.0


def test_segmented_pressure_fitness_uses_segment_weights():
    t_exp = np.array([0.04, 0.08, 0.20, 0.80, 3.00])
    p_exp = np.array([5.0, 6.0, 5.2, 5.3, 2.0])
    result = {
        'time': np.array([0.0, 0.04, 0.08, 0.20, 0.80, 3.00]),
        'P_head': np.array([0.1, 5.0, 6.0, 6.2, 5.3, 2.0]) * 1e6,
        'summary': {},
    }
    weights = {
        'mse_spike': 0.0,
        'mse_post_spike': 1.0,
        'mse_plateau': 0.0,
        'mse_taildown': 0.0,
    }

    fitness = segmented_pressure_fitness(
        t_exp, p_exp, t_min=0.01, weights=weights,
    )

    assert fitness(result) == pytest.approx(1.0)
