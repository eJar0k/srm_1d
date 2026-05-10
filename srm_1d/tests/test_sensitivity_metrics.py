import numpy as np
import pytest

from srm_1d.tools.sensitivity import (
    pressure_trace_metrics,
    segmented_pressure_fitness,
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
