import math

import numpy as np
import pytest

from srm_1d.solid_thermal import (
    _compute_T_surf,
    _step_goodman_ode,
    _surface_has_ignited,
)


def test_compute_T_surf_limits():
    T_initial = 293.0
    T_gas = 2400.0
    k_solid = 0.3

    assert _compute_T_surf(
        1.0e-15, 1000.0, T_gas, T_initial, k_solid
    ) == pytest.approx(T_initial, abs=1.0e-4)
    assert _compute_T_surf(
        100.0, 1000.0, T_gas, T_initial, k_solid
    ) == pytest.approx(T_gas, rel=1.0e-4)
    assert _compute_T_surf(
        1.0e-3, 0.0, T_gas, T_initial, k_solid
    ) == pytest.approx(T_initial)


def test_goodman_matches_constant_flux_early_limit():
    alpha = 1.0e-7
    k_solid = 0.3
    T_initial = 293.0
    q_flux = 1.0e6

    # Low h keeps h*delta << 3*k through 100 ms. With T_gas fixed so
    # h*(T_gas - T_initial) = q, this exercises the Goodman constant-q
    # early-time limit from equations_goodman_integral.md.
    h_c = 50.0
    T_gas = T_initial + q_flux / h_c

    delta = 1.0e-6
    T_surf = T_initial
    dt = 1.0e-4
    sample_times = (0.001, 0.010, 0.050, 0.100)
    samples = {}
    t = 0.0

    for _ in range(int(sample_times[-1] / dt)):
        delta, T_surf = _step_goodman_ode(
            delta, T_surf, h_c, T_gas, T_initial, alpha, k_solid, dt
        )
        t += dt
        rounded_t = round(t, 4)
        if rounded_t in sample_times:
            samples[rounded_t] = T_surf

    assert set(samples) == set(sample_times)
    for t_sample, actual_T in samples.items():
        expected_rise = q_flux * 2.0 * math.sqrt(alpha * t_sample / 3.0) / k_solid
        actual_rise = actual_T - T_initial
        assert actual_rise == pytest.approx(expected_rise, rel=0.03)


def test_surface_has_ignited_uses_strict_threshold():
    assert not _surface_has_ignited(849.9, 850.0)
    assert not _surface_has_ignited(850.0, 850.0)
    assert _surface_has_ignited(850.1, 850.0)


def test_goodman_stable_and_monotonic_near_small_delta():
    alpha = 1.0e-7
    k_solid = 0.3
    T_initial = 293.0
    T_gas = 2400.0
    h_c = 500.0

    delta = 1.0e-6
    T_surf = T_initial
    prev_delta = delta
    prev_T = T_surf

    for _ in range(1000):
        delta, T_surf = _step_goodman_ode(
            delta, T_surf, h_c, T_gas, T_initial, alpha, k_solid, 1.0e-5
        )
        assert np.isfinite(delta)
        assert np.isfinite(T_surf)
        assert delta > 0.0
        assert delta >= prev_delta
        assert T_surf >= prev_T
        prev_delta = delta
        prev_T = T_surf


@pytest.mark.parametrize(
    "delta,h_c,T_gas,T_initial,alpha,k_solid,dt",
    [
        (1.0e-6, 0.0, 2400.0, 293.0, 1.0e-7, 0.3, 1.0e-5),
        (1.0e-6, -1.0, 2400.0, 293.0, 1.0e-7, 0.3, 1.0e-5),
        (1.0e-6, 500.0, 2400.0, 293.0, 0.0, 0.3, 1.0e-5),
        (1.0e-6, 500.0, 2400.0, 293.0, 1.0e-7, 0.0, 1.0e-5),
        (1.0e-6, 500.0, 2400.0, 293.0, 1.0e-7, 0.3, 0.0),
        (0.0, 500.0, 2400.0, 293.0, 1.0e-7, 0.3, 0.0),
    ],
)
def test_degenerate_inputs_remain_finite(
    delta, h_c, T_gas, T_initial, alpha, k_solid, dt
):
    new_delta, new_T = _step_goodman_ode(
        delta, T_initial, h_c, T_gas, T_initial, alpha, k_solid, dt
    )
    assert np.isfinite(new_delta)
    assert np.isfinite(new_T)
    assert new_delta > 0.0
