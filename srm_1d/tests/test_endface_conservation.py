"""End-face mass-injection conservation (v0.6.0 gating test).

The hat-function distribution kernel introduced in v0.6.0 distributes a
single end-face's mass across two adjacent cells using weights that form
a partition of unity. This test verifies that the total instantaneous
end-face mass rate matches the analytic value `2 * rho * r * A_face`
for both:

  - **Snapped grids** (segment ends land exactly on cell edges)
  - **Unsnapped grids** (segment ends fall inside cells)

If snapped passes but unsnapped fails (or vice versa), the kernel is
correct only for one regime — that's a real bug. If both fail, the
boundary clamp logic at update_cell_geometry's top-of-loop branch is
likely culprit.
"""

import numpy as np
import pytest

from srm_1d.grain_geometry import (
    build_snapped_geometry, GrainSegment, MotorGeometry,
    update_cell_geometry,
)
from srm_1d.propellant import Propellant, PropellantTab


def _make_test_propellant():
    """Hasegawa-Propellant-1-equivalent constructed inline (no factory dep)."""
    return Propellant(
        name="endface-conservation-test",
        tabs=[PropellantTab(
            min_pressure=0.0, max_pressure=20e6,
            a=4.821e-5, n=0.3,
            gamma=1.19, T_flame=3041.0, molecular_weight=0.0254,
        )],
        rho_propellant=1700.0, Cps=1500.0,
        T_surface=1000.0, T_initial=293.0,
        mu_gas=8.842e-5, k_gas=0.3685, Cp_gas=2060.0,
    )


def _endface_mass_rate(geo, prop, P_ref=3e6):
    """Run a single update_cell_geometry call; return Σ endface_msource·dx [kg/s]."""
    ga = geo.compile_geometry_arrays()
    N, dx = geo.N_cells, geo.dx
    D_port = ga['D_port'].copy()
    regress = ga['regress'].copy()
    A_port = np.zeros(N)
    C_burn = np.zeros(N)
    D_hyd = np.zeros(N)
    is_grain = np.zeros(N, dtype=np.bool_)
    endface_msource = np.zeros(N)
    P = np.full(N, P_ref)

    update_cell_geometry(
        regress, D_port, ga['x_centers'], dx, N, ga['N_seg'], ga['D_outer'],
        ga['seg_x_start'], ga['seg_length'],
        ga['seg_fwd_regression'], ga['seg_aft_regression'],
        ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'],
        ga['cell_segment_id'], P, prop.rho_propellant,
        *prop.tab_arrays(), len(prop.tabs),
        A_port, C_burn, D_hyd, is_grain, endface_msource,
        ga['cell_D_bore_init'], ga['cell_wall_web'],
        ga['cell_segment_type'], ga['cell_fmm_idx'],
        ga['fmm_offset'], ga['fmm_reg_flat'],
        ga['fmm_perim_flat'], ga['fmm_port_flat'],
    )
    return float(np.sum(endface_msource) * dx)


def _analytical_mass_rate(D_outer, D_bore, prop, P_ref=3e6):
    """Both faces free; each face area = (π/4)(D_outer² − D_bore²)."""
    A_face = np.pi / 4.0 * (D_outer**2 - D_bore**2)
    r_normal = prop.burn_rate_normal(P_ref)
    return 2.0 * prop.rho_propellant * r_normal * A_face


D_OUTER = 0.080
D_BORE = 0.040
L_SEG = 0.500
TOLERANCE_PCT = 0.1


@pytest.mark.parametrize("target_cells", [50, 100, 200, 500])
def test_endface_snapped(target_cells):
    """Snapped grid: forward face lands on a cell edge by construction.

    Old interval-containment kernel double-counts here because both
    adjacent cells satisfy `x_lo <= x_face <= x_hi`. New hat-function
    kernel must split mass with weights summing to 1.
    """
    prop = _make_test_propellant()
    geo = build_snapped_geometry(
        [{'D_bore_fwd': D_BORE, 'length': L_SEG}],
        D_outer=D_OUTER,
        target_propellant_cells=target_cells,
    )
    sim = _endface_mass_rate(geo, prop)
    expected = _analytical_mass_rate(D_OUTER, D_BORE, prop)
    error_pct = abs(sim - expected) / expected * 100
    assert error_pct < TOLERANCE_PCT, (
        f"Snapped (target_cells={target_cells}, N_cells={geo.N_cells}): "
        f"error {error_pct:.4f}% (sim={sim:.6e} kg/s, expected={expected:.6e} kg/s)"
    )


@pytest.mark.parametrize("N_cells", [51, 103, 199, 503])
def test_endface_unsnapped(N_cells):
    """Unsnapped grid: face position deliberately non-aligned with cell boundaries.

    Selected N_cells values produce non-integer L_seg/dx ratios so the
    hat-function kernel exercises the typical (sub-cell) split path
    rather than the cell-edge degenerate case.
    """
    prop = _make_test_propellant()
    leading_gap = 0.020
    trailing_gap = 0.020
    L_motor = leading_gap + L_SEG + trailing_gap
    seg = GrainSegment(
        x_start=leading_gap, length=L_SEG,
        D_bore_fwd=D_BORE, D_outer=D_OUTER,
    )
    geo = MotorGeometry(
        L_motor=L_motor, D_outer=D_OUTER,
        segments=[seg], N_cells=N_cells,
    )
    sim = _endface_mass_rate(geo, prop)
    expected = _analytical_mass_rate(D_OUTER, D_BORE, prop)
    error_pct = abs(sim - expected) / expected * 100
    assert error_pct < TOLERANCE_PCT, (
        f"Unsnapped (N_cells={N_cells}, dx={geo.dx*1000:.4f}mm): "
        f"error {error_pct:.4f}% (sim={sim:.6e} kg/s, expected={expected:.6e} kg/s)"
    )
