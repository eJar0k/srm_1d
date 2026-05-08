"""Inline test motor/propellant fixtures.

Replaces the v0.5.x ``make_hasegawa_propellant_1`` /
``make_hasegawa_motor_A_geo`` / etc. factory functions, which were
deleted in v0.6.0. Tests construct propellants and geometries inline
or via ``build_snapped_geometry`` / ``run_from_ric``.

These helpers exist so the test suite has a single source of truth
for the validated Hasegawa Propellant 1 numbers without resurrecting
production-side factory functions.
"""

from pathlib import Path

from srm_1d.propellant import Propellant, PropellantTab
from srm_1d.grain_geometry import build_snapped_geometry


MOTORS_DIR = Path(__file__).resolve().parents[1] / 'motors'


def hasegawa_propellant_1():
    """RPA-validated 69AP/17HTPB/14Al composite — same numbers the old
    ``make_hasegawa_propellant_1`` factory produced."""
    return Propellant(
        name="Hasegawa Propellant 1 (69AP/17HTPB/14Al)",
        tabs=[PropellantTab(
            min_pressure=0.0, max_pressure=20e6,
            a=4.821e-5, n=0.3,
            gamma=1.19, T_flame=3041.0, molecular_weight=0.0254,
        )],
        rho_propellant=1700.0, Cps=1500.0,
        T_surface=1000.0, T_initial=293.0,
        mu_gas=8.842e-5, k_gas=0.3685, Cp_gas=2060.0,
    )


def king_propellant_4525():
    """73AP/27HTPB no-metal propellant — old ``make_king_propellant_4525``."""
    return Propellant(
        name="King 4525 (73AP(20um)/27HTPB)",
        tabs=[PropellantTab(
            min_pressure=0.0, max_pressure=20e6,
            a=0.83e-2 / 6e6**0.3, n=0.3,
            gamma=1.25, T_flame=1667.0, molecular_weight=0.025,
        )],
        rho_propellant=1500.0, Cps=1200.0,
        T_surface=750.0, T_initial=293.0,
        mu_gas=6.0e-5, k_gas=0.25, Cp_gas=1800.0,
    )


def single_cylinder_geo(D_bore=0.040, D_outer=0.080, length=0.500,
                        target_propellant_cells=50):
    """Single inhibited-end cylinder. Replaces ``make_single_cylinder``."""
    return build_snapped_geometry(
        [{
            'D_bore_fwd': D_bore,
            'length': length,
            'inhibit_fwd': True,
            'inhibit_aft': True,
        }],
        D_outer=D_outer,
        target_propellant_cells=target_propellant_cells,
    )


def example_bates_geo(target_propellant_cells=120):
    """4-segment BATES with both ends free. Replaces ``make_example_bates``."""
    seg_spec = {
        'D_bore_fwd': 0.040,
        'length': 0.120,
    }
    return build_snapped_geometry(
        [{**seg_spec, 'gap_after': 0.005} for _ in range(3)] + [seg_spec],
        D_outer=0.070,
        target_propellant_cells=target_propellant_cells,
    )


def bates_motor_geo(D_bore, D_outer, L_segment, N_segments, spacing,
                    target_propellant_cells=None):
    """Parametric BATES builder. Replaces ``make_bates_motor``."""
    if target_propellant_cells is None:
        target_propellant_cells = max(50, int(N_segments * L_segment / 0.005))
    base = {
        'D_bore_fwd': D_bore,
        'length': L_segment,
    }
    specs = [{**base, 'gap_after': spacing} for _ in range(N_segments - 1)]
    specs.append(base)
    return build_snapped_geometry(specs, D_outer,
                                  target_propellant_cells=target_propellant_cells)
