"""
test_fmm.py — Tests for FMM grain support.

Skips if openMotor + scikit-fmm aren't installed/reachable.
"""
import numpy as np
import pytest


# Skip the whole module if openMotor isn't reachable.
try:
    from srm_1d.fmm_grain import (
        FmmTable, from_ric_grain, fmm_table_lookup, _setup_openmotor_path,
    )
    _setup_openmotor_path()
    HAS_OPENMOTOR = True
except (ImportError, Exception):
    HAS_OPENMOTOR = False

pytestmark = pytest.mark.skipif(
    not HAS_OPENMOTOR,
    reason="openMotor checkout / scikit-fmm not available",
)


# Use a smaller mapDim for fast tests (still accurate enough).
TEST_MAP_DIM = 301


# ================================================================
# FmmTable extraction
# ================================================================

class TestFmmTableExtraction:

    def test_finocyl_basic(self):
        """A 6-fin finocyl extracts a sensible table."""
        ric = {
            'type': 'Finocyl',
            'properties': {
                'diameter':     0.080,
                'length':       0.200,
                'inhibitedEnds':'Neither',
                'coreDiameter': 0.020,
                'numFins':      6,
                'finWidth':     0.005,
                'finLength':    0.015,
                'invertedFins': False,
            },
        }
        table = from_ric_grain(ric, map_dim=TEST_MAP_DIM)

        assert table.geom_name == 'Finocyl'
        assert table.grain_outer_diameter == pytest.approx(0.080)
        assert table.grain_length == pytest.approx(0.200)
        assert table.inhibited_fwd is False
        assert table.inhibited_aft is False

        # Wall web should be roughly the casing radius minus the core
        # radius (10mm), once you account for fins extending to r=25mm.
        # FMM gives ~20mm regardless of fin count.
        assert 0.005 < table.wall_web < 0.040

        # Initial perimeter: bigger than just the core circle (62.8mm).
        # Initial port area: bigger than just the core (314mm²).
        assert table.initial_perimeter > np.pi * 0.020
        assert table.initial_port_area > np.pi / 4 * 0.020 ** 2

        # At burnout, port area equals casting tube area.
        casting_area = np.pi / 4 * 0.080 ** 2
        assert table.port_area[-1] == pytest.approx(casting_area, rel=1e-4)
        assert table.perimeter[-1] == 0.0

    def test_table_lookup_clamped_at_boundaries(self):
        """fmm_table_lookup clamps to first/last sample at the edges."""
        ric = {
            'type': 'Finocyl',
            'properties': {
                'diameter':     0.080, 'length': 0.100,
                'inhibitedEnds':'Neither', 'coreDiameter': 0.020,
                'numFins': 4, 'finWidth': 0.004, 'finLength': 0.010,
                'invertedFins': False,
            },
        }
        table = from_ric_grain(ric, map_dim=TEST_MAP_DIM)

        # At reg=0 → first sample
        v0 = fmm_table_lookup(0.0, table.reg_depth, table.perimeter, table.n_samples)
        assert v0 == pytest.approx(table.perimeter[0], rel=1e-12)

        # At reg=2×wall_web (well past burnout) → clamped to last sample
        vp = fmm_table_lookup(
            2 * table.wall_web, table.reg_depth, table.perimeter, table.n_samples
        )
        assert vp == pytest.approx(table.perimeter[-1], rel=1e-12)

    def test_supported_fmm_types(self):
        """All 7 FMM grain types load. Each has different properties so
        we just check they don't error out. Type names match openMotor's
        `geomName` strings exactly (they include spaces)."""
        configs = [
            ('Finocyl', {
                'diameter': 0.080, 'length': 0.100, 'inhibitedEnds': 'Neither',
                'coreDiameter': 0.020, 'numFins': 4, 'finWidth': 0.004,
                'finLength': 0.010, 'invertedFins': False,
            }),
            ('Star Grain', {
                'diameter': 0.080, 'length': 0.100, 'inhibitedEnds': 'Neither',
                'numPoints': 5, 'pointLength': 0.020, 'pointWidth': 0.008,
            }),
            ('Moon Burner', {
                'diameter': 0.080, 'length': 0.100, 'inhibitedEnds': 'Neither',
                'coreDiameter': 0.020, 'coreOffset': 0.005,
            }),
            ('C Grain', {
                'diameter': 0.080, 'length': 0.100, 'inhibitedEnds': 'Neither',
                'slotWidth': 0.005, 'slotOffset': 0.020,
            }),
            ('D Grain', {
                'diameter': 0.080, 'length': 0.100, 'inhibitedEnds': 'Neither',
                'slotOffset': 0.020,
            }),
            ('X Core', {
                'diameter': 0.080, 'length': 0.100, 'inhibitedEnds': 'Neither',
                'slotWidth': 0.005, 'slotLength': 0.025,
            }),
        ]
        for gtype, props in configs:
            ric = {'type': gtype, 'properties': props}
            table = from_ric_grain(ric, map_dim=TEST_MAP_DIM)
            assert table.geom_name == gtype, f"Wrong geom_name for {gtype}"
            assert table.wall_web > 0.0, f"{gtype} produced zero wall_web"
            assert table.initial_perimeter > 0.0, f"{gtype} produced zero perimeter"

    def test_unsupported_type_raises(self):
        ric = {'type': 'Bogus', 'properties': {}}
        with pytest.raises(ValueError, match="not a supported FMM grain"):
            from_ric_grain(ric)


# ================================================================
# End-to-end FMM simulation
# ================================================================

class TestFmmSimulation:

    def _build_finocyl_motor(self):
        """Construct a small Finocyl motor for sim tests."""
        from srm_1d.grain_geometry import MotorGeometry, GrainSegment

        ric = {
            'type': 'Finocyl',
            'properties': {
                'diameter':     0.080, 'length': 0.200,
                'inhibitedEnds':'Neither', 'coreDiameter': 0.020,
                'numFins':      6, 'finWidth': 0.005, 'finLength': 0.015,
                'invertedFins': False,
            },
        }
        table = from_ric_grain(ric, map_dim=TEST_MAP_DIM)
        spacing = 0.001
        seg = GrainSegment(
            x_start=spacing, length=0.200,
            D_bore_fwd=0.080, D_outer=0.080,
            inhibit_fwd=False, inhibit_aft=False,
            fmm_table=table,
        )
        geo = MotorGeometry(
            L_motor=0.200 + 2*spacing, D_outer=0.080,
            segments=[seg], N_cells=60,
        )
        return geo, table

    def test_total_propellant_volume_uses_fmm(self):
        """For an FMM segment, total_propellant_volume should compute
        (casting_area − initial_port_area) × length, not the cylindrical
        annulus formula (which would give zero since D_bore_fwd is
        a placeholder)."""
        geo, table = self._build_finocyl_motor()
        casting_area = np.pi / 4 * 0.080 ** 2
        expected = (casting_area - table.initial_port_area) * 0.200
        assert geo.total_propellant_volume() == pytest.approx(expected, rel=1e-6)

    def test_compile_geometry_arrays_tags_fmm_cells(self):
        """compile_geometry_arrays should mark all segment cells with
        cell_segment_type=1 and the FMM-table index."""
        geo, _ = self._build_finocyl_motor()
        ga = geo.compile_geometry_arrays()
        seg_cells = ga['cell_segment_id'] >= 0
        assert np.all(ga['cell_segment_type'][seg_cells] == 1)
        assert np.all(ga['cell_fmm_idx'][seg_cells] == 0)
        assert ga['n_fmm_segs'] == 1
        # cell_wall_web for FMM cells should match the table's wall_web.
        assert np.all(np.isclose(
            ga['cell_wall_web'][seg_cells], geo.segments[0].fmm_table.wall_web
        ))

    def test_finocyl_end_to_end_mass_balance(self):
        """Full simulation should produce mass balance < 1% (mass produced
        through bore vs. mass through nozzle)."""
        from srm_1d import run_simulation
        from srm_1d.propellant import make_hasegawa_propellant_1
        from srm_1d.nozzle import Nozzle

        geo, _ = self._build_finocyl_motor()
        prop = make_hasegawa_propellant_1()
        nozzle = Nozzle(D_throat=0.020, D_exit=0.035, efficiency=0.95)

        result = run_simulation(
            geo, prop, nozzle,
            roughness=20e-6, P_ignition=0.05e6, ignition_ramp_tau=0.010,
            P_cutoff=0.5e6, snapshot_interval=2.0, print_interval=10.0,
            igniter_mass=0.005, t_max=10.0,
        )

        # Mass produced through bore should match mass through nozzle
        # within numerical tolerance.
        s = result['summary']
        mass_balance_err = abs(s['mass_produced'] - s['mass_nozzle']) / max(
            s['mass_produced'], 1e-6
        )
        assert mass_balance_err < 0.01, (
            f"Mass balance error {mass_balance_err*100:.2f}% exceeds 1%"
        )
        # Some plausibility checks.
        assert s['P_peak'] > 1e6, f"P_peak={s['P_peak']/1e6:.2f}MPa too low"
        assert s['t_burn'] > 1.0, f"t_burn={s['t_burn']:.2f}s too short"
