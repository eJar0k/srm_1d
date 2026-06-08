"""
test_taper.py — Parametric axial-taper FMM grain support.

Skips if openMotor + scikit-fmm aren't installed/reachable.
"""
import numpy as np
import pytest


try:
    from srm_1d.fmm_grain import (
        from_ric_grain, linear_taper, taper_profile, resolve_taper,
        TaperSpec, _interpolate_props, _setup_openmotor_path,
    )
    _setup_openmotor_path()
    HAS_OPENMOTOR = True
except (ImportError, Exception):
    HAS_OPENMOTOR = False

pytestmark = pytest.mark.skipif(
    not HAS_OPENMOTOR,
    reason="openMotor checkout / scikit-fmm not available",
)

# Imported unconditionally — pure-numpy, no openMotor needed.
from srm_1d.grain_geometry import build_snapped_geometry, _nearest_station


MAP_DIM = 201  # small for fast tests


def _finocyl_props(fin_length, diameter=0.080, length=0.200):
    return {
        'diameter': diameter,
        'length': length,
        'inhibitedEnds': 'Neither',
        'coreDiameter': 0.020,
        'numFins': 6,
        'finWidth': 0.005,
        'finLength': fin_length,
        'invertedFins': False,
    }


# ================================================================
# Property interpolation + validation
# ================================================================

class TestInterpolateProps:

    def test_float_props_blend_linearly(self):
        a = _finocyl_props(0.010)
        b = _finocyl_props(0.020)
        mid = _interpolate_props(a, b, 0.5)
        assert mid['finLength'] == pytest.approx(0.015)
        # Constant float held constant.
        assert mid['coreDiameter'] == pytest.approx(0.020)
        # Endpoints reproduce inputs exactly.
        assert _interpolate_props(a, b, 0.0)['finLength'] == pytest.approx(0.010)
        assert _interpolate_props(a, b, 1.0)['finLength'] == pytest.approx(0.020)

    def test_numfins_mismatch_raises(self):
        a = _finocyl_props(0.010)
        b = _finocyl_props(0.010)
        b['numFins'] = 8
        with pytest.raises(ValueError):
            _interpolate_props(a, b, 0.5)

    def test_bool_mismatch_raises(self):
        a = _finocyl_props(0.010)
        b = _finocyl_props(0.010)
        b['invertedFins'] = True
        with pytest.raises(ValueError):
            _interpolate_props(a, b, 0.5)

    def test_key_mismatch_raises(self):
        a = _finocyl_props(0.010)
        b = _finocyl_props(0.010)
        del b['finWidth']
        with pytest.raises(ValueError):
            _interpolate_props(a, b, 0.5)


# ================================================================
# TaperSpec construction / validation
# ================================================================

class TestTaperSpec:

    def test_linear_taper_two_points_sorted(self):
        t = linear_taper('Finocyl', _finocyl_props(0.010),
                         _finocyl_props(0.020), map_dim=MAP_DIM)
        fracs = [f for f, _ in t.control_stations]
        assert fracs == [0.0, 1.0]

    def test_profile_sorts_control_points(self):
        t = taper_profile('Finocyl', [
            (1.0, _finocyl_props(0.020)),
            (0.0, _finocyl_props(0.010)),
            (0.5, _finocyl_props(0.030)),
        ], map_dim=MAP_DIM)
        fracs = [f for f, _ in t.control_stations]
        assert fracs == [0.0, 0.5, 1.0]

    def test_out_of_range_fraction_raises(self):
        with pytest.raises(ValueError):
            taper_profile('Finocyl', [(0.0, _finocyl_props(0.01)),
                                      (1.5, _finocyl_props(0.02))])


# ================================================================
# resolve_taper — real per-station FMM tables
# ================================================================

class TestResolveTaper:

    def test_endpoints_reproduce_exact_tables(self):
        fwd = _finocyl_props(0.010)
        aft = _finocyl_props(0.020)
        t = linear_taper('Finocyl', fwd, aft, map_dim=MAP_DIM)
        tables, fracs = resolve_taper(t, 5)

        assert fracs[0] == pytest.approx(0.0)
        assert fracs[-1] == pytest.approx(1.0)
        assert len(tables) == 5

        ref_fwd = from_ric_grain({'type': 'Finocyl', 'properties': fwd},
                                 map_dim=MAP_DIM)
        ref_aft = from_ric_grain({'type': 'Finocyl', 'properties': aft},
                                 map_dim=MAP_DIM)
        assert np.allclose(tables[0].perimeter, ref_fwd.perimeter)
        assert np.allclose(tables[0].port_area, ref_fwd.port_area)
        assert np.allclose(tables[-1].perimeter, ref_aft.perimeter)
        assert np.allclose(tables[-1].port_area, ref_aft.port_area)

    def test_wall_web_varies_monotonically(self):
        t = linear_taper('Finocyl', _finocyl_props(0.008),
                         _finocyl_props(0.022), map_dim=MAP_DIM)
        tables, _ = resolve_taper(t, 6)
        webs = np.array([tab.wall_web for tab in tables])
        diffs = np.diff(webs)
        # Geometry genuinely changes along the taper...
        assert webs[0] != pytest.approx(webs[-1])
        # ...and does so monotonically (one sign throughout).
        assert np.all(diffs >= -1e-9) or np.all(diffs <= 1e-9)

    def test_degenerate_taper_dedupes_solves(self):
        same = _finocyl_props(0.015)
        t = linear_taper('Finocyl', same, dict(same), map_dim=MAP_DIM)
        tables, _ = resolve_taper(t, 4)
        # Identical cross-sections share one cached table object.
        assert all(tab is tables[0] for tab in tables)


# ================================================================
# nearest-station mapping
# ================================================================

class TestNearestStation:

    def test_picks_closest(self):
        fr = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        assert _nearest_station(0.0, fr) == 0
        assert _nearest_station(0.6, fr) == 2
        assert _nearest_station(0.7, fr) == 3
        assert _nearest_station(1.0, fr) == 4


# ================================================================
# Mesh-dependent station count via build_snapped_geometry
# ================================================================

class TestMeshDependentResolution:

    def _spec(self, max_stations, length=0.200):
        t = linear_taper('Finocyl', _finocyl_props(0.010, length=length),
                         _finocyl_props(0.020, length=length),
                         map_dim=MAP_DIM, max_stations=max_stations)
        return [{'length': length, 'taper': t}]

    def test_fine_mesh_caps_at_max_stations(self):
        # Many cells -> station count clamped to max_stations.
        geo = build_snapped_geometry(self._spec(max_stations=8),
                                     D_outer=0.080,
                                     target_propellant_cells=100)
        assert len(geo.segments[0].fmm_tables) == 8

    def test_coarse_mesh_tracks_cell_count(self):
        # Few cells -> fewer stations than the cap (tracks the mesh).
        geo = build_snapped_geometry(self._spec(max_stations=32),
                                     D_outer=0.080,
                                     target_propellant_cells=6)
        n_tables = len(geo.segments[0].fmm_tables)
        assert 2 <= n_tables <= 8  # ~6 cells, never the 32 cap

    def test_taper_and_fmm_table_mutually_exclusive(self):
        t = linear_taper('Finocyl', _finocyl_props(0.010),
                         _finocyl_props(0.020), map_dim=MAP_DIM)
        ref = from_ric_grain({'type': 'Finocyl',
                              'properties': _finocyl_props(0.010)},
                             map_dim=MAP_DIM)
        with pytest.raises(ValueError):
            build_snapped_geometry(
                [{'length': 0.2, 'taper': t, 'fmm_table': ref}],
                D_outer=0.080, target_propellant_cells=50)


# ================================================================
# compile_geometry_arrays — per-cell station assignment
# ================================================================

class TestCompileTaper:

    def _tapered_geo(self):
        t = linear_taper('Finocyl', _finocyl_props(0.008),
                         _finocyl_props(0.022), map_dim=MAP_DIM,
                         max_stations=16)
        return build_snapped_geometry([{'length': 0.200, 'taper': t}],
                                      D_outer=0.080,
                                      target_propellant_cells=40)

    def test_cells_span_multiple_tables(self):
        geo = self._tapered_geo()
        ga = geo.compile_geometry_arrays()
        grain = ga['cell_segment_id'] >= 0
        used = set(ga['cell_fmm_idx'][grain].tolist())
        used.discard(-1)
        assert len(used) > 1  # cells point at different station tables

    def test_wall_web_varies_along_axis(self):
        geo = self._tapered_geo()
        ga = geo.compile_geometry_arrays()
        grain = ga['cell_segment_id'] >= 0
        webs = ga['cell_wall_web'][grain]
        assert webs.std() > 0.0  # axially-varying burnout depth
        assert np.all(ga['cell_segment_type'][grain] == 1)


# ================================================================
# total_propellant_volume — per-cell Riemann sum (any taper shape)
# ================================================================

class TestTaperVolume:

    def test_volume_equals_per_cell_riemann_sum(self):
        # The returned volume must equal the exact cell->station Riemann sum
        # the solver burns (valid for arbitrary, incl. nonlinear, tapers).
        t = taper_profile('Finocyl', [
            (0.0, _finocyl_props(0.008)),
            (0.5, _finocyl_props(0.020)),   # nonlinear: peak in the middle
            (1.0, _finocyl_props(0.012)),
        ], map_dim=MAP_DIM, max_stations=12)
        geo = build_snapped_geometry([{'length': 0.200, 'taper': t}],
                                     D_outer=0.080,
                                     target_propellant_cells=40)
        V = geo.total_propellant_volume()

        ga = geo.compile_geometry_arrays()
        casting = np.pi / 4.0 * ga['D_outer'] ** 2
        grain = ga['cell_segment_id'] >= 0
        V_expected = float(
            np.sum(casting - ga['cell_A_port_init'][grain]) * ga['dx']
        )
        assert V == pytest.approx(V_expected, rel=1e-9)
        assert V > 0.0

    def test_degenerate_taper_matches_uniform(self):
        # fwd == aft taper reduces to a uniform FMM segment's volume.
        props = _finocyl_props(0.015)
        t = linear_taper('Finocyl', props, dict(props), map_dim=MAP_DIM)
        geo = build_snapped_geometry([{'length': 0.200, 'taper': t}],
                                     D_outer=0.080,
                                     target_propellant_cells=40)
        V = geo.total_propellant_volume()

        tab = geo.segments[0].fmm_tables[0]
        casting = np.pi / 4.0 * 0.080 ** 2
        seg_len = geo.segments[0].length  # snapped
        V_uniform = (casting - tab.initial_port_area) * seg_len
        assert V == pytest.approx(V_uniform, rel=1e-9)
