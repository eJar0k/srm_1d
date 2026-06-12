"""
test_taper.py — Parametric axial-taper FMM grain support.

Skips if openMotor + scikit-fmm aren't installed/reachable.
"""
import numpy as np
import pytest


try:
    from srm_1d.fmm_grain import (
        from_ric_grain, linear_taper, taper_profile, resolve_taper,
        taper_spec_from_props, TaperSpec, _interpolate_props,
        _setup_openmotor_path,
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


# ================================================================
# .ric taper block -> TaperSpec (cross-solver read path)
# ================================================================

def _ric_taper_grain(fin_fwd=0.010, fin_aft=0.020):
    """An openMotor .ric Finocyl grain dict carrying a bore-taper block."""
    props = _finocyl_props(fin_fwd)
    props['taper'] = {
        'enabled': True,
        'bore': {'profile': 'linear',
                 'controlStations': [{'frac': 1.0, 'props': {'finLength': fin_aft}}]},
    }
    return {'type': 'Finocyl', 'properties': props}


class TestRicTaperRead:

    def test_convert_geometry_builds_tapered_segment(self):
        from srm_1d.openmotor_adapter import convert_geometry
        geo = convert_geometry([_ric_taper_grain()],
                               target_propellant_cells=40, fmm_map_dim=MAP_DIM)
        seg = geo.segments[0]
        assert seg.is_tapered
        assert len(seg.fmm_tables) > 1
        perims = [t.initial_perimeter for t in seg.fmm_tables]
        assert perims[-1] != perims[0]            # the taper took effect

    def test_ric_read_matches_python_authored(self):
        # The .ric read path must produce the same tables as the Python API.
        from srm_1d.openmotor_adapter import convert_geometry
        geo_ric = convert_geometry([_ric_taper_grain()],
                                   target_propellant_cells=40, fmm_map_dim=MAP_DIM)
        t = linear_taper('Finocyl', _finocyl_props(0.010),
                         _finocyl_props(0.020), map_dim=MAP_DIM)
        geo_py = build_snapped_geometry([{'length': 0.200, 'taper': t}],
                                        D_outer=0.080, target_propellant_cells=40)
        ric_tables = geo_ric.segments[0].fmm_tables
        py_tables = geo_py.segments[0].fmm_tables
        assert len(ric_tables) == len(py_tables)
        for a, b in zip(ric_tables, py_tables):
            assert np.allclose(a.perimeter, b.perimeter)
            assert np.allclose(a.port_area, b.port_area)

    def test_taper_spec_from_props_helper(self):
        spec = taper_spec_from_props(
            'Finocyl', _finocyl_props(0.010),
            {'enabled': True, 'bore': {'controlStations':
                [{'frac': 1.0, 'props': {'finLength': 0.020}}]}},
            map_dim=MAP_DIM)
        assert isinstance(spec, TaperSpec)
        # frac 0 = base (0.010), frac 1 = override (0.020)
        assert spec.control_stations[0][1]['finLength'] == pytest.approx(0.010)
        assert spec.control_stations[-1][1]['finLength'] == pytest.approx(0.020)

    def test_non_fmm_taper_raises(self):
        from srm_1d.openmotor_adapter import convert_geometry
        bates = {'type': 'BATES', 'properties': {
            'diameter': 0.080, 'length': 0.200, 'coreDiameter': 0.020,
            'inhibitedEnds': 'Neither',
            'taper': {'enabled': True, 'bore': {'controlStations':
                [{'frac': 1.0, 'props': {'coreDiameter': 0.040}}]}}}}
        with pytest.raises(NotImplementedError):
            convert_geometry([bates], target_propellant_cells=40, fmm_map_dim=MAP_DIM)


# ================================================================
# OD / end taper — transient cell_D_outer (the casing tapers)
# ================================================================

_OD_AFT_CONE = [{'end': 'aft', 'length': 0.06, 'endDiameter': 0.050,
                 'profile': 'linear'}]
_OD_FWD_DOME = [{'end': 'fwd', 'length': 0.05, 'endDiameter': 0.040,
                 'profile': 'elliptical'}]


class TestTransientOdTaper:
    """OD / end taper drives a per-cell casing diameter (cell_D_outer) for both
    FMM (per-station OD-clipped tables) and analytic (BATES) grains."""

    def _fmm_od_geo(self, od_ends=_OD_AFT_CONE, inhibit_aft=True):
        t = taper_profile('Finocyl', [(0.0, _finocyl_props(0.012))],
                          map_dim=MAP_DIM, max_stations=12,
                          od_ends=od_ends, grain_length=0.200)
        return build_snapped_geometry(
            [{'length': 0.200, 'taper': t, 'od_ends': od_ends,
              'inhibit_aft': inhibit_aft}],
            D_outer=0.080, target_propellant_cells=40)

    def test_fmm_od_cone_varies_cell_d_outer(self):
        geo = self._fmm_od_geo()
        seg = geo.segments[0]
        assert seg.has_od and seg.is_tapered
        ga = geo.compile_geometry_arrays()
        grain = ga['cell_segment_id'] >= 0
        cD = ga['cell_D_outer'][grain]
        # Casing is full-OD in the uniform region and shrinks over the cone.
        assert cD.max() == pytest.approx(0.080, abs=1e-6)
        assert cD.min() < 0.080 - 1e-3
        assert cD.std() > 1e-4
        # The aft cone clips the FMM tables -> smaller wall_web toward the aft.
        ww = ga['cell_wall_web'][grain]
        assert ww[-1] < ww[0]

    def test_fmm_od_mass_is_per_cell_riemann_sum(self):
        geo = self._fmm_od_geo()
        ga = geo.compile_geometry_arrays()
        grain = ga['cell_segment_id'] >= 0
        cell_casting = np.pi / 4.0 * ga['cell_D_outer'] ** 2
        V_expected = float(
            np.sum(cell_casting[grain] - ga['cell_A_port_init'][grain]) * ga['dx']
        )
        assert geo.total_propellant_volume() == pytest.approx(V_expected, rel=1e-9)
        assert geo.total_propellant_volume() > 0.0

    def test_analytic_bates_od_dome(self):
        geo = build_snapped_geometry(
            [{'D_bore_fwd': 0.030, 'length': 0.200, 'od_ends': _OD_FWD_DOME,
              'inhibit_fwd': True}],
            D_outer=0.080, target_propellant_cells=40)
        seg = geo.segments[0]
        assert seg.has_od and not seg.is_tapered
        ga = geo.compile_geometry_arrays()
        grain = ga['cell_segment_id'] >= 0
        cD = ga['cell_D_outer'][grain]
        ww = ga['cell_wall_web'][grain]
        assert cD.min() < 0.080 - 1e-3 and cD.std() > 1e-4
        # Forward (domed) end has the thinner web; bore (port) is unchanged.
        assert ww[0] < ww[-1]
        assert np.allclose(ga['cell_A_port_init'][grain],
                           np.pi / 4.0 * ga['cell_D_bore_init'][grain] ** 2)
        # wall_web is exactly (cell_D_outer - bore)/2 for analytic cells.
        assert np.allclose(ww, (cD - ga['cell_D_bore_init'][grain]) / 2.0)

    def test_ric_od_block_converts_and_inhibits_end(self):
        from srm_1d.openmotor_adapter import convert_geometry
        props = _finocyl_props(0.012)
        props['taper'] = {'enabled': False, 'od': {'enabled': True,
            'ends': [{'end': 'aft', 'length': 0.06, 'endDiameter': 0.050,
                      'profile': 'linear'}]}}
        geo = convert_geometry([{'type': 'Finocyl', 'properties': props}],
                               target_propellant_cells=40, fmm_map_dim=MAP_DIM)
        seg = geo.segments[0]
        # OD-only on an FMM grain still forces the per-station path, and the
        # coned end is auto-inhibited (bonded to the closure).
        assert seg.has_od and seg.is_tapered and seg.inhibit_aft
        ga = geo.compile_geometry_arrays()
        cD = ga['cell_D_outer'][ga['cell_segment_id'] >= 0]
        assert cD.std() > 1e-4

    def test_bates_od_only_does_not_raise(self):
        # A *bore* taper of BATES raises; an OD-only taper is supported.
        from srm_1d.openmotor_adapter import convert_geometry
        bates = {'type': 'BATES', 'properties': {
            'diameter': 0.080, 'length': 0.200, 'coreDiameter': 0.030,
            'inhibitedEnds': 'Neither',
            'taper': {'enabled': False, 'od': {'enabled': True,
                'ends': [{'end': 'aft', 'length': 0.05, 'endDiameter': 0.050,
                          'profile': 'linear'}]}}}}
        geo = convert_geometry([bates], target_propellant_cells=40,
                               fmm_map_dim=MAP_DIM)
        seg = geo.segments[0]
        assert seg.has_od and not seg.is_tapered and seg.inhibit_aft

    def test_od_degenerate_no_od_is_flat(self):
        # A bore-only taper (no OD) keeps a flat casing == D_outer.
        t = linear_taper('Finocyl', _finocyl_props(0.012),
                         _finocyl_props(0.020), map_dim=MAP_DIM)
        geo = build_snapped_geometry([{'length': 0.200, 'taper': t}],
                                     D_outer=0.080, target_propellant_cells=40)
        ga = geo.compile_geometry_arrays()
        assert np.allclose(ga['cell_D_outer'], 0.080)
        assert not geo.segments[0].has_od

    def test_transient_od_run_mass_balance(self):
        # A full transient run on the OD-coned finocyl conserves mass.
        from srm_1d import run_simulation
        from srm_1d.nozzle import Nozzle
        from srm_1d.propellant import Pyrogen
        from srm_1d.igniter_plenum import PyrogenChamber
        from tests._motor_fixtures import hasegawa_propellant_1

        geo = self._fmm_od_geo()
        prop = hasegawa_propellant_1()
        nozzle = Nozzle(D_throat=0.020, D_exit=0.035, efficiency=0.95)
        pyro = Pyrogen("t", 3.0e-5, 0.5, 1700.0, 2800.0, 0.030, 1.25,
                       heat_flux_cal_cm2_s=69.4)
        chamber = PyrogenChamber(pyro, 0.005, 5.0e-4, 2.0e-5, 3.0e-6,
                                 "end_burning")
        res = run_simulation(geo, prop, nozzle, chamber, roughness=20e-6,
                             T_ignition=294.0, P_cutoff=0.5e6,
                             snapshot_interval=2.0, print_interval=1e9,
                             t_max=3.0)
        s = res['summary']
        mb = abs(s['mass_produced'] - s['mass_nozzle']) / max(s['mass_produced'], 1e-6)
        assert mb < 0.01, f"mass balance {mb*100:.2f}% > 1%"
        assert np.isfinite(s['P_peak']) and s['P_peak'] > 1e6
        # The result carries the per-cell casing for the slice viewer.
        assert 'cell_D_outer' in res
        assert np.asarray(res['cell_D_outer']).std() > 1e-4
