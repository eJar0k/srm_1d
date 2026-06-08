"""
test_station_viz.py — headless tests for the per-station axial-viz backend
(``srm_1d.station_viz``; design ``docs/v0_8_0/STATION_VIZ_DESIGN.md``).

The station model + payload builder are pure-logic / Qt-free, so most tests
run on a synthetic ``cell_segment_id`` or a geometry's compiled arrays — no
Numba simulation needed. One end-to-end test confirms a real
``run_simulation`` result carries the v0.8.x data contract.
"""
import numpy as np
import pytest

from srm_1d.station_viz import (
    GAP_SENTINEL,
    AxialPayload,
    Station,
    build_axial_payload,
    default_stations,
    make_station,
    grain_cell_spans,
    gap_cell_indices,
    cell_categories,
    grain_role,
    classify_cell,
    station_full_label,
    _decimate_frame_indices,
)


# ================================================================
# Synthetic cell→grain maps
# ================================================================
# Two grains of 4 cells each, separated by a 2-cell gap:
#   cells:  0 1 2 3 | 4 5 | 6 7 8 9
#   grain:  0 0 0 0 |-1 -1| 1 1 1 1
SEG_TWO_GRAIN = np.array([0, 0, 0, 0, -1, -1, 1, 1, 1, 1], dtype=np.int64)
X_TWO_GRAIN = np.arange(10, dtype=float) * 0.01  # 10 mm spacing


class TestSpansAndGaps:

    def test_grain_cell_spans(self):
        spans = grain_cell_spans(SEG_TWO_GRAIN)
        assert set(spans.keys()) == {0, 1}
        np.testing.assert_array_equal(spans[0], [0, 1, 2, 3])
        np.testing.assert_array_equal(spans[1], [6, 7, 8, 9])

    def test_gap_cell_indices(self):
        np.testing.assert_array_equal(gap_cell_indices(SEG_TWO_GRAIN), [4, 5])

    def test_no_gap_motor(self):
        seg = np.array([0, 0, 0, 1, 1, 1], dtype=np.int64)
        assert gap_cell_indices(seg).size == 0
        assert set(grain_cell_spans(seg)) == {0, 1}


class TestDefaultStations:

    def test_three_per_grain_fore_on(self):
        st = default_stations(SEG_TWO_GRAIN, X_TWO_GRAIN)
        # 3 stations per grain × 2 grains.
        assert len(st) == 6
        for grain in (0, 1):
            grp = [s for s in st if s.grain == grain]
            assert [s.role for s in grp] == ['fore', 'mid', 'aft']
            assert [s.active for s in grp] == [True, False, False]
            # all cells belong to the grain's span
            span = grain_cell_spans(SEG_TWO_GRAIN)[grain]
            assert all(s.cell_index in span for s in grp)

    def test_fore_is_foremost_cell(self):
        st = default_stations(SEG_TWO_GRAIN, X_TWO_GRAIN)
        g1_fore = next(s for s in st if s.grain == 1 and s.role == 'fore')
        assert g1_fore.cell_index == 6  # foremost cell of grain 1
        g1_aft = next(s for s in st if s.grain == 1 and s.role == 'aft')
        assert g1_aft.cell_index == 9

    def test_no_stations_in_gaps(self):
        st = default_stations(SEG_TWO_GRAIN, X_TWO_GRAIN)
        assert all(s.grain != GAP_SENTINEL for s in st)
        assert all(s.cell_index not in (4, 5) for s in st)

    def test_position_label_in_mm(self):
        st = default_stations(SEG_TWO_GRAIN, X_TWO_GRAIN)
        g0_aft = next(s for s in st if s.grain == 0 and s.role == 'aft')
        assert g0_aft.position_m == pytest.approx(0.03)  # cell 3 at 30 mm
        assert '30 mm' in g0_aft.label
        assert 'Grain 1' in g0_aft.label  # grain index 0 -> "Grain 1"

    def test_short_grain_dedupes(self):
        # A 2-cell grain: fore=0, mid=cells[1]=1, aft=1 → mid and aft collide.
        seg = np.array([0, 0], dtype=np.int64)
        x = np.array([0.0, 0.01])
        st = default_stations(seg, x)
        cells = [s.cell_index for s in st]
        assert len(cells) == len(set(cells))  # no duplicate cells
        assert 0 in cells and 1 in cells

    def test_single_cell_grain(self):
        seg = np.array([0], dtype=np.int64)
        st = default_stations(seg, np.array([0.0]))
        assert len(st) == 1
        assert st[0].cell_index == 0
        assert st[0].active is True


# head 0-1 | Grain1 2-5 | Gap1 6-7 | Grain2 8-11 | Aft 12-13
SEG_FULL = np.array([-1, -1, 0, 0, 0, 0, -1, -1, 1, 1, 1, 1, -1, -1], dtype=np.int64)
X_FULL = np.arange(14, dtype=float) * 0.01


class TestCellCategories:

    def test_ordered_head_grain_gap_aft(self):
        cats = cell_categories(SEG_FULL)
        kinds = [(c['kind'], c['lo'], c['hi']) for c in cats]
        assert kinds == [
            ('head', 0, 1), ('grain', 2, 5), ('gap', 6, 7),
            ('grain', 8, 11), ('aft', 12, 13),
        ]
        # gaps numbered, grains short-labeled
        gap = next(c for c in cats if c['kind'] == 'gap')
        assert gap['gap'] == 1 and gap['label'] == 'Gap 1'
        g2 = [c for c in cats if c['kind'] == 'grain'][1]
        assert g2['grain'] == 1 and g2['short'] == 'G2'

    def test_no_head_when_grain_at_zero(self):
        seg = np.array([0, 0, -1, 1, 1], dtype=np.int64)
        kinds = [c['kind'] for c in cell_categories(seg)]
        assert kinds == ['grain', 'gap', 'grain']   # no head, no aft

    def test_multiple_gaps_numbered(self):
        seg = np.array([0, -1, 1, -1, 2], dtype=np.int64)
        gaps = [c for c in cell_categories(seg) if c['kind'] == 'gap']
        assert [g['gap'] for g in gaps] == [1, 2]


class TestGrainRole:

    def test_fore_mid_aft(self):
        assert grain_role(2, 0, SEG_FULL) == 'fore'
        assert grain_role(5, 0, SEG_FULL) == 'aft'
        assert grain_role(4, 0, SEG_FULL) == 'mid'   # cells[len//2] of [2,3,4,5]
        assert grain_role(3, 0, SEG_FULL) == ''      # interior, no role

    def test_non_grain_has_no_role(self):
        assert grain_role(0, -1, SEG_FULL) == ''


class TestClassifyCell:

    def test_categories(self):
        assert classify_cell(0, SEG_FULL)['kind'] == 'head'
        g = classify_cell(2, SEG_FULL, X_FULL)
        assert g['kind'] == 'grain' and g['grain'] == 0 and g['role'] == 'fore'
        assert g['position_m'] == pytest.approx(0.02)
        assert classify_cell(7, SEG_FULL)['kind'] == 'gap'
        assert classify_cell(7, SEG_FULL)['gap'] == 1
        assert classify_cell(12, SEG_FULL)['kind'] == 'aft'

    def test_reclassifies_on_index_change(self):
        # The same station 'moves' categories purely by its cell index.
        assert classify_cell(5, SEG_FULL)['role'] == 'aft'      # grain aft
        assert classify_cell(6, SEG_FULL)['kind'] == 'gap'      # next cell -> gap

    def test_full_labels(self):
        assert station_full_label(classify_cell(2, SEG_FULL)) == 'G1 fore (c2)'
        assert station_full_label(classify_cell(3, SEG_FULL)) == 'G1 (c3)'
        assert station_full_label(classify_cell(0, SEG_FULL)) == 'Head (c0)'
        assert station_full_label(classify_cell(7, SEG_FULL)) == 'Gap1 (c7)'
        assert station_full_label(classify_cell(12, SEG_FULL)) == 'Aft (c12)'


class TestMakeStation:

    def test_classifies_grain(self):
        s = make_station(7, SEG_TWO_GRAIN, X_TWO_GRAIN)
        assert s.grain == 1
        assert s.role == 'custom'
        assert s.cell_index == 7

    def test_gap_cell_forces_gap_role(self):
        s = make_station(4, SEG_TWO_GRAIN, X_TWO_GRAIN)
        assert s.grain == GAP_SENTINEL
        assert s.role == 'gap'
        assert 'Gap' in s.label

    def test_out_of_range_raises(self):
        with pytest.raises(IndexError):
            make_station(99, SEG_TWO_GRAIN, X_TWO_GRAIN)


class TestDecimation:

    def test_keeps_all_below_budget(self):
        np.testing.assert_array_equal(_decimate_frame_indices(5, 240), [0, 1, 2, 3, 4])

    def test_keeps_first_and_last(self):
        idx = _decimate_frame_indices(1000, 50)
        assert idx[0] == 0 and idx[-1] == 999
        assert len(idx) <= 50

    def test_zero_budget_keeps_all(self):
        np.testing.assert_array_equal(_decimate_frame_indices(3, 0), [0, 1, 2])

    def test_empty(self):
        assert _decimate_frame_indices(0, 10).size == 0


class TestBuildAxialPayloadSynthetic:

    def _fake_result(self, n_snaps, n_cells):
        x = np.arange(n_cells, dtype=float) * 0.01
        snaps = []
        for s in range(n_snaps):
            snaps.append({
                't': float(s),
                'x': x,
                'P': np.full(n_cells, float(s)),
                'u': np.arange(n_cells, dtype=float) + s,
                'Mach': np.zeros(n_cells),
                'T': np.zeros(n_cells),
                'r_total': np.zeros(n_cells),
                'r_erosive': np.zeros(n_cells),
                'D_port': np.zeros(n_cells),
                'regress': np.zeros(n_cells),
                'rho': np.full(n_cells, 2.0),
            })
        return {
            'snapshots': snaps,
            'cell_segment_id': SEG_TWO_GRAIN[:n_cells],
            'x_cell': x,
            'dx': 0.01,
            'D_outer': 0.06,
            'cell_wall_web': np.full(n_cells, 0.02),
        }

    def test_shapes_and_contents(self):
        res = self._fake_result(n_snaps=10, n_cells=10)
        pl = build_axial_payload(res, max_frames=0)
        assert isinstance(pl, AxialPayload)
        assert pl.n_frames == 10
        assert pl.n_cells == 10
        # P field at frame s, any cell == s (we set P = s everywhere).
        assert pl.fields['P'][3, 5] == 3.0
        # series slice
        np.testing.assert_array_equal(pl.series('P', 0), np.arange(10.0))

    def test_decimation_applied(self):
        res = self._fake_result(n_snaps=100, n_cells=10)
        pl = build_axial_payload(res, max_frames=20)
        assert pl.n_frames <= 20
        assert pl.snap_times[0] == 0.0
        assert pl.snap_times[-1] == 99.0
        # decimated field rows align with kept frame times (P == t).
        np.testing.assert_allclose(pl.fields['P'][:, 0], pl.snap_times)

    def test_mass_flux_G_derived(self):
        # rho = 2.0 everywhere; u = arange + s. G = rho * u must hold cellwise.
        res = self._fake_result(n_snaps=5, n_cells=10)
        pl = build_axial_payload(res, max_frames=0)
        assert 'rho' in pl.fields
        assert 'G' in pl.fields
        np.testing.assert_allclose(pl.fields['G'], pl.fields['rho'] * pl.fields['u'])
        assert pl.fields['G'].shape == (pl.n_frames, pl.n_cells)

    def test_slice_geometry_carried(self):
        # Roadmap #2: build_axial_payload carries dx / D_outer / cell_wall_web.
        res = self._fake_result(n_snaps=4, n_cells=10)
        pl = build_axial_payload(res, max_frames=0)
        assert pl.dx == 0.01
        assert pl.D_outer == 0.06
        assert pl.cell_wall_web.shape == (10,)
        np.testing.assert_allclose(pl.cell_wall_web, 0.02)

    def test_slice_geometry_defaults_when_absent(self):
        # Pre-v0.8.x results lack the geometry → graceful defaults, no error.
        res = self._fake_result(n_snaps=3, n_cells=10)
        for k in ('dx', 'D_outer', 'cell_wall_web'):
            del res[k]
        pl = build_axial_payload(res)
        assert pl.dx == 0.0
        assert pl.D_outer == 0.0
        assert pl.cell_wall_web.size == 0

    def test_G_absent_without_rho(self):
        # Old-style result without a rho snapshot: G degrades to absent, not faked.
        res = self._fake_result(n_snaps=3, n_cells=10)
        for snap in res['snapshots']:
            del snap['rho']
        pl = build_axial_payload(res)
        assert 'rho' not in pl.fields
        assert 'G' not in pl.fields

    def test_missing_field_skipped(self):
        res = self._fake_result(n_snaps=3, n_cells=10)
        pl = build_axial_payload(res, fields=('P', 'nonexistent'))
        assert 'P' in pl.fields
        assert 'nonexistent' not in pl.fields

    def test_no_snapshots_returns_none(self):
        assert build_axial_payload({'snapshots': []}) is None
        assert build_axial_payload({}) is None

    def test_missing_cell_segment_id_raises(self):
        res = self._fake_result(n_snaps=3, n_cells=10)
        del res['cell_segment_id']
        with pytest.raises(ValueError, match="cell_segment_id"):
            build_axial_payload(res)

    def test_cell_segment_id_override(self):
        res = self._fake_result(n_snaps=3, n_cells=10)
        del res['cell_segment_id']
        pl = build_axial_payload(res, cell_segment_id=SEG_TWO_GRAIN)
        np.testing.assert_array_equal(pl.cell_segment_id, SEG_TWO_GRAIN)


# ================================================================
# End-to-end: a real run_simulation result carries the contract
# ================================================================

class TestEndToEndContract:

    def test_multi_grain_result_carries_payload(self):
        from srm_1d import run_simulation
        from srm_1d.nozzle import Nozzle
        from srm_1d.propellant import Pyrogen
        from srm_1d.igniter_plenum import PyrogenChamber
        from tests._motor_fixtures import bates_motor_geo, hasegawa_propellant_1

        geo = bates_motor_geo(
            D_bore=0.020, D_outer=0.060, L_segment=0.080,
            N_segments=2, spacing=0.010, target_propellant_cells=40,
        )
        prop = hasegawa_propellant_1()
        nozzle = Nozzle(D_throat=0.012, D_exit=0.024, efficiency=0.95)
        pyro = Pyrogen("t", 3.0e-5, 0.5, 1700.0, 2800.0, 0.030, 1.25,
                       heat_flux_cal_cm2_s=69.4)
        chamber = PyrogenChamber(pyro, 0.004, 5.0e-4, 2.0e-5, 3.0e-6,
                                 "end_burning")
        result = run_simulation(
            geo, prop, nozzle, chamber,
            roughness=20e-6, T_ignition=294.0, P_cutoff=0.5e6,
            snapshot_interval=0.5, print_interval=1e9, t_max=6.0,
        )

        # Contract fields present.
        assert 'cell_segment_id' in result
        assert 'x_cell' in result
        n_cells = result['x_cell'].shape[0]
        assert result['cell_segment_id'].shape[0] == n_cells

        # Two grains + a gap present.
        spans = grain_cell_spans(result['cell_segment_id'])
        assert set(spans.keys()) == {0, 1}
        assert gap_cell_indices(result['cell_segment_id']).size > 0

        # Payload + stations build cleanly off the real result.
        pl = build_axial_payload(result)
        assert pl.n_cells == n_cells
        assert pl.n_frames == len(result['snapshots'])
        for name in ('P', 'u', 'Mach', 'T', 'r_total', 'regress', 'D_port'):
            assert pl.fields[name].shape == (pl.n_frames, n_cells)

        # v0.8.x: per-cell density snapshot + derived mass flux G = rho * u.
        assert 'rho' in pl.fields and 'G' in pl.fields
        assert np.all(pl.fields['rho'] > 0.0)          # physical density
        np.testing.assert_allclose(pl.fields['G'], pl.fields['rho'] * pl.fields['u'])

        # Roadmap #2 slice geometry carried + physical.
        assert pl.dx > 0.0
        assert pl.D_outer > 0.0
        assert pl.cell_wall_web.shape == (n_cells,)
        assert np.all(pl.cell_wall_web >= 0.0)
        # Bore never exceeds the casing: D_port/2 <= D_outer/2 on grain cells.
        grain = pl.cell_segment_id >= 0
        assert np.all(pl.fields['D_port'][:, grain] <= pl.D_outer + 1e-9)

        # Roadmap #2 face burnback: per-segment geometry carried; faces recede.
        sg = pl.seg_geom
        assert sg and {'seg_x_start', 'seg_length', 'seg_fwd_reg', 'seg_aft_reg'} <= set(sg)
        n_seg = sg['seg_x_start'].size
        assert sg['seg_fwd_reg'].shape == (pl.n_frames, n_seg)
        assert sg['seg_aft_reg'].shape == (pl.n_frames, n_seg)
        # Faces start unburned and recede (uninhibited faces present this fixture).
        assert np.all(sg['seg_fwd_reg'] >= -1e-12)
        assert (sg['seg_fwd_reg'][-1].max() + sg['seg_aft_reg'][-1].max()) > 0.0

        st = default_stations(result['cell_segment_id'], result['x_cell'])
        assert len([s for s in st if s.role == 'fore' and s.active]) == 2
        # every default station samples a real grain cell
        assert all(result['cell_segment_id'][s.cell_index] == s.grain for s in st)
