"""G2 step 1 (memo 2.2): mass-preserving correction for the length rounding.

Segment lengths are rounded to whole cells on CUMULATIVE position, and each
segment's burning perimeter is scaled by L_true / L_snapped, so the burning area
and the burned propellant mass equal the .ric values at every grid.
"""
import numpy as np
import pytest

from srm_1d.grain_geometry import build_snapped_geometry, update_cell_geometry
from srm_1d.propellant import Propellant, PropellantTab

D_OUT = 0.080
D_BORE = 0.030
LENS = (0.1234, 0.2345, 0.0987)
GAP = 0.004


def _spec(gap=GAP, lens=LENS, inhibit=False):
    out = []
    for i, L in enumerate(lens):
        d = {'D_bore_fwd': D_BORE, 'length': L, 'inhibit_fwd': inhibit, 'inhibit_aft': inhibit}
        if i < len(lens) - 1:
            d['gap_after'] = gap
        out.append(d)
    return out


def _prop():
    from tests.test_endface_conservation import _make_test_propellant
    return _make_test_propellant()


def _init(geo):
    ga = geo.compile_geometry_arrays()
    N, dx = geo.N_cells, geo.dx
    prop = _prop()
    A_port, C, Dh = np.zeros(N), np.zeros(N), np.zeros(N)
    ig, ems = np.zeros(N, dtype=np.bool_), np.zeros(N)
    update_cell_geometry(
        ga['regress'].copy(), ga['D_port'].copy(), ga['x_centers'], dx, N, ga['N_seg'],
        ga['cell_D_outer'], ga['seg_x_start'], ga['seg_length'],
        ga['seg_fwd_regression'], ga['seg_aft_regression'],
        ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'], ga['cell_segment_id'],
        np.full(N, 3e6), prop.rho_propellant, *prop.tab_arrays(), len(prop.tabs),
        A_port, C, Dh, ig, ems, ga['cell_D_bore_init'], ga['cell_wall_web'],
        ga['cell_segment_type'], ga['cell_fmm_idx'], ga['fmm_offset'],
        ga['fmm_reg_flat'], ga['fmm_perim_flat'], ga['fmm_port_flat'],
        ga['cell_burn_scale'])
    return ga, C, ems


@pytest.mark.parametrize('target', [17, 23, 40, 51, 77, 100, 137, 200, 331])
def test_burning_area_and_mass_are_the_ric_values_at_every_grid(target):
    geo = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=target)
    ga, C, _ = _init(geo)
    L_true = sum(LENS)
    # burning (bore) area: sum over cells of C_burn * dx = pi D L_true
    assert float(np.sum(C) * geo.dx) == pytest.approx(np.pi * D_BORE * L_true, rel=1e-12)
    # propellant volume = annulus * true length
    assert geo.total_propellant_volume() == pytest.approx(
        np.pi / 4 * (D_OUT ** 2 - D_BORE ** 2) * L_true, rel=1e-12)


@pytest.mark.parametrize('target', [200, 331, 500])
def test_rounding_is_on_cumulative_position(target):
    # lumped gaps (cells = propellant only): every segment end lies within dx/2 of its TRUE
    # cumulative propellant position, not up to k*dx/2; the legacy gap-cell build, when the
    # 1-cell minimum gap does not bind (dx < gap), is held to the same bound including its gaps
    geo = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=target)
    dx = geo.dx
    pos = 0.0
    for seg, L in zip(geo.segments, LENS):
        pos += L
        assert abs(seg.x_start + seg.length - pos) <= 0.5 * dx + 1e-12
        assert seg.length_true == L and seg.burn_scale == pytest.approx(L / seg.length)
    leg = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=target, lump_gaps=False)
    assert leg.dx < GAP
    pos = 0.0
    for seg, L in zip(leg.segments, LENS):
        pos += L
        assert abs(seg.x_start + seg.length - pos) <= 0.5 * leg.dx + 1e-12
        pos += GAP


def test_lengths_that_are_whole_cells_are_unchanged():
    dx = 0.01
    spec = [{'D_bore_fwd': D_BORE, 'length': 0.10, 'inhibit_fwd': True, 'inhibit_aft': True},
            {'D_bore_fwd': D_BORE, 'length': 0.20, 'inhibit_fwd': True, 'inhibit_aft': True}]
    geo = build_snapped_geometry(spec, D_OUT, target_propellant_cells=30)
    assert geo.dx == pytest.approx(dx)
    assert [s.burn_scale for s in geo.segments] == [1.0, 1.0]
    ga, C, _ = _init(geo)
    assert np.all(ga['cell_burn_scale'] == 1.0)


def test_conical_bore_volume_uses_the_true_length():
    spec = [{'D_bore_fwd': 0.02, 'D_bore_aft': 0.04, 'length': 0.3333}]
    for n in (13, 29, 64):
        geo = build_snapped_geometry(spec, D_OUT, target_propellant_cells=n)
        expect = np.pi / 4 * (D_OUT ** 2 - (0.02 ** 2 + 0.02 * 0.04 + 0.04 ** 2) / 3.0) * 0.3333
        assert geo.total_propellant_volume() == pytest.approx(expect, rel=1e-12)


def test_endface_mass_rate_is_unchanged_by_the_scale():
    from tests.test_endface_conservation import _analytical_mass_rate, _make_test_propellant
    prop = _make_test_propellant()
    geo = build_snapped_geometry([{'D_bore_fwd': 0.04, 'length': 0.3011, 'gap_after': 0.02},
                                  {'D_bore_fwd': 0.04, 'length': 0.2026}], 0.08, target_propellant_cells=47)
    assert any(sg.burn_scale != 1.0 for sg in geo.segments)
    ga, C, ems = _init(geo)
    # two segments, four burning faces: 4 rho r A_face (= 2 x the 2-face value), whatever the rounding
    assert float(np.sum(ems) * geo.dx) == pytest.approx(2.0 * _analytical_mass_rate(0.08, 0.04, prop), rel=1e-3)


# ---------------------------------------------------------------------------
# G2 step 2 (memo 2.1): inter-grain gaps are lumped volume, not cells
# ---------------------------------------------------------------------------

def test_lumped_gaps_have_no_cells_and_no_dx_clamp():
    spec = _spec(gap=0.003)
    lumped = build_snapped_geometry(spec, D_OUT, target_propellant_cells=100)
    legacy = build_snapped_geometry(spec, D_OUT, target_propellant_cells=100, lump_gaps=False)
    assert lumped.gaps_lumped and not legacy.gaps_lumped
    # no clamp: dx = L_prop / N; the legacy build clamps dx to the 3 mm gap
    assert lumped.dx == pytest.approx(sum(LENS) / 100, rel=0.02)
    assert legacy.dx <= 0.003 + 1e-12
    assert lumped.N_cells < 0.7 * legacy.N_cells
    ga = lumped.compile_geometry_arrays()
    assert np.all(ga['cell_segment_id'] >= 0)             # every cell is a grain cell


def test_gap_volume_is_split_half_and_half_into_the_two_neighbours():
    geo = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=100)
    ga = geo.compile_geometry_arrays()
    A_case = np.pi / 4 * D_OUT ** 2
    ve = ga['vol_extra']
    edges = [int(round((s.x_start + s.length) / geo.dx)) for s in geo.segments[:-1]]
    for e in edges:
        assert ve[e - 1] == pytest.approx(0.5 * A_case * GAP)
        assert ve[e] == pytest.approx(0.5 * A_case * GAP)
    # the default _spec head face is free -> the head gap is in cell 0; nothing else is lumped
    assert ve[0] == pytest.approx(A_case * geo.head_gap)
    assert np.count_nonzero(ve) == 2 * len(edges) + 1
    assert ve.sum() == pytest.approx(len(edges) * A_case * GAP + A_case * geo.head_gap)


def test_total_gas_volume_equals_the_legacy_gap_cells_volume():
    # same physical motor, legacy cells vs lumped gaps: the gap volume is the same
    # (fine grid so the legacy build's whole-cell gap rounding is small)
    lumped = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=3000)
    legacy = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=3000, lump_gaps=False)
    assert lumped.initial_gas_volume() == pytest.approx(legacy.initial_gas_volume(), rel=5e-3)


def test_x_offset_is_the_gaps_upstream():
    geo = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=100)
    ga = geo.compile_geometry_arrays()
    xo = ga['x_offset']
    edges = [int(round((s.x_start + s.length) / geo.dx)) for s in geo.segments[:-1]]
    assert xo[0] == 0.0 and xo[edges[0] - 1] == 0.0
    assert xo[edges[0]] == pytest.approx(GAP) and xo[edges[1]] == pytest.approx(2 * GAP)
    assert xo[-1] == pytest.approx(2 * GAP)
    # physical axis length = propellant + gaps
    assert ga['x_centers'][-1] + xo[-1] + 0.5 * geo.dx == pytest.approx(sum(LENS) + 2 * GAP, abs=geo.dx)


def test_lumped_gap_faces_stay_free_not_auto_inhibited():
    geo = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=100)
    ga = geo.compile_geometry_arrays()
    assert not ga['seg_inhibit_fwd'][1] and not ga['seg_inhibit_aft'][0]
    # ... while a bonded interface (gap 0) still auto-inhibits
    b = build_snapped_geometry(_spec(gap=0.0), D_OUT, target_propellant_cells=100)
    gb = b.compile_geometry_arrays()
    assert gb['seg_inhibit_fwd'][1] and gb['seg_inhibit_aft'][0]


@pytest.mark.parametrize('frac', [0.0, 0.3, 1.0, 1.5, 3.7])
def test_a_face_regressing_through_a_lump_conserves_mass_and_volume(frac):
    """All four free faces recede by delta: the end-face mass source stays 4 rho r A_face (partition of
    unity, also at the wall cells and the lumped interface), and the gas volume grows by exactly
    4 A_annulus delta (the vacated propellant), whatever cell the faces are in."""
    from tests.test_endface_conservation import _analytical_mass_rate, _make_test_propellant
    prop = _make_test_propellant()
    geo = build_snapped_geometry(_spec(inhibit=False), D_OUT, target_propellant_cells=100)
    ga0, C0, ems0 = _init(geo)
    N, dx = geo.N_cells, geo.dx
    ga = geo.compile_geometry_arrays()
    seg_id0 = ga['cell_segment_id'].copy()         # update_cell_geometry rewrites it in place
    delta = frac * dx
    fwd = ga['seg_fwd_regression'].copy(); fwd[:] = delta
    aft = ga['seg_aft_regression'].copy(); aft[:] = delta
    A_port, C, Dh = np.zeros(N), np.zeros(N), np.zeros(N)
    ig, ems = np.zeros(N, dtype=np.bool_), np.zeros(N)
    update_cell_geometry(
        ga['regress'].copy(), ga['D_port'].copy(), ga['x_centers'], dx, N, ga['N_seg'],
        ga['cell_D_outer'], ga['seg_x_start'], ga['seg_length'], fwd, aft,
        ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'], ga['cell_segment_id'],
        np.full(N, 3e6), prop.rho_propellant, *prop.tab_arrays(), len(prop.tabs),
        A_port, C, Dh, ig, ems, ga['cell_D_bore_init'], ga['cell_wall_web'],
        ga['cell_segment_type'], ga['cell_fmm_idx'], ga['fmm_offset'],
        ga['fmm_reg_flat'], ga['fmm_perim_flat'], ga['fmm_port_flat'], ga['cell_burn_scale'])
    n_faces = 2 * len(geo.segments)
    expect = n_faces / 2.0 * _analytical_mass_rate(D_OUT, D_BORE, prop)
    assert float(np.sum(ems) * dx) == pytest.approx(expect, rel=1e-3)
    A_init = np.where(seg_id0 >= 0, np.pi / 4 * D_BORE ** 2, np.pi / 4 * D_OUT ** 2)
    dV = float(np.sum(A_port - A_init) * dx)
    annulus = np.pi / 4 * (D_OUT ** 2 - D_BORE ** 2)
    assert dV == pytest.approx(n_faces * annulus * delta, rel=1e-9, abs=1e-15)
