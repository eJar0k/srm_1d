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
    # fine enough that the 1-cell minimum gap does not bind (dx < gap): every segment end lies
    # within dx/2 of its TRUE cumulative position (propellant + gaps), not up to k*dx/2
    geo = build_snapped_geometry(_spec(), D_OUT, target_propellant_cells=target)
    dx = geo.dx
    assert dx < GAP
    pos = 0.0
    for k, (seg, L) in enumerate(zip(geo.segments, LENS)):
        pos += L
        assert abs(seg.x_start + seg.length - pos) <= 0.5 * dx + 1e-12
        assert seg.length_true == L and seg.burn_scale == pytest.approx(L / seg.length)
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
