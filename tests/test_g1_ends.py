"""G1: physical head and aft ends (no 3-cell full-diameter gas buffers).

* head: the first grain cell is cell 0; a closure-to-grain gap exists only when
  the head face is uninhibited and is LUMPED into cell 0's storage (no cells);
* aft: no full-diameter buffer; optional nozzle-convergent cone cells + a lumped
  remainder of the convergent volume;
* the solver stores gas in ``A_vol * dx`` (``A_vol = A_port + vol_extra / dx``)
  while the flow areas stay ``A_port``.
"""
import numpy as np
import pytest

from srm_1d.grain_geometry import (
    build_snapped_geometry, nozzle_convergent_cells, update_cell_geometry,
)
from srm_1d.nozzle import Nozzle
from srm_1d.solver import _piso_step_with_energy_diagnostics

D_OUT = 0.080
D_BORE = 0.040
D_T = 0.025


def _geo(inhibit_fwd=False, **kw):
    return build_snapped_geometry(
        [{'D_bore_fwd': D_BORE, 'length': 0.5, 'inhibit_fwd': inhibit_fwd,
          'inhibit_aft': True}],
        D_outer=D_OUT, target_propellant_cells=50, **kw)


class TestHead:
    def test_no_head_cells(self):
        geo = _geo(inhibit_fwd=True)
        ga = geo.compile_geometry_arrays()
        assert geo.N_cells == 50
        assert ga['cell_segment_id'][0] == 0           # cell 0 is a grain cell
        assert ga['x_centers'][0] == pytest.approx(0.5 * geo.dx)

    def test_bonded_head_has_no_gap(self):
        geo = _geo(inhibit_fwd=True)
        assert geo.head_gap == 0.0
        assert np.all(geo.compile_geometry_arrays()['vol_extra'] == 0.0)

    def test_uninhibited_head_gets_default_gap_lumped_in_cell0(self):
        geo = _geo(inhibit_fwd=False)
        assert geo.head_gap == pytest.approx(max(0.003, 0.05 * D_OUT))
        ga = geo.compile_geometry_arrays()
        assert geo.N_cells == 50                        # no cells added
        assert ga['vol_extra'][0] == pytest.approx(
            np.pi / 4 * D_OUT ** 2 * geo.head_gap)
        assert np.all(ga['vol_extra'][1:] == 0.0)

    def test_grid_independent_gap_volume(self):
        v = []
        for n in (25, 50, 100):
            g = build_snapped_geometry(
                [{'D_bore_fwd': D_BORE, 'length': 0.5}], D_outer=D_OUT,
                target_propellant_cells=n)
            v.append(float(np.sum(g.compile_geometry_arrays()['vol_extra'])))
        assert v[0] == pytest.approx(v[1]) == pytest.approx(v[2])

    def test_explicit_gap_overrides(self):
        geo = _geo(inhibit_fwd=True, head_gap=0.002)
        assert geo.head_gap == 0.002

    def test_free_volume_includes_gap(self):
        a = _geo(inhibit_fwd=True).initial_gas_volume()
        b = _geo(inhibit_fwd=True, head_gap=0.01).initial_gas_volume()
        assert b - a == pytest.approx(np.pi / 4 * D_OUT ** 2 * 0.01)


class TestAft:
    def test_none_without_nozzle(self):
        geo = _geo()
        assert geo.n_cone_cells == 0 and geo.aft_lump_volume == 0.0

    def test_cone_cells_follow_the_cone(self):
        noz = Nozzle(D_throat=D_T, D_exit=0.05, conv_angle=30.0)
        geo = _geo(nozzle=noz, aft_end='cone')
        ga = geo.compile_geometry_arrays()
        n = geo.n_cone_cells
        # 30 deg half angle: L_cone = (0.08-0.025)/2/tan30 = 47.6 mm, dx = 10 mm
        assert n >= 1
        D = ga['cell_D_outer'][-n:]
        assert D[0] < D_OUT and np.all(np.diff(D) < 0.0)
        for m in range(n):
            expect = D_OUT - 2 * (m + 0.5) * geo.dx * np.tan(np.radians(30.0))
            assert D[m] == pytest.approx(expect)
        assert np.all(ga['cell_segment_id'][-n:] == -1)

    def test_cut_keeps_isentropic_mach_below_grain_exit(self):
        # every kept cone cell has area >= the grain's initial aft port area
        noz = Nozzle(D_throat=D_T, D_exit=0.05, conv_angle=30.0)
        geo = _geo(nozzle=noz, aft_end='cone')
        ga = geo.compile_geometry_arrays()
        n = geo.n_cone_cells
        A_cone = np.pi / 4 * ga['cell_D_outer'][-n:] ** 2
        assert np.all(A_cone >= np.pi / 4 * D_BORE ** 2 - 1e-12)

    def test_convergent_volume_is_conserved(self):
        # resolved cone cells + lumped remainder = frustum volume
        L = (D_OUT - D_T) / 2 / np.tan(np.radians(30.0))
        V = np.pi / 12 * L * (D_OUT ** 2 + D_OUT * D_T + D_T ** 2)
        for mode in ('cone', 'lump'):
            n, lump, _ = nozzle_convergent_cells(
                D_OUT, D_T, 30.0, np.pi / 4 * D_BORE ** 2, 0.01, mode)
            Vc = sum(np.pi / 4 * (D_OUT - 2 * (m + .5) * 0.01 * np.tan(np.radians(30))) ** 2
                     * 0.01 for m in range(n))
            assert Vc + lump == pytest.approx(V)

    def test_modes(self):
        n, lump, _ = nozzle_convergent_cells(D_OUT, D_T, 30.0, 0.0, 0.01, 'none')
        assert (n, lump) == (0, 0.0)
        n, lump, _ = nozzle_convergent_cells(D_OUT, D_T, 30.0, 0.0, 0.01, 'cone_nolump')
        assert n > 0 and lump == 0.0
        with pytest.raises(ValueError):
            nozzle_convergent_cells(D_OUT, D_T, 30.0, 0.0, 0.01, 'bogus')

    def test_lump_goes_to_last_cell(self):
        noz = Nozzle(D_throat=D_T, D_exit=0.05, conv_angle=30.0)
        geo = _geo(nozzle=noz, aft_end='lump')
        ga = geo.compile_geometry_arrays()
        assert geo.n_cone_cells == 0
        assert ga['vol_extra'][-1] == pytest.approx(geo.aft_lump_volume)
        assert geo.aft_lump_volume > 0.0


class TestEndFaceAtWalls:
    """With no buffer cells the burning end faces sit on the wall: the hat
    function must still return the full face mass (partition of unity)."""

    @pytest.mark.parametrize('reg', [0.0, 0.0003, 0.004, 0.0137])
    def test_both_faces_conserved_when_regressed(self, reg):
        from tests.test_endface_conservation import (
            _make_test_propellant, _analytical_mass_rate)
        prop = _make_test_propellant()
        geo = build_snapped_geometry(
            [{'D_bore_fwd': D_BORE, 'length': 0.5}], D_outer=D_OUT,
            target_propellant_cells=50, head_gap=0.0)
        ga = geo.compile_geometry_arrays()
        N, dx = geo.N_cells, geo.dx
        fwd = ga['seg_fwd_regression'].copy(); fwd[:] = reg
        aft = ga['seg_aft_regression'].copy(); aft[:] = reg
        A_port = np.zeros(N); C = np.zeros(N); Dh = np.zeros(N)
        ig = np.zeros(N, dtype=np.bool_); ems = np.zeros(N)
        update_cell_geometry(
            ga['regress'].copy(), ga['D_port'].copy(), ga['x_centers'], dx, N,
            ga['N_seg'], ga['cell_D_outer'], ga['seg_x_start'], ga['seg_length'],
            fwd, aft, ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'],
            ga['cell_segment_id'], np.full(N, 3e6), prop.rho_propellant,
            *prop.tab_arrays(), len(prop.tabs), A_port, C, Dh, ig, ems,
            ga['cell_D_bore_init'], ga['cell_wall_web'], ga['cell_segment_type'],
            ga['cell_fmm_idx'], ga['fmm_offset'], ga['fmm_reg_flat'],
            ga['fmm_perim_flat'], ga['fmm_port_flat'])
        # head face is uninhibited here (seg_inhibit_fwd False); aft forced free
        expected = _analytical_mass_rate(D_OUT, D_BORE, prop)
        assert float(np.sum(ems) * dx) == pytest.approx(expected, rel=1e-3)


def _piso_args(N, A_port, mass_source):
    gamma, R, Cp = 1.2, 300.0, 2000.0
    rho = np.full(N, 101325.0 / (R * 300.0))
    u = np.zeros(N + 1)
    P = np.full(N, 101325.0)
    T = np.full(N, 300.0)
    return dict(
        rho=rho, u=u, P=P, T=T, A_port=A_port, D_hyd=np.full(N, 0.035),
        mass_source=mass_source, thermal_source=np.zeros(N),
        momentum_source=np.zeros(N + 1), f_darcy=np.zeros(N),
        dx=0.01, dt=1.0e-6, gamma_arr=np.full(N, gamma), R_arr=np.full(N, R),
        Cp_arr=np.full(N, Cp), T_ceiling_arr=np.full(N, 3000.0),
        A_throat=1.0e-9, P_ambient=101325.0, T_ambient=300.0, N=N)


def _step(args, vol_extra=None):
    a = dict(args)
    for k in ('rho', 'u', 'P', 'T'):
        a[k] = a[k].copy()
    return _piso_step_with_energy_diagnostics(
        a['rho'], a['u'], a['P'], a['T'], a['A_port'], a['D_hyd'],
        a['mass_source'], a['thermal_source'], a['momentum_source'],
        a['f_darcy'], a['dx'], a['dt'], a['gamma_arr'], a['R_arr'],
        a['Cp_arr'], a['T_ceiling_arr'], a['A_throat'], a['P_ambient'],
        a['T_ambient'], a['N'], vol_extra=vol_extra)


class TestLumpedStorage:
    N = 5

    def _args(self):
        A = np.full(self.N, 1.0e-3)
        src = np.zeros(self.N)
        src[0] = 50.0            # kg/(s m) injected into cell 0
        return _piso_args(self.N, A, src)

    def test_zero_extra_is_bitwise_default(self):
        args = self._args()
        base = _step(args)
        zero = _step(args, np.zeros(self.N))
        for b, z in zip(base[:4], zero[:4]):
            assert np.array_equal(b, z)

    def test_gas_mass_ledger_counts_the_lumped_volume(self):
        args = self._args()
        extra = np.zeros(self.N); extra[0] = 2.0e-5     # = 2 cells' worth
        out = _step(args, extra)
        # gas_mass_before is index -? find by value: sum(rho * (A dx + extra))
        expected = float(np.sum(args['rho'] * (args['A_port'] * args['dx'] + extra)))
        assert any(abs(x - expected) < 1e-12 * expected for x in out
                   if np.isscalar(x))

    def test_extra_volume_lowers_the_pressure_rise_of_the_same_injection(self):
        args = self._args()
        extra = np.zeros(self.N); extra[0] = 3.0e-5     # V0: 1e-5 -> 4e-5
        P0 = _step(args)[2]
        P1 = _step(args, extra)[2]
        dP0 = P0[0] - args['P'][0]
        dP1 = P1[0] - args['P'][0]
        assert dP0 > 0.0 and dP1 > 0.0
        # acoustic compliance ~ V: rise ratio ~ V_without / V_with = 1/4
        assert dP1 / dP0 == pytest.approx(0.25, rel=0.25)
