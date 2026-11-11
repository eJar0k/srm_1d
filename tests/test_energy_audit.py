"""
test_energy_audit.py -- solver-verification T6: sensible-enthalpy energy audit.

Pins (a) the PISO accounting identity (ledger storage + EOS jump + floor jump
== carried storage), (b) the per-step ledgers/jump decomposition on a short
real Chunc run, and (c) the audit post-processing on real and synthetic data.
All of it is bookkeeping: gate B0 enforces that the trajectory is unchanged.
"""
import numpy as np
import pytest

pytest.importorskip("numba")

from srm_1d.solver import _piso_step_with_energy_diagnostics  # noqa: E402
from srm_1d.tools.ignition_diagnostics import (  # noqa: E402
    audit_phase_masks, omitted_term_magnitude, sensible_enthalpy_audit,
)

R_GAS, CP_GAS, GAMMA = 300.0, 2000.0, 1.2


def _duct_step(N=4, dt=2.0e-5, floor=False):
    dx, A = 0.01, 1.0e-3
    P = np.linspace(2.0e5, 1.5e5, N) if not floor else np.full(N, 1500.0)
    T = np.linspace(900.0, 600.0, N)
    rho = P / (R_GAS * T)
    u = np.zeros(N + 1)
    u[1:N] = np.linspace(20.0, 60.0, N - 1)
    A_port = np.full(N, A)
    zN = np.zeros(N)
    src = np.full(N, 1.0e-2)
    out = _piso_step_with_energy_diagnostics(
        rho, u, P, T, A_port, np.full(N, 0.035),
        src, src * CP_GAS * 2500.0, np.zeros(N + 1), zN + 0.02,
        dx, dt, np.full(N, GAMMA), np.full(N, R_GAS), np.full(N, CP_GAS),
        np.full(N, 3000.0), A, 1.0 if floor else 101325.0, 300.0, N)
    return out, rho, A_port, dx


class TestPisoAccountingIdentity:
    @pytest.mark.parametrize("floor", [False, True])
    def test_ledger_plus_eos_plus_floor_equals_carried_storage(self, floor):
        out, rho, A_port, dx = _duct_step(floor=floor, dt=1.0e-4 if floor else 2e-5)
        rho_new, T_new = out[0], out[3]
        e_after = out[5]
        m_floor, e_floor = out[13], out[14]
        m_before, m_after, m_eos, e_eos = out[16], out[17], out[18], out[19]
        carried_m = float(np.sum(rho_new * A_port * dx))
        carried_e = float(np.sum(rho_new * A_port * dx * CP_GAS * T_new))
        assert m_before == pytest.approx(float(np.sum(rho * A_port * dx)), rel=1e-13)
        assert m_after + m_eos + m_floor == pytest.approx(carried_m, rel=1e-12)
        assert e_after + e_eos + e_floor == pytest.approx(carried_e, rel=1e-12)
        if floor:
            assert out[12] > 0 and m_floor > 0
        assert np.isfinite(out[20])  # nozzle_mdot_piso

    def test_tuple_length(self):
        out, *_ = _duct_step()
        assert len(out) == 24


@pytest.fixture(scope="module")
def chunc_short():
    from srm_1d.openmotor_adapter import run_from_ric
    return run_from_ric(
        "motors/machbusterNew.ric", roughness=32e-6, kappa=0.44,
        pyrogen="bpnv", pyrogen_mass=None, T_ignition=756.0, t_max=0.05,
        print_interval=0.0, verbose=False, snapshot_interval=5e-4,
        # D1 (2026-10-31): pinned to the pre-D1 sensible-enthalpy step. These tests pin the
        # Cp*T convention (energy_convention 'sensible_enthalpy_CpT') and the T10 identity
        # 'regression fill debit cancels the geometry jump', whose debit is an ENTHALPY sink; under
        # the total-energy default the geometry jump is internal+kinetic energy, so the two differ
        # by the p dV work of the opened volume (about 16 % of the debit in this short run, ~0.1 % of
        # the source power). That difference is the physical boundary work, not a leak: QF1 pins it in
        # TestRegressionFillDebitTotalEnergy below (total-energy default).
        energy_formulation='sensible', momentum_formulation='legacy',
        ignition_driver_temperature='static', head_end_jet_h='off')[0]


class TestRunLedgers:
    def test_summary_keys(self, chunc_short):
        s = chunc_short["summary"]
        assert s["energy_convention"] == "sensible_enthalpy_CpT"
        assert s["energy_omitted_terms"] == [
            "kinetic_energy", "pressure_work",
            "formation_energy_implicit_via_T_flame", "solid_energy_storage",
            "nozzle_kinetic_energy_static_outflow"]
        assert "energy_residual_convention" in s

    def test_history_keys_and_lengths(self, chunc_short):
        n = len(chunc_short["time"])
        keys = ["convective_wall_sink_power", "pressure_floor_energy",
                "source_ledger_mismatch_W", "nozzle_ledger_mismatch_W",
                "nozzle_enthalpy_power_post_piso", "nozzle_massflow_piso"]
        for c in ("eos", "floor", "mixture", "geometry", "remainder", "total"):
            keys += [f"storage_jump_{c}_J", f"mass_jump_{c}_kg"]
        for k in keys:
            a = np.asarray(chunc_short[k])
            assert a.shape == (n,), k
            assert np.all(np.isfinite(a)), k

    def test_source_ledger_closes(self, chunc_short):
        mis = np.abs(np.asarray(chunc_short["source_ledger_mismatch_W"]))
        peak = np.max(np.abs(np.asarray(chunc_short["thermal_source_power"])))
        assert peak > 0
        assert mis.max() <= 1e-6 * peak

    def test_conv_wall_sink_is_nonnegative_and_active(self, chunc_short):
        a = np.asarray(chunc_short["convective_wall_sink_power"])
        assert np.all(a >= 0.0) and a.max() > 0.0

    def test_jump_decomposition_is_exhaustive(self, chunc_short):
        r = chunc_short
        E = np.asarray(r["gas_sensible_energy"])
        Eb = np.asarray(r["gas_sensible_energy_before"])
        tot = np.asarray(r["storage_jump_total_J"])
        np.testing.assert_allclose(tot[:-1], Eb[1:] - E[:-1], rtol=0, atol=0)
        scale = np.max(np.abs(tot))
        assert scale > 0
        # remainder is round-off: the four named operations are exhaustive
        assert np.max(np.abs(r["storage_jump_remainder_J"])) <= 1e-9 * scale
        mscale = np.max(np.abs(r["mass_jump_total_kg"]))
        assert np.max(np.abs(r["mass_jump_remainder_kg"])) <= 1e-9 * mscale
        # the mixture refresh never changes mass (rho untouched)
        assert np.all(np.asarray(r["mass_jump_mixture_kg"]) == 0.0)
        # last slot is incomplete by construction
        assert tot[-1] == 0.0
        # floor alias = T3's mass + the new per-step energy
        np.testing.assert_array_equal(r["storage_jump_floor_J"], r["pressure_floor_energy"])
        np.testing.assert_array_equal(r["mass_jump_floor_kg"], r["pressure_floor_mass"])

    def test_geometry_jump_only_on_refresh_steps(self, chunc_short):
        g = np.asarray(chunc_short["storage_jump_geometry_J"])
        interval = chunc_short["summary"]["geometry_refresh_interval_steps"]
        # slot n holds the refresh that happens at the start of step n+1
        idx = np.flatnonzero(g[:-1] != 0.0)
        assert idx.size > 0
        assert np.all((idx + 1) % interval == 0)


class TestAuditFunction:
    def test_keys_and_closure(self, chunc_short):
        a = sensible_enthalpy_audit(chunc_short)
        assert a["convention"] == "sensible_enthalpy_CpT"
        for blk in [a["run"]] + list(a["phases"].values()):
            for k in ("energy_residual_max_abs_W", "energy_residual_rms_W",
                      "source_ledger_mismatch_max_abs_W",
                      "nozzle_ledger_mismatch_max_abs_W", "jumps_J",
                      "jumps_kg", "budget_J"):
                assert k in blk
            for c in ("eos", "floor", "mixture", "geometry", "remainder", "total"):
                for s in ("max_abs_single", "sum_abs", "sum_signed"):
                    assert np.isfinite(blk["jumps_J"][c][s])
                    assert np.isfinite(blk["jumps_kg"][c][s])
        run = a["run"]
        assert abs(run["budget_J"]["closure_error_with_jumps_frac_of_source"]) < 1e-9
        # phases partition the steps
        assert sum(p["n_steps"] for p in a["phases"].values()) == run["n_steps"]

    def test_phase_masks_partition(self, chunc_short):
        m = audit_phase_masks(chunc_short)["masks"]
        total = m["fill"].astype(int) + m["spike"] + m["plateau"]
        assert np.all(total == 1)

    def test_omitted_terms_from_snapshots(self, chunc_short):
        o = omitted_term_magnitude(chunc_short)
        assert 1.0 < o["gamma"] < 1.5
        n = sum(p["n_snapshots"] for p in o["phases"].values())
        assert n == len(chunc_short["snapshots"])


class TestOmittedTermSynthetic:
    def test_ratio_is_gamma_minus_one_over_two_mach_squared(self):
        N = 5
        mach = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
        snap = {"t": 0.0, "Mach": mach, "rho": np.ones(N),
                "D_port": np.ones(N)}
        t = np.linspace(0.0, 0.1, 11)
        res = {"time": t, "P_head": np.linspace(1e5, 1e6, 11),
               "snapshots": [snap]}
        o = omitted_term_magnitude(res, gamma=1.2)
        ph = [p for p in o["phases"].values() if p["n_snapshots"] == 1]
        assert len(ph) == 1
        p = ph[0]
        assert p["bore_max"]["max"] == pytest.approx(0.1 * 0.16)
        assert p["nozzle_inlet"]["max"] == pytest.approx(0.1 * 0.16)
        assert p["mass_weighted_mean"]["max"] == pytest.approx(
            0.1 * np.mean(mach ** 2))


# T9: the left-neighbour radiation sink is now debited (was a strict xfail, T6).
def test_radiation_sink_is_fully_debited_from_gas():
    from srm_1d.openmotor_adapter import run_from_ric
    r = run_from_ric(
        "motors/machbusterNew.ric", roughness=32e-6, kappa=0.44,
        pyrogen="bpnv", pyrogen_mass=None, T_ignition=756.0, t_max=0.02,
        radiation_emissivity=0.05, print_interval=0.0, verbose=False)[0]
    peak = np.max(np.abs(np.asarray(r["thermal_source_power"])))
    assert np.max(np.abs(np.asarray(r["source_ledger_mismatch_W"]))) <= 1e-6 * peak


class TestRegressionFillDebit:
    """T10: the gas filling newly opened port volume is debited at the source,
    so it cancels the geometry-refresh storage jump (which adds rho_g*dA*dx)."""

    def test_debit_cancels_geometry_jump(self, chunc_short):
        r = chunc_short
        dt = np.asarray(r["dt"])[:-1]
        geo_kg = np.sum(np.asarray(r["mass_jump_geometry_kg"])[:-1])
        geo_J = np.sum(np.asarray(r["storage_jump_geometry_J"])[:-1])
        deb_kg = np.sum(np.asarray(r["regression_fill_debit_massflow"])[:-1] * dt)
        deb_J = np.sum(np.asarray(r["regression_fill_debit_power"])[:-1] * dt)
        assert deb_kg > 0.0 and deb_J > 0.0
        assert geo_kg == pytest.approx(deb_kg, rel=1e-2)
        assert geo_J == pytest.approx(deb_J, rel=1e-2)

    def test_source_ledger_includes_debit(self, chunc_short):
        r = chunc_short
        peak = np.max(np.abs(np.asarray(r["thermal_source_power"])))
        assert np.max(np.abs(np.asarray(r["source_ledger_mismatch_W"]))) <= 1e-6 * peak


@pytest.fixture(scope="module")
def chunc_short_total():
    """Same run under the shipped default (total energy, conservative momentum, recovery driver, jet)."""
    from srm_1d.openmotor_adapter import run_from_ric
    return run_from_ric(
        "motors/machbusterNew.ric", roughness=32e-6, kappa=0.44,
        pyrogen="bpnv", pyrogen_mass=None, T_ignition=756.0, t_max=0.05,
        print_interval=0.0, verbose=False, snapshot_interval=5e-4)[0]


class TestRegressionFillDebitTotalEnergy:
    """QF1: control volume with a moving boundary, d(rho E V)/dt = mdot h_f - p dV/dt.

    The solver advances a fixed-area step and the geometry refresh then grows A_port at
    fixed (rho, T, u), adding rho dV (Cv T + KE). The fill debit removes rho dV (Cp T + KE),
    the TOTAL enthalpy of the gas that fills the opened volume, so debit + jump = p dV: the
    work the gas does on the receding wall. (Before QF1 the debit carried Cp T only, so the
    residual was the kinetic energy rho dV KE, ~0.4 % of the debit.)
    """

    @staticmethod
    def _sums(r):
        dt = np.asarray(r["dt"])[:-1]
        f = lambda k: float(np.sum(np.asarray(r[k])[:-1] * dt))  # noqa: E731
        g = lambda k: float(np.sum(np.asarray(r[k])[:-1]))  # noqa: E731
        return (f("regression_fill_debit_massflow"), f("regression_fill_debit_power"),
                f("regression_fill_pdv_power"),
                g("mass_jump_geometry_kg"), g("storage_jump_geometry_J"))

    def test_convention_is_total_energy(self, chunc_short_total):
        assert chunc_short_total["summary"]["energy_convention"] == "total_energy_CvT_plus_KE"

    def test_debit_minus_geometry_jump_is_the_boundary_work(self, chunc_short_total):
        deb_kg, deb_J, pdv_J, geo_kg, geo_J = self._sums(chunc_short_total)
        assert deb_J > 0.0 and pdv_J > 0.0
        assert geo_kg == pytest.approx(deb_kg, rel=1e-2)
        # the physical work is the ~16 % of the debit the total-energy jump does not cancel
        assert 0.10 < pdv_J / deb_J < 0.25
        # closure: debit - jump = p dV to 1 % of the work (the 0.3 % left is the step-vs-refresh
        # sampling of rho, T, u: the same lag the mass identity shows at 6e-4)
        assert deb_J - geo_J == pytest.approx(pdv_J, rel=1e-2)

    def test_ke_is_in_the_debit(self, chunc_short_total):
        """The kinetic-energy part of the jump is cancelled too: with a Cp T-only debit the
        residual would be rho dV KE, which this run's velocities make ~3.6e-3 of the debit."""
        deb_kg, deb_J, pdv_J, geo_kg, geo_J = self._sums(chunc_short_total)
        assert abs(deb_J - geo_J - pdv_J) / deb_J < 1.5e-3

    def test_energy_ledger_closes(self, chunc_short_total):
        r = chunc_short_total
        peak = np.max(np.abs(np.asarray(r["thermal_source_power"])))
        assert np.max(np.abs(np.asarray(r["source_ledger_mismatch_W"]))) <= 1e-6 * peak
        a = sensible_enthalpy_audit(r)
        assert abs(a["run"]["budget_J"]["closure_error_with_jumps_frac_of_source"]) < 1e-9

    def test_sensible_form_is_unchanged(self, chunc_short):
        """Sensible keeps debit == jump (Cp T both sides, no work term) to the sampling lag."""
        _, deb_J, _, _, geo_J = self._sums(chunc_short)
        assert geo_J == pytest.approx(deb_J, rel=1e-2)
