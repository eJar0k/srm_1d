"""
v0.7.1 Phase 1f — N-species mass-fraction transport tests.

Tests cover:
- ``_advect_species`` kernel invariants (sum=1, bounds)
- Pure-flow advection without sources
- Source-only update without flow
- Integration: short Hasegawa A run + species accounting from the
  returned ``Y_species_final``.
"""
import numpy as np
import pytest

from srm_1d.simulation import _advect_species, run_simulation
from srm_1d.openmotor_adapter import load_ric, load_transport, load_pyrogen, convert_propellant
from srm_1d.igniter_plenum import PyrogenChamber


# ================================================================
# Direct kernel tests
# ================================================================

def test_advect_species_pure_advection_pyrogen_pulse_moves_right():
    """A 100% pyrogen blob at cell 0 with rightward flow should migrate
    into cell 1 over time, preserving Y invariants."""
    N = 5
    S = 3
    dx = 0.01
    dt = 1.0e-5
    A_port = np.full(N, 1.0e-4)
    rho_old = np.full(N, 1.0)
    rho_new = np.full(N, 1.0)  # steady flow
    # Uniform rightward face velocity
    u = np.full(N + 1, 50.0)
    u[0] = 0.0  # wall
    Y = np.zeros((N, S))
    Y[0, 0] = 1.0  # cell 0: 100% pyrogen
    for i in range(1, N):
        Y[i, 2] = 1.0  # rest: 100% ambient

    mass_source = np.zeros((N, S))
    nozzle_mdot = 0.0  # nothing leaves yet

    Y_before = Y.copy()
    _advect_species(Y, rho_old, rho_new, u, A_port, nozzle_mdot, dx, dt,
                    mass_source, N, S)

    # After one step, some pyrogen should have moved from cell 0 to cell 1
    assert Y[0, 0] < Y_before[0, 0], "Pyrogen should leave cell 0"
    assert Y[1, 0] > Y_before[1, 0], "Pyrogen should enter cell 1"

    # Y invariants
    for i in range(N):
        s_total = Y[i, :].sum()
        assert abs(s_total - 1.0) < 1.0e-10, f"sum(Y[{i}, :]) = {s_total}"
        for s in range(S):
            assert 0.0 <= Y[i, s] <= 1.0


def test_advect_species_source_only_no_flow():
    """Source-only update (no advection): pyrogen mass added at cell 0
    should increase Y[0, 0] proportionally."""
    N = 5
    S = 3
    dx = 0.01
    dt = 1.0e-3
    A_port = np.full(N, 1.0e-4)
    rho_old = np.full(N, 1.0)
    # rho_new reflects added mass at cell 0: rho_new = rho_old + source/V * dt
    mdot_per_dx = 0.1  # kg/s/m at cell 0
    rho_new = rho_old.copy()
    V0 = A_port[0] * dx
    rho_new[0] = rho_old[0] + (mdot_per_dx * dx * dt) / V0  # V_i = A_port * dx
    # Actually: added_mass = mdot_per_dx * dx * dt; rho_new[0] = (m_old + added)/V
    # m_old[0] = rho_old[0] * V_i = 1.0 * 1e-6 = 1e-6 kg
    # added = 0.1 * 0.01 * 1e-3 = 1e-6 kg → rho_new[0] = 2.0

    u = np.zeros(N + 1)
    Y = np.zeros((N, S))
    Y[:, 2] = 1.0  # all ambient initially

    mass_source = np.zeros((N, S))
    mass_source[0, 0] = mdot_per_dx  # pyrogen source at cell 0

    _advect_species(Y, rho_old, rho_new, u, A_port, 0.0, dx, dt,
                    mass_source, N, S)

    # Cell 0 mass doubled (50% old ambient + 50% new pyrogen)
    assert abs(Y[0, 0] - 0.5) < 1.0e-6, f"Y[0, pyrogen]={Y[0, 0]}"
    assert abs(Y[0, 2] - 0.5) < 1.0e-6, f"Y[0, ambient]={Y[0, 2]}"
    # Other cells untouched
    for i in range(1, N):
        assert abs(Y[i, 2] - 1.0) < 1.0e-10
        assert Y[i, 0] == 0.0

    # Y invariants
    for i in range(N):
        assert abs(Y[i, :].sum() - 1.0) < 1.0e-10


def test_advect_species_renormalization_repairs_fp_drift():
    """Tiny FP drift in Y[i, :] should be repaired by the renormalization
    step inside the kernel."""
    N = 3
    S = 3
    dx = 0.01
    dt = 1.0e-5
    A_port = np.full(N, 1.0e-4)
    rho_old = np.full(N, 1.0)
    rho_new = np.full(N, 1.0)
    u = np.zeros(N + 1)
    Y = np.zeros((N, S))
    # Inject deliberate drift
    Y[0, 0] = 0.5 + 1e-9
    Y[0, 1] = 0.5
    Y[0, 2] = 0.0
    # Cells 1, 2 valid
    Y[1, 2] = 1.0
    Y[2, 2] = 1.0

    mass_source = np.zeros((N, S))
    _advect_species(Y, rho_old, rho_new, u, A_port, 0.0, dx, dt,
                    mass_source, N, S)
    # Renormalization should clean up cell 0
    assert abs(Y[0, :].sum() - 1.0) < 1.0e-12


def test_advect_species_nozzle_outflow_drains_last_cell():
    """Nozzle outflow at face N should remove mass from cell N-1, carrying
    its current Y composition."""
    N = 3
    S = 3
    dx = 0.01
    dt = 1.0e-4
    A_port = np.full(N, 1.0e-4)
    rho_old = np.full(N, 1.0)
    # Simulate a nozzle drain: rho_new[N-1] < rho_old[N-1]
    rho_new = rho_old.copy()
    nozzle_mdot = 0.005  # kg/s
    # mass loss = nozzle_mdot * dt = 5e-7 kg; V_{N-1} = 1e-6 m^3
    rho_new[N - 1] = rho_old[N - 1] - (nozzle_mdot * dt) / (A_port[N - 1] * dx)

    u = np.zeros(N + 1)
    Y = np.zeros((N, S))
    Y[:, 1] = 1.0  # all propellant species

    mass_source = np.zeros((N, S))
    Y_before = Y.copy()
    _advect_species(Y, rho_old, rho_new, u, A_port, nozzle_mdot, dx, dt,
                    mass_source, N, S)
    # Y composition at cell N-1 should remain 100% propellant (mass left
    # carries the same composition, so fraction is unchanged).
    assert abs(Y[N - 1, 1] - 1.0) < 1.0e-9
    assert abs(Y[N - 1, :].sum() - 1.0) < 1.0e-10


# ================================================================
# Integration test against real Hasegawa A motor
# ================================================================

def _build_hasegawa_a_inputs():
    """Build the kwargs for a short Hasegawa A run with v0.7.0 calibrated
    parameters (rank-1 LHS). Mirrors examples/hasegawa_motor_a.py."""
    motor = load_ric('srm_1d/motors/hasegawa_a.ric')
    gp = load_transport('srm_1d/motors/hasegawa_a.transport.yaml')
    prop = convert_propellant(motor['propellant'], gp)
    pyr = load_pyrogen('bpnv')
    pyc = PyrogenChamber(
        pyrogen=pyr,
        m_pyrogen_initial=12.3e-3,
        A_burn_initial=1e-4,
        A_throat=38.5e-6,
        V_plenum=3.2e-6,
    )
    return {'motor': motor, 'propellant': prop, 'pyrogen_chamber': pyc}


def test_yns_hasegawa_a_short_run_invariants():
    """Short Hasegawa A run; verify final Y invariants and that the
    species composition is sensible (mostly propellant after a burn)."""
    pytest.importorskip("numba")
    from srm_1d.openmotor_adapter import run_from_ric
    result, _perf, _nz, _geo, _prop = run_from_ric(
        'srm_1d/motors/hasegawa_a.ric',
        pyrogen='bpnv',
        pyrogen_mass=12.3e-3,
        pyrogen_throat_area=38.5e-6,
        pyrogen_volume=3.2e-6,
        T_ignition=927.0,
        roughness=37.5e-6,
        kappa=0.429,
        t_max=0.05,   # 50 ms — enough for ignition + early QSS
        cfl_target=0.3,
        snapshot_interval=0.01,
        verbose=False,
    )

    Y = result['Y_species_final']
    species_names = result['species_names']
    assert Y.shape[1] == 3, "Expected 3 species"
    assert species_names == ['BPNV_gas',
                             'Hasegawa A Prop_gas',
                             'ambient_air']

    # Y invariants per cell
    for i in range(Y.shape[0]):
        s_total = Y[i, :].sum()
        assert abs(s_total - 1.0) < 1.0e-6, (
            f"sum(Y[{i}, :])={s_total} after Hasegawa A run")
        for s in range(Y.shape[1]):
            assert -1.0e-9 <= Y[i, s] <= 1.0 + 1.0e-9, (
                f"Y[{i}, {s}]={Y[i, s]} outside [0, 1]")

    # After 50 ms: ambient should be ~purged from most cells; the
    # average ambient fraction should be small.
    avg_ambient = Y[:, 2].mean()
    assert avg_ambient < 0.5, (
        f"Average ambient Y={avg_ambient:.3f} after 50 ms — pre-fill "
        f"should be largely purged"
    )

    # Sanity: at least one cell should have non-trivial propellant mass
    # fraction (grain has ignited and burned).
    max_propellant = Y[:, 1].max()
    assert max_propellant > 0.1, (
        f"max Y_propellant={max_propellant:.3f} — grain should have "
        f"contributed mass during 50 ms"
    )
