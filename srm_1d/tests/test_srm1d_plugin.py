"""
test_srm1d_plugin.py — v0.8.0 Phase 5: srm_1d openMotor solver plugin (D1).

Verifies that importing the plugin registers srm_1d's transient solver with
openMotor's solver registry, that the built-in quasi-steady solver remains
registered, and that the registry can select and run the srm_1d solver
headlessly on a canonical openMotor Motor — returning a populated openMotor
SimulationResult.
"""

import os

import pytest

pytest.importorskip("skfmm")
pytest.importorskip("numba")

from srm_1d.fmm_grain import _setup_openmotor_path  # noqa: E402

MOTORS = os.path.join(os.path.dirname(__file__), '..', 'motors')


@pytest.fixture(scope="module")
def om():
    _setup_openmotor_path()
    import srm_1d.srm1d_plugin  # noqa: F401  (registers on import)
    from motorlib import solvers, motor as om_motor  # type: ignore
    return solvers, om_motor


def _canonical_motor(om_motor):
    from srm_1d.openmotor_adapter import load_ric
    motor_dict = load_ric(os.path.join(MOTORS, 'hasegawa_a.ric'))
    return om_motor.Motor(motor_dict)


def test_solver_registered_alongside_quasi_steady(om):
    solvers, _ = om
    names = solvers.list_solvers()
    assert 'srm_1d-transient' in names
    assert solvers.QUASI_STEADY in names
    plugin = solvers.get_solver('srm_1d-transient')
    assert plugin.capabilities['transient'] is True
    assert plugin.capabilities['axial_fields'] is True


def test_quasi_steady_solver_preserved(om):
    """The QS solver wraps Motor.runSimulation — coexistence (D6)."""
    solvers, om_motor = om
    qs = solvers.get_solver(solvers.QUASI_STEADY)
    assert qs.capabilities['transient'] is False
    motor = _canonical_motor(om_motor)
    # QS solver runs the built-in path and returns a SimulationResult.
    sr = qs.simulate(motor)
    assert sr is not None


def test_srm1d_solver_runs_via_registry(om):
    solvers, om_motor = om
    motor = _canonical_motor(om_motor)
    solver = solvers.get_solver('srm_1d-transient')
    sr = solver.simulate(motor, config={'t_max': 0.03, 'P_cutoff': 1.0})

    # A populated, successful openMotor SimulationResult.
    assert sr.success
    t = sr.channels['time'].getData()
    assert len(t) > 0
    assert sr.getMaxPressure() > 1e5          # produced real pressure
    assert sr.channels['force'].getMax() >= 0.0
    # dThroat channel starts at zero (relative to initial throat).
    assert sr.channels['dThroat'].getData()[0] == pytest.approx(0.0)


def test_per_grain_channels_populated(om):
    """Phase 6 task 3: the per-grain multi-value channels are populated,
    share the per-step time base, and carry one value per grain."""
    solvers, om_motor = om
    motor = _canonical_motor(om_motor)
    solver = solvers.get_solver('srm_1d-transient')
    sr = solver.simulate(motor, config={'t_max': 0.03, 'P_cutoff': 1.0})

    n_steps = len(sr.channels['time'].getData())
    n_grains = len(motor.grains)
    for name in ('mass', 'massFlow', 'massFlux', 'regression', 'web',
                 'machNumber'):
        data = sr.channels[name].getData()
        assert len(data) == n_steps, f"{name} length != time length"
        assert all(len(frame) == n_grains for frame in data), \
            f"{name} not one value per grain"

    # Sanity on physical content: initial propellant mass is positive, and
    # the motor consumes mass over the burn (mass channel is non-increasing
    # in total, web regresses, mass flow stays non-negative).
    assert sr.getPropellantMass() > 0.0
    assert sr.channels['web'].getData()[0][0] >= sr.channels['web'].getLast()[0]
    assert sr.channels['massFlow'].getMin() >= 0.0
    assert sr.getISP() > 0.0                   # needs the 'mass' channel

    # volumeLoading is a per-step scalar the GUI reads at index 0 (motor
    # stats panel) — it must be populated, in [0, 100], and start positive.
    vol = sr.channels['volumeLoading'].getData()
    assert len(vol) == n_steps
    assert 0.0 < sr.getVolumeLoading() <= 100.0


def test_decimate_indices_preserves_peaks_and_endpoints():
    """The GUI decimation keeps the sample count under the cap while always
    retaining the first, last, peak-pressure and peak-thrust samples."""
    import numpy as np
    from srm_1d.srm1d_plugin import _decimate_indices

    n = 50000
    p_head = np.zeros(n); p_head[12345] = 99.0      # unique pressure peak
    thrust = np.zeros(n); thrust[45678] = 88.0      # unique thrust peak
    idx = _decimate_indices(n, p_head, thrust, max_points=5000)

    assert len(idx) <= 5000 + 4                      # cap + the forced samples
    assert np.all(np.diff(idx) > 0)                  # sorted, unique
    assert idx[0] == 0 and idx[-1] == n - 1
    assert 12345 in idx and 45678 in idx
    # Decimated peak equals the true peak (the spike sample is retained).
    assert p_head[idx].max() == 99.0
    assert thrust[idx].max() == 88.0

    # Under the cap, every sample is kept (no decimation).
    assert len(_decimate_indices(100, p_head[:100], thrust[:100], 5000)) == 100


def test_progress_callback_driven_to_completion(om):
    """Phase 6: the @njit loop publishes live progress; the plugin's poller
    forwards it to the callback, which lands at 1.0 on a successful finish."""
    from srm_1d.srm1d_plugin import simulate_motor
    solvers, om_motor = om
    motor = _canonical_motor(om_motor)

    seen = []
    def cb(progress):
        seen.append(progress)
        return False  # don't cancel

    sr = simulate_motor(motor, callback=cb, t_max=0.03, P_cutoff=1.0)
    assert sr.success
    assert len(seen) > 0
    assert all(0.0 <= p <= 1.0 for p in seen)
    assert seen[-1] == pytest.approx(1.0)  # final tick after the loop ends
    # The bar never regresses (the loop enforces monotonic progress).
    assert all(b >= a for a, b in zip(seen, seen[1:]))


def test_progress_callback_cancels_run(om):
    """A truthy callback return sets the cancel flag the loop reads each step;
    the run aborts cooperatively and is reported as not-successful."""
    from srm_1d.srm1d_plugin import simulate_motor
    solvers, om_motor = om
    motor = _canonical_motor(om_motor)

    # Always request cancel. t_max is large so the loop is still running when
    # the first poll fires (~50 ms), making the cancel deterministic.
    sr = simulate_motor(motor, callback=lambda p: True, t_max=3.0, P_cutoff=1.0)
    assert not sr.success
    assert sr.motor is motor  # a real (if partial) SimulationResult


def test_motor_round_trips_igniter_block(om):
    """Phase 6 task 4: openMotor's Motor carries the data.igniter block
    through getDict/applyDict (closes the Phase 4 caveat)."""
    solvers, om_motor = om
    motor = _canonical_motor(om_motor)

    block = motor.getDict()['igniter']
    assert 'pyrogen' in block
    assert 'injection_topology' in block

    # Mutate, round-trip through a fresh Motor, and confirm it survives.
    block['pyrogen']['name'] = 'TestPyro'
    block['injection_topology'] = 'head_basket'
    d = motor.getDict()
    d['igniter'] = block
    motor2 = om_motor.Motor(d)
    out = motor2.getDict()['igniter']
    assert out['pyrogen']['name'] == 'TestPyro'
    assert out['injection_topology'] == 'head_basket'


def test_motor_without_igniter_keeps_defaults(om):
    """A motor dict lacking the igniter block loads with default igniter."""
    solvers, om_motor = om
    motor = _canonical_motor(om_motor)
    d = motor.getDict()
    d.pop('igniter', None)
    motor2 = om_motor.Motor(d)
    # Default igniter is present and self-describing (BPNV / forward_plenum).
    block = motor2.getDict()['igniter']
    assert block['injection_topology'] == 'forward_plenum'
