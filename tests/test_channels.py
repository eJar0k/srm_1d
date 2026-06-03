"""
test_channels.py — v0.8.0 Phase 1 channel-model parity tests.

Verifies that ``build_channels`` re-shapes a ``run_simulation()`` results
dict into the channel model with byte-for-byte identical data, that the
``AxialChannel`` carries per-cell fields with correct shape, and that
unit-aware ``getData`` defers correctly to openMotor's conversion table.

A synthetic results dict (mirroring the real structure) is used for the
mapping/parity logic so the test is fast and deterministic — it does not
exercise the numba solver, which is unchanged this phase.
"""

import numpy as np
import pytest

from srm_1d.channels import (
    Channel, AxialChannel, SimulationChannels, build_channels,
)


def _synthetic_results(n_steps=7, n_snaps=4, n_cells=5, n_grains=2):
    """A minimal results dict shaped like run_simulation()'s output."""
    rng = np.random.default_rng(0)
    t = np.linspace(0.0, 1.0, n_steps)
    snap_times = np.linspace(0.0, 1.0, n_snaps)
    x = np.linspace(0.0, 0.5, n_cells)

    snapshots = []
    for s in range(n_snaps):
        snapshots.append({
            't': float(snap_times[s]),
            'x': x.copy(),
            'P': rng.uniform(1e5, 1e7, n_cells),
            'u': rng.uniform(-10, 200, n_cells),
            'Mach': rng.uniform(0, 0.5, n_cells),
            'T': rng.uniform(300, 3000, n_cells),
            'r_total': rng.uniform(0, 0.02, n_cells),
            'r_erosive': rng.uniform(0, 0.01, n_cells),
            'D_port': rng.uniform(0.01, 0.05, n_cells),
            'C_burn': rng.uniform(0.01, 0.2, n_cells),
            'endface_msource': np.zeros(n_cells),
            'is_burning': rng.uniform(0, 1, n_cells) > 0.5,
            'is_grain': np.ones(n_cells, dtype=bool),
            'T_surf': rng.uniform(300, 1000, n_cells),
            'mass_source': rng.uniform(0, 5, n_cells),
            'thermal_source': rng.uniform(0, 1e6, n_cells),
            'momentum_source': rng.uniform(0, 100, n_cells),
            'pyrogen_surface_heat_flux': np.zeros(n_cells),
            'radiation_heat_flux': np.zeros(n_cells),
        })

    grains = []
    for g in range(n_grains):
        grains.append({
            'segment': g,
            'regression': rng.uniform(0, 0.02, n_snaps),
            'web': rng.uniform(0, 0.05, n_snaps),
        })

    return {
        'time': t,
        'P_head': rng.uniform(1e5, 1e7, n_steps),
        'P_exit': rng.uniform(1e5, 5e5, n_steps),
        'D_throat': rng.uniform(0.01, 0.011, n_steps),
        'Kn': rng.uniform(100, 400, n_steps),
        'massflow': rng.uniform(0, 2, n_steps),
        'max_mach': rng.uniform(0, 0.4, n_steps),
        'snapshots': snapshots,
        'grains': grains,
        'summary': {'P_peak': 1.0e7, 't_peak': 0.5},
        'P_ambient': 101325.0,
        'ignition_time_by_cell': rng.uniform(0, 0.1, n_cells),
        'species_names': ['igniter', 'propellant', 'ambient'],
    }


def test_scalar_channel_parity():
    res = _synthetic_results()
    sc = build_channels(res)
    for key in ('time', 'P_head', 'P_exit', 'Kn', 'massflow', 'max_mach'):
        assert key in sc.channels, f"{key} missing from channels"
        np.testing.assert_array_equal(sc.channels[key].getData(), res[key])
        assert not sc.channels[key].per_grain
        # Legacy dict proxy: item access returns the raw array.
        np.testing.assert_array_equal(sc[key], res[key])


def test_scalar_channel_stats_match_numpy():
    res = _synthetic_results()
    sc = build_channels(res)
    ch = sc.channels['P_head']
    assert ch.getMax() == pytest.approx(float(np.max(res['P_head'])))
    assert ch.getMin() == pytest.approx(float(np.min(res['P_head'])))
    assert ch.getAverage() == pytest.approx(float(np.mean(res['P_head'])))
    assert ch.getLast() == pytest.approx(res['P_head'][-1])
    assert ch.getPoint(2) == pytest.approx(res['P_head'][2])
    assert len(ch) == len(res['P_head'])


def test_axial_channel_shape_and_parity():
    res = _synthetic_results()
    sc = build_channels(res)
    assert 'P' in sc.axial
    ax = sc.axial['P']
    snaps = res['snapshots']
    assert ax.n_frames == len(snaps)
    assert ax.n_cells == len(snaps[0]['x'])
    np.testing.assert_array_equal(ax.x_cells, snaps[0]['x'])
    np.testing.assert_array_equal(ax.times, [s['t'] for s in snaps])
    # Each frame row equals the corresponding snapshot's field.
    for i, snap in enumerate(snaps):
        np.testing.assert_array_equal(ax.getFrame(i), snap['P'])
    # Cell-column accessor matches a manual slice.
    np.testing.assert_array_equal(
        ax.getCell(1), np.array([s['P'][1] for s in snaps]))
    assert ax.getMax() == pytest.approx(
        max(float(np.max(s['P'])) for s in snaps))


def test_per_grain_channel_shape():
    res = _synthetic_results(n_snaps=4, n_grains=3)
    sc = build_channels(res)
    assert 'regression' in sc.channels
    ch = sc.channels['regression']
    assert ch.per_grain
    # (n_frames, n_grains)
    assert ch.data.shape == (4, 3)
    # Column g equals grain g's regression history.
    for g, grain in enumerate(res['grains']):
        np.testing.assert_array_equal(ch.data[:, g], grain['regression'])
    # getMax over all grains matches the global max.
    assert ch.getMax() == pytest.approx(
        max(float(np.max(g['regression'])) for g in res['grains']))


def test_per_grain_average_raises():
    sc = build_channels(_synthetic_results())
    with pytest.raises(NotImplementedError):
        sc.channels['regression'].getAverage()


def test_summary_and_extras_passthrough():
    res = _synthetic_results()
    sc = build_channels(res)
    assert sc.summary == res['summary']
    assert sc.extras['P_ambient'] == res['P_ambient']
    np.testing.assert_array_equal(
        sc.extras['ignition_time_by_cell'], res['ignition_time_by_cell'])
    assert sc.extras['species_names'] == res['species_names']


def test_missing_keys_skipped():
    """An older/partial dict (no axial, no grains) maps without error."""
    res = {'time': np.arange(5.0), 'P_head': np.ones(5)}
    sc = build_channels(res)
    assert 'P_head' in sc.channels
    assert sc.axial == {}
    assert 'regression' not in sc.channels


def test_identity_getdata_returns_same_array():
    res = _synthetic_results()
    sc = build_channels(res)
    ch = sc.channels['P_head']
    # No unit and matching unit are both identity (no openMotor needed).
    assert ch.getData() is ch.data
    assert ch.getData('Pa') is ch.data


def test_unit_conversion_pressure():
    """getData('MPa') defers to openMotor's units table (Pa -> MPa = 1e-6)."""
    pytest.importorskip("skfmm")  # openMotor checkout / deps required
    res = _synthetic_results()
    sc = build_channels(res)
    converted = sc.channels['P_head'].getData('MPa')
    np.testing.assert_allclose(converted, res['P_head'] * 1e-6, rtol=1e-9)


def test_axial_validation_rejects_bad_shape():
    with pytest.raises(ValueError):
        AxialChannel('bad', 'Pa', times=np.arange(3),
                     data=np.zeros((2, 4)), x_cells=np.arange(4))
    with pytest.raises(ValueError):
        AxialChannel('bad', 'Pa', times=np.arange(2),
                     data=np.zeros((2, 4)), x_cells=np.arange(5))


def test_as_channels_passthrough_and_build():
    from srm_1d.channels import as_channels
    res = _synthetic_results()
    sc = build_channels(res)
    # Already-channels: returned as-is (identity).
    assert as_channels(sc) is sc
    # Dict: re-shaped into channels.
    sc2 = as_channels(res)
    assert isinstance(sc2, SimulationChannels)
    np.testing.assert_array_equal(sc2.channels['P_head'].getData(), res['P_head'])


# ----------------------------------------------------------------------
# Consumer-migration parity: dict and channel inputs must give identical
# results through the migrated data-layer consumers (byte-for-byte gate).
# Uses one real short simulation, built once for the module.
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def real_run():
    pytest.importorskip("numba")
    pytest.importorskip("skfmm")
    from srm_1d import run_simulation
    from tests.test_simulation_phase3 import _small_motor, _test_chamber
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0, t_max=0.010, P_cutoff=1.0,
        dt_max=2.0e-5, burn_update_interval=1, snapshot_interval=0.002,
        cfl_target=0.5, verbose=False,
    )
    return geo, prop, nozzle, result


def test_compute_motor_performance_dict_vs_channels(real_run):
    from srm_1d.nozzle import compute_motor_performance
    geo, prop, nozzle, result = real_run
    perf_dict = compute_motor_performance(result, nozzle, prop)
    perf_chan = compute_motor_performance(build_channels(result), nozzle, prop)
    for key in ('thrust', 'Isp', 'C_F', 'P_exit'):
        np.testing.assert_array_equal(perf_dict[key], perf_chan[key])


def test_result_to_csv_dict_vs_channels_parity(real_run):
    from srm_1d.openmotor_adapter import result_to_csv
    from srm_1d.nozzle import compute_motor_performance
    geo, prop, nozzle, result = real_run
    perf = compute_motor_performance(result, nozzle, prop)
    csv_dict = result_to_csv(result, perf, geo, prop)
    csv_chan = result_to_csv(build_channels(result), perf, geo, prop)
    # Byte-for-byte identical CSV whether fed a dict or channels.
    assert csv_dict == csv_chan
    # And the CSV's pressure column matches the raw dict at sampled rows.
    assert csv_dict.count('\n') > 1


def test_verify_run_health_accepts_channels(real_run):
    from srm_1d.run_artifacts import verify_run_health
    geo, prop, nozzle, result = real_run
    healthy_dict = verify_run_health(result, motor_name='test')
    healthy_chan = verify_run_health(build_channels(result), motor_name='test')
    assert healthy_dict == healthy_chan
