"""
test_flame_front_gate.py — v0.7.4 Phase F flame-spread front gate
=================================================================

Covers the `_advance_flame_front` kernel (constant literature flame-spread
velocity, cartridge + behind-front exposure mask, propagation direction,
grid-independence) and a full-sim integration check that enabling the gate
spreads ignition in time rather than lighting the whole grain at once.

The front advances at a single bounded physical velocity
(Propellant.flame_front_velocity, ~1-10 m/s for AP/HTPB; Peretz-Kuo-Caveny-
Summerfield 1973, Kumar & Kuo 1984), decoupled from the acoustic fill.

GOTCHA #1: after editing the @njit kernels in simulation.py, delete
srm_1d/__pycache__/ (+ .nbi/.nbc) before running.
"""

import numpy as np
import pytest

from srm_1d.simulation import _advance_flame_front


DT = 1.0e-5
V_FLAME = 3.0   # m/s


def _kernel(N, x_centers, is_burning, x_front, front_direction, seed,
            cart_lo, cart_hi, v_flame=V_FLAME, dt=DT):
    ignitable = np.zeros(N, dtype=np.bool_)
    xf = _advance_flame_front(
        is_burning, x_centers, v_flame, x_front, front_direction, seed,
        cart_lo, cart_hi, dt, N, ignitable,
    )
    return xf, ignitable


def test_no_advance_before_burning_cartridge_exposed():
    """No burning cell → front does not move; only cartridge/seed exposed."""
    N = 4
    x = np.array([0.005, 0.015, 0.025, 0.035])
    is_burning = np.zeros(N, dtype=np.bool_)
    xf, ignitable = _kernel(N, x, is_burning, x[0], 1, 0, 0, 0)
    assert xf == pytest.approx(x[0])             # no advance
    assert ignitable[0]                          # cartridge/seed exposed
    assert not ignitable[1:].any()               # nothing downstream


def test_constant_velocity_advance():
    """Front advances by front_direction * v_flame * dt once a cell burns."""
    N = 4
    x = np.array([0.005, 0.015, 0.025, 0.035])
    is_burning = np.array([True, False, False, False])
    xf, _ = _kernel(N, x, is_burning, x[0], 1, 0, 0, 0)
    assert xf == pytest.approx(x[0] + V_FLAME * DT, rel=1e-12)


def test_velocity_is_grid_independent():
    """Same v_flame, dt ⇒ identical advance regardless of dx/N (it's a velocity)."""
    def advance(N, dx):
        x = (np.arange(N) + 0.5) * dx
        is_burning = np.zeros(N, dtype=np.bool_); is_burning[0] = True
        xf, _ = _kernel(N, x, is_burning, x[0], 1, 0, 0, 0)
        return xf - x[0]
    assert advance(4, 0.01) == pytest.approx(advance(8, 0.005), rel=1e-12)
    assert advance(4, 0.01) == pytest.approx(V_FLAME * DT, rel=1e-12)


def test_aft_basket_front_propagates_toward_head():
    """front_direction=-1 advances x_front toward index 0 and exposes aft cells."""
    N = 4
    x = np.array([0.005, 0.015, 0.025, 0.035])
    is_burning = np.array([False, False, False, True])
    xf, ignitable = _kernel(N, x, is_burning, x[N - 1], -1, N - 1, N - 1, N - 1)
    assert xf < x[N - 1]
    assert ignitable[N - 1]
    assert not ignitable[0]


def test_cartridge_cells_always_exposed():
    """All cells in [cart_lo, cart_hi] are ignitable regardless of the front."""
    N = 5
    x = (np.arange(N) + 0.5) * 0.01
    is_burning = np.zeros(N, dtype=np.bool_)
    # head_basket-style cartridge spanning cells 1..2, seed at 1.
    xf, ignitable = _kernel(N, x, is_burning, x[1], 1, 1, 1, 2)
    assert ignitable[1] and ignitable[2]         # cartridge exposed
    assert not ignitable[3] and not ignitable[4]  # downstream gated


# ---------------------------------------------------------------------------
# Integration: enabling the front gate must NOT collapse / starve the run and
# must spread ignition in time (a genuine fore→aft front).
# ---------------------------------------------------------------------------

def test_flame_front_spreads_ignition_in_full_sim():
    from srm_1d.openmotor_adapter import run_from_ric
    from srm_1d.tools.ignition_diagnostics import ignition_spread_metrics

    result, _perf, _nz, _geo, _prop = run_from_ric(
        'srm_1d/motors/hasegawa_a.ric',
        roughness=37.1e-6, kappa=0.45, pyrogen='bpnv', pyrogen_mass=None,
        T_ignition=850.0, P_cutoff=0.05e6,
        snapshot_interval=0.05, print_interval=20.0, verbose=False,
        flame_front_enabled=True, flame_front_velocity=3.0,
    )

    assert result['summary']['termination_code'] != 4, "run collapsed"
    ign = np.asarray(result['ignition_time_by_cell'], dtype=float)
    finite = ign[ign < 1.0e9]
    assert finite.size >= 2, "fewer than two cells ignited (possible stall)"

    spread = ignition_spread_metrics(result)
    order = np.asarray(spread['axial_ignition_order'], dtype=int)
    assert order.size >= 2
    # Genuine fore→aft front ignites cells in increasing axial index.
    assert np.all(np.diff(order) > 0), \
        f"ignition not a monotonic fore-to-aft front: {order[:12]}"
    assert spread['spread_10_90_s'] > 5.0e-4, \
        f"ignition spread implausibly tight: {spread['spread_10_90_s']}"
