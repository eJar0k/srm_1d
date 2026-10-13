"""
test_sm2_replay.py -- SM2: ignition-event diagnostics, the open-loop forcing recorder and the replay tool.

The recorder / diagnostics are bookkeeping only (the trajectory must not change), and replaying a run's own recorded
forcing through the same solid model must reproduce that run's per-cell ignition times exactly.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from srm_1d.tools import solid_replay as rp  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def _run(model, rec=0):
    from srm_1d.openmotor_adapter import run_from_ric
    return run_from_ric(os.path.join(ROOT, "motors", "machbusterNew.ric"), pyrogen="bpnv", pyrogen_mass=4.5e-3,
                        injection_topology="head_basket", roughness=32e-6, kappa=0.44, T_ignition=756.0,
                        t_max=0.012, print_interval=0.0, verbose=False, solid_heating_model=model,
                        record_forcing_steps=rec)


def test_recorder_does_not_change_the_trajectory():
    a = _run("goodman", 0)[0]
    b = _run("goodman", 20000)[0]
    assert np.array_equal(np.asarray(a["P_head"]), np.asarray(b["P_head"]))
    assert np.array_equal(np.asarray(a["ignition_time_by_cell"]), np.asarray(b["ignition_time_by_cell"]))
    assert a["forcing"] is None and b["forcing"] is not None


def test_ignition_diag_filled_for_lit_cells_only():
    r = _run("goodman", 0)[0]
    ig = np.asarray(r["ignition_time_by_cell"])
    d = np.asarray(r["ignition_diag"])
    lit = ig < 1e9
    assert lit.any()
    assert np.all(np.isfinite(d[lit, 3]))
    assert np.all(d[lit, 3] > 756.0)                       # T_s at ignition above T_ignition
    assert np.allclose(d[lit, 8], ig[lit])                 # time column = ignition time
    assert np.all(d[lit, 2] > 0.0)                         # lit by a positive flux
    assert np.all(np.isnan(d[~lit]))


@pytest.mark.parametrize("model", ["goodman", "peretz", "conduction"])
def test_self_replay_reproduces_closed_loop_ignition_times(model):
    out = _run(model, 20000)
    r, prop = out[0], out[4]
    fo = r["forcing"]
    assert fo["h"].shape[0] > 1000
    ig = np.asarray(r["ignition_time_by_cell"])
    rho = float(prop.rho_propellant)
    o = rp.replay(fo, model, rho=rho)
    window = float(fo["t"][-1] + fo["dt"][-1])
    ok = (ig < window) & (o[:, 0] < 1e9)
    assert ok.sum() > 10
    assert np.abs(o[ok, 0] - ig[ok]).max() < 1e-9


def test_lighting_metrics():
    t = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, rp.NEVER])
    m = rp.lighting_metrics(t)
    assert m["n_never"] == 1 and m["first"] == 1.0 and m["t50"] == 5.0 and m["t90"] == 9.0
    assert m["dur90"] == 8.0
