"""History pre-allocation (solver-verification T8 finding, 2026-10-06).

run_simulation sizes its per-step history arrays from an estimated step
count. The former estimate assumed a fixed 1000 m/s wave speed, below the
~1100-1200 m/s hot-gas sound speed, so fine-grid / low-CFL runs filled the
history and stopped early (Chunc 806 cells @ CFL 0.15 at t ~ 0.83 s). The
estimate now uses the hottest species' sound speed times a velocity margin,
``history_capacity`` sets the capacity exactly, and a truncated run warns.
"""
import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.propellant import R_UNIVERSAL
from srm_1d.simulation import _estimate_history_steps, _HISTORY_WAVE_SPEED_MARGIN

ROOT = Path(__file__).resolve().parents[1]
CANON = dict(roughness=32e-6, kappa=0.44, T_ignition=756.0, pyrogen='bpnv',
             pyrogen_mass=None, print_interval=0.0, verbose=False)


def _committed_chunc():
    """The committed Chunc .ric (the working copy may carry local edits)."""
    blob = subprocess.run(
        ['git', 'show', 'HEAD:motors/machbusterNew.ric'], cwd=ROOT,
        capture_output=True, check=True).stdout
    f = tempfile.NamedTemporaryFile(suffix='.ric', delete=False)
    f.write(blob)
    f.close()
    return f.name


def test_estimate_uses_hottest_species_sound_speed():
    # (gamma, Cp, M, T_flame): a cold ambient species and a hot gas.
    species = np.array([[1.40, 1005.0, 0.02897, 300.0],
                        [1.20, 2000.0, 0.025, 3000.0]])
    a_hot = (1.20 * R_UNIVERSAL / 0.025 * 3000.0) ** 0.5
    dt = 0.3 * 1e-3 / (_HISTORY_WAVE_SPEED_MARGIN * a_hot)
    assert _estimate_history_steps(1.0, 0.3, 1e-3, 0.002, species) == \
        int(1.0 / dt) + 1000


def test_estimate_respects_dt_max():
    species = np.array([[1.20, 2000.0, 0.025, 3000.0]])
    # CFL step here is 0.3*1e-3/(1.5*~1095 m/s) ~ 1.8e-7 s; a dt_max below
    # it (1e-8 s) must take over the estimate.
    assert _estimate_history_steps(1.0, 0.3, 1e-3, 1e-8, species) == \
        int(1.0 / 1e-8) + 1000


@pytest.fixture(scope='module')
def chunc_short():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return run_from_ric(_committed_chunc(), t_max=0.3, **CANON)[0]


def test_estimate_covers_actual_steps_and_old_one_did_not(chunc_short):
    s = chunc_short['summary']
    assert s['termination'] == 't_max reached'
    assert s['steps'] <= s['history_estimated_steps']
    # The former fixed-1000 m/s estimate under-counts this same run (the
    # 5M-row floor used to hide that at 100 cells).
    dx = float(chunc_short['dx'])
    old_estimate = int(0.3 / (0.3 * dx / 1000.0)) + 1000
    assert old_estimate < s['steps']


def test_explicit_capacity_truncates_with_warning():
    with pytest.warns(RuntimeWarning, match='history_capacity='):
        r = run_from_ric(_committed_chunc(), t_max=0.05, history_capacity=500,
                         **CANON)[0]
    assert r['summary']['termination'] == 'history array full'
    assert r['summary']['history_capacity'] == 500
    assert len(r['time']) == 500


@pytest.mark.parametrize('bad', [0, -5, 1.5])
def test_invalid_capacity_raises(bad):
    with pytest.raises(ValueError, match='history_capacity'):
        run_from_ric(_committed_chunc(), t_max=0.01, history_capacity=bad,
                     **CANON)
