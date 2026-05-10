"""
sensitivity.py — Latin Hypercube parameter sweeps with parallel execution.
==========================================================================

Lightweight, motor-agnostic LHS driver promoted from Gemini's v0.6.0
``haseOptimizer.py`` scratch script. Runs N independent simulations in
parallel via ``ProcessPoolExecutor``, computes a user-supplied scalar
fitness (default: full-trace MSE vs experimental), and persists every
run to CSV.

Usage:

    from srm_1d.tools.sensitivity import run_lhs, mse_fitness

    bounds = {
        'roughness':         (5e-6, 50e-6),
        'pyrogen_mass':      (0.001, 0.050),
        'pyrogen_throat_area': (1e-6, 5e-5),
        'T_ignition':        (700.0, 950.0),
        'kappa':             (0.30, 0.60),
    }

    results = run_lhs(
        motor_path='srm_1d/motors/hasegawa_a.ric',
        bounds=bounds,
        n_samples=500,
        fitness_fn=mse_fitness(experimental_t, experimental_p, t_min=0.01),
        csv_path='hasegawa_a_lhs.csv',
        seed=42,
        # Locked sim kwargs (passed verbatim to run_simulation per call):
        kappa=0.45, t_max=6.0, P_cutoff=0.05e6,
        snapshot_interval=2.0, print_interval=20.0,
    )

    # Sort by fitness, take top-K, ...
    top = sorted(results, key=lambda r: r['fitness'])[:5]

The default ``mse_fitness`` factory builds a closure over the
experimental trace; alternative fitness functions (impulse error,
spike-region MSE, peak-P error) live below as factories with the
same ``(result) -> float`` shape.
"""

import csv
import concurrent.futures
from typing import Callable, Dict, Tuple

import numpy as np
from scipy.stats import qmc

from srm_1d.openmotor_adapter import run_from_ric


# ================================================================
# Fitness functors (picklable for ProcessPoolExecutor)
# ================================================================
#
# Note: Windows uses spawn for multiprocessing, which pickles the
# fitness function before sending it to workers. Closures over numpy
# arrays don't pickle. Instead, fitness is implemented as top-level
# callable classes — instances pickle cleanly.

class MSEFitness:
    """Mean-squared error in MPa² between simulated and experimental traces.

    The simulated head-end pressure (``result['P_head']``, Pa) is
    interpolated onto the experimental time grid; samples with
    ``t_exp < t_min`` are dropped (typically used to skip the ignition
    transient that the LHS is *trying* to fit).
    """
    def __init__(self, t_exp, p_exp_mpa, t_min=0.0):
        self.t_exp = np.asarray(t_exp)
        self.p_exp_mpa = np.asarray(p_exp_mpa)
        self.t_min = float(t_min)

    def __call__(self, result):
        t_sim = result['time']
        p_sim_mpa = result['P_head'] / 1e6
        if len(t_sim) < 100 or t_sim[-1] < self.t_min + 1e-3:
            return 1e6
        p_sim_at_exp = np.interp(self.t_exp, t_sim, p_sim_mpa)
        mask = self.t_exp >= self.t_min
        return float(np.mean((p_sim_at_exp[mask] - self.p_exp_mpa[mask]) ** 2))


class ImpulseErrorFitness:
    """Absolute relative error in total impulse vs target [N·s].

    Note: ``result['summary']`` doesn't include impulse by default —
    pair this with a wrapper that injects it from
    ``compute_motor_performance``.
    """
    def __init__(self, impulse_target_n_s):
        self.target = float(impulse_target_n_s)

    def __call__(self, result):
        s = result.get('summary', {})
        I = s.get('impulse', None)
        if I is None:
            return 1e6
        return abs(I - self.target) / self.target


class PeakPressureErrorFitness:
    """Absolute relative error in peak head-end pressure [Pa]."""
    def __init__(self, p_peak_target_pa):
        self.target = float(p_peak_target_pa)

    def __call__(self, result):
        s = result.get('summary', {})
        Pp = s.get('P_peak', None)
        if Pp is None:
            return 1e6
        return abs(Pp - self.target) / self.target


# Backwards-friendly factory aliases (instantiate the class directly)
def mse_fitness(t_exp, p_exp_mpa, t_min=0.0):
    return MSEFitness(t_exp, p_exp_mpa, t_min=t_min)


def impulse_error_fitness(impulse_target_n_s):
    return ImpulseErrorFitness(impulse_target_n_s)


def peak_pressure_error_fitness(p_peak_target_pa):
    return PeakPressureErrorFitness(p_peak_target_pa)


# ================================================================
# Worker (must be top-level for ProcessPoolExecutor)
# ================================================================

def _run_one(args: tuple):
    """Worker entrypoint. Importable by pickling."""
    idx, params, motor_path, sim_kwargs, fitness_fn = args
    try:
        result, _perf, _noz, _geo, _prop = run_from_ric(
            motor_path,
            **{**sim_kwargs, **params},
        )
    except Exception as exc:
        return idx, params, 1e6, str(exc)
    fitness = float(fitness_fn(result))
    return idx, params, fitness, None


# ================================================================
# Public LHS driver
# ================================================================

def run_lhs(
    motor_path: str,
    bounds: Dict[str, Tuple[float, float]],
    n_samples: int,
    fitness_fn: Callable,
    n_workers: int = None,
    seed: int = 42,
    csv_path: str = None,
    progress_every: int = 25,
    **sim_kwargs,
) -> list:
    """
    Run a Latin Hypercube sweep of ``run_from_ric`` against a fitness
    function, in parallel.

    Parameters
    ----------
    motor_path : str
        Path to a .ric motor file (transport.yaml auto-resolved).
    bounds : dict[str, (low, high)]
        Maps parameter names to ``(min, max)`` ranges. Names must be
        kwargs accepted by ``run_from_ric`` (e.g. ``roughness``,
        ``pyrogen_mass``, ``pyrogen_throat_area``, ``T_ignition``,
        ``kappa``).
    n_samples : int
        Number of LHS samples to draw and evaluate.
    fitness_fn : callable
        ``(result_dict) -> float`` (lower is better).
    n_workers : int or None
        Worker process count; defaults to ``os.cpu_count()``.
    seed : int
        RNG seed for the LHS sampler — locks reproducibility.
    csv_path : str or None
        If given, write all runs to a CSV (one row per sample, columns =
        bounds keys + ``fitness`` + ``error`` for any worker exceptions).
    progress_every : int
        Print a one-line progress message every N completed samples.
    **sim_kwargs
        Locked simulation kwargs passed verbatim to ``run_from_ric``
        (e.g. ``kappa=0.45``, ``t_max=6.0``).

    Returns
    -------
    list of dicts, one per sample, each with the parameter values plus
    ``fitness`` (and ``error`` if the worker raised).
    """
    keys = list(bounds.keys())
    l_bounds = [bounds[k][0] for k in keys]
    u_bounds = [bounds[k][1] for k in keys]

    sampler = qmc.LatinHypercube(d=len(keys), seed=seed)
    raw = sampler.random(n=n_samples)
    scaled = qmc.scale(raw, l_bounds, u_bounds)

    sim_kwargs = dict(sim_kwargs)
    sim_kwargs.setdefault('pyrogen', 'bpnv')

    work = [
        (i, dict(zip(keys, row.tolist())), motor_path, sim_kwargs, fitness_fn)
        for i, row in enumerate(scaled)
    ]

    print(f"sensitivity.run_lhs: {n_samples} samples, "
          f"{len(keys)} dims, motor={motor_path}")

    rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as exe:
        for i, (idx, params, fitness, err) in enumerate(exe.map(_run_one, work)):
            rec = dict(params)
            rec['fitness'] = fitness
            if err is not None:
                rec['error'] = err
            rows.append(rec)
            if (i + 1) % progress_every == 0:
                print(f"  {i+1}/{n_samples} done")

    if csv_path is not None:
        fieldnames = keys + ['fitness', 'error']
        with open(csv_path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for rec in rows:
                # Ensure 'error' column exists for every row
                rec.setdefault('error', '')
                w.writerow(rec)
        print(f"sensitivity.run_lhs: wrote {csv_path}")

    return rows
