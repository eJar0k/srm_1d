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
import time
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

def _windowed_peak_time(t, p, window):
    """Return the time of the maximum-pressure sample inside ``window``.

    Returns ``np.nan`` if the window contains no samples or the input
    is empty. Used by the peak-alignment hook in ``MSEFitness`` and
    ``PressureTraceMetrics`` to compute a per-sample ignition-timing
    offset before MSE evaluation.

    Parameters
    ----------
    t, p : np.ndarray
        Time and pressure arrays of equal length. Pressure can be in
        any consistent unit (Pa or MPa) — only ``argmax`` is used.
    window : tuple(float, float)
        ``(t_lo, t_hi)`` inclusive bounds in the same time units as t.
    """
    t = np.asarray(t)
    p = np.asarray(p)
    if t.size == 0 or window is None:
        return float('nan')
    lo, hi = float(window[0]), float(window[1])
    mask = (t >= lo) & (t <= hi)
    if not np.any(mask):
        return float('nan')
    t_win = t[mask]
    p_win = p[mask]
    return float(t_win[int(np.argmax(p_win))])


def _peak_alignment_offset(t_sim, p_sim, t_exp, p_exp, peak_align_window):
    """Compute ``t_offset`` so that shifting t_sim by it aligns the sim
    and experimental peak-pressure samples inside ``peak_align_window``.

    Returns 0.0 (no shift) if alignment cannot be computed — empty
    arrays, no samples in window on either trace, or NaN peak times.
    """
    if peak_align_window is None:
        return 0.0
    t_sim_peak = _windowed_peak_time(t_sim, p_sim, peak_align_window)
    t_exp_peak = _windowed_peak_time(t_exp, p_exp, peak_align_window)
    if not (np.isfinite(t_sim_peak) and np.isfinite(t_exp_peak)):
        return 0.0
    return t_exp_peak - t_sim_peak


class MSEFitness:
    """Mean-squared error in MPa² between simulated and experimental traces.

    The simulated head-end pressure (``result['P_head']``, Pa) is
    interpolated onto the experimental time grid; samples with
    ``t_exp < t_min`` are dropped (typically used to skip the ignition
    transient that the LHS is *trying* to fit).

    Parameters
    ----------
    t_exp, p_exp_mpa : np.ndarray
        Experimental head-end pressure trace.
    t_min : float
        Drop samples before this time from the MSE sum.
    peak_align_window : tuple(float, float) or None, optional
        When set, find the sim and experimental peak times inside this
        window each evaluation, then shift the sim time array by
        ``t_exp_peak - t_sim_peak`` BEFORE interpolating. This removes
        ignition-timing residual from the MSE so the optimizer fits
        trace SHAPE rather than ignition PHASING — useful when the
        igniter model has known timing drift (e.g. post-Phase-3.5
        Hasegawa A, where the pyrogen-to-surface sensible-power cap
        changed the ignition delay). Default ``None`` preserves the
        original behavior.
    """
    def __init__(self, t_exp, p_exp_mpa, t_min=0.0, peak_align_window=None):
        self.t_exp = np.asarray(t_exp)
        self.p_exp_mpa = np.asarray(p_exp_mpa)
        self.t_min = float(t_min)
        if peak_align_window is None:
            self.peak_align_window = None
        else:
            self.peak_align_window = (
                float(peak_align_window[0]), float(peak_align_window[1])
            )

    def __call__(self, result):
        t_sim = np.asarray(result['time'])
        p_sim_mpa = np.asarray(result['P_head']) / 1e6
        if len(t_sim) < 100 or t_sim[-1] < self.t_min + 1e-3:
            return 1e6
        t_offset = _peak_alignment_offset(
            t_sim, p_sim_mpa, self.t_exp, self.p_exp_mpa,
            self.peak_align_window,
        )
        t_sim_aligned = t_sim + t_offset
        p_sim_at_exp = np.interp(self.t_exp, t_sim_aligned, p_sim_mpa)
        mask = self.t_exp >= self.t_min
        if not np.any(mask):
            return 1e6
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


DEFAULT_PRESSURE_SEGMENTS = (
    ('spike', 0.03, 0.12),
    ('post_spike', 0.12, 0.60),
    ('plateau', 0.60, 2.45),
    ('taildown', 2.45, 4.75),
)


class PressureTraceMetrics:
    """Named pressure-trace metrics for calibration diagnostics.

    Pressures are reported in MPa-based units. Segment MSE columns are
    ``mse_<segment>`` in MPa^2 and are computed after interpolating the
    simulation onto the experimental timestamps.
    """
    def __init__(self, t_exp, p_exp_mpa, t_min=0.0,
                 segments=DEFAULT_PRESSURE_SEGMENTS,
                 peak_window=(0.03, 0.18),
                 trough_window=(0.12, 0.60),
                 peak_align_window=None):
        self.t_exp = np.asarray(t_exp, dtype=float)
        self.p_exp_mpa = np.asarray(p_exp_mpa, dtype=float)
        self.t_min = float(t_min)
        self.segments = tuple((str(name), float(lo), float(hi))
                              for name, lo, hi in segments)
        self.peak_window = tuple(float(v) for v in peak_window)
        self.trough_window = tuple(float(v) for v in trough_window)
        # v0.7.1 Phase 5: optional per-sample ignition-time alignment.
        # When set, find sim/exp peaks in this window and shift t_sim
        # by t_exp_peak - t_sim_peak BEFORE interpolation + segment MSE
        # evaluation. None preserves original (no-shift) behavior.
        if peak_align_window is None:
            self.peak_align_window = None
        else:
            self.peak_align_window = (
                float(peak_align_window[0]), float(peak_align_window[1])
            )

    def _window_peak(self, t, p, window):
        mask = (t >= window[0]) & (t <= window[1])
        if not np.any(mask):
            return np.nan, np.nan
        t_win = t[mask]
        p_win = p[mask]
        idx = int(np.argmax(p_win))
        return float(t_win[idx]), float(p_win[idx])

    def _window_trough(self, t, p, window):
        mask = (t >= window[0]) & (t <= window[1])
        if not np.any(mask):
            return np.nan, np.nan
        t_win = t[mask]
        p_win = p[mask]
        idx = int(np.argmin(p_win))
        return float(t_win[idx]), float(p_win[idx])

    def __call__(self, result):
        t_sim = np.asarray(result.get('time', []), dtype=float)
        p_sim = np.asarray(result.get('P_head', []), dtype=float) / 1e6
        metrics = {}

        if len(t_sim) < 2:
            return {'mse_all': 1e6}

        # v0.7.1 Phase 5: optional ignition-timing alignment. The shift
        # is applied to t_sim BEFORE the interpolation step so all
        # downstream segment MSEs see the aligned trace. The applied
        # offset is exposed in metrics for CSV / diagnostic capture.
        t_offset = _peak_alignment_offset(
            t_sim, p_sim, self.t_exp, self.p_exp_mpa,
            self.peak_align_window,
        )
        t_sim_aligned = t_sim + t_offset
        metrics['t_offset_applied_s'] = t_offset

        p_sim_at_exp = np.interp(self.t_exp, t_sim_aligned, p_sim)
        all_mask = self.t_exp >= self.t_min
        if np.any(all_mask):
            residual = p_sim_at_exp[all_mask] - self.p_exp_mpa[all_mask]
            metrics['mse_all'] = float(np.mean(residual ** 2))
            metrics['mae_all'] = float(np.mean(np.abs(residual)))
        else:
            metrics['mse_all'] = 1e6
            metrics['mae_all'] = 1e3

        for name, lo, hi in self.segments:
            mask = (self.t_exp >= lo) & (self.t_exp <= hi)
            if np.any(mask):
                err = p_sim_at_exp[mask] - self.p_exp_mpa[mask]
                metrics[f'mse_{name}'] = float(np.mean(err ** 2))
                metrics[f'mae_{name}'] = float(np.mean(np.abs(err)))
                metrics[f'bias_{name}'] = float(np.mean(err))
            else:
                metrics[f'mse_{name}'] = np.nan
                metrics[f'mae_{name}'] = np.nan
                metrics[f'bias_{name}'] = np.nan

        # Peak/trough diagnostics use the ALIGNED sim trace so reported
        # peak_error_pct / trough_error_pct compare like-for-like.
        t_peak_sim, p_peak_sim = self._window_peak(
            t_sim_aligned, p_sim, self.peak_window,
        )
        t_peak_exp, p_peak_exp = self._window_peak(
            self.t_exp, self.p_exp_mpa, self.peak_window,
        )
        metrics['t_peak_sim'] = t_peak_sim
        metrics['P_peak_sim_MPa'] = p_peak_sim
        metrics['t_peak_exp'] = t_peak_exp
        metrics['P_peak_exp_MPa'] = p_peak_exp
        if np.isfinite(p_peak_exp) and abs(p_peak_exp) > 1e-12:
            metrics['peak_error_pct'] = float(
                100.0 * (p_peak_sim - p_peak_exp) / p_peak_exp
            )
        else:
            metrics['peak_error_pct'] = np.nan

        t_trough_sim, p_trough_sim = self._window_trough(
            t_sim_aligned, p_sim, self.trough_window,
        )
        t_trough_exp, p_trough_exp = self._window_trough(
            self.t_exp, self.p_exp_mpa, self.trough_window,
        )
        metrics['t_trough_sim'] = t_trough_sim
        metrics['P_trough_sim_MPa'] = p_trough_sim
        metrics['t_trough_exp'] = t_trough_exp
        metrics['P_trough_exp_MPa'] = p_trough_exp
        if np.isfinite(p_trough_exp) and abs(p_trough_exp) > 1e-12:
            metrics['trough_error_pct'] = float(
                100.0 * (p_trough_sim - p_trough_exp) / p_trough_exp
            )
        else:
            metrics['trough_error_pct'] = np.nan

        summary = result.get('summary', {})
        metrics['t_burn_sim'] = float(summary.get(
            't_burn', t_sim[-1] if len(t_sim) else np.nan,
        ))
        metrics['pyrogen_duration_ms'] = float(
            summary.get('pyrogen_duration', np.nan)
        ) * 1000.0
        metrics['pyrogen_peak_P_MPa'] = float(
            summary.get('pyrogen_peak_P', np.nan)
        ) / 1e6
        metrics['pyrogen_mass_burned_g'] = float(
            summary.get('pyrogen_mass_burned', np.nan)
        ) * 1000.0
        return metrics


class SegmentedPressureFitness:
    """Weighted average of segment MSE metrics.

    The score remains in MPa^2 when weights sum to one. The default
    weighting deliberately gives the post-spike shoulder a separate vote
    instead of allowing plateau/tail duration to dominate the scalar MSE.
    """
    def __init__(self, t_exp, p_exp_mpa, t_min=0.0,
                 segments=DEFAULT_PRESSURE_SEGMENTS, weights=None,
                 peak_align_window=None):
        self.metrics = PressureTraceMetrics(
            t_exp, p_exp_mpa, t_min=t_min, segments=segments,
            peak_align_window=peak_align_window,
        )
        if weights is None:
            weights = {
                'mse_spike': 0.25,
                'mse_post_spike': 0.35,
                'mse_plateau': 0.20,
                'mse_taildown': 0.20,
            }
        self.weights = dict(weights)

    def __call__(self, result):
        metrics = self.metrics(result)
        score = 0.0
        weight_sum = 0.0
        for key, weight in self.weights.items():
            value = metrics.get(key, np.nan)
            if np.isfinite(value):
                score += float(weight) * float(value)
                weight_sum += float(weight)
        if weight_sum <= 0.0:
            return 1e6
        return float(score / weight_sum)


# Backwards-friendly factory aliases (instantiate the class directly)
def mse_fitness(t_exp, p_exp_mpa, t_min=0.0, peak_align_window=None):
    return MSEFitness(
        t_exp, p_exp_mpa, t_min=t_min, peak_align_window=peak_align_window,
    )


def impulse_error_fitness(impulse_target_n_s):
    return ImpulseErrorFitness(impulse_target_n_s)


def peak_pressure_error_fitness(p_peak_target_pa):
    return PeakPressureErrorFitness(p_peak_target_pa)


def pressure_trace_metrics(t_exp, p_exp_mpa, t_min=0.0,
                           segments=DEFAULT_PRESSURE_SEGMENTS,
                           peak_window=(0.03, 0.18),
                           trough_window=(0.12, 0.60),
                           peak_align_window=None):
    return PressureTraceMetrics(
        t_exp, p_exp_mpa, t_min=t_min, segments=segments,
        peak_window=peak_window, trough_window=trough_window,
        peak_align_window=peak_align_window,
    )


def segmented_pressure_fitness(t_exp, p_exp_mpa, t_min=0.0,
                               segments=DEFAULT_PRESSURE_SEGMENTS,
                               weights=None,
                               peak_align_window=None):
    return SegmentedPressureFitness(
        t_exp, p_exp_mpa, t_min=t_min, segments=segments, weights=weights,
        peak_align_window=peak_align_window,
    )


# ================================================================
# Worker (must be top-level for ProcessPoolExecutor)
# ================================================================

def _run_one(args: tuple):
    """Worker entrypoint. Importable by pickling."""
    idx, params, motor_path, sim_kwargs, fitness_fn, metrics_fn = args
    try:
        result, _perf, _noz, _geo, _prop = run_from_ric(
            motor_path,
            **{**sim_kwargs, **params},
        )
    except Exception as exc:
        return idx, params, 1e6, {}, str(exc)
    fitness = float(fitness_fn(result))
    metrics = {}
    if metrics_fn is not None:
        metrics = dict(metrics_fn(result))
    return idx, params, fitness, metrics, None


# ================================================================
# Public LHS driver
# ================================================================

def run_lhs(
    motor_path: str,
    bounds: Dict[str, Tuple[float, float]],
    n_samples: int,
    fitness_fn: Callable,
    metrics_fn: Callable = None,
    n_workers: int = None,
    seed: int = 42,
    csv_path: str = None,
    progress_every: int = 25,
    progress_mode: str = 'brief',
    sim_verbose: bool = False,
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
    metrics_fn : callable or None
        Optional ``(result_dict) -> dict``. Returned keys are added to
        each result row and CSV output.
    n_workers : int or None
        Worker process count; defaults to ``os.cpu_count()``.
    seed : int
        RNG seed for the LHS sampler — locks reproducibility.
    csv_path : str or None
        If given, write all runs to a CSV (one row per sample, columns =
        bounds keys + ``fitness`` + ``error`` for any worker exceptions).
    progress_every : int
        Print a one-line progress message every N completed samples.
    progress_mode : {'brief', 'verbose', 'none'}
        Controls LHS progress output. ``brief`` updates a compact status
        line, ``verbose`` prints normal progress lines, and ``none``
        suppresses LHS progress output.
    sim_verbose : bool
        Passed to ``run_from_ric`` as ``verbose`` unless already supplied
        in ``sim_kwargs``. Defaults to False so worker simulations do not
        print setup/summary blocks.
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
    sim_kwargs.setdefault('verbose', sim_verbose)

    work = [
        (
            i,
            dict(zip(keys, row.tolist())),
            motor_path,
            sim_kwargs,
            fitness_fn,
            metrics_fn,
        )
        for i, row in enumerate(scaled)
    ]

    progress_mode = str(progress_mode).lower()
    if progress_mode not in {'brief', 'verbose', 'none'}:
        raise ValueError("progress_mode must be 'brief', 'verbose', or 'none'")
    progress_every = max(1, int(progress_every))
    start = time.time()
    best = float('inf')
    n_errors = 0

    def emit_progress(done, final=False):
        if progress_mode == 'none':
            return
        elapsed = time.time() - start
        msg = (
            f"LHS {done}/{n_samples}  best={best:.4g}  "
            f"errors={n_errors}  elapsed={elapsed:.0f}s"
        )
        if progress_mode == 'brief':
            print('\r' + msg, end='' if not final else '\n', flush=True)
        else:
            print('  ' + msg)

    if progress_mode != 'none':
        workers = n_workers if n_workers is not None else 'auto'
        print(
            f"sensitivity.run_lhs: {n_samples} samples, {len(keys)} dims, "
            f"workers={workers}, sim_verbose={sim_kwargs.get('verbose', True)}"
        )

    # Incremental CSV checkpointing: a small warmup buffer establishes the
    # fieldnames (capturing any metrics_fn keys), then every completed row is
    # written and flushed immediately — so an interruption (crash, reboot)
    # preserves all completed samples instead of losing the whole motor.
    rows = []
    _ckpt = {'file': None, 'writer': None, 'warmup': []}
    _warmup_target = min(n_samples, max(progress_every, 16)) if csv_path else 0

    def _open_csv(sample_recs):
        metric_keys = sorted({
            k for rec in sample_recs for k in rec.keys()
            if k not in set(keys + ['idx', 'fitness', 'error'])
        })
        fieldnames = ['idx'] + keys + ['fitness'] + metric_keys + ['error']
        f = open(csv_path, 'w', newline='')
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore',
                           restval='')
        w.writeheader()
        for rec in sample_recs:
            rec.setdefault('error', '')
            w.writerow(rec)
        f.flush()
        _ckpt['file'], _ckpt['writer'] = f, w

    def _checkpoint(rec):
        if csv_path is None:
            return
        if _ckpt['writer'] is None:
            _ckpt['warmup'].append(rec)
            if len(_ckpt['warmup']) >= _warmup_target:
                _open_csv(_ckpt['warmup'])
                _ckpt['warmup'] = []
        else:
            rec.setdefault('error', '')
            _ckpt['writer'].writerow(rec)
            _ckpt['file'].flush()

    def _consume(idx, params, fitness, metrics, err):
        nonlocal best, n_errors
        rec = dict(params)
        rec['idx'] = idx
        rec['fitness'] = fitness
        rec.update(metrics)
        if err is not None:
            rec['error'] = err
            n_errors += 1
        rows.append(rec)
        if fitness < best:
            best = fitness
        _checkpoint(rec)

    if n_workers == 1:
        for done, result in enumerate((_run_one(item) for item in work), start=1):
            _consume(*result)
            if done % progress_every == 0 or done == n_samples:
                emit_progress(done, final=done == n_samples)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as exe:
            futures = [exe.submit(_run_one, item) for item in work]
            for done, fut in enumerate(concurrent.futures.as_completed(futures), start=1):
                _consume(*fut.result())
                if done % progress_every == 0 or done == n_samples:
                    emit_progress(done, final=done == n_samples)

    if csv_path is not None:
        if _ckpt['writer'] is None:        # fewer rows than the warmup target
            _open_csv(_ckpt['warmup'])
        if _ckpt['file'] is not None:
            _ckpt['file'].close()
        if progress_mode != 'none':
            print(f"sensitivity.run_lhs: wrote {csv_path}")

    return rows
