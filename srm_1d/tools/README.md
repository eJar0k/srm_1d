# srm_1d.tools — analysis & calibration tooling

Dev/analysis helpers that ship with the package (`import srm_1d.tools`) but
are **not** part of the simulation core — nothing here runs inside the time
loop. Two modules:

| Module | Purpose |
|---|---|
| [`sensitivity.py`](sensitivity.py) | Latin-Hypercube parameter sweeps with parallel execution + trace-fitting fitness functions |
| [`ignition_diagnostics.py`](ignition_diagnostics.py) | Post-processing of a completed run to dissect the ignition/startup transient |

Both consume the plain results **dict** returned by `run_simulation` /
`run_from_ric` (the source of truth — see the contributor guide
[`01_SIM_CORE.md`](../../docs/contributor_guide/01_SIM_CORE.md) §10), so they
never touch the numerics.

---

## `sensitivity.py` — Latin-Hypercube sweeps

Runs `n_samples` independent simulations in parallel
(`concurrent.futures.ProcessPoolExecutor`), scores each with a user-supplied
scalar **fitness**, and persists every run to CSV. Motor-agnostic — used for
all the calibration LHS work (`docs/v0_7_5/`, the cross-motor re-LHS, etc.).

```python
from srm_1d.tools.sensitivity import run_lhs, mse_fitness

bounds = {                       # each knob -> (lo, hi), sampled by LHS
    'roughness':   (5e-6, 50e-6),
    'kappa':       (0.30, 0.60),
    'T_ignition':  (700.0, 950.0),
}
results = run_lhs(
    motor_path='motors/hasegawa_a.ric',
    bounds=bounds,
    n_samples=500,
    fitness_fn=mse_fitness(exp_t, exp_p_mpa, t_min=0.01),
    n_workers=None,              # None = all cores; 1 = serial (debuggable)
    seed=42,
    csv_path='hasegawa_a_lhs.csv',
    pyrogen='bpnv',              # extra **sim_kwargs passed verbatim per run
)
```

**Public API:**
- `run_lhs(motor_path, bounds, n_samples, fitness_fn, metrics_fn=None,
  n_workers=None, seed=42, csv_path=None, progress_mode='brief',
  sim_verbose=False, **sim_kwargs)` — the driver. Returns the per-sample
  rows (params + fitness + optional metrics) and writes `csv_path`.
- **Fitness factories** (return a picklable callable `fn(result) -> float`,
  lower = better):
  - `mse_fitness(t_exp, p_exp_mpa, t_min=0.0, peak_align_window=None)` — full
    head-pressure-trace MSE vs experimental (the default for trace fitting;
    optional peak-time alignment).
  - `segmented_pressure_fitness(...)` — weighted per-segment MSE so the spike,
    post-spike shoulder, plateau, and tail-down can be fit with separate
    weights.
  - `impulse_error_fitness(impulse_target_n_s)` /
    `peak_pressure_error_fitness(p_peak_target_pa)` — single-scalar targets.
- `pressure_trace_metrics(t_exp, p_exp_mpa, t_min=0.0, ...)` — a `metrics_fn`
  factory that records named all/segment MSE, MAE, bias, peak/trough error,
  and pyrogen summary values alongside each run (for post-hoc analysis).

**Parallelism notes:** the worker `_run_one` is module-level so it pickles
cleanly into `ProcessPoolExecutor` children; the solver's `@njit(cache=True)`
kernels let workers reuse the compiled artifacts. `n_workers=1` runs serially
for easier debugging.

**Example drivers:** `examples/hasegawa_a_lhs.py`, `examples/cross_motor_lhs*.py`,
`examples/zerox_lhs.py`.

---

## `ignition_diagnostics.py` — startup-transient dissection

Post-processing only: given a completed `result` (plus optional `geo` /
`propellant` for absolute source magnitudes), it reconstructs what happened
during the ignition transient and which source family drives the head-end
pressure spike. Built for the ignition-spike investigations
(`docs/v0_7_4/`); the spike is now a documented quasi-steady-erosive
limitation, but these tools remain the way to inspect any startup transient.

**By task:**
- **Landmarks** — `pressure_landmarks(result, ...)` (global peak + a separately
  labeled startup-window peak); `pyrogen_landmarks(result, ...)` (igniter
  active window).
- **Ignition spread** — `ignition_spread_metrics(result)` (first ignition,
  t10/t50/t90/t100, burning fraction, axial ignition order).
- **Source / energy time-series** — `source_timeseries(result, geo, propellant)`
  (normal/erosive/end-face/pyrogen/pyrogen-surface/adjacent-radiation source
  estimates at snapshot times); `energy_momentum_timeseries(result)` and
  `step_diagnostics_timeseries(result)` (per-step energy + pyrogen-momentum
  audit histories).
- **Fine-grained transient** — `early_time_diagnostics(result, window_s=...)`
  and `collapse_event_trace(result, early, ...)` for the first few ms.
- **Classification** — `classify_driver(...)` / `analyze_ignition_spike(result,
  geo, propellant)` reduce all of the above to a labeled "what drove the
  spike" verdict; `classification_report(...)` and
  `literature_evaluation_report(...)` render text reports.
- **Outputs** — `write_diagnostic_outputs(diagnostics, output_dir, ...)` (CSV
  artifacts) and `plot_diagnostic_figures(result, diagnostics, ...)` (overview
  + x–t plots; needs matplotlib).
- `stack_snapshots(result, key)` — helper returning an
  `(n_snapshots, n_cells)` array for any per-cell snapshot field.

**Example drivers:** `examples/chunc_ignition_2x2.py`,
`examples/hasegawa_a_lhs_mode_transport_2x2.py`.

---

*See also:* [`../ARCHITECTURE.md`](../ARCHITECTURE.md) (function-level map),
[`../DEVNOTES.md`](../DEVNOTES.md) (gotchas + calibration state), and the
[`contributor guide`](../../docs/contributor_guide/) (how the core works).
