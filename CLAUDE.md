# srm_1d — Claude Code orientation

A 1D transient finite-volume solid rocket motor internal ballistics
simulator with the Ma et al. (2020) erosive burning model. Numba-JIT
compiled time loop hits ~45–90k steps/s. Validated against Hasegawa
Motor A (P_peak within 2% of experiment).

This file is loaded on every session — keep it tight. Pointers to
deeper docs at the bottom.

## Quick start

```bash
# Tests (pyenv 3.10.5 — has numba, pytest, scikit-fmm installed)
"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" -m pytest srm_1d/tests/

# Hasegawa A example
"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" -m srm_1d.examples.hasegawa_motor_a
```

System Python (3.14 on PATH) does NOT have these deps; always use the
pyenv 3.10.5 path explicitly.

## Module map (one-liners)

```
srm_1d/
├── solver.py            PISO + TDMA + adaptive CFL (pure numerics, no project deps)
├── burn_rate.py         Ma 2020: Haaland → Gnielinski → bisection. Multi-tab Saint-Robert lookup.
├── propellant.py        Propellant + tabs (PropellantTab) + GasProperties + thermo utilities
├── grain_geometry.py    GrainSegment / MotorGeometry; per-cell regress[i] state; analytic + FMM paths
├── nozzle.py            openMotor-aligned Nozzle: thrust, Isp, CF, throat erosion. Adjusted-CF formula.
├── fmm_grain.py         Bridge to local openMotor checkout; FmmTable extraction + Numba lookup
├── simulation.py        run_simulation wrapper + 60-arg @njit _run_time_loop
├── plotting.py          matplotlib plots (pressure, thrust, flow snapshots, summary)
├── openmotor_adapter.py .ric file reader; convert_propellant/_geometry/_nozzle; CSV export
└── tests/               5 files, 98 tests
```

## Dev workflow

- **Versioning is git tags**, not folder names. Current: `v0.5.0`.
  Bump on hard API breaks; document each break in DEVNOTES "API
  Breaking Changes Log."
- **Hard API breaks are fine** — refactor cleanly, no backward-compat
  shims. (See `feedback_api_breaks` memory.)
- **Defer to openMotor's architecture** when adding data structures
  (field names, semantics). UNITS are the documented exception:
  srm_1d keeps human-readable engineering units (μm/(s·MPa) for
  erosion_coeff, etc.); adapter converts at the boundary. (See
  `feedback_openmotor_alignment` memory.)
- **Repo: github.com/eJar0k/srm_1d** (private).
- **Pytest before commit**; clear `__pycache__/` and `.nbi`/`.nbc`
  after edits to @njit functions (Numba cache persistence is the #1
  source of "the fix didn't work" bugs — see DEVNOTES gotchas).

## Critical gotchas (concentrated, full list in DEVNOTES)

1. **Numba cache** — delete `srm_1d/__pycache__/` after any @njit
   edit, or you'll run the old compiled code.
2. **Mass conservation** — the burnout ramp `f_active` MUST multiply
   both `C_burn` AND the regression rate. Either alone causes 3–40%
   mass error.
3. **End-face injection** uses interval containment
   (`x_lo ≤ x_face ≤ x_hi`), not distance. We tried distance — it
   double-counts at boundaries.
4. **Frozen vs effective gas transport** — tunable knob. Frozen
   (k=0.37, Cp=2060) under-predicts erosive spike; effective
   (k~0.65, Cp~1800) over-predicts plateau. Hasegawa A is sensitive
   to this — see calibration memory before tuning.

## External deps

- **openMotor checkout** sits as sibling: `Erosive Burning Solver/openMotor/openMotor/`.
  `srm_1d.fmm_grain` walks upward to find it; override with
  `SRM1D_OPENMOTOR_PATH` env var if needed.
- pip: `numba`, `numpy`, `scipy`, `pyyaml`, `matplotlib`, `pytest`,
  `scikit-fmm` (for FMM grains), `scikit-image` (for openMotor's
  Custom grain). The Cython `mathlib._find_perimeter_cy` is replaced
  at import-time with a Numba marching-squares shim, so MSVC build
  tools aren't needed.

## Where to look for more

- `srm_1d/README.md` — public API, motor designation, validated parameters
- `srm_1d/ARCHITECTURE.md` — function-level map of every module
- `srm_1d/DEVNOTES.md` — full gotchas, calibration state, complete
  API breaking-change log per minor version, performance profile
- Memory directory (`~/.claude/projects/.../memory/`):
  - `project_validation_targets.md` — figure of merit for the model
    (full pressure-trace match vs experimental is the gold standard)
  - `project_hasegawa_calibration_state.md` — known trace-fit issues
    on Hasegawa A; canonical sensitivity-tooling target
  - `feedback_*` — user preferences (defer to openMotor, hard breaks
    OK, terse responses preferred)
  - `reference_openmotor_source.md` — pointer to local openMotor checkout
  - `reference_validation_papers.md` — Ma 2020 + Hasegawa 2006 PDFs
    in repo root

## Open roadmap (priority order)

1. **Sensitivity tooling** — Latin-hypercube / OAT sweeps over
   propellant + transport + roughness + igniter knobs, with batched
   `run_simulation` calls and aggregate plotting. Figure of merit:
   RMS pressure-trace error vs experimental.
2. **RodTube grain support** — small extension (PerforatedGrain in
   addition to FmmGrain in `from_openmotor`).
3. **Per-step gas thermo for multi-tab** (deferred) — γ, T_flame, MW
   varying inside the hot loop. Documented in DEVNOTES; hold off
   until calibration shows it helps.
