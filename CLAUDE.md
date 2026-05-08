# srm_1d — Claude Code orientation

A 1D transient finite-volume solid rocket motor internal ballistics
simulator with the Ma et al. (2020) erosive burning model. Numba-JIT
compiled time loop hits ~45–90k steps/s. Validated against Hasegawa
Motor A.

This file is loaded on every session — keep it tight. Pointers to
deeper docs at the bottom.

## Quick start

```bash
# Tests (pyenv 3.10.5 — has numba, pytest, scikit-fmm installed; 107 tests)
"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" -m pytest srm_1d/tests/

# Hasegawa A example (loads srm_1d/motors/hasegawa_a.ric)
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
├── grain_geometry.py    GrainSegment / MotorGeometry / build_snapped_geometry; per-cell regress[i]
├── nozzle.py            openMotor-aligned Nozzle: thrust, Isp, CF, throat erosion. Adjusted-CF formula.
├── fmm_grain.py         Bridge to local openMotor checkout; FmmTable extraction + Numba lookup
├── simulation.py        run_simulation wrapper + @njit _run_time_loop (exp-decay igniter)
├── plotting.py          matplotlib plots (pressure, thrust, flow snapshots, summary)
├── openmotor_adapter.py .ric reader, transport YAML loader, convert_propellant/_geometry/_nozzle, CSV export
├── motors/              Canonical motor data: <motor>.ric + <motor>.transport.yaml pairs
├── tools/sensitivity.py Latin Hypercube parameter sweeps with parallel execution
├── examples/            hasegawa_motor_a, bates_4seg, hasegawa_a_lhs, Zerox_test, ZeroxOptimizer
└── tests/               7 files, 107 tests (incl. test_endface_conservation gating)
```

## Dev workflow

- **Versioning is git tags**, not folder names. Current: `v0.6.0`.
  Bump on hard API breaks; document each break in DEVNOTES "API
  Breaking Changes Log."
- **Hard API breaks are fine** — refactor cleanly, no backward-compat
  shims. (See `feedback_api_breaks` memory.)
- **Defer to openMotor's architecture** when adding data structures
  (field names, semantics). UNITS are the documented exception:
  srm_1d keeps human-readable engineering units (μm/(s·MPa) for
  erosion_coeff, etc.); adapter converts at the boundary. (See
  `feedback_openmotor_alignment` memory.)
- **Named motors live as data**, not Python factories: add a
  `<motor>.ric` + sibling `<motor>.transport.yaml` to `srm_1d/motors/`,
  load via `run_from_ric`. Parametric geometry uses
  `build_snapped_geometry` directly.
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
3. **End-face injection (v0.6.0+)** uses a partition-of-unity hat
   function: each face's mass splits over 2 adjacent cells with weights
   summing to 1.0. Coupled to snapping (which puts faces on cell edges).
   Gated by `tests/test_endface_conservation.py`.
4. **Igniter is a v0.6.0 placeholder** — exponential decay
   (`mdot = (m/τ)·exp(-t/τ)`) with no pressure feedback.
   `igniter_tau≈127ms` for Hasegawa A is acting as a numerical
   FSI-cushioning proxy, not a physical timescale. v0.7.0+ replaces
   with a hot-gas plenum.
5. **Frozen vs effective gas transport** — tunable knob. Frozen
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
- `gemini summary.md` (repo root) — historical record of the v0.6.0
  development cycle that originated build_snapped_geometry, the new
  end-face kernel, and the exponential-decay igniter
- Memory directory (`~/.claude/projects/.../memory/`):
  - `project_validation_targets.md` — figure of merit for the model
    (full pressure-trace match vs experimental is the gold standard)
  - `project_hasegawa_calibration_state.md` — v0.6.0 LHS-derived
    parameter set + FSI fat-spike caveat
  - `feedback_*` — user preferences (defer to openMotor, hard breaks
    OK, terse responses preferred)
  - `reference_openmotor_source.md` — pointer to local openMotor checkout
  - `reference_validation_papers.md` — Ma 2020 + Hasegawa 2006 PDFs
    in repo root

## Open roadmap (priority order)

1. **Hot-gas plenum igniter model** (v0.7.0) — replace exponential
   decay with prescribed P0(t), T0(t) and choked-orifice mass injection
   into cell 0. Removes `igniter_tau` as a tuning knob; restores
   pressure coupling. Pending literature review.
2. **Per-step gas thermo for multi-tab** (deferred) — γ, T_flame, MW
   varying inside the hot loop. Documented in DEVNOTES; hold off
   until calibration shows it helps.
3. **RodTube grain support** — small extension (PerforatedGrain in
   addition to FmmGrain in `from_openmotor`).
