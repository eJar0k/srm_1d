# srm_1d — Claude Code orientation

A 1D transient finite-volume solid rocket motor internal ballistics
simulator with the Ma et al. (2020) erosive burning model. Numba-JIT
compiled time loop hits ~45-90k steps/s. v0.7.0 ships a pyrogen-plenum
igniter + Goodman solid-heating ignition model and is calibrated
against Hasegawa A at `mse_all = 0.0968 MPa²` (was 0.24 in v0.6.0).

**v0.7.1 in progress on branch `v0.7.0-phase4`** (no tag yet): N-species
bore-gas refactor (SPINBALL-style "infinite-gases mixture"). Phases
1+2+3+4 complete — `Y[N, 3]` advected per step; per-cell `(γ, Cp, R, M)`
derived and consumed by PISO; sensible-enthalpy advection; cell-N-1
nozzle BC; T_ceiling per cell. Phase 4 validation suite (6 tests in
`tests/test_yns_phase4_validation.py`) confirms mixture arrays collapse
to the correct single-species thermo in pure-pyrogen and pure-propellant
limits. Hasegawa A baseline: 1.4M steps, P_peak 6.26 MPa, mass-balance
err 0.1%, c* 1543 m/s — within Phase 3's ±10% target. Two follow-ups
queued before Phase 5 LHS: strict T_ceiling formula (memory
`project_v0_7_1_t_ceiling_strict_form_pending`) and per-species Cp at
source sites (Phase 3.5). See `srm_1d/docs/v0_7_1/`.

This file is loaded on every session — keep it tight. Pointers to
deeper docs at the bottom.

## Quick start

```bash
# Tests (pyenv 3.10.5 -- has numba, pytest, scikit-fmm installed; 199 tests)
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
├── propellant.py        Propellant + tabs (PropellantTab) + GasProperties + GasSpecies + Pyrogen + thermo utilities
├── grain_geometry.py    GrainSegment / MotorGeometry / build_snapped_geometry; per-cell regress[i]
├── nozzle.py            openMotor-aligned Nozzle: thrust, Isp, CF, throat erosion. Adjusted-CF formula.
├── fmm_grain.py         Bridge to local openMotor checkout; FmmTable extraction + Numba lookup
├── igniter_plenum.py    Pyrogen chamber, choked/subsonic venting, Sutton sizing defaults
├── solid_thermal.py     Goodman integral solid-heating ignition subsolver
├── simulation.py        run_simulation wrapper + @njit _run_time_loop (pyrogen + Goodman); v0.7.1: _advect_species, _compute_mixture_cell, _refresh_mixture_arrays, _compute_T_ceiling_arr
├── plotting.py          matplotlib plots (pressure, thrust, flow snapshots, summary)
├── openmotor_adapter.py .ric reader, transport YAML loader, convert_propellant/_geometry/_nozzle, CSV export
├── motors/              Canonical motor data: <motor>.ric + <motor>.transport.yaml pairs
├── tools/sensitivity.py Latin Hypercube parameter sweeps with parallel execution
├── examples/            hasegawa_motor_a, bates_4seg, hasegawa_a_lhs, Zerox_test, ZeroxOptimizer
└── tests/               15 files, 199 tests
```

## Dev workflow

- **Versioning is git tags**, not folder names. Latest tag: `v0.7.0`
  (2026-05-21), still on branch `v0.7.0-phase4` (no merge to main yet).
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
4. **Igniter API hard-broke in v0.7.0** -- the exponential knobs
   (`igniter_mass`, `igniter_tau`, `ignition_ramp_tau`, `P_ignition`)
   are gone. Use `pyrogen_chamber` directly or `run_from_ric(...,
   pyrogen='bpnv')`. The plenum injects mass + enthalpy + axial
   momentum (via `_orifice_exit_velocity`) into cell 0 / face 1.
   The `igniter_axial_momentum_fraction` kwarg scales the momentum
   contribution (default 1.0 = pure axial head-end jet).
5. **Frozen vs effective gas transport** — tunable knob in v0.7.0,
   addressed structurally in v0.7.1. Frozen (k=0.37, Cp=2060)
   under-predicts erosive spike; effective (k~0.65, Cp~1800)
   over-predicts plateau. v0.7.1 Phase 3 now uses per-cell γ/Cp/R
   derived from the Y[N, 3] mixture; `k_thermal` and `mu_gas` remain
   scalar (deferred to v0.7.2 if Phase 5 calibration shows benefit).
   Until Phase 5 re-LHS lands, the v0.7.0 calibrated `k_solid = 0.482`
   compromise stands.
6. **v0.7.1 thermal_source units changed to W/m (Phase 3 step 1)**.
   Previously kg·K/(s·m) — multiply by Cp_gas to convert legacy
   external builds. Each source site multiplies its mdot·T contribution
   by Cp_gas; the PISO energy equation now treats the input as direct
   enthalpy injection per unit length. `_pyrogen_surface_thermal_sink`
   and `_thermal_source_power` signatures changed accordingly.
7. **v0.7.1 PISO + post-PISO take per-cell arrays (Phase 3 step 2)**.
   `_piso_step_with_energy_diagnostics` and `piso_step` now take
   `gamma_arr / R_arr / Cp_arr / T_ceiling_arr` instead of scalar gas
   thermo. Energy advection is sensible-enthalpy (Cp·T) — face fluxes
   carry upwind Cp·T to conserve energy across cells with different Cp.
   Nozzle BC uses cell-N-1 mixture. T_ceiling is **relaxed** from
   DESIGN §5: max-over-all-species T_flame * 1.01 (not the strict
   Y > 0.05 filter), to avoid clipping the v0.7.0 IC (T = T_flame_prop
   while Y = 100% ambient) on step 0. See
   `_compute_T_ceiling_arr` docstring.

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

- `srm_1d/README.md` -- public API, motor designation, validated parameters
- `srm_1d/ARCHITECTURE.md` -- function-level map of every module
- `srm_1d/DEVNOTES.md` -- full gotchas, calibration state, complete
  API breaking-change log per minor version, performance profile
- `srm_1d/docs/v0_7_0/` -- v0.7.0 hot-gas plenum design package:
  `DESIGN.md` (implemented architecture), `TASKS.md` (phase status),
  `references/` (extracted papers + Goodman integral derivation +
  Sutton/DeMar summaries).
- `srm_1d/docs/v0_7_1/` -- v0.7.1 N-species design package:
  `DESIGN.md` (mixture architecture, Phase plan, ambient species
  decision), `TASKS.md` (Phase 1+2 complete; Phase 3-5 pending).
- `srm_1d/docs/post_v0_7_0/references/` -- SPINBALL research that
  motivated v0.7.1: Cavallini 2009 + DiGiacinto 2008 extractions plus
  the `spinball_walkthrough.md` decision document (recommends Z-N
  dynamic burn rate as the spike-taildown candidate; N-species is the
  prerequisite infrastructure being delivered first).
- `gemini summary.md` (repo root) -- historical record of the v0.6.0
  development cycle that originated build_snapped_geometry, the new
  end-face kernel, and the exponential-decay igniter
- `generic agent instructions.md` -- short current handoff for future
  coding agents. Older external agent-memory references are historical;
  this repo's committed Markdown is the source of truth.

## Open roadmap (priority order)

1. **v0.7.1 Phase 3** -- thread per-cell γ/R/Cp arrays through the
   PISO step + energy advection (T → Cp·T sensible enthalpy) + nozzle
   BC + source-CFL cap + T_flame·1.01 clip ceiling. Estimated 30-40%
   of a fresh session. Phase 4-5 (validation + Hasegawa A re-LHS)
   follow. v0.7.1 tags when re-LHS rank-1 mse_all ≤ v0.7.0 (0.0968).
2. **Z-N dynamic burn rate** -- the spike-taildown candidate identified
   in the SPINBALL walkthrough. One state per cell + a relaxation ODE
   on the steady r_b. No fitted constants if `τ_ZN = κ·α/r²`. Slotted
   as v0.7.2 after v0.7.1 multi-species lands.
3. **Zerox v0.7.0+ re-calibration** -- the v0.6.0 Zerox LHS calibrated
   roughness, `erosionCoeff`, and `a` against the old exponential-decay
   igniter. Re-run after v0.7.1 ships so calibration uses the new
   per-cell mixture thermo.
4. **ε = 0.05 single-cell ignition spike** (radiation-collapse residual)
   -- the only remaining outlier in the 27-variant radiation matrix
   trips on the last-grain-cell ignition transition. Would require
   source sub-stepping (split-operator within one PISO step).
5. **RodTube grain support** -- small extension (PerforatedGrain in
   addition to FmmGrain in `from_openmotor`).
6. **Al2O3 two-phase thermal lag** -- Pardue 1992 form, the secondary
   spike-taildown candidate from the SPINBALL walkthrough. Higher
   implementation cost than Z-N (full re-cal needed). v0.7.3+ if Z-N
   alone doesn't close the residual.
