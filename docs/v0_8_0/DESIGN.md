# v0.8.0 — openMotor Frontend Integration (DESIGN)

> **Status: scoped.** Phase-0 open questions resolved with the user
> (2026-06-01). This is the durable design doc; [`TASKS.md`](TASKS.md)
> turns it into phases. Branch `openmotor-frontend`, parallel to the
> v0.7.x re-LHS line on `v0.7.0-phase4`. Backend left off at
> [`docs/v0_7_4/IGNITION_SPIKE_CLOSEOUT.md`](../v0_7_4/IGNITION_SPIKE_CLOSEOUT.md).

## North star

srm_1d's PISO transient internal-ballistics backend, **driven from and
displayed in openMotor's GUI** — a user designs a grain/nozzle in
openMotor and gets srm_1d's transient pressure/thrust traces back in the
same tool, alongside openMotor's own quasi-steady solver. Source:
`project_srm_1d_long_term_openmotor_integration` memory.

## Version tier & cross-line sync

This is the **API-breaking tier** (results become channel objects;
transport moves into the propellant; igniters become first-class data) →
**v0.8.0** by the project's tag-on-hard-break convention. **Cut v0.8.0
from a base that already contains v0.7.5 (cross-motor re-LHS)** so the
frontend inherits the recalibrated model, not the pre-recal one. When
v0.7.5 tags on `v0.7.0-phase4`/main, merge/rebase it into
`openmotor-frontend` before tagging v0.8.0. The frontend's API churn must
not block the calibration line, and must not drift from the validated
physics.

**Two version spaces, kept distinct:** srm_1d's git tag (`v0.8.0`) is
separate from openMotor's *file-format* version (`appVersion`, currently
`(0, 6, 1)` in `openMotor/uilib/fileIO.py`). The schema additions below
bump the openMotor file-format `appVersion` (proposed `(0, 7, 0)`) and
register migrations; that bump is independent of the srm_1d tag.

---

## Resolved decisions (Phase 0)

### D1 — Invocation: **motorlib plugin API**

openMotor calls srm_1d through a **formal extension point** in
`motorlib`, not an ad-hoc import or subprocess. srm_1d registers itself
as a *transient solver* that openMotor's simulation manager can select
alongside its built-in quasi-steady solver.

- **Contract (to be specified in Phase 1):** a solver plugin exposes
  `name`, `simulate(motor, config) -> SimulationResult`, and a
  capabilities descriptor (e.g. "produces axial fields", "requires
  transport properties"). openMotor discovers registered solvers and
  routes the active one.
- **Why plugin over in-process import / subprocess:** cleanest long-term
  architecture; the GUI stays solver-agnostic, transient/quasi-steady
  coexistence (D6) falls out of solver selection, and a future
  out-of-process or remote solver can implement the same contract without
  touching the GUI. Cost is the upfront contract design — accepted.
- **Dependency reality:** srm_1d already bridges *up* into the openMotor
  checkout (`fmm_grain.py`, `SRM1D_OPENMOTOR_PATH`). The plugin inverts
  the *call* direction (openMotor → srm_1d) while srm_1d keeps importing
  openMotor's grain/geometry classes. Numba/scikit-fmm stay srm_1d-side;
  the plugin boundary is where heavy deps are isolated.

### D2 — v0.8.0 scope: **data-model + channels backbone**

In v0.8.0:

1. **Channel-object results + unit-aware/generic plotting** (the original
   core). srm_1d emits openMotor-aligned channels; `plot_channels`
   replaces bespoke per-figure functions.
2. **Transport-in-propellant + migration** (D4). Retire the sidecar
   `.transport.yaml`; fold transport into the propellant, version-locked.
3. **Igniter as data/library** (D3).

**Deferred** (consume the backbone; see Roadmap): GUI render hook +
fore/mid/aft cross-section display → **v0.8.x**; tapering finocyl →
**v0.8.x / v0.9.0**; RocketCEA → **v0.9.0**; rich field/vector viz →
**v0.8.x+**.

Rationale: every deferred item depends on the channel + data-model
backbone existing first. Landing the backbone alone keeps v0.8.0
independently testable (byte-for-byte numerics unchanged) and lets the
GUI work start from a stable contract.

### D3 — Igniter: **library + motor-block** (mirror propellant)

openMotor today has **no igniter concept** — `igniterPressure` was
deliberately removed in the 0.4→0.5 motor migration. We reintroduce it as
first-class data modeled on the propellant pattern: a reusable **igniter
library** file (new `fileTypes.IGNITERS`) + an **embedded igniter block**
in each motor file. Captures srm_1d's `PyrogenChamber` parameters
(mass, form archetype, injection topology, Sutton sizing, transport).

- Self-describing motor files; a future GUI igniter picker is natural.
- New library ⇒ new default-library seeding + migration that injects a
  default igniter into pre-v0.7.0-format motors.

### D4 — Transport: **per-PropellantTab** (frozen + effective)

Transport (`mu`, `k`, `Cp`) attaches to **each `PropellantTab`**, not the
propellant as a whole — pressure-varying transport that parallels the
existing per-tab `a/n/k/t/m` combustion block and srm_1d's own multi-tab
Saint-Robert lookup in [`burn_rate.py`](../../burn_rate.py).

- **Frozen + effective both retained.** Each tab carries a named transport
  variant pair; the calibration line still A/Bs frozen vs effective
  (`project_v0_7_3_post_phaseB_state` reversed to frozen for Hasegawa A).
  Schema must hold both; the active variant is a motor/solver setting.
- **Implemented as (Phase 3, 2026-06-02):** a single per-tab `mu` (viscosity
  is invariant of the frozen/effective equilibrium shift — user tweak) +
  `kThermalFrozen/cpFrozen` + `kThermalEffective/cpEffective`, flat
  `FloatProperty`s mirroring oM's `a/n/k/t/m`. `transportVariant` on the
  `Propellant`, **default `frozen`** (user decision; matches the v0.7.3
  post-phaseB finding that frozen beats effective). `0.0` = D7 sentinel.
- **Shared by propellant (user decision):** the migrator builds a transport
  table keyed by propellant name and shares it across all motors using that
  propellant, so a sidecar on any one motor fills the rest (fills the 3
  sidecar-less repo motors). Migrated `.ric` keeps oM's YAML tags so it
  stays oM-loadable.
- **Adapter behavior today:** srm_1d's solver uses scalar `mu/k/Cp`, so
  `openmotor_adapter` collapses per-tab transport to the operating-tab
  value (or single tab). Per-pressure transport lookup in the solver is a
  later option the schema already permits — no second migration needed.
- **CEA-shaped:** a future RocketCEA solver (v0.9) populates exactly this
  block from a formula, so design the field set (Cp, viscosity,
  conductivity, Pr, plus the combustion `t/k/m`) to match CEA output now.

### D5 — Channel model: **`AxialChannel` subclass** for per-cell fields

openMotor's `LogChannel` (`simResult.py`) handles scalar and per-*grain*
list channels. srm_1d's diagnostics are per-*cell* axial fields
(time × N_cells). Add an **`AxialChannel`** alongside `LogChannel`:

- Stores a `time × N_cells` array + the axial coordinate (`x_cells`),
  with unit-aware `getData(unit)` matching `LogChannel`'s conversion
  contract.
- openMotor's existing per-grain channels are **untouched** (GUI parity
  preserved); srm_1d adds axial channels for its fields. The GUI ignores
  axial channels until the visualization work (v0.8.x) consumes them.
- srm_1d emits openMotor's `SimulationResult` populated with both the
  GUI-native scalar/per-grain channels **and** the new axial channels, so
  the GUI renders the standard graphs immediately and richer field views
  later read the same result object.

### D7 — Migrated transport defaults: **hard-fault sentinels, not fabricated values**

A consequence of the universal-schema decision (D6): when migration adds
transport (and any other physically required) properties to existing
propellants, what value goes in?

- **Rejected (a): fabricated "sensible" defaults.** Transport
  (`mu/k/Cp`) varies too much from propellant type to type for a generic
  default to be anything but misleading — a plausible-looking wrong number
  silently corrupts the erosive-burning chain (μ feeds Re).
- **Chosen (b): hard-faulting sentinel.** Migration injects a sentinel
  ("not provided") that **blocks simulation with the srm_1d transient
  solver** until the user supplies real numbers (RPA/CEA/measurement),
  with a clear alert pointing at the missing field. Robust over
  convenient — consistent with the no-fabrication standard
  (`feedback_defer_to_thermochem_solvers`, `feedback_no_unfounded_smoothing`).

Scope: applies to properties whose value is type-specific and unknowable
without data (transport). It does **not** apply where a physically
motivated default exists (igniter Sutton sizing, form archetype). The QS
solver, which never reads transport, is unaffected by the sentinel —
only the srm_1d solver hard-faults, preserving D6's "basic users keep
iterating on QS" path.

### D6 — Quasi-steady vs transient coexistence (QS solver is PRESERVED)

openMotor's built-in **0D quasi-steady solver is kept, not replaced** —
it is genuinely useful for rapid, cheap iteration. The plugin model (D1)
exists precisely so a basic user iterates on the QS solver while an
advanced user opts into the heavier srm_1d transient model. Both solvers
implement the plugin contract and return a `SimulationResult`; transient
results render in the same unit-tagged, plot-generic widgets.

**Capability-gated exposure (instead of a standalone "advanced mode"
flag).** srm_1d-specific affordances — axial/field viz, the igniter
library editor, per-tab transport editing — are surfaced by **which solver
is active**, driven by the plugin's **capabilities descriptor** (already
in D1), rather than a parallel global mode system. The active solver
declares the panels/affordances it needs; the GUI reveals them. This is
more discoverable and less modal than a separate switch, and it reuses the
plugin contract we're already building.

A thin **"expert/advanced" preference** still has a role: it gates
*visibility of the solver picker itself* (and the srm_1d transient option)
so beginners aren't confronted with solver choice — but it does not
duplicate per-feature toggles. Net: one preference unlocks "advanced
mode" = the solver picker; everything else follows from solver selection.

**Data vs UI separation (important for D3/D4).** The schema changes are
**universal** — after migration, every motor file carries an igniter block
and per-tab transport regardless of which solver runs (files stay
self-describing; the QS solver simply ignores fields it doesn't use). Only
the **editing UI** for those fields is capability/expert-gated. Don't
conflate "advanced feature" (UI exposure) with "advanced data" (the schema
is shared).

Side-by-side overlay of QS vs transient traces is a v0.8.x display
feature, not a data-model change.

---

## Architecture summary

```
openMotor GUI
  └─ simulationManager ── selects solver via motorlib plugin registry (D1)
                            ├─ built-in quasi-steady solver
                            └─ srm_1d transient solver  ◄── registers here
                                 │  imports openMotor grain/geometry (existing bridge)
                                 │  runs PISO (numba, srm_1d-side)
                                 └─ returns SimulationResult with:
                                      • LogChannel  (scalar + per-grain)  → GUI graphs
                                      • AxialChannel (time × N_cells)     → fields/viz (D5)

Motor / propellant / igniter files  (openMotor YAML, {version,type,data})
  • Propellant.tabs[i] += transport {frozen,effective}{mu,k,Cp}  (D4)
  • new fileTypes.IGNITERS  + embedded motor igniter block       (D3)
  • appVersion (0,6,1) → (0,7,0); chained migrations seed defaults
```

## Migration plan (openMotor file format)

Follow the established `migrations` chain in `uilib/fileIO.py`:

- Bump `appVersion` → `(0, 7, 0)`; add a `(0,6,1) → (0,7,0)` entry.
- `migrateProp_0_6_1_to_0_7_0`: inject transport (frozen+effective) into
  every tab of every library propellant — existing `.transport.yaml`
  values where known, else a **hard-fault sentinel** (D7), never a
  fabricated number.
- `migrateMotor_0_6_1_to_0_7_0`: inject transport into the embedded
  propellant's tabs (same sentinel rule) **and** a default igniter block
  (igniter defaults ARE physically motivated — Sutton sizing — so no
  sentinel needed there).
- New `IGNITERS` library: seed a default igniter library on first run
  (mirrors `DEFAULT_PROPELLANTS` seeding).
- Gate with a round-trip test (load old fixture → migrated → re-save →
  re-load) per the project's per-phase discipline.

## Remaining sub-questions (decide during phase work, recommendations noted)

- **`appVersion` coordination** — confirm `(0, 7, 0)` (vs `(0, 6, 2)`).
  *Rec:* `(0, 7, 0)` — schema additions are a minor-version concern.
- **Transport default source** — RESOLVED (D7): use existing
  `.transport.yaml` values where present; otherwise inject a
  **hard-faulting sentinel** that blocks the srm_1d solver until the user
  supplies real numbers. No fabricated defaults.
- **Igniter block fields** — exact `PyrogenChamber` subset to serialize
  (topology, form archetype, mass, Sutton knobs, transport). *Rec:*
  serialize the full constructor signature; defaults match
  `build_pyrogen_chamber`.
- **Plugin discovery mechanism** — entry-points vs explicit registry call
  in the openMotor checkout. *Rec:* explicit registry (the checkout is
  vendored; entry-points add packaging ceremony for no gain here).

## Deferred roadmap (post-backbone)

- **v0.8.x — GUI render hook + fore/mid/aft cross-section.** Render
  srm_1d traces in openMotor widgets (D6). Fix the misleading single
  cross-section under axial burnback: query one grain's `regressionMap`
  at fore/mid/aft cell `regDist` (configurable station count, default 3).
  The regression map already supports this — it's three contour queries on
  one map (`grain.getRegressionData` mechanics).
- **v0.8.x / v0.9.0 — tapering geometries (finocyl).** New `Grain`
  subclass: interpolate fin *parameters* axially → coreMap + FMM
  regression map **per axial station** → srm_1d consumes per-cell
  FmmTables. BALLSStick (real tapering finocyl, currently stepped) is the
  validation target. Interpolating endpoint *results* is the cheaper but
  geometrically approximate fallback.
- **v0.8.x+ — rich visualization.** Burnback animation, flow scalar/vector
  fields over time. Builds on the channel refactor + existing srm_1d
  `plot_flow_snapshots` / `plot_field_heatmap` helpers.
- **v0.9.0 — RocketCEA.** PropPep-style formula entry → CEA equilibrium →
  populate the per-tab transport + combustion block (D4), deprecating RPA
  inputs. Schema is designed CEA-shaped now; implementation deferred.

## Constraints / standing guidance

- **Hard API breaks are fine** (v0.8.0) — clean refactor, no compat shims;
  log each in DEVNOTES "API Breaking Changes Log" (`feedback_api_breaks`).
- **Defer to openMotor data-structure conventions**; **units are the
  exception** — keep μm/(s·MPa) etc. internal, convert at the boundary /
  display (`feedback_openmotor_alignment`).
- **Never hand-edit `.ric`** — regenerate from openMotor; format changes
  go through the migration system (`feedback_ric_files_openmotor_owned`).
- Delete `srm_1d/__pycache__/` (+ `.nbi/.nbc`) after any `@njit` edit.
- `pytest srm_1d/tests/` green before each phase close.

## References

- `project_srm_1d_long_term_openmotor_integration` memory (north star).
- openMotor checkout: `Erosive Burning Solver/openMotor/openMotor/`
  (`reference_openmotor_source`; prefer Read/Grep). Key files:
  `uilib/fileIO.py` (version+migration), `motorlib/propellant.py`,
  `motorlib/simResult.py` (LogChannel), `motorlib/grain.py`
  (regression display).
- srm_1d boundary today: [`openmotor_adapter.py`](../../openmotor_adapter.py),
  [`fmm_grain.py`](../../fmm_grain.py), [`plotting.py`](../../plotting.py).
- v0.7.4 backend state: [`docs/v0_7_4/IGNITION_SPIKE_CLOSEOUT.md`](../v0_7_4/IGNITION_SPIKE_CLOSEOUT.md).
