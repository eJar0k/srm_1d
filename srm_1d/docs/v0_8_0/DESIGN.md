# v0.8.0 — openMotor Frontend Integration (SCOPING SKELETON)

> **Status: scoping stub (not yet specified).** This is the first commit on
> the `openmotor-frontend` branch, created so the architectural work is
> front-loaded with a design package (matching v0.7.0/7.1/7.2 practice)
> before any code. Fill in the open questions below, THEN write `TASKS.md`
> phases as real specs. Branch is parallel to the v0.7.x re-LHS work on
> `v0.7.0-phase4`; see `IGNITION_SPIKE_CLOSEOUT.md` (v0.7.4) for where the
> backend left off.

## North star

srm_1d's PISO transient internal-ballistics backend, driven from and
displayed in **openMotor's GUI** — so a user designs a grain/nozzle in
openMotor and gets srm_1d's transient pressure/thrust traces back in the
same tool (alongside openMotor's own quasi-steady solver).

Source: `project_srm_1d_long_term_openmotor_integration` memory.

## Why v0.8.0 (version tier)

This is the **API-breaking tier**. The deferred refactors below change
public data structures (results become channel objects, plotting becomes
unit-aware/generic), which is a hard break vs v0.7.x. Per the project
convention (version = git tags, bump on hard API breaks) this lands as
**v0.8.0**, cut from a base that already contains v0.7.5 (re-LHS) — sync
`v0.7.0-phase4`/main in once v0.7.5 tags, so the frontend inherits the
recalibrated model rather than freezing the pre-recal one.

## Current boundary (what exists today)

- [`openmotor_adapter.py`](../../openmotor_adapter.py) — `.ric` reader,
  transport-YAML loader, `convert_propellant/_geometry/_nozzle`,
  `run_from_ric`, CSV export. This is the existing srm_1d ↔ openMotor data
  boundary.
- [`fmm_grain.py`](../../fmm_grain.py) — bridges to the local openMotor
  checkout (`Erosive Burning Solver/openMotor/openMotor/`) for FMM grain
  regression; walks upward to find it, `SRM1D_OPENMOTOR_PATH` override.
- [`plotting.py`](../../plotting.py) — matplotlib plots (pressure, thrust,
  flow snapshots, summary). Currently srm_1d-internal, not unit-aware, not
  channel-based.
- Results today: a `dict` of numpy arrays from `run_simulation`
  (`time`, `P_head`, `massflow`, … plus the many `*_hist` diagnostics).

## Deferred refactors now IN SCOPE (from the integration memory)

1. **Channel-object refactor** — results as named, unit-tagged channels
   (a `Channel` with name, unit, data, maybe metadata) instead of a bare
   array dict. openMotor's own results model is the alignment target
   (`feedback_openmotor_alignment`: defer to openMotor field names/structure).
2. **Unit-aware plotting** — plots carry units; convert at the boundary
   (srm_1d keeps engineering units internally per `feedback_openmotor_alignment`;
   openMotor uses its unit-preference system → convert at display).
3. **Generic `plot_channels`** — plot any set of channels by name, replacing
   the bespoke per-figure plot functions.

## Open questions to scope (resolve before TASKS phases)

- **Invocation model:** how does openMotor call srm_1d? (a) in-process import
  (both are Python; openMotor is the sibling checkout), (b) a `motorlib`
  plugin/extension point, (c) subprocess + serialized results. Trade: import
  is simplest but couples dependency/versions; subprocess isolates.
- **Data interchange:** what crosses the boundary in each direction?
  openMotor motor definition → srm_1d sim inputs (already partly via
  `.ric`/transport YAML); srm_1d results → openMotor display (the channel
  objects). Define the schema.
- **Channel-object shape:** name, unit, data array, sample-rate/time-base,
  and how diagnostics (the `*_hist` family) map to channels (expose all? a
  curated subset?).
- **Units handling:** which unit system is canonical at the boundary, and
  where conversion happens (srm_1d internal engineering units → openMotor
  user-preference units at display only).
- **What stays in srm_1d vs moves to openMotor:** the PISO solver + physics
  stay in srm_1d; the question is where the channel model and plotting live
  (srm_1d emits channels; openMotor renders them?).
- **Quasi-steady vs transient coexistence:** how srm_1d's transient result
  sits alongside openMotor's existing quasi-steady output in the GUI.

## Constraints / standing guidance

- **Hard API breaks are fine** here (v0.8.0) — refactor cleanly, no
  backward-compat shims (`feedback_api_breaks`). Document each break in the
  DEVNOTES "API Breaking Changes Log."
- **Defer to openMotor's architecture** for new data-structure shapes/field
  names/semantics; UNITS are the documented exception (keep μm/(s·MPa) etc.
  internal, convert at the boundary) — `feedback_openmotor_alignment`.
- **Never hand-edit `.ric`** — regenerate from openMotor
  (`feedback_ric_files_openmotor_owned`).
- Lightweight `save_figure` / `artifact_dir` helpers already exist and are
  fine to build on; the channel/plotting refactor is the deferred piece this
  package scopes.

## References

- `project_srm_1d_long_term_openmotor_integration` memory (north star).
- openMotor local checkout: `Erosive Burning Solver/openMotor/openMotor/`
  (`reference_openmotor_source` memory — prefer Read/Grep over WebFetch).
- v0.7.4 backend state: `docs/v0_7_4/IGNITION_SPIKE_CLOSEOUT.md`.
