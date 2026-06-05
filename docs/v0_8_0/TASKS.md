# v0.8.0 — openMotor Frontend Integration: Phase Plan

> Companion to [`DESIGN.md`](DESIGN.md) (scoped 2026-06-01). Branch
> `openmotor-frontend`. v0.8.0 = the **data-model + channels backbone**
> (decision D2); GUI render, fore/mid/aft display, tapering grains,
> RocketCEA, and rich viz are deferred (see DESIGN Roadmap).

## Sequencing intent

Land the **data model** (transport-in-propellant, igniter, channels)
before the **plugin contract**, because the plugin returns the channel
objects and consumes the new propellant/igniter schema. Each phase
independently testable; additive/default-OFF where possible; numerics
byte-for-byte unchanged through the refactor.

## Phase 0 — Scoping ✅ (complete)

DESIGN.md decisions D1–D6 resolved with the user. Output: this plan.

## Phase 1 — Channel-object results + `AxialChannel` (D5) — IN PROGRESS

- **DONE (2026-06-01) — channel model:** `srm_1d/channels.py` — `Channel`
  (scalar + per-grain, aligned to openMotor `LogChannel`), `AxialChannel`
  (time × N_cells + `x_cells`, unit-aware `getData` via `motorlib.units`),
  `SimulationChannels` container, `build_channels(results)` (pure re-shape
  of the dict: scalars + snapshot axial fields + per-grain regression/web,
  no recompute), and `as_channels(result)` normalizer (dict-or-channels).
  Exported from `srm_1d/__init__`.
- **DONE (2026-06-01) — data-layer consumer migration:** the
  post-processing consumers now accept dict **or** channels via
  `as_channels`: `nozzle.compute_motor_performance`,
  `openmotor_adapter.compute_grain_metrics` / `result_to_csv` / `save_csv`,
  `run_artifacts.verify_run_health`. Byte-for-byte parity tests prove
  dict-vs-channels identical output (perf arrays, **CSV export**, health
  flag) on a real run. Full suite green (305 passed).
- **DONE (2026-06-02) — capstone (return-type flip):** `run_simulation`
  (and therefore `run_from_ric`) now **returns a `SimulationChannels`**.
  The container **proxies item access / iteration to the raw results dict**
  (`self.raw`), so every legacy `result['P_head']` / `result['summary']`
  consumer — all 24 examples, the sensitivity tool, diagnostics — is
  unchanged, while `.channels[name]` (unit-aware `Channel`) and
  `.axial[name]` (`AxialChannel`) expose the channel API. No example
  rewrite needed; the proxy makes the flip transparent. Internal channel
  consumers updated to `.channels[name].getData(unit)` for conversion vs
  `sc[name]` for raw SI arrays.
- **Gate ✅:** results are channels (the return type); legacy reads work via
  the proxy; numerics byte-for-byte unchanged (channels wrap, never
  recompute); CSV export parity. Full suite green.

> **Why a proxy instead of a hard break:** flipping the return to a pure
> channel object would have broken ~320 tests + 24 examples' `result[key]`
> access. The drop-in mapping proxy (D-style, not a compat shim on the
> numerics) keeps the v0.8.0 API break confined to the *type* while
> preserving call-site behavior — the pragmatic cap that finishes the
> end-to-end channel path without a destabilizing sweep.

## Phase 2 — Unit-aware + generic plotting — DONE (2026-06-01)

- **`plot_channels(result, names, ...)`** — generic unit-aware engine in
  `plotting.py`: plots any scalar channel(s) vs a time channel, overlay or
  stacked, converting to display units at draw time via
  `Channel.getData(unit)`. Accepts a `SimulationChannels`, a results dict
  (normalized via `as_channels`), or a `{name: Channel}` mapping.
- **pressure / thrust / summary routed through the generic path:**
  `plot_pressure` draws its primary trace via `plot_channels` (Pa→MPa on
  the channel); `plot_thrust` wraps the performance series as channels and
  renders via `plot_channels`; `plot_summary` composes the generic
  pressure panel + channel reads.
- **flow-snapshot / snapshots / heatmap re-expressed over `AxialChannel`**
  via `_frames_from_channels(sc)` (rebuilds per-frame field dicts from the
  axial channels — rendering reads the channel model, not the raw dict).
- All plotters accept **dict or channels** (`as_channels`). `plot_channels`
  + `plot_flow_snapshots` + `plot_field_heatmap` exported from
  `srm_1d/__init__`. `plot_grain_regression` already takes
  `compute_grain_metrics` output (channel-agnostic) — unchanged.
- **Gate ✅:** pressure/thrust/summary reproduced via the generic path;
  flow-snapshot/heatmap re-expressed over channels. Parity tests
  (`tests/test_plotting_channels.py`, 9: dict-vs-channels identical plotted
  data, unit conversion, all figures render on Agg). Full suite 314 passed.

## Phase 3 — Transport-in-propellant + migration (D4) — DONE (2026-06-02)

- **openMotor checkout (now git-tracked):** `PropellantTab` gains a single
  `mu` (invariant of frozen/effective — viscosity doesn't shift with
  equilibrium chemistry) + `kThermalFrozen/cpFrozen` +
  `kThermalEffective/cpEffective` (`FloatProperty`, default `0.0` = D7
  sentinel). `Propellant` gains `transportVariant` (`EnumProperty`,
  **default `frozen`**). `fileIO appVersion (0,6,1)→(0,7,0)` + chained
  `(0,6,1)→(0,7,0)` migration seeding sentinels. Backward-compatible
  (oM ignores unknown keys; old files load with sentinel transport).
- **srm_1d adapter:** `convert_propellant` reads per-tab transport from the
  `.ric` by `transportVariant` (collapsed to the representative tab for the
  scalar solver), **D7 hard-fault** on the sentinel; explicit `gas_props`
  still overrides. `load_ric` runs a load-time sentinel migration mirroring
  oM's `doMigration`. **Sidecar retired** in `run_from_ric` (reads `.ric`;
  `transport_path`/`gas_props` are explicit overrides only).
- **Migrator:** `build_transport_library` groups `.ric` by propellant name
  and **shares transport across same-propellant motors** (a sidecar on any
  one fills all — e.g. BATES Test←Hasegawa A Prop, L1843←machbusterNew,
  hlbl38←hlbltest). `migrate_ric_transport` writes both slots + variant and
  **preserves oM's YAML tags** (`!!python/object/apply:uilib.fileIO.fileTypes`,
  `!!python/tuple`) so files stay oM-loadable. `migrate_all_motors` ran on
  all 19 — **0 sentinels** (every propellant had ≥1 sidecar).
- **Default-variant change is intentional:** frozen is now the active
  default (per the v0.7.3 post-phaseB finding that frozen beats effective).
  Motors with explicit frozen data (Hasegawa A Prop, Chase-Energy, 84-CHAR,
  Risky Batman) now solve with frozen k/Cp; lone-primary propellants get
  frozen=effective (behavior-preserving). Both slots are faithfully stored;
  `transportVariant: effective` or `gas_props=` recovers the effective set.
- **Gate ✅:** `tests/test_ric_transport_migration.py` (round-trip, tag
  preservation, sharing, sentinel hard-fault, all 19 repo motors load) +
  `test_adapter` D7 hard-fault / per-tab-read tests.

## Phase 4 — Igniter as data/library (D3) — DONE (2026-06-02)

- **Format:** a motor `data.igniter` block mirrors how the motor embeds
  its propellant — an embedded **pyrogen material** (same field set as
  `srm_1d/motors/pyrogens/*.yaml`: `a/n/rho/T_flame/M/gamma` + kappa_jet,
  particle geometry, heat-delivery) under `igniter.pyrogen`, plus per-motor
  **chamber sizing/topology** (`mass/throat_area/volume/burn_area/burn_law/
  injection_topology/cartridge_length_m/basket_fill_fraction/
  pellet_packing_fraction`). `-1.0` = "auto" sentinel.
- **openMotor checkout:** `motorlib/igniter.py` — `Pyrogen` (material,
  the reusable library item, mirrors `Propellant`) + `Igniter` (chamber
  sizing). `fileTypes.IGNITERS` library type; `DEFAULT_IGNITER` (BPNV +
  forward_plenum auto-sizing); the `0.6.1→0.7.0` motor migration seeds
  `data.igniter`. (oM's Motor/GUI wiring of the block is deferred to
  Phase 5 — oM ignores the unknown key on load.)
- **srm_1d adapter:** `load_igniter(motor)` → `(Pyrogen, sizing_kwargs)`;
  `default_igniter_block` / `_seed_motor_igniter` for migration;
  `load_ric` surfaces `data.igniter`. **`run_from_ric` reads the igniter
  from the motor file when no `pyrogen` kwarg is given** (self-describing);
  an explicit `pyrogen=` takes the legacy kwarg-sizing path. The migrator
  (`migrate_all_motors`) seeds the block; all 19 repo motors carry one
  (BPNV default).
- **Gate ✅:** round-trip (`_pyrogen_to_block`↔`_block_to_pyrogen`),
  all 19 motors parse, chamber builds from the block, self-describing run
  matches the explicit-`bpnv` run (P_peak 6.83 MPa). Tests in
  `tests/test_ric_igniter.py`.

> **Caveat (resolved in Phase 5):** until openMotor's `Motor`/GUI consumes
> `data.igniter`, an openMotor *re-save* would drop the block (oM serializes
> only what its Motor knows). srm_1d reads/writes it fully; the oM GUI
> round-trip lands with the plugin hook.

## Phase 5 — motorlib plugin contract (D1) — CORE DONE (2026-06-02)

- **openMotor checkout:** `motorlib/solvers.py` — a solver registry
  (`register_solver`/`get_solver`/`list_solvers`) + `SolverPlugin` base
  (`name`, `capabilities`, `simulate(motor, config, callback)`). The
  built-in quasi-steady solver is wrapped as `QuasiSteadySolver` and
  registered by default (D6 coexistence — `Motor.runSimulation` stays the
  QS entry point).
- **srm_1d:** `srm_1d/srm1d_plugin.py` — `Srm1dTransientSolver`
  (`capabilities={transient, axial_fields, needs_transport}`) registers on
  import. `simulate_motor(motor)` consumes `Motor.getDict()` (the
  transport-carrying `{nozzle,propellant,grains,config}` shape), runs
  `run_simulation`, and maps the channels into an openMotor
  `SimulationResult` (scalar per-step channels: time/pressure/force/
  exitPressure/dThroat/kn). Igniter defaults to BPNV until Motor carries
  `data.igniter` (GUI wiring).
- **Gate ✅:** `tests/test_srm1d_plugin.py` — the registry selects and runs
  the srm_1d solver headlessly on a canonical `Motor`, returning a
  populated, successful `SimulationResult`; QS solver still registered.
- **Deferred to Phase 5 GUI / Phase 6:** per-grain channel mapping
  (mass/massFlux/regression/web/machNumber), the GUI solver picker +
  capability-gated panels (D6), and Motor wiring of `data.igniter`/transport
  round-trip.

## Phase 6 — GUI integration + validation + docs (IN PROGRESS)

The headless plugin path works (Phase 5). Phase 6 wires it into openMotor's
GUI and finishes the deferred round-trips. Status by task:

1. **GUI solver picker + routing (D6) — DONE (2026-06-02, GUI-verified by
   user).** `motorlib/solvers.py` gains `discover_external_solvers()` (best-
   effort sibling-checkout import of `srm_1d.srm1d_plugin`; honors
   `OPENMOTOR_SRM1D_PATH`; swallows missing-dep failures). `app.py` calls it
   at startup. `SimulationManager` carries an `activeSolverName`
   (default `quasi-steady`) + `setActiveSolver(name)`; `_simThread` routes
   through `solvers.get_solver(name).simulate(motor, callback=...)` and wraps
   it in try/except so a solver fault (e.g. the D7 transport hard-fault)
   surfaces as an alert instead of hanging the progress dialog. `mainWindow`
   adds a **Solver** submenu under Simulate (exclusive checkable actions,
   shown only when >1 solver is registered). *Capability-gated panels
   (axial viz / igniter / transport editors) deferred — see task 2 below.*
2. **Render srm_1d results in openMotor's plot widgets — coexistence path
   DONE (GUI-verified).** The plugin returns a populated openMotor
   `SimulationResult` (now including per-grain channels, task 3), so the
   existing graph / motor-stat widgets render the transient trace with no
   widget changes. **Fix (2026-06-02):** the motor-stats panel reads
   `volumeLoading.getPoint(0)`, which the plugin hadn't populated → `IndexError`
   on result render (Hasegawa C). `_per_grain_series` now also returns the
   per-step `vol_loading` series (100·solid-volume / grain-bounding-volume),
   filled into the `volumeLoading` channel. *Capability-gated srm_1d-only
   panels remain deferred to a v0.8.x GUI pass.*

   **Progress bar + Stop button (2026-06-02).** The transient bar was inert
   (the whole advance is one `@njit` `while` loop). Now `_run_time_loop` takes
   a shared `progress_state` `float64[2]` — it **writes** a composite metric
   `max(web-consumed-fraction, t/t_max)` to `[0]` each step and **reads** a
   cancel flag from `[1]` (>0.5 ⇒ break, termination_code 5 "canceled"). The
   loop is now `@njit(cache=True, nogil=True)` so it runs GIL-free; the plugin
   runs it in a worker thread and a poller forwards `[0]` to openMotor's
   `callback`, setting `[1]` when the callback returns truthy (Stop). Canceled
   runs report `success=False`. `run_simulation` gained `progress_state=None`
   (auto-allocates when run headlessly). Gated by
   `test_progress_callback_driven_to_completion` / `_cancels_run`.

   **Graph-switch decimation (2026-06-02).** The transient solver emits one
   sample per solver step (10⁵+); openMotor's graph widget re-converts
   (per-point, in Python) and re-plots the whole series on every channel
   switch → multi-second hangs. `_result_to_om_simresult` now decimates the
   per-step channels to `GUI_MAX_POINTS=5000` for the GUI `SimulationResult`,
   always keeping the first/last/peak-pressure/peak-thrust samples
   (`_decimate_indices`) so peak P, burn time and impulse stay faithful.
   srm_1d's own full-resolution channels are untouched (analysis / CSV use
   those). Gated by `test_decimate_indices_preserves_peaks_and_endpoints`.
3. **Per-grain channel mapping — DONE (2026-06-02).** `_per_grain_series`
   in `srm1d_plugin.py` aggregates per-cell axial snapshots to per-grain via
   the `cell_segment_id` map (`mass` from ρ·dx·Σ(A_outer−A_port); `machNumber`
   = peak core Mach per grain; `regression`/`web` from `result['grains']`;
   `massFlow` = cumulative −d(mass)/dt forward→aft; `massFlux` = massFlow /
   aft port area), interpolated from the snapshot time base onto the per-step
   time grid so every channel shares `time`'s length.
   `_result_to_om_simresult` now takes `geo`/`prop` and fills all six
   multi-value channels. Gated by `test_per_grain_channels_populated`.
4. **Motor wiring of `data.igniter` — DONE (2026-06-02).** openMotor's
   `Motor` carries `self.igniter` (`Igniter`) + `self.igniterPyrogen`
   (`Pyrogen`); `getDict` emits `data.igniter` (sizing + nested `pyrogen`),
   `applyDict` reads it (backward-compatible — absent block keeps defaults).
   The plugin's `simulate_motor` now resolves the igniter from the motor's
   own block via `load_igniter` (explicit `igniter_pyrogen=` overrides; BPNV
   fallback for pre-block motors) — closes the Phase 4 caveat. Gated by
   `test_motor_round_trips_igniter_block` / `test_motor_without_igniter_keeps_defaults`.
   Transport already round-trips via the per-tab propellant schema (Phase 3).
5. **Validation + docs — PENDING.** End-to-end: a motor designed in
   openMotor → srm_1d transient trace in-GUI (user verifying). Remaining:
   DEVNOTES "API Breaking Changes Log" entries (results→channels,
   transport→propellant, igniter data, plugin) + v0.8.0 close-out doc.

**Test status:** full srm_1d suite green — **337 passed** (10 plugin tests:
solver registry, per-grain channels, volumeLoading, igniter round-trip,
live progress, cancel, decimation). openMotor GUI/registry edits
syntax-checked (`py_compile`); the app builds offscreen, the Solver picker
lists both solvers, and switching/cancel/progress are GUI-verified by the
user. (The progress bar's taildown fill is heuristic — flagged for a
standards-based burn-time revisit, see Deferred.)

## Phase 7 — capability-gated GUI panels (v0.8.x, IN PROGRESS)

Porting QS-only GUI surfaces to be solver-aware. First slice landed
2026-06-02 (offscreen-verified; user verifies from-source):

1. **Per-solver config schema + screen.** `SolverPlugin.get_config_schema()`
   (default `None` = uses the global `MotorConfig`); `Srm1dTransientSolver`
   returns a `PropertyCollection` of its run params (`t_max`, `P_cutoff`,
   `cfl_target`, `roughness`, `kappa`, `T_ignition`, `snapshot_interval` —
   keys match `run_simulation` kwargs; `ambPressure` stays in the shared
   global config; transport variant stays a per-propellant property). The
   **Preferences → General** tab gains a **Solver** dropdown that swaps the
   entire field set between registered solvers (`preferencesMenu` builds an
   independent editable collection per solver; quasi-steady edits `general`,
   plugin solvers edit their schema). Per-solver values persist in
   `Preferences.solverConfigs` (tolerant load; no version bump).
   `SimulationManager._simThread` passes the active solver's saved config to
   `solver.simulate(config=...)` (→ `simulate_motor` overrides → run). Gated
   by `tests/test_srm1d_plugin.py::test_solver_config_schema[_drives_run]`.
   **Refinements (2026-06-02):** kappa relabeled **"Gnielinski Temp-Ratio
   Exponent"** (Ma Eq.9 `(T_gas/T_surface)^κ`, NOT erosive); `snapshot_interval`
   dropped (pending deprecation); **roughness in µm** with its own `units.py`
   'Surface Roughness' category (selectable in the Units tab; converted µm→m
   at the run boundary). **Shared-values design (user-chosen option 2):** each
   solver pane shows a bold **"Shared settings"** section (`maxPressure/
   maxMassFlux/maxMachNumber/minPortThroat/ambPressure/mapDim`, values synced
   across panes via `CollectionEditor.loadGrouped`) plus a solver-specific
   section.
2. **Active-solver coupling (DONE).** Active solver is a persisted preference
   (`Preferences.activeSolver`); `PreferencesManager.setActiveSolver` is the
   single setter (saves + emits `activeSolverChanged`). The Sim→Solver menu,
   the Preferences dropdown, and the per-motor config dropdown all read/write
   it and stay in sync; `SimulationManager` resolves it from preferences.
3. **Per-motor config — true per-motor override (DONE).** `Motor` carries
   `solverConfigs` (dict keyed by solver), serialized in `getDict`/`applyDict`
   and seeded by the `0.6.1→0.7.0` migration (mirrors the igniter block — oM
   pattern). The grain-table **Config** row opens the solver-aware editor
   (`MotorEditor.loadMotorConfig`): shared + QS fields edit the motor's own
   `MotorConfig`; the srm_1d section edits this motor's per-motor override
   (persists in the `.ric` + flows through `getCurrentMotor`/history like all
   motor edits). `SimulationManager` prefers per-motor over the global default.
   Shared logic in `widgets/solverConfigController.py` (used by both screens —
   DRY). Gated by `test_motor_round_trips_solver_configs`.
4. **Basic-box verification (2026-06-03).** Drove openMotor's export +
   editor surfaces against a real srm_1d-transient result offscreen:
   **Eng export ✅** (valid `.eng`), **CSV export ✅** (13 cols, all rows on a
   clean result), **Image export ✅** (same `GraphWidget.saveImage` path as the
   GUI-verified on-screen graph), **per-tab transport editing ✅** (transport
   fields are `Property` objects on `PropellantTab`, rendered by oM's existing
   tabular propellant editor). **Bug found (pre-existing, openMotor-side,
   affects QS too):** `EngExporter.doConversion` mutates the shared result
   in place — `getData()` returns the live channel list and Eng `.append()`s a
   0-thrust point to `time`/`force`, so *Eng-then-CSV on the same result* →
   ragged channels → `getCSV` IndexError. One-line fix (copy the lists before
   appending). PENDING.
5. **Igniter + pyrogen library — DONE + user-verified (2026-06-03; committed
   openMotor fork `fe5c271`, srm_1d this commit).** Final design after
   iteration (full record: `IGNITER_LIBRARY_DESIGN.md`):
   - **Pyrogen library** (Edit → Pyrogen Library): `PyrogenManager` +
     code-built `PyrogenMenu` mirror the propellant subsystem; reusable
     materials persisted via the `IGNITERS` file section (`DEFAULT_PYROGENS`
     seed; duplicate-name guard). `form` is informational (since v0.7.4) so it
     was dropped from the editor; particle dims tooltipped with pellet defaults.
   - **Igniter consolidated into the per-motor Config screen** (a separate
     grain-table "Igniter" row was built, then removed per user request): for
     igniter-capable solvers (`capabilities['igniter']`),
     `MotorEditor.loadMotorConfig` adds a **Pyrogen picker** (copies a library
     material into `motor.igniterPyrogen`) + an **Igniter-chamber** section via
     `SolverConfigController` (igniterCopy group). Topology drives conditional
     fields; the `-1` auto-size sentinel shows as **"auto"** (clear / ≤0
     resets). `motorConfigApplied(general, solverConfigs, pyrogenName,
     igniterProps)` → `mainWindow.applyMotorConfig`. Solver dropdown is exposed
     only inside the open Config screen (`close()` drops state so it can't
     linger) + Simulate→Solver. Editor pane wrapped in a `QScrollArea`.
   - **Run-health watchdog** (`srm1d_plugin`): collapsed / non-igniting runs
     (term code 4, or code 0 with `P_peak < 0.3 MPa`) → `success=False` + an
     ERROR alert the GUI surfaces (no more silent no-plot, e.g. machbusterFAIL).
     `SimulationAlertsDialog` made resizable / word-wrapped / content-fit.
   - Also: `burn_law` → enum; `EngExporter` copies channel lists before
     appending (was mutating the shared result → Eng-then-CSV crash).
6. **Per-station axial viz — SCOPED, design-note-first.** Replace/augment the
   per-grain plot selector with an arbitrary-station selector reading srm_1d's
   native per-cell axial fields, redrawing dynamically. QS keeps its existing
   grain-channel driver untouched; srm_1d gets a capability-gated station panel
   over a carried per-cell field. Defaults = 3 stations/grain (fore/mid/aft;
   fore on, mid/aft off); add/remove/toggle; gap/next-grain boundary
   reassignment; parametric-grain auto-placement deferred. Full design +
   data contract + phases in [`STATION_VIZ_DESIGN.md`](STATION_VIZ_DESIGN.md).
   Phases 1–2 (data payload + station model) are headless/testable in the
   canonical repo and good first lands.

**Tag gate:** cut **v0.8.0** only from a base containing **v0.7.5** (the
cross-motor re-LHS) — see Cross-line sync below. **v0.7.5 re-LHS RAN +
COMPLETE (2026-06-03):** N=3000/motor, 16 workers, ~13.8h (worktree
`../srm_1d-v075-lhs` on `v0.7.0-phase4`). **Rank-1 cross-motor knobs
(combined 1.633, all physical): roughness 32.2 µm, kappa 0.439,
T_ignition 756 K, k_solid 0.271 W/(m·K)** — top-5 tightly clustered, no
unphysical pegging; per-motor fit hasegawa 0.37 / zerox 1.48 / chunc 6.56
(chunc worst, the documented high-L/D QS-erosive limitation). LHS
script + per-sample checkpoint + a result doc committed on `v0.7.0-phase4`.

## NEXT SESSION (viz) — handoff

The igniter/pyrogen GUI feature is DONE + committed (task 5). The clean next
thread is **visualization**, in two queued pieces:
1. **FMM regression-value fix — DONE (2026-06-04).** The `regression` channel
   used to carry **`-web`** for FMM grains — `(avg_D - D_bore_init)/2` is a
   BATES-only formula and Finocyl has `cell_D_bore_init = D_outer`. Fix
   landed: a per-cell `regress` snapshot field (`simulation.py` `_SNAP_REGRESS`
   = 17, `N_SNAP_CHANNELS` 17→18, written in the njit loop, reconstructed in
   the snapshot dict as `snap['regress']`); per-grain `regression` =
   **fore-cell** regress, `web` = **min** over cells of `(wall_web - regress)`.
   Verified positive / monotonic / bounded for a Finocyl (was negative before);
   value assertions added to `tests/test_fmm.py`. Sole consumer is the burnback
   cross-section (→ long-term station-driven). The plugin already stacks
   `grains[k]['regression'|'web']` directly, so no plugin change was needed.
   Detail: `STATION_VIZ_DESIGN.md` §8a. **341/341 pytest green.**
2. **Per-station axial viz** (task 6 / `STATION_VIZ_DESIGN.md`):
   - **Phases 1–2 — DONE (2026-06-04), headless.** `srm_1d/station_viz.py`
     (Qt-free): `build_axial_payload` → decimated `AxialPayload`
     (`[n_frames × n_cells]` field matrices, first+last frame kept) and
     `default_stations` / `make_station` / `grain_cell_spans` /
     `gap_cell_indices` (fore/mid/aft per grain, fore-ON, short-grain
     dedupe, gap classification). `run_simulation` now exports
     `result['cell_segment_id']` + `result['x_cell']` (the data contract);
     the per-cell `regress` from (1) is one of the carried fields. Exported
     from `srm_1d/__init__.py`. 23 tests in `tests/test_station_viz.py`
     (incl. a real multi-grain-result end-to-end contract test). Mass flux
     `G` deferred (needs a per-cell ρ/R snapshot field — not fabricated).
     **364/364 pytest green.**
   - **Phases 3–5 — DONE (2026-06-04).** Plugin: `_axial_payload_for_gui`
     attaches `sr.srm1d_axial` (arrays + embedded default station model, plain
     GUI structures); test `test_srm1d_axial_payload_attached`. oM fork (4
     widgets): `GrainSelector.setupStations`/`getSelectedStations` (the
     **Grains** selector becomes a **Stations** selector — checkboxes grouped
     by grain, fore-ON); `ChannelSelector.appendStationFields` (APPENDS the
     `Axial:` per-cell fields after the normal Y channels — the X/Y channel
     selectors are otherwise UNCHANGED, so pressure/force/Kn/exit pressure/
     dThroat + per-grain channels still plot vs time, default kn/pressure/force);
     `GraphWidget.plotData` (one backward-compatible path: scalar channels vs X,
     per-grain channels for the stations' grains, `Axial:` fields sliced per
     station vs time); `resultsWidget` `rebuildGrainColumns` generates **one
     grain-tab cross-section COLUMN per active station** (the interim per-grain
     dropdown was AXED), each slicing per-cell `regress` (burnout gate = grain
     t=0 `web`). Gated on `_stationMode`/`_yMode`; QS↔srm_1d both directions
     verified. **Correction:** a first cut wrongly replaced the channel
     selectors (lost the standard transient traces) — reverted to this additive
     design. Offscreen smokes: default kn/pressure/force restored (3 lines),
     +Axial field = line/station, +per-grain channel = line for stations' grain,
     columns track stations, `saveImage` signature intact, real native **QS
     unchanged**.
   - **Rich station selector — DONE (2026-06-04).** Replaced the grain
     checkbox list (srm_1d only) with `uilib/widgets/stationSelector.py`
     (`StationSelector`): a global cell-index **slider + spinbox + distance +
     Add** editor and a scrollable list of stations grouped into auto-classified
     categories **Head / Grain N / Gap N / Aft** (each header with its `(c{lo}-
     {hi})` span). A station is just `{cell_index, active}`; category + role are
     DERIVED via new headless `station_viz` helpers (`cell_categories`,
     `grain_role`, `classify_cell`, `station_full_label` — 8 tests). Rows show a
     visibility checkbox + `fore (c3/106) · 59 mm` label and reveal edit/delete
     on hover; double-click/Edit loads the station into the editor and redraws
     live as the slider drags; Add disabled on duplicate; defaults = fore cell
     per grain (visible). Only grain-owned stations get a burnback column;
     head/gap/aft plot only. resultsWidget swaps grainSelector↔stationSelector
     by mode. Offscreen-verified (add/edit/delete/toggle/reclassify/columns/QS
     fallback). 373 pytest green.
   - **Station viz core COMPLETE + COMMITTED (2026-06-05):** srm_1d `0d24cde`
     (backend) + oM fork `dc9cdf3` (GUI) + `eb182cf` (results-column proportions
     fix — main `horizontalLayout` had no stretch → editor ate ~half the screen;
     `setStretch(0,0)+(2,1)` gives the plot the slack; v0.6.1 diligence-checked).
     Nozzle "bars" = upstream casing section (`09155b8`, not a bug).
   - **Post-completion roadmap (user-prioritized 2026-06-05), in order:**
     1. **Mass-flux `G` field + solver-scraping plots** — carry per-cell `G=ρ·u`
        (needs a per-cell ρ/R snapshot field; do NOT fabricate) + other per-cell
        quantities worth exposing as axial fields.
     2. **Axial-profile-at-a-time plotting in the grain view** — time-scrubbed
        field-vs-x view (vector plots / heatmaps) animating alongside the
        station/regression cross-section.
     3. **Parametric tapering geometry for arbitrary FMM grains** — tapered/
        finocyl lack fore/mid/aft anchors; station placement may need
        web-fraction/geometric anchoring. (Geometry effort.)
     - Beyond: open, but the standing high-value target is the **high-L/D
       ignition-transient overshoot** (QS-erosive limitation,
       `docs/v0_7_4/IGNITION_SPIKE_CLOSEOUT.md`) — impairs validating aggressive
       designs without static fires.
     - **Unit-aware station distance — DONE (2026-06-05):** `StationSelector`
       distances (editor readout + row labels) now convert to the user's length
       unit via `motorlib.units` (`setLengthUnit`/`_fmtDist`); resultsWidget
       passes `preferences.getUnit('m')`. Verified mm→`58.8 mm`, in→`2.31 in`.
     - **Pocketed (low priority):** per-station CSV/image export.
   - **Deferred to core-solver sessions (user, 2026-06-05):** #4 sub-cell grain
     eaten by snapping (geometry; affects the sim).

## Cross-line sync (important)

When **v0.7.5 (cross-motor re-LHS)** tags on `v0.7.0-phase4`/main,
merge/rebase it into `openmotor-frontend` so the frontend inherits the
recalibrated model. Cut **v0.8.0** only from a base that already contains
v0.7.5.

## Deferred (post-v0.8.0 — see DESIGN Roadmap)

- v0.8.x: GUI render hook + fore/mid/aft cross-section display.
- v0.8.x/v0.9.0: tapering finocyl grains (BALLSStick target).
- v0.9.0: RocketCEA formula-driven transport (deprecate RPA).
- v0.8.x+: rich burnback/flow-field visualization.
- **v0.8.x: principled burn-time / progress metric (user-flagged
  2026-06-02).** The current progress bar is an ad-hoc composite
  (web-consumed fraction in the burn, steady wall-clock fill through the
  pressure taildown — see Phase 6 task 2). It works but the taildown fill
  is heuristic (fills to ~0.96 then snaps to 100% at termination). Revisit
  using an **industry burn-time definition** so both the progress bar *and*
  the reported burn time/impulse window are standards-based — e.g. the
  **tangent-bisector** method or a **% -of-max-pressure (or % -of-max-
  thrust) threshold** (the 5%–5% / 10%–90% conventions from NFPA 1125 /
  ISO / amateur-rocketry practice). Define burn time on the full-resolution
  trace, then normalize progress against it (so the bar reaches ~100% at the
  standards-defined end-of-burn rather than at the numerical P_cutoff). Also
  consider reporting the chosen burn-time window in the motor-stats panel.

## openMotor fork integration — future work (Phases 1 & 2)

The srm_1d ↔ openMotor integration follows a **one-core-two-frontends**
model: the canonical srm_1d repo is the single source of truth; the
openMotor fork *vendors* the core; upstream gets only the generic,
dependency-free pieces. Full strategy + remote topology are in the
`project_openmotor_fork_integration_strategy` memory.

**Done (2026-06-03):** fork created (`eJar0k/openMotor`; `origin`=fork,
`upstream`=reilleya); srm_1d restructured into a clean installable **flat
package** (`pyproject.toml`; package = `srm_1d/` with `tools/` + `pyrogens/`;
`tests/ examples/ docs/ motors/ static_fire_data/` at repo root; run from
repo root via `python -m pytest tests/` / `python -m examples.<name>`; 340
pytest green; GUI + transient run user-verified). Both repos pushed.

**Ordering decision (2026-06-03):** finish the basic GUI feature set →
**srm_1d code cleanup in the canonical repo** → *then* vendor. Vendoring
does NOT change the plugin-facing API, so deferring it creates no rework;
doing the cleanup first means we vendor a finished, clean core exactly once
(no subtree re-sync). If a GUI feature needs a new srm_1d capability, add it
to the canonical repo + the plugin adapter — never hack it into the fork.

- **Phase 1 — subtree-vendor the core into the fork.** `git subtree split
  --prefix=srm_1d` on the canonical repo → add as a subtree in the oM fork
  (location TBD, e.g. `openMotor/openMotor/srm_1d/`); rewire
  `discover_external_solvers` from the sibling-path hack to the in-tree
  copy; declare srm_1d's heavy deps (numba, scikit-fmm) in the fork's
  requirements; verify the GUI runs self-contained. One-way sync
  (canonical → fork).
- **Phase 2 — upstream PR series to reilleya.** Carve the generic,
  dependency-free pieces off the clean commits (DROP the local `.gitignore`
  rewrite): the **solver-plugin registry** (`motorlib/solvers.py`) is the
  best first PR, then the **GUI solver-picker hooks**. The vendoring itself
  is NOT an upstream PR (upstream won't take numba/scikit-fmm as hard deps).

## Gates / discipline (carried)

- Delete `srm_1d/__pycache__/` (+ `.nbi/.nbc`) after any `@njit` edit.
- Hard API breaks OK — log each in DEVNOTES "API Breaking Changes Log."
- Defer to openMotor data-structure conventions; units convert at the boundary.
- Never hand-edit `.ric`; format changes go through the migration system.
- `pytest tests/` green before each phase close (run from repo root,
  flat-package layout).
