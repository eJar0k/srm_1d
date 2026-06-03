# Per-station axial visualization — design note (v0.8.x)

> Scoped 2026-06-03. A srm_1d-mode GUI feature: replace/augment the
> per-grain plot selector with an **arbitrary-station** selector that reads
> srm_1d's native per-cell axial fields and redraws dynamically. Companion
> to [`TASKS.md`](TASKS.md) (Phase 7 — capability-gated GUI panels).

## 1. Motivation

srm_1d is a **1D axial** solver — per-cell data is its native resolution.
The current GUI per-grain channels (`massFlux`, `regression`, `web`,
`machNumber`, `mass`, `massFlow`) are an **aggregation** of per-cell axial
snapshots (`_per_grain_series` maps cells→grains via `cell_segment_id`).
A station selector exposes the resolution that already exists, and lets the
user inspect **axial variation within a grain** — e.g. the fore-to-aft mass-
flux gradient that drives erosive burning. "Per-grain" becomes the special
case "one station per grain."

## 2. Architecture decision (dynamic over a carried field)

**Two independent plot drivers, gated by solver:**

- **Quasi-steady keeps the existing openMotor path unchanged** — fixed
  multi-value channels + the grain-selector checkboxes. No regression risk
  to oM's native solver.
- **srm_1d gets a new per-station panel** whose backend **carries the full
  per-cell axial field in the result and redraws dynamically** from it. It
  does NOT route through openMotor's static channel/grain-selector model
  (that model fixes its columns at result-build time and can't add a station
  post-run). The panel is **capability-gated**: shown only for results that
  carry the srm_1d axial payload.

Rationale: the user's requirements (slider, add/remove stations, per-grain
defaults, post-run flexibility) need on-demand slicing of a carried field —
which openMotor's fixed-channel plumbing can't express. Keeping QS on its own
driver avoids destabilizing the native solver.

## 3. Data contract (what the result carries)

srm_1d's `result['snapshots']` already holds, per snapshot time `t`, these
per-cell arrays (length = n_cells): `P` (pressure), `u` (velocity), `Mach`,
`T` (gas temp), `r_total` / `r_erosive` (burn rate, incl. erosive
component), `D_port`, `x` (cell-center positions), `T_surf`,
`pyrogen_surface_heat_flux`, `mass_source`, `is_burning`, `is_grain`, plus
`grains[k]['regression'|'web']` per grain.

The plugin will attach a compact **axial payload** to the openMotor
`SimulationResult` (a side attribute, e.g. `sr.srm1d_axial`, NOT new
channels) containing:

- `snap_times` — `[n_frames]` time base.
- `x_cell` — `[n_cells]` cell-center positions (for labels + slider).
- `cell_segment_id` — `[n_cells]` cell→grain index (and gap sentinel) for
  classification + default placement.
- per-quantity field matrices `[n_frames × n_cells]` for the plottable set
  (start with: mass flux `G`, `r_total`, `r_erosive`, `Mach`, `P`, `u`,
  `T`, `D_port`; mass flux per cell `G = ρ·u` or from `mass_source`).
- decimated to a sane `n_frames × n_cells` budget (reuse the GUI-decimation
  idea) to bound payload size.

Open: whether to carry the raw snapshots or a curated subset (subset
preferred for payload size + a stable GUI contract).

## 4. Station model

- **A station = a cell index** (selection primitive). Distance-from-head is
  ambiguous under integer cell-snapping, so selection is by cell; the cell's
  **center position is shown as a derived label** (mm from head). Slider and
  numeric entry both map to a cell index.
- Each station carries: owning grain (or gap), cell index, an **active
  (displayed) checkbox**, and a label.
- Plotted series for a station = the chosen field sliced at that cell,
  interpolated onto the common time base.

## 5. Defaults (mirror oM, avoid clutter)

- Pre-populate **three stations per grain — fore / mid / aft cells** (from
  each grain's `cell_segment_id` span).
- **Fore station default-ON; mid + aft default-OFF.** Mirrors openMotor's
  "first grain shown" default and keeps the initial plot uncluttered.
- The selector groups stations **under their owning grain object**.

## 6. Interactions

- Toggle any station active/inactive (checkbox) → live redraw.
- **Add** a station to a grain (pick a cell via slider/entry).
- **Remove** a user-added station. (Default fore/mid/aft are restorable.)
- Switch the plotted quantity (mass flux / burn rate / Mach / P / u / T …)
  → redraw all active stations.

## 7. Boundary handling

- When auto-placing fore/mid/aft, if a nominal cell falls **outside a
  grain's cell span** (very short grains, or it lands in an inter-segment
  **gap**), the station **reclassifies to the gap category** (if the gap
  exists) or **shifts to the next relevant grain**. No station should point
  at an undefined/empty cell.
- Gaps between segments are a first-class category in the selector (a cell
  with `cell_segment_id` = gap sentinel), so gap stations are explicit.

## 8. Coexistence with QS

The results widget chooses the plot driver by the active solver / result
capability: QS → existing grain-selector + channels; srm_1d → station panel.
Shared chrome (axes, export, units) is reused where practical; the selector
column is what swaps. Export (CSV/Image) for srm_1d stations is a follow-up
(the per-station series can feed the existing exporters or a dedicated path).

## 8a. Regression channel & the grain cross-section (decided 2026-06-03)

Tracing consumers: the openMotor `regression` channel has **exactly one
consumer** — `resultsWidget.updateGrainTab`, which draws the **grain burnback
cross-section** (`mapDist = regDist/(0.5·dia)`; `image = grainImages[gid] >
mapDist`). The `web` channel is consumed by the grain-table "web remaining"
number and the `hasWebLeft` burnout-threshold check. Neither feeds burn time,
port ratio, or exports.

**Key semantic point:** that cross-section assumes **one `regDist` for the
whole grain** (uniform burnback). srm_1d's regression is **axially varying**,
so any per-grain scalar (average/max/…) renders a slice that is correct at
*no actual axial location*. Averaging is therefore misleading, not merely
lossy.

**Decisions:**
- **The cross-section is a per-axial-station concept.** Long-term it becomes
  **station-driven**: it renders at the *selected station's* cell `regDist`
  (the same station selection that drives the time-series plots). "Regression
  for grain *k*" → "regression at the displayed station."
- **`web` (grain-table + burnout) → `min` over the grain's cells** — burnout
  is governed by the first cell to break through, so min-web is the
  physically meaningful "web remaining" / burnout trigger (not average).
- **Interim, before the station UI lands:** per-grain `regression` = the
  **fore (foremost) cell** of each grain's cell span. Rationale: the fore
  cell is typically the **least-regressed** (lowest mass flux at the head),
  so it is **last to burn out** — the legacy cross-section stays non-blank as
  long as the head retains web, avoiding a confusing prematurely-blank
  display when the time slider scrolls past the burnout of more-regressed aft
  cells. (Chosen over max-regression for exactly this UX reason.)
- **Crash prerequisite (separate):** the per-grain cross-section only renders
  once `grainImages`/`regressionMap` exist, which requires the plugin to call
  `grain.simulationSetup(config)` on the motor's grains (see TASKS Phase 7
  FMM-crash fix). That is independent of the regression *value* question.

## 9. Open decisions / future work

- **Parametric grains (tapered finocyl, BALLSStick target):** auto fore/mid/
  aft assumes clean axial anchors; tapered/parametric grains may lack them.
  Future auto-population may anchor by **web fraction** or geometric features
  rather than fixed fore/mid/aft cells. Deferred until tapering lands.
- Which quantities to expose by default vs. an "advanced" set.
- Time handling: the main plot is series-vs-time; a future **axial-profile-
  at-a-time** mode (scrub time, plot field vs. x) is a natural sibling but
  out of scope here.
- Payload size / decimation budget for the carried field (perf).
- Whether srm_1d station series should also populate hidden openMotor
  channels for export reuse, or use a dedicated exporter.
- Naming/labeling convention for stations (grain + position + cell #).

## 10. Implementation phases (proposed)

1. **Data payload** — plugin attaches `sr.srm1d_axial` (carried field +
   geometry + cell map), decimated; unit-tested for shape/consistency.
2. **Station model + default population** — pure-logic station list
   (fore/mid/aft per grain, boundary reassignment), unit-tested headless.
3. **Station panel widget** — selector (grouped by grain, checkboxes, slider/
   entry, add/remove) replacing the grain selector for srm_1d results.
4. **Dynamic plot backend** — redraw active stations from the carried field;
   quantity switcher.
5. **Coexistence + gating** — driver selection by solver; QS untouched.
6. **(Later)** export, axial-profile-at-time mode, parametric auto-placement.

Steps 1–2 are headless + testable in the canonical srm_1d repo (good
candidates to land first, before any GUI work); 3–5 are openMotor-fork side.
