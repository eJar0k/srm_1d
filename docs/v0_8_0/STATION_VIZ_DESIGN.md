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

## 8a. Regression channel & the grain cross-section (decided 2026-06-03; IMPLEMENTED 2026-06-04)

> **Status (2026-06-04): the per-cell `regress` snapshot + the interim
> per-grain regression/web rule below are IMPLEMENTED** in `simulation.py`
> (`_SNAP_REGRESS`; `snap['regress']`; per-grain `regression` = fore-cell
> regress, `web` = min over cells of `wall_web - regress`). Verified
> positive/monotonic/bounded for a Finocyl FMM grain (negative before the
> fix); value assertions in `tests/test_fmm.py`.
>
> **The station-driven cross-section is now IMPLEMENTED (2026-06-04):**
> `resultsWidget` renders **one grain-tab COLUMN per active station** (the
> graph-view station selection), each slicing the carried per-cell `regress`
> at its cell (grain wall web = its t=0 `web` channel value → burnout gate).
> The earlier interim per-grain dropdown was AXED per the user's direction.
> Default active = fore (matches the per-grain `regression` channel).
> Remaining: add/remove arbitrary (non fore/mid/aft) cells in the selector.


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

> **Status (2026-06-04): phases 1–2 IMPLEMENTED headless** in
> `srm_1d/station_viz.py` (Qt-free) + `tests/test_station_viz.py` (23 tests).
> `run_simulation` now also exports `result['cell_segment_id']` and
> `result['x_cell']` (the v0.8.x data contract) so the payload + station
> model build from `result` alone. Phases 3–5 (the GUI panel) are next, on
> the openMotor-fork side.

1. **Data payload — DONE.** `build_axial_payload(result, fields, max_frames)`
   → `AxialPayload` (decimated `[n_frames × n_cells]` field matrices +
   `snap_times` + `x_cell` + `cell_segment_id`). First+last frame always kept.
   Mass flux `G` deferred (needs a per-cell ρ/R snapshot field — not
   fabricated). The plugin will attach this as `sr.srm1d_axial` in phase 3.
2. **Station model + default population — DONE.** `default_stations(...)`
   (fore/mid/aft per grain; fore-ON, mid/aft-OFF; short-grain dedupe),
   `make_station(...)` (user add, gap classification), `grain_cell_spans` /
   `gap_cell_indices`. Boundary reassignment is implicit for auto-placement
   (cells chosen from each grain's actual span, never a position nominal).
3. **Station selector — DONE (2026-06-04; corrected).** ONLY the graph view's
   **Grains** selector becomes a **Stations** selector for srm_1d results
   (`GrainSelector.setupStations`/`getSelectedStations` — fore/mid/aft
   checkboxes grouped under per-grain headers, fore-ON). The X/Y channel
   selectors are left intact (see 4). This selection is the source of truth
   for the plot AND the grain-tab columns. (The earlier per-grain dropdown was
   AXED.) **Correction note:** a first cut wrongly REPLACED the Y/X channel
   selectors with a field-only picker, losing pressure/thrust/Kn/exit pressure/
   dThroat over time — reverted to the additive design below.
4. **Plot backend — DONE (2026-06-04; additive, non-destructive).** The X/Y
   selectors keep ALL openMotor channels (scalars + per-grain), so the standard
   transient traces still plot vs time; default Y stays kn/pressure/force.
   `ChannelSelector.appendStationFields` APPENDS the carried per-cell fields
   (`Axial: Burn Rate/Erosive Burn Rate/Regression/Mach/Pressure/Velocity/Gas
   Temp/Port Diameter`) after the channels. `GraphWidget.plotData` (single
   path, backward-compatible signature) plots: scalar channels once vs X;
   per-grain channels for the grains owning the selected stations; and any
   selected `Axial:` field sliced at each selected station vs time (one line
   per field×station). The **grain burnback tab generates one cross-section
   COLUMN per active station** (`rebuildGrainColumns`), each rendering its
   cell's per-cell `regress` over the slider. Verified offscreen: default
   kn/pressure/force = 3 lines restored; +Axial Burn Rate adds a line per
   station; +per-grain Regression adds a line for the stations' grain; columns
   track active stations; `saveImage`/image-export signature unchanged.

   **Polish (2026-06-04, post-review):**
   - **No duplicate quantities.** In station mode the openMotor channels with an
     axial equivalent (`pressure`, `machNumber`, `regression`) are EXCLUDED
     (`stationExcludedChannels`); the axial fields take their canonical names
     (Pressure/Mach Number/Regression — no `Axial:` prefix). Default Y =
     kn/force + axial **P** (replaces the chamber-pressure channel) so the
     initial three traces still come up.
   - **Units.** Axial fields carry a unit category (`Pa`/`m/s`/`m`/`K`/'') and
     convert to the user's display unit (MPa, mm, …) via `motorlib.units`; the
     time base converts too. (Was raw SI.)
   - **Snapshot resolution.** `simulate_motor` sets `snapshot_interval =
     max(0.005, t_max/2000)` (was 0.2 s) and `_axial_payload_for_gui` carries up
     to 4000 frames — at 0.2 s (~5 frames) the axial traces + interpolated
     per-grain channels were jagged and ate the ignition spike; now ~200–2000.
   - **Legend.** `loc='upper right'`, small font, wrapped into columns (~10 rows
     max) so it no longer jumps mid-plot or overflows downward.
5. **Coexistence + gating — DONE (2026-06-04).** Gated on `getattr(sr,
   'srm1d_axial', None)` via `_stationMode`; Y is rebuilt only on a real mode
   change (`_yMode`), both switch directions (QS↔srm_1d) exercised. A real
   native **QS result is unchanged** (kn/pressure/force defaults, grain columns,
   X visible, no axial fields).
6. **Station viz core COMPLETE + COMMITTED (2026-06-05)** — canonical srm_1d
   `0d24cde` (backend) + oM fork `dc9cdf3` (GUI) + `eb182cf` (column-proportions
   fix). Phases 1–5 + §8a + the rich selector all shipped. QS untouched.

## 11. Post-completion roadmap (user-prioritized 2026-06-05)

Viz core is cleared. Next work, in order:

1. **Mass-flux `G` field + other solver-scraping plots.** Carry a per-cell
   mass flux (`G = ρ·u`) — needs a per-cell ρ or R snapshot field in
   `simulation.py` (do NOT fabricate from incomplete state); then expose `G`
   as an axial field, plus any other per-cell quantities worth scraping from
   the solver state.
2. **Axial-profile-at-a-time plotting in the grain view.** A time-scrubbed
   *field-vs-x* view (vector plots or heatmaps) that animates alongside the
   station/regression cross-section — i.e. scrub time and see the whole axial
   profile, not just per-station time series. Lives in the grain tab next to
   the burnback cross-section.
3. **Parametric tapering geometry for arbitrary FMM grains.** Tapered/finocyl
   grains lack clean fore/mid/aft anchors — station auto-placement may need
   web-fraction or geometric anchoring. Larger geometry effort.

Beyond that the field is open, but the standing high-value target is the
**high-L/D igniter / ignition-transient overshoot** (the QS-erosive limitation
documented in `docs/v0_7_4/IGNITION_SPIKE_CLOSEOUT.md`): it impairs a key use
case — validating extremely aggressive motor designs without expensive/
dangerous static fires. Any non-tuned transient closure here is high-impact.

**Done (2026-06-05):** unit-aware station distances — the selector readout +
row labels convert to the user's length unit (`StationSelector.setLengthUnit`/
`_fmtDist` via `motorlib.units`; resultsWidget passes `getUnit('m')`).

**Pocketed (low priority, most users won't need):** per-station CSV/image
export.

Steps 1–2 (data payload, station model) were headless + testable in the
canonical srm_1d repo; 3–5 are openMotor-fork side.

> **Roadmap #1 (mass-flux G) — DONE 2026-06-06.** Per-cell `rho` snapshot
> (`_SNAP_RHO`) + derived `G = rho*u` in `build_axial_payload`; oM-fork
> `resultsWidget.stationFields` plots `G` (Mass Flux) + `rho` (Gas Density).

## 12. Roadmap #2 — longitudinal motor-slice viewer (scoped 2026-06-06)

A side-on 2-D axial cut of the whole motor that animates **burnback and the
flow field together**. Each cell is a vertical strip: solid propellant is
drawn from the bore wall out to the casing (so the core widens as the web
burns), and the open bore is filled with a chosen flow field as a heatmap
(later + velocity/G vectors). A time slider scrubs all frames.

**User decisions (2026-06-06):**
- **Radial axis = auto-stretch** (independent X/Y scaling; radius
  exaggerated to fill the panel; axes labeled real units + "exaggerated"
  note). A true-1:1 toggle is deferred to polish.
- **Placement = inside the Grains tab** (alongside the existing per-station
  radial cross-section columns), not a separate tab.
- **v1 scope = Phases A–C** (data + animated solid burnback + bore heatmap
  + colorbar + scrub). Vectors (D) and polish (E) are a second pass.

**Data — already carried** (`sr.srm1d_axial`, per frame×cell): `D_port`
(hydraulic bore → R_bore), flow fields `P/u/G/T/Mach/rho`, `regress`,
`cell_segment_id`, `x_cell`, `snap_times`.

**Data — Phase A canonical additions** (genuine geometry, no fabrication;
all present in the solver's geometry arrays today, just not exported):
`dx` (uniform cell width → axial extents), `D_outer` (scalar → R_outer /
casing wall), per-cell `cell_wall_web` (→ %web-remaining shading). Export
from `run_simulation` → carry through `build_axial_payload` /
`_axial_payload_for_gui`. Add value-asserting tests.

**Rendering** (oM fork, matplotlib `FigureCanvas`): `fill_between` for the
mirrored solid (R_bore↔R_outer), `pcolormesh` over a non-uniform quad mesh
(per-cell y ∈ [−R_bore, +R_bore]) for the bore heatmap + colorbar; later a
`quiver` of axial u/G at cell centers. Scrub updates artist data for the
active frame (reuse the grain-tab time slider). Field dropdown selects the
bore quantity; unit-aware via `motorlib.units`.

**Geometry mapping:** R_bore(i,t) = `D_port`/2; R_outer = `D_outer`/2
(constant); %web(i,t) = 1 − `regress`/`cell_wall_web`. Non-grain cells
(gap/head/aft, `cell_segment_id` < 0) draw open chamber (no solid).

**Phases:** A) canonical data + tests · B) static slice (solid burnback +
bore heatmap + colorbar) · C) time-slider animation · D) u/G vector
overlay · E) polish (1:1-aspect toggle, gap/head/aft edge cases, perf via
artist-update/blit, unit-aware colorbar).

**Limitations to surface in-UI:** non-circular grains (finocyl/star)
render as an *equivalent hydraulic* radius, not literal geometry; the
deferred sub-cell-grain snapping gap means a grain thinner than one cell
won't appear; auto-stretched radial axis is not to scale (labeled).
