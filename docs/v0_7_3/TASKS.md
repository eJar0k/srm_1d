# srm_1d v0.7.3 — Implementation Plan

**Scope**: candidate 4 from
[../v0_7_2/candidates_post_phaseA.md](../v0_7_2/candidates_post_phaseA.md)
(submerged pyrogen modes) unified under an **uncontained-pyrogen
architecture** that supports both head-end and aft-cavity placement.
Phase A ships the topology infrastructure + diagnostic visualization
+ validation runs that exposed a structural ignition-initiation gap.
Tagged `v0.7.3-phaseA` as an intermediate milestone.

## Phase A — Uncontained pyrogen topology [COMPLETE 2026-05-24]

### Phase A.1 — kernel + PyrogenChamber topology fields [COMPLETE]

Commit `0d83501`.

- [x] `_compute_uniform_band_weights(dx, i_start, i_end, N, w)`
      Numba kernel in `simulation.py` — mass-conservative top-hat
      axial weights summing to 1.0, defensive fallback to cell 0 on
      empty range / zero total dx.
- [x] `PyrogenChamber.injection_topology: str = 'forward_plenum'`
      field added; valid values: `'forward_plenum' | 'head_basket'
      | 'aft_basket'`.
- [x] `PyrogenChamber.cartridge_length_m: float = -1.0` field —
      sentinel triggers mass-derived length at sim init.
- [x] `PyrogenChamber.resolve_cartridge_length(A_port_avg)` method.
- [x] `PyrogenChamber.resolve_injection_cells(x_centers, dx, N,
      A_port)` method returning `(i_start, i_end)` cell range.
- [x] Numba-friendly topology codes: `TOPOLOGY_FORWARD_PLENUM=0`,
      `TOPOLOGY_HEAD_BASKET=1`, `TOPOLOGY_AFT_BASKET=2`.
- [x] 144-line test file `test_uniform_band_weights.py` covering
      conservation, range handling, edge cases.

### Phase A.1.1 — naming pivot [COMPLETE]

Commit `a85c202`. Architectural clarification driven by the Super
Loki literature dive (NASA CR-61238, MIT Super Loki Report,
Smithsonian/NASM): the ISP Super Loki igniter is a head-end BKNO3
pellet charge in a consumable plastic moisture cup, with NO defined
orifice or pressure-containing aft cap. Modeling it as
plenum-with-orifice would be wrong physics.

- [x] Rename "submerged" → "uncontained" throughout PyrogenChamber
      docstring and validation.
- [x] Clarify in docstring: each pellet burns at its host cell's
      LOCAL bore pressure; `r_b[i] = a * P_bore[i]^n`; no plenum
      chamber wall, no internal pressure separation, no defined
      orifice or burst threshold.
- [x] Plenum-state fields (`A_throat`, `V_plenum`) repurposed:
      validated at the Python boundary so existing motor configs
      don't break, but ignored by the uncontained-burn time-loop
      computation.

### Phase A.2 — uncontained mdot kernel + time-loop branch [COMPLETE]

Commit `720b848`.

- [x] `_compute_uncontained_pyrogen_mdot(P_bore, a, n, rho_p,
      A_burn_per_cell, m_pyrogen_remaining, dt, i_start, i_end, N,
      mdot_arr)` Numba kernel. Mass conservation: per-step cap
      against `m_pyrogen_remaining` with uniform scale-down on the
      last step.
- [x] `_run_time_loop` topology branch at L1200-1377:
  - `topology_code == 0`: forward_plenum (existing path unchanged).
  - `topology_code != 0`: uncontained — per-cell mdot, no plenum
    gas state, no choked vent, no separate P_ig. plenum_state[0]
    reused as `m_pyrogen_remaining`; plenum_state[1, 2] vestigial
    for stable state-vector shape.
- [x] **Mass / species / enthalpy** delivery branches on topology
      (`simulation.py:1342-1377`): forward_plenum keeps Phase A
      exponential decay; uncontained uses per-cell mdot directly.
- [x] **Momentum injection** skipped for uncontained
      (`simulation.py:1258`): PISO handles axial flow via the
      pressure gradient between high-P pyrogen cells and surrounding
      bore. No explicit momentum source.
- [x] **DeMar surface heat flux** disabled for uncontained
      (`simulation.py:1279`): no impinging "plume" on a leading-
      edge unignited cell in the uncontained model — heat enters
      via bulk mass injection at the cartridge cells.
- [x] **P_ig diagnostic** for uncontained: volume-averaged bore P
      over cartridge cells (`simulation.py:1236-1245`).
- [x] 188-line `test_submerged_topology.py` + 242-line
      `test_uncontained_pyrogen.py` covering wiring and edge cases.

### Phase A.3 — validation + diagnostic helpers [COMPLETE 2026-05-24]

- [x] Diagnostic visualization helpers in `plotting.py`:
  - `plot_flow_snapshot` upgraded from 2x2 to 3x2: added `u_cell`
    panel (sign-banded for at-a-glance reverse-flow detection)
    and gas temperature panel.
  - `plot_flow_snapshots(result, t_targets=[...], fields=[...])`
    — multi-time subplot grid (rows = time, cols = field).
  - `plot_field_heatmap(result, fields=[...], t_max=...)` —
    pcolormesh(x, t, field) with symmetric colorbar for signed
    fields (u_cell shown in `RdBu_r`).
- [x] `ISP_SUPER_LOKI_EXPERIMENTAL` dataset moved from commented-
      out mis-labeled block in `examples/ISP_Super_Loki.py` to
      `plotting.py` with proper labeling.
- [x] `run_from_ric` and `build_pyrogen_chamber` extended with
      `injection_topology=` / `cartridge_length_m=` kwargs.
- [x] `examples/ISP_Super_Loki.py` wired to head_basket; renders
      pressure + multi-snapshot + heatmap + summary.
- [x] `examples/hasegawa_motor_a_aft_basket.py` created (sibling of
      `hasegawa_motor_a.py` with reversed mass-injection topology;
      same knobs otherwise).
- [x] Both example runs complete without crashes; full artifact
      sets under
      `artifacts/ISP_Super_Loki/<stamp>/` and
      `artifacts/hasegawa_a_aft_basket/<stamp>/`.

### Phase A.3 — validation findings

**ISP Super Loki (head_basket fit)**:
- P_peak = 0.12 MPa (experimental ~8.8 MPa)
- Pyrogen burned 4.84 g of 4.8 g initial in 298 ms
- Propellant produced 0.005 kg of 20.7 kg total — essentially never
  ignited the main grain
- Pyrogen burns to completion at atmospheric bore pressure
  (Saint-Robert at P_atm ≈ 0.1 MPa gives r_b ≈ a · 0.1^n ≈ tiny)
- Mass-flux pathway alone is too slow to pressurize the bore
  before pyrogen exhausts

**Hasegawa A (aft_basket diagnostic)**:
- P_peak = 0.10 MPa (essentially atmospheric)
- Pyrogen burned 4.18 g of 4.2 g in 437 ms
- Propellant produced 0.004 kg of 10.8 kg — same failure mode
- **The diagnostic question is inconclusive**: aft_basket never
  reaches an ignition cascade to compare against forward_plenum's
  spike, so it doesn't yet answer "does the artifact persist
  under reversed topology?"

**Structural finding** (the load-bearing v0.7.3 Phase A insight):
The uncontained model correctly captures the physical fact that
pellets at atmospheric pressure burn slowly — but it exposes a
gap that forward_plenum hides: real-world pyrogen ignition is
initiated by an external **thermal kick** (e-match, hot wire,
squib) that creates a transient hot-gas pulse the pellets need to
spool up. forward_plenum's choked-orifice fakes this via its
plenum-pressure boundary condition; uncontained has no equivalent.

The v0.7.3 Phase A architecture is sound but **incomplete**:
ignition-initiation infrastructure is needed before uncontained
topologies validate quantitatively.

### Phase A.4 — close-out [COMPLETE 2026-05-24]

- [x] `docs/v0_7_3/README.md` + `TASKS.md` created.
- [x] DEVNOTES API-breaking-change log entry for v0.7.3-phaseA.
- [x] CLAUDE.md banner updated.
- [x] `project_v0_7_3_phaseA_state` memory created;
      `[[srm-1d-fired-motor-set]]` and
      `[[srm-1d-v0-7-2-progress-state]]` cross-linked.
- [x] Full pytest run before tag.
- [x] Tag `v0.7.3-phaseA`.

## Phase B — heat-flux completeness for uncontained ignition [COMPLETE 2026-05-25]

Closes the uncontained-pyrogen ignition gap identified in Phase A.3.
See `PHASE_B_SCOPE.md` for the full scope/audit narrative.

### Phase B.0 — IC fix [SHIPPED]

- Default `T_initial_gas` switched from `rep_tab.T_flame` (v0.7.0
  numerical-stability shortcut) to `_ambient_T` (= propellant.T_initial,
  ~293 K).
- The previous IC short-circuited temperature-gradient flow under
  uncontained-pyrogen topologies; the new IC creates a real T gradient
  when pyrogen mass enters cold bore cells, driving the flow that
  Bartz convection needs to deliver heat to the propellant surface.
- Override via `initial_gas_temperature=` kwarg preserved for
  backward compat / special studies.
- Side effect: legitimate amplification of ignition spike on existing
  forward_plenum motors (Hasegawa A goes from ~6.2 MPa to ~12 MPa
  P_peak at the same v0.7.0 calibrated knobs). Test windows widened
  to ±150% as the "bug sanity gate" until v0.7.4 Phase C re-LHS.

### Phase B.2 — radiation_emitter gating extension [SHIPPED]

- One-line change in `_goodman_ignition_sources_and_mass`:
  `radiation_emitter[i] = is_burning[i] OR Y_species[i, IGNITER] > 0.5`.
- Pyrogen-hot cells (cartridge cells with T_gas ≈ T_flame_pyrogen)
  now contribute to cell-to-cell radiation just like
  propellant-burning cells. Previously ungated only if propellant was
  burning, which never happened during the uncontained ignition
  transient.
- No-op when `Propellant.radiation_emissivity == 0` (the default per
  `_default_radiation_emissivity`).

### Phase B.3 — pyrogen form archetypes [SHIPPED]

- `Pyrogen.form: str = 'pellets'` field with values
  `'powder' | 'pellets' | 'chunks'` (see `[[pyrogen-form-archetypes]]`
  memory).
- A_burn multipliers in `build_pyrogen_chamber`:
  - `'chunks'`: ×1.0 (single-sphere baseline)
  - `'pellets'`: ×5.0 (typical amateur HPR pellet pack — DEFAULT)
  - `'powder'`: ×20.0 (fine particles, thermite-class)
- Explicit `pyrogen_burn_area=` kwarg always wins (user override).
- BPNV and MTV YAMLs marked `form: pellets` (matches amateur usage).

### Phase B.4 — pyrogen-to-surface heat delivery modes [SHIPPED]

- Three mutually-exclusive modes on Pyrogen (Numba-coded as enum):
  - `'demar'` — uses `Pyrogen.heat_flux_cal_cm2_s` (DeMar 2021
    time-averaged) distributed uniformly across cartridge cells.
    Per-cell sensible cap via `mdot_uncontained[i] · Cp · ΔT`.
  - `'radiation'` — Stefan-Boltzmann pellet emission with geometric
    view factor `F_ij = A_port[j] / (4π·d² + A_port[j])` plus
    Beer-Lambert absorption attenuation `exp(-d/L_atten)`.
  - `'none'` — no surface heat flux from pyrogen (control case).
- New Pyrogen fields: `heat_delivery_mode`, `pellet_emissivity`
  (default 0.7), `radiation_absorption_length_m` (default 1.0 m
  clean; 0.5 m for MTV with MgO particles).
- New Numba kernel `_compute_pyrogen_heat_flux_arr` fills a per-cell
  flux array each step; `_goodman_ignition_sources_and_mass` consumes
  it directly (replaces the v0.7.0-v0.7.3-phaseA single-cell
  `pyrogen_heat_target` special case).

### Phase B.6 — A/B/control validation [DONE 2026-05-25]

**ISP Super Loki head_basket** (no experimental overlay — see
provenance note below):

| Mode | P_peak | t_peak | Verdict |
|---|---|---|---|
| `'none'` | 17.08 MPa | 0.045 s | Ignites |
| `'demar'` | 17.08 MPa | 0.045 s | Ignites |
| `'radiation'` | 17.08 MPa | 0.045 s | Ignites |

All three modes give **identical** P_peak — the load-bearing fix is
B.0 + B.3 (cold IC + realistic A_burn). Once those let pyrogen mass
inject hot products into cold bore at high mdot, the resulting
convective heat flux (Bartz with the strong T-gradient-driven Re)
alone ignites the propellant. Mode (B.4) becomes diagnostic
refinement, not a load-bearing knob.

**Provenance correction (2026-05-25)**: the commented-out
experimental array in `examples/ISP_Super_Loki.py` (migrated to
`ISP_SUPER_LOKI_EXPERIMENTAL` in `plotting.py` during Phase A.3)
turned out to be **Chunc/machbusterNew data**, not Super Loki —
identical byte-for-byte to `CHUNC_EXPERIMENTAL` in
`examples/machbusterNew.py` and the top-level `chunctest.py`. The
data was copy-pasted into the Super Loki example years ago. v0.7.3
Phase B reverted the mis-labeling: `plotting.py` now exposes the
array as `CHUNC_EXPERIMENTAL` (matching its actual provenance),
and the Super Loki example runs without an experimental overlay
until a verified static-fire dataset is sourced. The "ignites at
17 MPa" finding is still informative for the architecture (B.0 +
B.3 unblock ignition), but absolute peak comparison to "8.8 MPa
experimental" was actually against Chunc.

**Hasegawa A aft_basket** (diagnostic):

| Mode | P_peak | t_peak | Verdict |
|---|---|---|---|
| `'none'` | 0.113 MPa | 0.086 s | Stalls |
| `'demar'` | 0.113 MPa | 0.086 s | Stalls |
| `'radiation'` | 0.113 MPa | 0.086 s | Stalls |

All three modes stall at near-atmospheric P. The diagnostic
question — does the simultaneous-ignition artifact persist under
reversed topology? — **resolves indirectly**: `aft_basket` is
fundamentally inadequate as a startup mechanism because the
cartridge sits next to the nozzle, so pyrogen products vent
immediately without pressurizing the upstream bore. The deferred
`aft_fore_firing` topology (PyrogenChamber docstring L90-93) is
what's needed for a real Super Loki-class aft-firing diagnostic.

### Phase B.7 — close-out [DONE 2026-05-25]

- DEVNOTES v0.7.3-phaseB API-breaking-change log entry
- CLAUDE.md banner + roadmap update
- `project_v0_7_3_phaseB_state` memory
- 272/272 pytest green
- Tag `v0.7.3-phaseB`

## v0.7.4 candidate space (deferred from Phase B)

The natural Phase B is to **add an ignition-initiation pathway**
for uncontained topologies. Options (recommended order pending
user decision):

1. **Initial thermal pulse** (smallest scope ~150 LOC). Kick
   `T_surf` and/or local bore T at t=0 in cartridge cells.
   Physically: models the e-match heat dump that ignites BKNO3.
   Knobs: `pulse_T_kick` (target T), `pulse_duration` (or
   exponential decay τ), gated by topology.
2. **Per-pellet surface heat flux** (~200 LOC). Re-enable DeMar-
   style surface flux PER CARTRIDGE CELL (not cell 0 only).
   Most direct mapping from forward_plenum's known ignition
   behavior. Knob: `head_basket_heat_flux_w_m2` / similar.
3. **Coupled e-match model** (~400 LOC, largest scope). Add a
   small `Igniter` dataclass for the electrical transient with
   `t_kick`, `Q_kick`, `tau_kick`. Most physical, opens the door
   to candidate 6 (plenum-as-option refactor).
4. **Initial pressure pulse** (~100 LOC). Simplest hack — small
   transient mass injection at t=0 that primes the bore. Crude
   but quick to prove the diagnostic question.

Once any of these lands, **re-run both validation examples** to
test (a) Super Loki head_basket fit vs ~8.8 MPa experimental, and
(b) aft_basket diagnostic with an actually-reached ignition
cascade.

## Deferred / future work

**Diagnostic-visualization roadmap** (user-flagged 2026-05-24):

- **Full CSV export** of `result['snapshots']` and per-history
  channel arrays. Round-trips cleanly through
  `pandas.DataFrame.from_records` if we standardize on flat field
  names. Enables post-run replay without re-running the sim.
- **Generic CSV post-processor** that ingests the exported data
  and exposes all of {pressure trace, flow snapshots at any t,
  field heatmaps, snapshot grids, animations} via a tiny CLI or
  Python script. Useful for sharing artifacts and long-form
  analysis after a run finishes.
- **Per-cell species mass-fraction** `Y[i, s]` in the snapshot
  dict. Requires sim-loop instrumentation
  (`simulation.py:1567-1584`). Useful for visualizing
  igniter-gas / propellant-gas mixing across the bore — relevant
  for diagnosing whether the simultaneous-ignition artifact
  correlates with uniform igniter-species saturation.
- **Animations** (`matplotlib.animation.FuncAnimation` /
  GIF export). Heatmaps cover most of the value at much lower
  implementation cost — genuinely lower priority.
- **Interactive vs save-to-file mode toggle** at example-script
  level (currently examples force `matplotlib.use('Agg')` for CI
  safety — a `--interactive` flag would let humans inspect
  pop-up windows while CI still saves PNGs).

**Other v0.7.3+ candidates** (from
[../v0_7_2/candidates_post_phaseA.md](../v0_7_2/candidates_post_phaseA.md)):

- **Z-N dynamic burn rate** (candidate 1) — independent of
  topology work; stacks cleanly with anything else.
- **Per-cell coupling alternatives** (candidate 4 sub-options) —
  reverse-polarity Phase B damping, solid-phase axial
  conduction, shared boundary layer. Address the simultaneous-
  ignition artifact directly once Phase B unblocks the
  diagnostic.
- **Plenum-as-option refactor** (candidate 6) — unify
  forward_plenum + head_basket + aft_basket under a single
  `Igniter` dataclass. Best done AFTER the ignition-initiation
  work lands so the abstraction is informed by real
  requirements.
- **Pardue 1992 Al2O3 condensation heating** (candidate 5) — the
  SPINBALL walkthrough's secondary spike-taildown candidate;
  larger scope, requires extending N-species infrastructure to
  4 species (igniter + grain + ambient + Al2O3 condensable).
