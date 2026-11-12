# Changelog

Versions are git tags. Per-change API details (signatures, removed or renamed
keys, migration notes) are in `srm_1d/DEVNOTES.md`, "API Breaking Changes
Log". Each entry below names its log entry or worklog in brackets.

## Unreleased: 0.9.0-dev (branch `openmotor-frontend`, since v0.8.1)

This cycle verified the solver against exact solutions and rebuilt the
ignition and igniter models. It found and fixed several conservation and
pressure-solver bugs.

**Most defaults changed, so results differ from v0.8.1**, especially the
startup transient. The calibration constants (roughness 32 µm, κ 0.44,
T_ignition 756 K, k_solid 0.271) were fitted under the v0.7.5/v0.8.1 defaults
and have **not yet been re-fitted**. The startup pressure spike is still
over-predicted on the reference motors (see "Validation state" below); it is
the main open problem.

### Default changes

| Area | v0.8.1 | Now | Why / record |
|---|---|---|---|
| PISO second pressure corrector | omitted the storage term, so transient pressure increments were applied about twice | `piso_corr2_storage=True` | Bug fix; startup spike about −11 %. [T11] |
| Energy-pressure coupling | none | `piso_energy_coupling_iters=1` | [T11] |
| Energy equation | sensible enthalpy | `energy_formulation='total'` (internal + kinetic) | The only pairing verified against exact Sod, shock and vessel solutions. [EN1, D1] |
| Momentum predictor | `'legacy'` | `momentum_formulation='conservative'` | Pairs with total energy. [EN1, D1] |
| Ignition heat transfer | Ma / Gnielinski chain (3-5× high at a cold wall) | `ignition_heat_transfer='peretz'` (Peretz et al. 1973, local port area) | [IG3 stage 1] |
| Ignition driver temperature | gas static temperature | `ignition_driver_temperature='recovery'` (adiabatic-wall recovery temperature) | [RT1, D1] |
| Head-end igniter jet heating | not modelled | `head_end_jet_h='max'` (Peretz Eqs. III-28..30) | Removes the head-end dead zones. [HE1, D1] |
| Solid heating before ignition | Goodman integral | `solid_heating_model='conduction'` (exact 1-D transient conduction) | Goodman's surface temperature jumped whenever h changed. [SM1, SM2, conduction default] |
| Igniter surface heat flux (DeMar) | on | off | It is a time-averaged sizing figure, not a measured flux. [IG3 stage 2] |
| Auto igniter charge | Sutton m = 0.12·V^0.7 | DeMar M = V·P/W at 300 psi, raised ×1.25 per attempt until ≥ 95 % of the grain lights within 0.1 s | [IG3, IG4, AS1, AW1] |
| Pyrogen pellets | archetype area multipliers | MIL-P-46994B Type I-A cylinders with exact burn-back; head-basket area shrinks as pellets burn | [OB1] |
| Motor ends | 3-cell full-diameter gas buffers at head and aft | Physical ends: head gap set by head-face inhibition and lumped into the first cell; nozzle convergent lumped into the last cells | Startup is now grid-independent. [G1] |
| Grain-to-grid mapping | length snapping changed burn area; gaps had their own cells | Snapping preserves burn area and volume; gap volumes lumped; every lump spread so no cell's storage more than doubles | [G2] |
| Nozzle boundary | last cell's static P, T used as stagnation values | `nozzle_bc_stagnation='isentropic'`: true stagnation state from the port-exit Mach | Plateau −0.7 to −1 %, startup −4 to −10 %. [NB1] |
| Thrust / Isp / impulse | used head-end pressure | use the nozzle-entrance stagnation pressure `P0_nozzle`; head-end pressure stays the chamber-pressure trace | Simulated C* now matches the propellant's to 0.1 %. Impulse −3 to −6 %. [TP1] |

### Fixes

- **Mass and energy:**
  - The regression mass source is now (ρ_p − ρ_g)·r·A. [T10]
  - The left-neighbour radiation sink is debited from the gas. [T9]
  - The total-energy fill debit carries kinetic energy. [QF1]
- **Ma erosive burn rate:** the bisection has a guaranteed upper bracket.
  [T2]
- **Nozzle input:** the `slagCoeff` unit conversion was wrong. [slag
  coefficient]
- **.ric igniter blocks:** they round-trip losslessly, warn about stale
  igniter materials, and explicit sizing arguments to `run_from_ric` now
  override the block. [T1, run_from_ric precedence]
- **Pressure cutoff:** it no longer ends a slow start as a "burnout" before
  the grain has lit. [AS1]
- **Fine grids:** a fill collapse from a large lumped nozzle convergent is
  fixed. [G2]
- **History memory:** the history buffer is sized from the hottest sound
  speed and bounded (no 5M-row floor). [history capacity, MEM1]
- **GUI:** results frames are now 1 ms, and the station slider is linked to
  the frames.

### New features and options

- **Erosive-model registry:** `erosive_model` ∈ {`ma2020` (default),
  `mukunda_paul`, `saint_robert`}, plus an opt-in `erosive_mp_cap`.
- **Igniter supplied by the plugin:** srm_1d supplies the igniter library and
  the openMotor GUI solver-config defaults, and upgrades stale saved configs.
  [I1, I4]
- **Ignition diagnostics in `summary`:**
  - `prompt_start`, `prompt_lit_fraction`, `ignition_t50_s`,
    `ignition_t90_s`;
  - never-lit grain cells;
  - the charge attempts, and the charge's DeMar-equivalent pressure.
  - `PyrogenAutoSizeError` when no charge lights the motor; the plugin shows
    WARNING / ERROR alerts.
- **New result channels:** `P0_nozzle`, `T0_nozzle`, `gamma_nozzle`,
  `R_nozzle`; C* from the simulation. The CSV export gains a nozzle
  stagnation-pressure column.
- **Per-motor head-gap override:** the `run_from_ric(head_gap=...)` argument,
  the igniter-block field `head_gap_m`, and a plugin chamber field.
- **Solver diagnostics:**
  - per-cell numerical-limit diagnostics;
  - a `pressure_floor_pa` argument;
  - energy-audit ledger channels (including p·dV);
  - split burn and geometry refresh cadences.
- **`srm_1d.tools.guarded_run`:** runs a job under a whole-process-tree memory
  cap and a timeout. Sensitivity sweeps now default to 2 workers.

### Breaking API changes (summary)

See the DEVNOTES log for each one.
- `_run_time_loop` and several other `@njit` kernels changed signature
  repeatedly. Rebuild any external call against the current code.
- `run_simulation` gained the keyword arguments listed above;
  `ignition_retry_*` and the prompt-start arguments wrap the single run.
- `diagnostic_history_capacity` was renamed to `history_capacity`.
- For motors with inter-grain gaps, `L_motor` is now the propellant length;
  the physical length is `L_motor + sum(gap_after)`.

### Validation state at the current defaults (100 cells)

| Motor | Simulated | Experiment |
|---|---|---|
| Chunc | startup peak 11.09 MPa at 23 ms, plateau 8.16; peak/plateau **1.36** | 1.027 |
| BALLSstick | peak/trough 1.94, peak/mid 1.74 | ~1.50 / ~1.35 (2-inch subscale) |
| Hasegawa A | late peak 5.957 MPa (0.925×) | 6.44 MPa |
| pathfin54 | about 6 % low in thrust and impulse; peak pressure 8.6 % low | measured thrust and pressure |

On pathfin54, thrust divided by pressure agrees to 1.5 %, so the nozzle side
is consistent with the measurement.

### Repository

- **Reference PDFs are no longer distributed:** `docs/references/*.pdf` is
  untracked and gitignored. `docs/references/README.md` lists the papers;
  obtain them from the publisher.
- **History note:** to keep four PDFs out of the public history, the commits
  made since 2026-10-05 were rewritten before this push. Worklogs written
  earlier cite the old commit hashes. `docs/HISTORY_REWRITE_2026-11.md` maps
  old hashes to new.
- **New docs:**
  - `docs/solver_verification/`: verification specs, worklogs and the B0
    baseline gate;
  - `docs/pyrogen_realism/`: ignition and igniter reports;
  - `docs/fable/`: independent-audit prompts.

## v0.8.1 (2026-06-19)

Parametric grain taper (bore and OD/end, in both solvers and the GUI); the
longitudinal slice viewer; `fastmath` on the hot kernels (about +30 %); the
contributor guide. See `CLAUDE.md` and `docs/core_loop_opt/`.

## v0.8.0

openMotor frontend integration (solver-plugin contract, channels), per-station
axial visualisation, igniter as data, v0.7.5 cross-motor calibration defaults.
See `docs/v0_8_0/`.
