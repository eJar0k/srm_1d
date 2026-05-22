# srm_1d v0.7.1 — Multi-Species Bore Gas (N-species)

**Status**: Phases 1 + 2 + 3 + 4 complete (2026-05-23). 199/199 tests
pass. PISO consumes per-cell `(γ, R, Cp, T_ceiling)` and advects
sensible enthalpy; Phase 4 validation confirms the mixture arrays
collapse to the correct single-species thermo in the pure-pyrogen and
pure-propellant limits. Two follow-ups queued behind Phase 5: strict
T_ceiling formula (DESIGN §5 with IC guard) and per-species Cp at
source sites (Phase 3.5). Hasegawa A baseline smoke is within ±10% of
v0.7.0 without re-calibration.

**Target**: Replace the single calorically-perfect gas assumption (v0.7.0
"must-have" #5 from [v0_7_0/DESIGN.md](../v0_7_0/DESIGN.md) §3) with a
generalizable **N-species mixture** in the bore. Each cell carries
mass fractions `Y[i, s]` for `s = 0..S-1` species; per-cell γ, Cp, R
derive from these. The current scalar (γ, R, Cp) become per-cell arrays.

Forward-looking: v0.8.0's head-end primary motor architecture will need
at least 3 species (head-motor combustion products, sustainer combustion
products, igniter combustion products). Building the N-species
infrastructure now avoids a 2-species → N-species refactor later.

## Why this is in v0.7.1

The post-v0.7.0 [SPINBALL research](../post_v0_7_0/references/spinball_walkthrough.md)
established that Cavallini's "infinite-gases mixture" formulation is the
single architectural feature SPINBALL has that we lack. It is **not** a
spike-taildown solution (Z-N is, see [project_spinball_research_state](https://...)
memory) but it:

1. Addresses the frozen-vs-effective gas-transport tension principledly
   rather than via a single calibrated `k_solid` / `k_thermal` knob.
2. Provides the chamber-fill realism that SPINBALL highlights for Z23.
3. Sets up v0.8.0 head-end motor architecture.
4. Is conceptually cheap — the implementation cost is mostly
   "thread arrays through the solver where scalars are now," not new
   physics. No new PDE, no new fitted constants.

## Architectural decisions

### 1. Gas species as first-class objects

A new `GasSpecies` dataclass in [propellant.py](../../propellant.py)
captures the thermophysical state of one gas species at the bulk-flow
level:

```python
@dataclass
class GasSpecies:
    name: str                # e.g. "hasegawa_prop1_gas", "bpnv_pyrogen_gas"
    gamma: float             # ratio of specific heats [-]
    Cp: float                # specific heat at constant pressure [J/(kg·K)]
    molecular_weight: float  # mean MW [kg/mol]
    T_flame: float           # adiabatic flame temperature [K], source-only
```

`R_specific` is exposed as a `@property` derived from
`R_UNIVERSAL / molecular_weight`. Transport properties (`k_thermal`,
`mu_gas`, `Pr`) are intentionally omitted from `GasSpecies` for v0.7.1
because the solver keeps them scalar (§6). v0.7.2 may add optional
`k_thermal: float = 0.0` and `mu_gas: float = 0.0` fields when per-cell
transport is wired in.

A species **does not** carry burn-rate parameters; those stay on
`Pyrogen` (for the igniter chamber) and `PropellantTab` (for the main
grain). What changes:

- `Pyrogen.species` (property) returns a `GasSpecies` built from
  the pyrogen's (γ, M, T_flame). `Cp_ig` derives from
  `Cp = γ · R_universal / [M · (γ - 1)]` if not given explicitly.
- `Propellant.species` (property) returns a `GasSpecies` from a
  representative tab. Multi-tab Cp-vs-pressure remains a deferred
  improvement (see [DEVNOTES](../../DEVNOTES.md) gotchas).

### 2. Species registry built at sim start

`run_simulation` constructs a list `species_array` of `GasSpecies`
objects (length S). Each mass source declares which species index it
contributes to:

- **Pyrogen plenum** → species 0 (always; igniter is always first by
  convention so future code that special-cases it stays simple).
- **Main grain combustion** → species 1.
- **Future v0.8.0 head-end motor** → species 2.
- **Future ablation** → species 3.

The registry is held in a fixed-size Numba-friendly 2D array
`species_params[S, 4]` carrying `(γ, Cp, M, T_flame)` per species.
Transport properties stay scalar in the hot loop for v0.7.1 (see §6).

**Default species lineup for v0.7.1** (`S = 3`):

| s | name | role | source after t=0 |
|---|---|---|---|
| 0 | igniter | from `pyrogen.species` | pyrogen plenum injection at cell 0 |
| 1 | propellant | from `propellant.species()` | grain combustion (normal + erosive + endface) |
| 2 | ambient | air-like defaults (γ=1.40, Cp=1005, M=0.0289, T_flame=`T_initial`) | none — initial condition only, purged through the nozzle |

The ambient species exists so the pre-fill mass is properly accounted
for as it gets pushed out the nozzle during chamber fill. With this
arrangement, **per-species conservation is testable**: cumulative
species mass in the bore + cumulative species mass out the nozzle =
cumulative species mass injected (zero for ambient after the initial
fill).

### 3. Y[N, S] state array

A new state variable `Y_species` shape `(N, S)` tracks the mass fraction
of each species in each cell. Invariant: `sum(Y[i, :]) = 1.0` at all
times.

**Initial condition.** `Y[:, 2] = 1.0` — all bore cells are 100%
ambient species at sim start. The initial gas *temperature* keeps the
v0.7.0 default (`T_initial_gas = rep_tab.T_flame`, overridable via
`initial_gas_temperature`) for numerical stability. The shift is
**species composition**, not temperature: pre-fill thermo derives from
the ambient species (γ_air, Cp_air, M_air) rather than the propellant
species. This is mildly unphysical (cold air composition at hot
temperature) but the pre-fill mass is small and is purged within ~1 ms
of pyrogen flow. The benefit is per-species mass conservation testing.

If physically-correct ambient initialization is needed later
(`T_initial_gas = T_initial`, ambient species, ambient pressure),
that's a follow-up improvement; v0.7.0 numerical-stability behavior is
preserved here intentionally.

**Source update** (per timestep, before transport):
```
for each cell i with mass source mdot_s contributing to species s_src:
    m_old = rho[i] * V[i]
    m_added = mdot_s * dt
    Y[i, :] = Y[i, :] * m_old / (m_old + m_added)
    Y[i, s_src] += m_added / (m_old + m_added)
```

In words: existing mixture stays in proportion, new species adds its own
mass fraction.

**Transport update** (per timestep, after sources): standard
finite-volume update using the **same face mass fluxes** the PISO
continuity step uses. Upwinded on face velocity. No new PDE — the
existing fluxes determine which neighboring Y_row enters each cell.

```
for each cell i, for each species s:
    Y_new[i, s] = (Y_old[i, s] * m_old[i]
                  + flux_in_left * Y_upwind_left[s]
                  - flux_out_left * Y_old[i, s]
                  + flux_in_right * Y_upwind_right[s]
                  - flux_out_right * Y_old[i, s]) / m_new[i]
```

Cleaned-up after summing: a passive scalar advected by the existing
mass fluxes. Mass conservation guarantees `sum(Y_new[i, :]) = 1.0` to
roundoff if implemented carefully.

### 4. Per-cell mixture derivation

After Y is updated, compute per-cell (γ, Cp, R, M) arrays:

```
Cp_mix[i]   = sum_s ( Y[i, s] * Cp[s] )                       (mass-weighted)
1/M_mix[i]  = sum_s ( Y[i, s] / M[s] )                        (harmonic in molar)
R_mix[i]    = R_universal / M_mix[i]
gamma_mix[i]= Cp_mix[i] / (Cp_mix[i] - R_mix[i])              (ideal gas)
```

These are standard textbook ideal-gas mixing rules for mass fractions.
Computed once per timestep; cached as 1D arrays of length N.

### 5. Solver threading

Convert every scalar `(gamma, R_specific, Cp_gas)` callsite in
[solver.py](../../solver.py) and [simulation.py](../../simulation.py) to
take arrays:

- `_piso_step(... gamma_arr, R_arr, Cp_arr ...)` — indexed per cell.
- `_signed_choked_or_subsonic_flow(..., gamma_b, R_b, ...)` — uses the
  upstream cell's mixed thermo at the boundary face.
- `_critical_pressure_ratio(gamma_b)`, `_critical_flow_function(gamma_b)`
  — pure functions, take a scalar; called once per boundary face per
  step with the appropriate `gamma_b`.
- `_subsonic_throat_mdot_mag(..., gamma_b, R_b)` — similarly.
- EOS update: `rho_new[i] = P_new[i] / (R_arr[i] * T_new[i])`.
- Energy audit: `sensible[i] = rho[i] * Cp_arr[i] * T[i]`.
- `_pyrogen_surface_heat_power(..., Cp_arr[i_target])` — use the
  *cell's* Cp, since surface heat flows to/from the cell's gas mixture.
- `_goodman_ignition_sources_and_mass(..., Cp_arr ...)` — uses per-cell
  Cp wherever it currently uses scalar Cp_gas.
- T_ceiling clipping at [solver.py:612](../../solver.py#L612): switch
  from `T_flame * 1.01` to `T_ceiling_arr[i]` where:
  ```
  T_ceiling_arr[i] = max(T_flame[s] for s such that Y[i, s] > 0.05) * 1.01
  ```
  i.e., the ceiling reflects whichever species's combustion produced
  most of cell i's mass.

### 6. What stays scalar in v0.7.1

For v0.7.1, the following remain scalars (frozen at sim start from a
representative species):

- `k_thermal`, `mu_gas`, `Pr` — used in Ma 2020 erosive burning. These
  are transport properties; refactoring them per-cell doubles the work
  and the frozen-vs-effective sensitivity is mostly captured by γ/Cp
  variation. Deferred to v0.7.2 if calibration shows it helps.
- `roughness`, `kappa` — calibration constants, propellant-side, not
  gas-side.
- `T_flame` as a **source temperature** when grain combustion injects
  mass — uses propellant tab's T_flame (this is the temperature of the
  combustion *products*, which is intrinsic to the propellant, not a
  bore-mixture property). Same for pyrogen's T_ig.

### 7. Numerical stability considerations

- **Energy conservation across mixing interfaces.** When two cells with
  very different Cp share a face, advecting T directly can fail to
  conserve energy. **Fix**: advect `(Cp * T)` (sensible enthalpy) and
  back out T at the cell update. This is the standard fix in multi-Cp
  finite-volume codes.
- **Y bound enforcement.** `Y[i, s]` should stay in [0, 1] and sum to 1.
  Floating-point drift can violate this over thousands of steps. **Fix**:
  at the end of each step, clamp negative Y to 0, renormalize so
  `sum_s Y[i, s] = 1`. Cheap and robust.
- **Numba cache** — kernel signatures change substantially. Per
  [CLAUDE.md](../../CLAUDE.md) gotcha #1, clear `srm_1d/__pycache__/`
  after the first refactor commit lands.
- **Source-CFL cap** with variable Cp: cap currently uses scalar Cp.
  Switch to using the cell's Cp_arr[i]. The cap formula at
  [solver.py:762](../../solver.py#L762) remains structurally the same;
  just the input scalar becomes an indexed array.

## Phase plan

### Phase 1 — Species registry, Y[N, S] state, transport [COMPLETE 2026-05-23]

**Files touched**: `propellant.py` (GasSpecies, Pyrogen.species,
Propellant.species, ambient_air_species, species_array),
`simulation.py` (state allocation, source updates, transport update),
`tests/test_simulation_phase3.py` (signature update for the existing
`_goodman_ignition_sources_and_mass` test fixtures), new test file
`tests/test_yns_transport.py`.

**Delivered**:
- `GasSpecies` dataclass + `Pyrogen.species` property +
  `Propellant.species()` method + `ambient_air_species(T)` helper +
  `species_array(list)` packer.
- 3-species registry built in `run_simulation` (igniter, propellant,
  ambient).
- `Y_species[N, 3]` allocated, initialized to 100% ambient.
- `mass_source_by_species[N, 3]` parallel array; pyrogen injection
  writes s=0, `_goodman_ignition_sources_and_mass` writes s=1 for
  normal + erosive + endface grain mass.
- `rho_pre_step[N]` snapshot before PISO so the post-PISO advection
  has the old density.
- `_advect_species` Numba kernel: mass-fraction-conservative upwind,
  face density `0.5*(ρ_old[j-1] + ρ_old[j])` matching PISO, nozzle
  face carries `Y[N-1, :]` at `nozzle_mdot`, per-cell renormalize.
- Result dict gains `Y_species_final / species_params / species_names
  / rho_final / A_port_final`.
- 5 transport tests in `tests/test_yns_transport.py` cover pure
  advection, source-only, FP-drift renormalize, nozzle drain, and a
  50 ms Hasegawa A integration.

**Solver behavior unchanged this phase**: γ/R/Cp still scalar in the
solver. Y is tracked + advected but not consumed.

### Phase 2 — Per-cell mixture derivation [COMPLETE 2026-05-23]

**Files touched**: `simulation.py` (mixture helpers), new test
`tests/test_yns_mixture.py`.

**Delivered**:
- `_compute_mixture_cell(Y_row, species_params) → (γ, Cp, R, M)`
  per-cell Numba helper.
- `_refresh_mixture_arrays(Y, species_params, gamma_arr, Cp_arr,
  R_arr, M_arr, N)` per-step looper.
- Per-step recompute of `gamma_mix_arr[N]`, `Cp_mix_arr[N]`,
  `R_mix_arr[N]`, `M_mix_arr[N]` after each `_advect_species` call.
- Result dict gains `gamma_mix_final / Cp_mix_final / R_mix_final /
  M_mix_final`.
- 8 mixture tests: pure-species limits (all 3), 50/50 binary mix,
  γ=Cp/(Cp-R) consistency, degenerate Y=0 fallback,
  `_refresh_mixture_arrays` matches per-cell calls, Hasegawa A
  integration arrays in physical bounds.

### Phase 3 — Thread arrays through solver

**Files touched**: `solver.py` (every PISO callsite),
`simulation.py` (hot loop), energy diagnostics. Largest phase by LOC.

**Deliverables**:
- All scalar (γ, R, Cp) arguments converted to arrays.
- Sensible-enthalpy advection where Cp varies cell-to-cell.
- Nozzle BC uses cell-N-1 mixture.
- T_ceiling per cell.
- Source-CFL cap per cell.
- Existing Hasegawa A baseline still runs and produces a comparable
  trace (within ±10% MSE without re-calibration; full re-cal is
  Phase 5).

### Phase 4 — Validation tests

**Files touched**: `tests/test_yns_*.py`.

**Deliverables**:
1. `test_yns_pure_pyrogen_limit`: head-end pyrogen with no grain ignition
   → Y → 100% species-0, thermo arrays match pyrogen single-gas.
2. `test_yns_pure_propellant_limit`: long burn after pyrogen exhausted
   → Y → 100% species-1, thermo arrays match propellant single-gas
   (exact agreement with v0.7.0 baseline modulo machine epsilon).
3. `test_yns_mass_conservation`: cumulative species mass in bore +
   cumulative species mass out the nozzle = cumulative species mass
   injected, to <1e-3 relative.
4. `test_yns_hasegawa_a_baseline_runs`: existing calibrated parameters
   still produce a finite trace; MSE within ±50% of v0.7.0 baseline
   (we expect drift; full re-cal is Phase 5).
5. `test_yns_y_invariants`: `sum(Y[i, :]) ≈ 1.0` and
   `0 ≤ Y[i, s] ≤ 1` at every recorded step.

### Phase 5 — Hasegawa A re-LHS (separate work item)

After Phases 1-4 merge, re-run the 7-var Hasegawa A LHS. Hypothesis:
calibrated `k_solid` and possibly `k_thermal` shift because γ/Cp
variation captures some of what those parameters were absorbing as a
global compromise between frozen and effective gas regimes.

## API breaking-change log

v0.7.1 introduces hard breaks (per `feedback_api_breaks` memory; v0.7.x
allows this):

- `Pyrogen` gains a default-derived `Cp` (no signature change for
  existing callers).
- `Propellant`/`PropellantTab` callers that read `gas.Cp_gas`,
  `gas.gamma`, `gas.R_specific` as scalars on the GasProperties object
  see no change — the registry consumes the same fields and exposes
  per-cell arrays internally.
- `run_simulation` signature unchanged externally.
- Internal kernel signatures change substantially. Custom callers of
  `_piso_step` (none expected; this is a private function) must adapt.

## What's out of scope for v0.7.1

- **Per-cell transport properties** (k, μ, Pr) — deferred to v0.7.2 if
  calibration shows benefit.
- **Multi-tab per-species Cp(p)** — deferred. PropellantTab's pressure
  banding still applies only to (a, n); thermo stays at the
  representative tab value.
- **Z-N dynamic burn rate** — that's a separate v0.7.1 feature, not
  blocked by this work. They can be developed in parallel and merged
  independently. Order TBD; see TASKS.md.
- **N-species nozzle exit thermo for thrust** — for v0.7.1 the nozzle
  uses the cell-N-1 mixture for outflow. The thrust-side CF/Isp
  calculation may need its own per-mixture treatment in v0.7.2.

## Links

- Motivating research: [post_v0_7_0/references/spinball_walkthrough.md](../post_v0_7_0/references/spinball_walkthrough.md)
- v0.7.0 baseline: [v0_7_0/DESIGN.md](../v0_7_0/DESIGN.md)
- Current architecture map: [ARCHITECTURE.md](../../ARCHITECTURE.md)
- Calibration state: `project_hasegawa_calibration_state` memory
- Decision provenance: `project_spinball_research_state` memory
