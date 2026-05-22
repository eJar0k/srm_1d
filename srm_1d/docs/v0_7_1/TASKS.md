# srm_1d v0.7.1 — Phase Checklist

See [DESIGN.md](DESIGN.md) for the full architecture.

## Phase 1 — Species registry, Y[N,S] state, transport [COMPLETE 2026-05-23]

- [x] Add `GasSpecies` dataclass to `propellant.py` with (name, γ, Cp, M, T_flame).
- [x] Add `Pyrogen.species` property returning a `GasSpecies` (derive Cp from γ, M if not given).
- [x] Add `Propellant.species` method returning a `GasSpecies` from rep tab.
- [x] Add `ambient_air_species(T)` helper for the pre-fill species (s=2).
- [x] Add `species_array(list)` packer returning `np.ndarray[S, 4]`.
- [x] Build species registry in `run_simulation` (S=3: igniter, propellant, ambient).
- [x] Allocate `Y_species: np.ndarray[N, S]`; init pre-fill to species 2 (ambient).
- [x] Allocate `mass_source_by_species[N, S]` and `rho_pre_step[N]` working arrays.
- [x] Split mass-source writes per species in `_goodman_ignition_sources_and_mass` (grain → s=1) and the pyrogen injection block (pyrogen → s=0).
- [x] Thread `Y_species, species_params_arr, mass_source_by_species` into the `_run_time_loop` signature and call site.
- [x] Add `_advect_species` Numba kernel: mass-fraction-conservative upwind with face density matching PISO, nozzle drain, per-cell renormalize.
- [x] Call `_advect_species` after PISO nozzle bookkeeping using `rho_pre_step` + `rho` + `u`.
- [x] Expose `Y_species_final, species_params, species_names, rho_final, A_port_final` in result dict.
- [x] 5 new tests in `tests/test_yns_transport.py`:
  - [x] `test_advect_species_pure_advection_pyrogen_pulse_moves_right`
  - [x] `test_advect_species_source_only_no_flow`
  - [x] `test_advect_species_renormalization_repairs_fp_drift`
  - [x] `test_advect_species_nozzle_outflow_drains_last_cell`
  - [x] `test_yns_hasegawa_a_short_run_invariants` (50 ms integration; verifies Y invariants + species fractions sensible)

## Phase 2 — Per-cell mixture derivation [COMPLETE 2026-05-23]

- [x] Add `_compute_mixture_cell(Y_row, species_params) → (γ, Cp, R, M)` helper, Numba-decorated.
- [x] Add `_refresh_mixture_arrays(Y, species_params, gamma_arr, Cp_arr, R_arr, M_arr, N)` loop helper.
- [x] Allocate `gamma_mix_arr, Cp_mix_arr, R_mix_arr, M_mix_arr` of length N.
- [x] Recompute per cell after `_advect_species` each step.
- [x] Expose `gamma_mix_final / Cp_mix_final / R_mix_final / M_mix_final` in result dict.
- [x] 8 new tests in `tests/test_yns_mixture.py`:
  - [x] Pure-species limits for all 3 species
  - [x] 50/50 binary mix matches hand calc
  - [x] gamma = Cp/(Cp-R) consistency
  - [x] Degenerate Y=0 fallback
  - [x] `_refresh_mixture_arrays` matches per-cell calls
  - [x] Hasegawa A integration: arrays have sensible physical bounds

## Phase 3 — Thread arrays through solver [COMPLETE 2026-05-23]

Delivered in two commits on `v0.7.0-phase4`:
1. `thermal_source` units shifted from kg·K/(s·m) to W/m (behavior
   preserved by multiplying T_flame contributions by scalar Cp_gas
   at source sites and dividing back in PISO).
2. PISO + post-PISO + nozzle BC + source-CFL cap converted to take
   per-cell `gamma_arr / R_arr / Cp_arr / T_ceiling_arr`. Energy
   advection switched to sensible-enthalpy (Cp·T) form.

- [x] Convert `_piso_step` to take `gamma_arr, R_arr, Cp_arr` arrays; index per cell.
- [x] Nozzle BC uses cell-N-1 mixture (γ, R) for both pressure
      corrections, the velocity update, and the bookkeeping
      `_nozzle_boundary_flow` call.
- [x] EOS updates `rho = P / (R_arr[i] * T_new[i])` per cell.
- [x] Energy advection: switched from advecting T to advecting
      (Cp_upwind · T_upwind) — sensible enthalpy. Cell update solves
      for T_new = h_new / Cp_arr[i].
- [x] T_ceiling per cell from species_params — **relaxed** from
      DESIGN §5: uses `max(T_flame[s] for s in all species) * 1.01`,
      not the strict `Y > 0.05` filter, to avoid clipping the v0.7.0
      initial condition (T = T_flame_prop, Y = 100% ambient) to
      ~300 K on step 0. The relaxed ceiling equals v0.7.0's scalar
      `T_flame * 1.01` when propellant is the hottest species, so the
      baseline trace is preserved. See `_compute_T_ceiling_arr` doc.
- [x] `_goodman_ignition_sources_and_mass` writes thermal_source in
      W/m using scalar Cp_gas (Phase 3 step 1). True per-species
      Cp lookups are slotted for Phase 3.5 — see "Phase 3.5" section
      below.
- [x] `_pyrogen_surface_thermal_sink` returns W/m sink (Phase 3 step 1).
      `_pyrogen_surface_heat_power` retains its `Cp_gas` argument —
      this is the pyrogen species' Cp at the surface heat-transfer
      step, not the cell mixture Cp.
- [x] Energy diagnostics (`_gas_sensible_energy`, `_thermal_source_power`)
      updated for per-cell Cp_arr / W/m thermal_source.
- [x] Source-CFL cap (`compute_dt_source_cap`) takes per-cell Cp_arr.
- [x] `_post_piso_update` takes (gamma_arr, R_arr); a_max is now the
      max local sound speed across cells.
- [x] Clear `srm_1d/__pycache__/` + `.nbi` + `.nbc` after each kernel
      refactor commit.
- [x] Smoke: Hasegawa A baseline runs 1.4M steps end-to-end without
      NaN. P_peak 6.26 MPa, mass-balance err 0.1%, c* 1543 m/s,
      O5347 designation — within Phase 3's "±10% of v0.7.0 baseline"
      target.

### Phase 3.5 — True per-species Cp lookups [COMPLETE 2026-05-23]

Delivered in commit `6f0789e`:

- [x] `_goodman_ignition_sources_and_mass` multiplies grain T_flame
      contributions by `Cp_propellant` (= `species_params[_SPECIES_PROPELLANT, 1]`).
- [x] Pyrogen injection in `_run_time_loop` multiplies by
      `Cp_pyrogen` (= `species_params[_SPECIES_IGNITER, 1]`).
- [x] `_pyrogen_surface_heat_power` arg renamed `Cp_gas` → `Cp_pyrogen`;
      caller passes pyrogen Cp explicitly.

**Trace impact (Hasegawa A smoke)**:

| Metric        | Phase 3      | Phase 3.5    |
|---------------|--------------|--------------|
| P_peak        | 6.26 MPa     | 6.19 MPa     |
| t at P_peak   | 0.041 s      | 3.361 s      |
| pyrogen dur.  | 152 ms       | 576 ms       |
| c*            | 1543 m/s     | 1543 m/s     |
| mass balance  | 0.1 %        | 0.1 %        |

The pre-Phase-3.5 trace had a sharp ignition spike at ~40 ms followed
by a lower steady-state plateau. Phase 3.5 suppresses the spike
because the pyrogen-to-surface sensible-power cap (`mdot · Cp_pyrogen
· (T_ig - T_surf)`) is now ~33% lower than under the scalar `Cp_gas =
Cp_propellant` placeholder (BPNV Cp ≈ 1385 vs propellant 2060). This
delays grain ignition; pressure climbs more gradually to a mid-burn
peak. Phase 4 gates (P_peak ±50%, c*, mass balance) still pass.
Phase 5 LHS will calibrate against this physics-faithful build.

### Strict T_ceiling formula [COMPLETE 2026-05-23]

Delivered in commit `78209fb`. The relaxed Phase 3 ceiling
(max-over-all-species, global constant) is replaced by:

```
T_ceiling[i] = max(T_flame[s] for s with Y[i, s] > 0.05) * 1.01
T_ceiling[i] = max(T_ceiling[i], T_initial_gas * 1.01)
```

The IC guard preserves the v0.7.1 initial condition (T =
T_flame_propellant while Y = 100% ambient). Without it the strict
filter would clip the bore gas to T_ambient · 1.01 on step 0.

Tightens overshoot detection in:
- Pyrogen-only cells: cap at T_flame_pyrogen · 1.01 (~2828 K for
  BPNV) once pyrogen displaces ambient. Prior relaxed form allowed
  ~9% overshoot to T_flame_propellant · 1.01.
- v0.8.0 multi-grain configurations: pure-sustainer cells cap at
  sustainer T_flame, not the hottest species across the registry.

7 new direct-kernel tests in `tests/test_yns_mixture.py` cover
pyrogen / propellant / mixed cells, the Y_min filter, IC guard
activation and deactivation, and array length.

## Phase 4 — Validation tests [COMPLETE 2026-05-23]

Delivered as `tests/test_yns_phase4_validation.py` (6 tests, +6 to the
suite for 199/199 overall).

- [x] `test_yns_pure_pyrogen_limit_thermo_matches_pyrogen_species`:
      T_ignition = 20000 K suppresses grain ignition; after 80 ms,
      cells with Y[:, pyrogen] > 0.95 have (Cp_mix, R_mix) matching
      the pyrogen species exactly and γ_mix matching the
      ideal-gas-**derived** γ (= Cp/(Cp-R)).
- [x] `test_yns_pure_propellant_limit_thermo_matches_propellant_species`:
      1-second run; cells with Y[:, prop] > 0.99 have mixture thermo
      collapsing to the propellant species in the same way.
- [x] `test_yns_overall_mass_balance_closes`: 1-second run;
      `mass_balance_error` < 2% (v0.7.0 was ~0.1%, Phase 3 leaves
      headroom).
- [x] `test_yns_ambient_species_purges_through_nozzle`: ambient bore
      mass < 1% of total bore mass after 1 s — the per-species
      conservation guarantee the S=3 registry exists to enable.
- [x] `test_yns_hasegawa_a_baseline_within_phase3_tolerance`: P_peak
      and c* within ±50% of v0.7.0 baseline. Actual: P_peak +1% drift,
      c* unchanged. Phase 5 LHS will tighten this gate.
- [x] `test_yns_y_invariants_over_full_3s_hasegawa_run`: sum(Y) = 1 ±
      1e-6 and 0 ≤ Y ≤ 1 over a 3-second run (~420 k advection steps).

### Phase 4 findings (status)

- **Species γ is derived, not declared** — still open. The mixture rule
  computes γ_mix = Cp_mix / (Cp_mix - R_mix) and ignores
  `species_params[s, 0]`. Hasegawa A YAML declares γ=1.19 vs derived
  1.189 (~0.04% sound-speed drift). Open question: validate at YAML
  load (refuse inconsistent inputs) or document the derivation policy
  and drop the species γ column. No solver code change needed.
- **Strict T_ceiling formula** [LANDED 2026-05-23, commit `78209fb`]:
  per-cell, Y > 0.05 filter, with an IC guard at T_initial_gas · 1.01.
  See strict-form section above.

## Phase 5 — Hasegawa A re-LHS (post-merge)

- [ ] Re-run [hasegawa_a_lhs.py](../../examples/hasegawa_a_lhs.py) with the new model.
- [ ] Compare rank-1 mse_all to v0.7.0 baseline (0.0968 MPa²).
- [ ] Check whether k_solid and k_thermal calibrations shift meaningfully (hypothesis: they relax once γ/Cp variation absorbs some of their compromise role).
- [ ] Update `project_hasegawa_calibration_state` memory.

## Tag criteria

`v0.7.1` ships when:

- All Phase 4 tests green.
- Phase 5 rank-1 mse_all ≤ v0.7.0 baseline OR ≤ +10% (with documented physical justification for any regression).
- Existing pytest baseline (199 with current N-species infrastructure + Phase 4 validation suite) remains green.
- DEVNOTES API breaking-change log updated.
- CLAUDE.md "Critical gotchas" updated if new ones surface.

## Current state (2026-05-23)

- **Phases 1 + 2 + 3 + 3.5 + 4 complete + strict T_ceiling.** 206/206
  tests pass. PISO consumes per-cell (γ, R, Cp, T_ceiling) with
  sensible-enthalpy advection; cell-N-1 nozzle BC; each species
  injects its own Cp; strict per-cell ceiling with IC guard.
- **Phase 5 (Hasegawa A re-LHS) is the only remaining work item.**
  The species γ inconsistency (Phase 4 finding) is YAML-side cleanup
  that doesn't affect the solver; it can land alongside Phase 5 or
  independently.
