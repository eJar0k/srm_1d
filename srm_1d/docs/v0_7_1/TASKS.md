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

### Phase 3.5 — True per-species Cp lookups (not yet wired)

Phase 3 deliberately kept the source-side Cp multiplications using
scalar `Cp_gas` (the representative-tab gas Cp) instead of per-species
Cp from `species_params`. This is a CONSCIOUS deviation: at unit-shift
time the PISO step still expects enthalpy balanced against `Cp_gas`,
and threading per-species Cp into the source sites without simultaneous
per-species accounting in the energy equation would introduce a
mass-rate–Cp asymmetry that drifts the v0.7.0 calibration. With Phase 3
now consuming per-cell `Cp_arr` in PISO + energy advection, the next
incremental refinement is:

- [ ] `_goodman_ignition_sources_and_mass` multiplies grain T_flame
      contributions by `species_params[_SPECIES_PROPELLANT, 1]`
      (propellant Cp) instead of scalar Cp_gas.
- [ ] Pyrogen injection in `_run_time_loop` multiplies by
      `species_params[_SPECIES_IGNITER, 1]` (pyrogen Cp) instead of
      Cp_gas.
- [ ] `_pyrogen_surface_heat_power` caller passes the pyrogen
      species' Cp (it already takes a Cp arg; this is just a call-
      site change).

Order: Phase 3.5 lands before Phase 5 (Hasegawa A re-LHS) but after
Phase 4 validation tests. The motivation is faithfulness, not
calibration fit — the v0.7.0 fit already absorbs the Cp_gas
compromise, so Phase 5 LHS is the right place to absorb the
per-species refinement instead of trying to chase it pre-cal.

## Phase 4 — Validation tests

- [ ] `test_yns_pure_pyrogen_limit`: head-end pyrogen, no grain ignition → thermo matches pyrogen single-gas.
- [ ] `test_yns_pure_propellant_limit`: long burn after pyrogen depleted → thermo matches v0.7.0 single-gas baseline to machine epsilon.
- [ ] `test_yns_mass_conservation`: cumulative species accounting closes to <1e-3 relative.
- [ ] `test_yns_hasegawa_a_baseline_runs`: peak P, t_burn, mse_all within ±50% of v0.7.0 baseline; trace shape qualitatively similar.
- [ ] `test_yns_y_invariants_long_run`: sum≈1, 0≤Y≤1 over a full Hasegawa A run.

## Phase 5 — Hasegawa A re-LHS (post-merge)

- [ ] Re-run [hasegawa_a_lhs.py](../../examples/hasegawa_a_lhs.py) with the new model.
- [ ] Compare rank-1 mse_all to v0.7.0 baseline (0.0968 MPa²).
- [ ] Check whether k_solid and k_thermal calibrations shift meaningfully (hypothesis: they relax once γ/Cp variation absorbs some of their compromise role).
- [ ] Update `project_hasegawa_calibration_state` memory.

## Tag criteria

`v0.7.1` ships when:

- All Phase 4 tests green.
- Phase 5 rank-1 mse_all ≤ v0.7.0 baseline OR ≤ +10% (with documented physical justification for any regression).
- Existing pytest baseline (193 with current N-species infrastructure) remains green.
- DEVNOTES API breaking-change log updated.
- CLAUDE.md "Critical gotchas" updated if new ones surface.

## Current state (2026-05-23)

- **Phases 1 + 2 + 3 complete.** 193/193 tests pass. Solver behavior
  changes (sensible-enthalpy advection, per-cell EOS, cell-N-1 nozzle
  BC) but the Hasegawa A baseline smoke is well within Phase 3's
  ±10% target.
- **Phase 3.5 is next** (per-species Cp lookups at source sites),
  followed by **Phase 4** validation tests and **Phase 5** Hasegawa A
  re-LHS. See `Phase 3.5` section above + the `Phase 4` checklist below.
