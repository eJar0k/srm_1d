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

## Phase 3 — Thread arrays through solver

- [ ] Convert `_piso_step` to take `gamma_arr, R_arr, Cp_arr` arrays; index per cell.
- [ ] Update `_signed_choked_or_subsonic_flow` to accept boundary γ_b, R_b (from upstream cell).
- [ ] EOS updates `rho = P / (R_arr[i] * T_new[i])` per cell.
- [ ] Energy advection: switch from advecting T to advecting (Cp_arr·T) — sensible enthalpy.
- [ ] T_ceiling per cell from max(T_flame[s] for s such that Y[i,s] > 0.05) * 1.01.
- [ ] Update `_pyrogen_surface_heat_power` to take cell-target Cp_arr[i_target].
- [ ] Update `_goodman_ignition_sources_and_mass` to use per-cell Cp_arr.
- [ ] Update energy diagnostics (`_gas_sensible_energy`, `_thermal_source_power`) to use per-cell Cp_arr.
- [ ] Update source-CFL cap (`_source_cfl_cap_from_thermal_source`) to use per-cell Cp_arr.
- [ ] Nozzle BC uses cell-N-1 mixture.
- [ ] `Cp_arr` references in audit histories.
- [ ] Clear `srm_1d/__pycache__/` after first solver refactor commit (CLAUDE.md gotcha #1).
- [ ] Smoke: existing Hasegawa A baseline runs end-to-end without NaN.

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

- **Phases 1 + 2 complete.** 193/193 tests pass. Solver behavior unchanged (mixture arrays computed but unused by PISO).
- **Phase 3 is next.** This is where simulation behavior actually changes; estimated 30-40% of a focused session. See DESIGN.md §5 for the full callsite list.
