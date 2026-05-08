# Architecture Map — srm_1d v0.6.0

Function-level map of every module. Use as a reference when making
changes to verify nothing is lost or silently modified.

## solver.py (pure numerics, no project dependencies)

- `thomas_solve(a, b, c, d, N)` — TDMA tridiagonal solver. O(N).
- `piso_step(rho, u, P, T, ..., N)` — One PISO time step on staggered
  grid. Momentum predictor → 2 pressure corrections → energy → EOS.
  Returns updated (rho, u, P, T).
- `compute_dt_cfl(u, a_sound, dx, N, cfl, dt_max)` — Adaptive CFL step.

## burn_rate.py (Ma et al. 2020, no project dependencies)

- `haaland_friction(Re, roughness, D)` — Darcy factor: laminar/transition/turbulent.
- `gnielinski_nusselt(Re, Pr, D, L, f, T_gas, T_surf, kappa)` — Nu with entrance + temp correction.
- `transpiration_correction(beta)` — h/h0 = β/(exp(β)-1) blowing correction.
- `select_tab_idx(P, tab_min_p, tab_max_p, n_tabs)` — Hard-switchover
  tab lookup (strict containment first, closest-boundary fallback).
- `saint_robert_from_tabs(P, tab_min_p, tab_max_p, tab_a, tab_n, n_tabs)` —
  r₀ = a(P)·P^n(P) using tab lookup.
- `burn_rate_cell(P, Re, D_hyd, x, ..., tab_min_p, tab_max_p, tab_a, tab_n, n_tabs, kappa)` —
  Single-cell total burn rate via bisection; uses tab lookup for r₀.
- `compute_burn_rates(P, Re, D_hyd, x, is_burning, ..., tab_arrays, kappa, N)` —
  Vectorized wrapper.

## propellant.py (leaf module, no project dependencies)

- `class PropellantTab` — dataclass mirroring openMotor's `PropellantTab`
  schema: (min_pressure, max_pressure, a, n, gamma, T_flame, molecular_weight).
- `class Propellant` — dataclass: name, tabs (list of PropellantTab),
  rho_propellant, Cps, T_surface, T_initial, mu_gas, k_gas, Cp_gas.
  Methods: `select_tab(P)` (hard-switchover, closest-boundary fallback,
  matches openMotor's `getCombustionProperties`), `representative_tab(P)`,
  `burn_rate_normal(P)` (scalar, uses select_tab),
  `tab_arrays()` → (min_p, max_p, a, n) numpy arrays for Numba.
- `class GasProperties` — dataclass: gamma, MW, R_specific, mu, k, Pr, Cp.
- `create_gas_properties(gamma, MW, T, mu, k, Cp)` — From CEA/RPA.
- `create_gas_properties_estimated(gamma, MW, T)` — Sutherland/Eucken fallback.
- `speed_of_sound(gamma, R, T)` — a = √(γRT).
- `density_from_ideal_gas(P, R, T)` — ρ = P/(RT).
- `critical_flow_function(gamma)` — Γ for choked nozzle flow.
- `characteristic_velocity(gamma, R, T)` — c* = √(RT)/Γ.

(Named propellants live as data in `srm_1d/motors/<motor>.ric` plus
sibling `<motor>.transport.yaml`. The v0.5.x `make_hasegawa_propellant_1`
/ `make_king_propellant_4525` factories were deleted in v0.6.0.)

## grain_geometry.py (no project dependencies)

### Config classes (setup only)
- `class GrainSegment` — D_bore_fwd, D_bore_aft, length, inhibition,
  optional `fmm_table` (FmmTable from srm_1d.fmm_grain). When
  fmm_table is set, the segment uses FMM table lookup; otherwise
  it's analytic cylindrical/conical.
- `class MotorGeometry` — segments list, D_outer, N_cells.
  Throat lives on the separate `Nozzle` object — see nozzle.py.
  - `compile_geometry_arrays()` — exports to Numba-compatible arrays.
    Produces both analytic (cell_D_bore_init) and FMM-packed
    (fmm_offset, fmm_reg_flat, fmm_perim_flat, fmm_port_flat,
    cell_fmm_idx) data, plus the per-cell switch (cell_segment_type)
    and per-cell regress[] state. Auto-inhibits touching faces.
  - `total_propellant_volume()` — exact for conical (integrated annular
    formula); for FMM segments uses (casting_area − port_area_init) × L.

### Numba-compiled per-step functions
- `update_cell_geometry(regress, D_port, ..., cell_D_bore_init,
  cell_wall_web, cell_segment_type, cell_fmm_idx, fmm_offset,
  fmm_reg_flat, fmm_perim_flat, fmm_port_flat)` — Full geometry
  recomputation. Branches per cell on cell_segment_type:
  type 0 (analytic): D_port = D_bore_init + 2·regress, A_port = π/4·D²,
  base_perimeter = π·D, D_hyd = D_port.
  type 1 (FMM): A_port and base_perimeter from CSR-packed FMM tables;
  D_hyd = 4·A_port/perimeter (correct for non-circular ports).
  Then axial overlap fraction (grain_frac), radial burnout ramp
  (f_active in regression-depth space), C_burn = base_perimeter·grain_frac·f_active.
  v0.6.0+ also: a **volumetric overlap accumulator** (cells straddling
  two segments sum C_burn from both) and a **partition-of-unity
  end-face kernel** (each face's mass split over 2 adjacent cells with
  weights summing to 1.0).
- `advance_bore_regression(regress, r_total, dt, N, cell_wall_web,
  cell_segment_id)` — regress[i] += r_total[i]·f_active·dt. Same
  burnout ramp logic for analytic and FMM cells; primary state is
  always regression depth.
- `advance_endface_regression(...)` — Saint-Robert rate only (no erosion).
- `_saint_robert_local`, `_fmm_lookup_flat` — private helpers.

### `build_snapped_geometry(segments_spec, D_outer, target_propellant_cells=100)`
The canonical builder. Computes `dx = L_propellant / target_cells`,
applies a Nyquist-CFD clamp so the smallest gap gets ≥1 cell, then
integer-snaps every segment length and inter-segment gap to multiples
of `dx`. Returns a `MotorGeometry` whose cell boundaries align with
segment edges by construction. `segments_spec` keys: `D_bore_fwd`
(required), `length` (required), `D_bore_aft`, `gap_after`,
`inhibit_fwd`/`inhibit_aft`, `fmm_table` (all optional).

(The v0.5.x parametric factories — `make_bates_motor`,
`make_single_cylinder`, `make_conical_grain`, `make_stepped_motor`,
`make_example_bates`, `make_hasegawa_motor_A/B/C_geo/_nozzle` — were
deleted in v0.6.0. Build geometry directly via `build_snapped_geometry`
or load named motors from `srm_1d/motors/*.ric`.)

## nozzle.py (imports propellant.critical_flow_function; openMotor-aligned)

- `class Nozzle` — D_throat, D_exit, efficiency, div_angle, conv_angle,
  throat_length, erosion_coeff, slag_coeff. Field names mirror
  openMotor's `motorlib.nozzle.Nozzle` (snake_case'd; units kept
  human-readable internally).
  Methods: `throat_area(d_throat=0)`, `exit_area()`, `expansion_ratio`,
  `divergence_losses()`, `throat_losses(d_throat)` (RasAero aspect-ratio),
  `skin_losses()` (constant 0.99), `exit_pressure(gamma, P_c)`,
  `ideal_thrust_coeff(...)`, `adjusted_thrust_coeff(...)`.
- `exit_pressure_from_expansion_ratio(gamma, eps)` — Newton iteration.
- `ideal_thrust_coefficient(gamma, eps, Pc, Pa)` — Sutton eq. 3-30.
- `compute_thrust_isp(Pc, gamma, eps, div_loss, efficiency, throat_loss,
  skin_loss, ..., c_star)` — Single instant; applies openMotor's
  adjusted-CF formula:
      CF_adj = divLoss × throatLoss × efficiency
             × (skinLoss × CF_ideal + (1 − skinLoss))
- `compute_thrust_history(t, P, ..., erosion, slag, div_loss, efficiency,
  throat_length, skin_loss, ...)` — Integrates throat diameter over the
  pressure history. throat_loss is recomputed each step from the current
  D_throat. Returns thrust, CF, Isp, Pe, D_throat arrays.
- `compute_motor_performance(result, nozzle, propellant, P_ambient=...)` —
  Post-processing. Uses simulation's D_throat history when available
  (in-loop coupling).
- `print_performance_summary(perf, nozzle)` — Formatted output with
  divergence/throat/skin loss factors and throat change.

## simulation.py (imports solver, burn_rate, propellant, grain_geometry)

### Fused compiled helpers
- `_post_piso_update(rho, u, P, T, ..., gamma_R, roughness)` — Fused:
  velocity interpolation + Re + Mach + friction + T_max. One pass over N cells.
- `_ignition_source_and_mass(P, ..., rho_propellant, N)` — Fused:
  ignition check + ramp + source assembly + mass sum. Returns n_burning,
  n_ignited, mass_sum.

### Compiled time loop
- `_run_time_loop` — Single @njit function containing the entire
  while loop. Steps: geometry → burn rates → ignition+source → throat
  evolution → PISO → post-PISO → bookkeeping+snapshots.
  Returns: n_steps, n_snaps, mass_produced, mass_nozzle, burnthrough_time,
  D_throat_final, termination_code.

### Public API
- `run_simulation(geo, propellant, nozzle, P_ambient=101325, roughness,
  kappa, ...)` — Extracts D_throat / erosion_coeff / slag_coeff from
  the Nozzle, calls _run_time_loop, wraps results into dict with time
  histories, snapshots, per-grain data, and structured summary.

### Igniter model (v0.6.0 placeholder)
Single-knob exponential decay distributed full-grain:
`mdot_igniter(t) = (igniter_mass/igniter_tau) · exp(-t/igniter_tau)`,
injected uniformly across all N cells (`n_ign_cells = N`). No pressure
feedback. The v0.5.x Saint-Robert pyrogen model and `v_flame_*`
flame-front tracking were deleted. This model is documented as a
placeholder pending a hot-gas plenum replacement (v0.7.0+); see
DEVNOTES "Igniter".

## plotting.py (imports matplotlib)

- `plot_pressure(result, experimental, time_offset=0.0, ...)` —
  Head-end pressure with experimental overlay. `time_offset` (seconds)
  is applied uniformly to all experimental datasets to align ignition
  events; pre-shift the time arrays for per-dataset offsets.
- `plot_thrust(result, perf, save_path)` — Thrust + Isp two-panel.
- `plot_flow_snapshot(result, t_target)` — 2×2: P, Mach, burn rate
  (with endface orange bars), port diameter.
- `plot_summary(result, perf, experimental, time_offset=0.0)` — 2×2 combined.
- `plot_comparison(result, perf, reference)` — Overlay with openMotor CSV.
- `plot_grain_regression(grain_metrics, geo)` — Per-grain regression/web.
- `load_experimental_csv(filepath, ...)` — CSV with unit conversion.
- `HASEGAWA_MOTOR_A_EXPERIMENTAL` — Embedded digitized data (36 points).
  v0.6.0 removed the in-dict `time_offset` key; pass it as a kwarg
  to `plot_pressure`/`plot_summary`.

## openmotor_adapter.py (imports propellant, grain_geometry, nozzle, simulation)

- `load_ric(filepath)` — YAML parser with Python-tag handling.
- `load_transport(transport_path)` — sibling YAML loader returning a
  `gas_props` dict (`mu`, `k`, `Cp`).
- `convert_propellant(ric_prop, gas_props)` — MW g/mol→kg/mol.
- `convert_geometry(ric_grains, target_propellant_cells=100,
  fmm_map_dim=1001)` — Routes through `build_snapped_geometry`.
  Auto-applies inter-segment gap of `max(3mm, 5%·D_outer)` via
  per-segment `gap_after`. BATES + Conical analytic, FMM via
  `from_ric_grain`. v0.6.0 removed `N_cells` and `spacing` kwargs.
- `convert_nozzle(ric_nozzle)` — Returns a Nozzle. Maps
  `throat/exit/efficiency/divAngle/convAngle/throatLength/erosionCoeff/slagCoeff`
  to our snake_case fields. erosionCoeff converted m/(s·Pa)→μm/(s·MPa).
- `ric_to_sim_args(motor, gas_props=None, target_propellant_cells=100,
  **sim_overrides)` — Returns a kwargs dict including `geo`,
  `propellant`, `nozzle`, `P_ambient`, `P_cutoff`.
- `run_from_ric(filepath, gas_props=None, transport_path=None,
  **sim_overrides)` — Full pipeline. If `gas_props` is None,
  auto-resolves a sibling `<stem>.transport.yaml`. Returns
  `result, perf, nozzle, geo, prop`.
- `compute_grain_metrics(result, geo, propellant)` — Per-grain regression,
  web, mass remaining at each snapshot time.
- `result_to_csv(result, perf, geo, propellant)` — openMotor-compatible CSV.
- `save_csv(filepath, ...)` — Write to file.
- `load_openmotor_csv(filepath)` — Read openMotor CSV for comparison.
- `print_ric_summary(filepath)` — Human-readable .ric summary.

## motors/ (data-driven motor library, v0.6.0+)

Each named motor lives as a `.ric` file (openMotor schema) plus an
optional sibling `.transport.yaml` (srm_1d-specific extension supplying
`mu`, `k`, `Cp` — combustion gas transport). Loaded via
`run_from_ric('srm_1d/motors/<motor>.ric')` which auto-discovers the
sibling YAML.

- `hasegawa_a.ric`, `hasegawa_a.transport.yaml` — Hasegawa Motor A
  (L/D=42 single-segment BATES; canonical validation target).
- `hasegawa_b.ric`, `hasegawa_b.transport.yaml` — Motor B (half-length).
- `hasegawa_c.ric`, `hasegawa_c.transport.yaml` — Motor C (wider bore).
- `example_bates.ric`, `example_bates.transport.yaml` — 4×120mm BATES.

## tools/sensitivity.py (Latin Hypercube parameter sweeps)

- `run_lhs(motor_path, bounds, n_samples, fitness_fn, n_workers,
  seed, csv_path, **sim_kwargs)` — Parallel LHS driver. Uses
  `scipy.stats.qmc.LatinHypercube` + `concurrent.futures.ProcessPoolExecutor`.
- `mse_fitness(t_exp, p_exp, t_min)` — Factory: full-trace MSE
  fitness function (default for trace-fitting).
- `impulse_error_fitness(I_target)`,
  `peak_pressure_error_fitness(P_target)` — Alternate fitness factories.
- Worker (`_run_one`) is module-level so it pickles cleanly into
  ProcessPoolExecutor children. All `@njit` functions in the solver
  use `cache=True` so workers reuse compiled artifacts.

## fmm_grain.py (openMotor bridge for FMM grain types)

- `class FmmTable` — dataclass holding sampled (reg_depth, perimeter,
  port_area) arrays plus wall_web, grain_outer_diameter, grain_length,
  inhibited_fwd/aft, geom_name. Built by `from_openmotor`/`from_ric_grain`.
- `_setup_openmotor_path()` — lazy adds local openMotor checkout to
  sys.path. Injects a Numba-JIT `_get_perimeter` shim into
  `mathlib._find_perimeter_cy` so the Cython build isn't needed.
- `_marching_squares_perimeter(arr, level)` — @njit perimeter
  computation, verbatim port of openMotor's Cython algorithm.
- `from_openmotor(om_grain, map_dim=1001)` — runs initGeometry +
  generateCoreMap + generateRegressionMap on the supplied openMotor
  FmmGrain instance and samples the regression map.
- `from_ric_grain(ric_grain, map_dim=1001)` — convenience: dispatches
  on `ric_grain['type']` to instantiate the right openMotor class
  (Finocyl, StarGrain, MoonBurner, CGrain, DGrain, XCore, CustomGrain).
- `fmm_table_lookup(regress, reg_arr, val_arr, n_samples)` — @njit
  O(1) linear interp on a uniform-grid FMM table (1D version).

## Dependency Graph

```
propellant.py        ← leaf (no project deps)
burn_rate.py         ← leaf (no project deps; propellant params passed as scalars)
solver.py            ← leaf (no project deps)
grain_geometry.py    ← leaf (no project deps)
nozzle.py            ← imports propellant (critical_flow_function, R_UNIVERSAL)
fmm_grain.py         ← imports motorlib + mathlib (lazy, optional)
simulation.py        ← imports solver, burn_rate, propellant, grain_geometry
plotting.py          ← imports matplotlib only (result dicts are plain data)
openmotor_adapter.py ← imports propellant, grain_geometry, nozzle,
                       simulation; lazily imports fmm_grain for FMM grain types
```
