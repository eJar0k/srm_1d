# Development Notes — srm_1d

Lessons learned, gotchas, and non-obvious decisions. Read this before
making changes to avoid repeating costly mistakes.

## Critical Gotchas

### Mass conservation requires matched factors
The burnout ramp `f_active` MUST multiply both `C_burn` (in
update_cell_geometry) AND `dD/dt` (in advance_bore_regression).
Applying it to only one causes 3-40% mass errors:
- C_burn only → 3.3% deficit (old bug)
- dD/dt only → 40% surplus (attempted fix that failed)
- Both → 0.01% error (correct)

### Numba cache must be cleared after edits
Delete `srm_1d/__pycache__/` after any change to `@njit` functions.
Numba's `.nbi`/`.nbc` cache files persist even when `.pyc` files are
removed. Stale compiled code is the #1 cause of "the fix didn't work."

### End-face injection uses interval containment
`x_lo <= x_face <= x_hi` (geometric containment), NOT distance-based.
We tried three approaches:
- `dist < dx` → double-counting at segment boundaries
- `dist < dx/2` → missed faces when face drifts between cell centers
- `dist < dx/2 + 1e-6` → missed faces after regression
- `x_lo <= x_face <= x_hi` → correct, geometrically exact

### compile_geometry_arrays must match update_cell_geometry
Both must use the same overlap matching for cell-segment assignment.
If compile uses point-in-range and update uses overlap, boundary cells
get wrong initial D_port, causing silent conservation errors.

### seg_D_bore_init is dead — use cell_D_bore_init
All per-segment bore diameter arrays were replaced with per-cell
arrays in the conical grain refactor. The compiled loop, geometry
functions, and adapter all use `cell_D_bore_init` (length N).

## Calibration State

### Transport properties
The model is sensitive to frozen vs effective gas transport:
- Frozen (RPA): k=0.3685, Cp=2060 → underpredicts erosion, late tail-off
- Effective: k~0.55, Cp~1900 → better tail-off but higher plateau
- Recommendation: use effective values; frozen are physically wrong

### Roughness
Higher roughness preferentially boosts nozzle-end erosion (high Re).
- 20 μm: conservative, validated for Motor A with frozen transport
- 35-50 μm: physically reasonable for AP/HTPB/Al (aluminum agglomeration)
- >50 μm: overpredicts; ignition spike becomes too large

### Igniter
- igniter_A_burn is the primary calibration knob for spike height
- Auto-default (~11 cm² for 10g charge) works for Motor A
- The 2/3-power taper controls decay rate — approximate, not tunable
- Spike height matches well; decay back to plateau is slightly slow

## Performance Profile (BATES 120 cells, 2.3M steps)

```
Before optimization:    140s (20k steps/s)
Pre-allocated arrays:   113s (20k steps/s)  — negligible
Compiled reductions:     92s (25k steps/s)  — numpy dispatch eliminated
Fused per-step calls:    88s (26k steps/s)  — fewer Numba transitions
Compiled time loop:      33s (60k steps/s)  — Python loop eliminated
```

The remaining 33s is ~20s in piso_step (irreducible numerics) + ~13s
in geometry/burn rate/ignition (called every step or every N steps).

## .ric File Format Notes

- YAML with Python-specific tags (!!python/object/apply, !!python/tuple)
- SafeLoader rejects these; use custom loader with multi_constructor
- erosionCoeff stored in m/(s·Pa), multiply by 1e12 for μm/(s·MPa)
- propellant.m is g/mol, divide by 1000 for kg/mol
- propellant.a is m/s per Pa^n (same as ours, no conversion needed)
- propellant.k is gamma (not thermal conductivity)
- No inter-segment spacing in file; openMotor is 0-D and doesn't model gaps
- No gas transport properties in file; must be supplied separately
- No Cps or T_surface; we default to 1500 J/(kg·K) and 1000 K

## Grid Resolution Guidelines

- Minimum ~25 cells per grain segment for adequate spatial resolution
- Minimum 2-3 cells per inter-segment gap to avoid boundary artifacts
- The gap anomaly at 120 cells (dx ≈ gap width) resolves at 240 cells
- Conical grains need sufficient cells to resolve the bore gradient
- Micro-steps in pressure trace converge to zero with increasing N

## API Breaking Changes Log

- v0.2.0: GrainSegment.D_bore_initial → D_bore_fwd + D_bore_aft
- v0.2.0: seg_D_bore_init (N_seg) → cell_D_bore_init (N)
- v0.2.0: igniter_duration removed → igniter_a, igniter_n, igniter_rho, igniter_A_burn
- v0.2.0: run_from_ric returns 5 values (added geo, prop)
- v0.2.0: compute_burnout_ramp deleted (logic moved into update_cell_geometry)
- v0.3.0: openMotor-aligned Nozzle refactor:
    - `Nozzle.divergence_half_angle` → `div_angle`;
      `Nozzle.discharge_coefficient` → `efficiency`
    - Added `Nozzle.conv_angle`, `Nozzle.throat_length` fields
    - Removed `Nozzle.P_ambient` (now a `run_simulation` kwarg, matching
      openMotor's `motor.config.ambPressure`)
    - Removed `Nozzle.divergence_loss_factor` property →
      `Nozzle.divergence_losses()` instance method
    - Added `Nozzle.throat_losses()`, `Nozzle.skin_losses()`,
      `Nozzle.ideal_thrust_coeff()`, `Nozzle.adjusted_thrust_coeff()`
- v0.3.0: `MotorGeometry.D_throat` removed — throat lives on `Nozzle`.
    - All `make_*_motor`/`make_*_grain` factories drop the `D_throat` arg.
    - Hasegawa A/B/C: split into `make_hasegawa_motor_X_geo()` (geometry)
      + `make_hasegawa_motor_X_nozzle()` (nozzle).
- v0.3.0: `run_simulation(geo, propellant, nozzle, P_ambient=..., ...)` —
    requires `nozzle: Nozzle`. `erosion_coeff` and `slag_coeff` scalars
    removed (now read from the Nozzle).
- v0.3.0: Adjusted CF formula matches openMotor exactly:
      CF_adj = divLoss × throatLoss × efficiency
             × (skinLoss × CF_ideal + (1 − skinLoss))
    Default-case thrust drops ~1.4% vs v0.2.0 (geometric losses now
    applied, where v0.2.0 used only `Cd`).
- v0.3.0: `convert_nozzle(ric_nozzle)` — single arg, returns a Nozzle
    directly (was 3-tuple). `convert_geometry(ric_grains, ...)` — no
    longer takes `ric_nozzle`. `ric_to_sim_args` returns one dict only
    (was a 2-tuple); `nozzle` and `P_ambient` are inside the dict.
- v0.5.0: FMM grain support (openMotor-bridge).
    - New module `srm_1d/fmm_grain.py`. Lazy `import motorlib` from
      local checkout (`Erosive Burning Solver/openMotor/openMotor/`).
      Numba-JIT marching-squares perimeter shim replaces openMotor's
      Cython `mathlib._find_perimeter_cy` so MSVC build tools aren't
      required. Optional deps: scikit-fmm, scikit-image.
    - New `PropellantTab`-style `FmmTable` dataclass with sampled
      (reg_depth, perimeter, port_area) arrays + wall_web.
    - `from_openmotor(om_grain)` and `from_ric_grain(ric_dict)`
      wrappers for converting openMotor's FmmGrain output.
    - All 7 openMotor FMM types supported: Finocyl, Star Grain,
      Moon Burner, X Core, C Grain, D Grain, Custom Grain. (RodTube
      is a PerforatedGrain not FmmGrain — could be added later by
      extending from_openmotor to accept PerforatedGrain.)
    - `GrainSegment.fmm_table` field — when set, segment uses FMM
      table lookup. `MotorGeometry.compile_geometry_arrays` packs all
      FMM tables in CSR layout (`fmm_offset`, `fmm_reg_flat`,
      `fmm_perim_flat`, `fmm_port_flat`) plus per-cell `cell_fmm_idx`,
      `cell_segment_type`, `cell_wall_web`.
    - **Hot loop primary state changed**: `D_port[i]` → `regress[i]`
      (per-cell radial regression depth). For analytic cells: D_port
      derived from `D_bore_init + 2·regress`. For FMM: A_port and
      perimeter from table lookup; D_port = √(4·A_port/π); D_hyd =
      4·A_port/perimeter. `advance_bore_regression` advances
      `regress[i]` instead of `D_port[i]`. Burnout ramp now in
      regression-depth space.
    - **mass conservation**: `total_propellant_volume()` handles FMM
      segments via `(casting_area − port_area_init) × length`.
    - **Mixed-grain motors** (BATES + FMM in one MotorGeometry) work
      via `cell_segment_type[i]` per-cell switching.
    - Adapter `convert_geometry` dispatches on grain type: BATES stays
      analytic; FMM types route through `fmm_grain.from_ric_grain`.
    - Hasegawa Motor A reproduces v0.4.0 numbers exactly (22974.9 N·s,
      4.28s, 205.2s Isp) — analytic path is bit-for-bit consistent
      with the prior release. Finocyl smoke test: 0.08% mass balance
      error, 93k steps/s.
    - 7 new tests in `tests/test_fmm.py` (98 total).
- v0.4.0: Multi-tab propellant refactor (openMotor-aligned).
    - New `PropellantTab` dataclass with all 7 openMotor fields:
      `min_pressure`, `max_pressure`, `a`, `n`, `gamma`, `T_flame`,
      `molecular_weight` (snake_case mirror of openMotor's
      `motorlib.PropellantTab`).
    - `Propellant` no longer holds top-level `a`, `n`, `gamma`,
      `T_flame`, `molecular_weight` — those moved to `Propellant.tabs`
      (a list of `PropellantTab`). Solid props (rho_propellant, Cps,
      T_surface, T_initial) and gas transport (mu_gas, k_gas, Cp_gas)
      stay at the propellant level.
    - New methods on `Propellant`: `select_tab(P)`, `representative_tab(P)`,
      `tab_arrays()` (parallel numpy arrays for Numba consumption).
    - `burn_rate_normal(P)` now uses tab lookup.
    - Adapter `convert_propellant` preserves all .ric tabs 1:1; the
      single-tab "selecting one tab" warning is gone.
    - `_select_propellant_tab` removed from openmotor_adapter (use
      `Propellant.select_tab` instead).
    - **v0.4.0 scope: only a/n vary in the hot loop.** gamma, T_flame,
      and molecular_weight are evaluated ONCE at simulation start from
      `propellant.representative_tab()`, not per-step. This freezes gas
      thermo across the burn. The future "(b)" upgrade — per-step gas
      thermo lookup — would require:
        * `_run_time_loop` taking tab arrays for γ/T_flame/MW too
        * Recomputing R_specific, c\*, Γ_crit, gamma_R, nozzle_denom each
          step (or per-tab, with selection inside the loop)
        * Updating the choked-flow BC and EOS to use current-step values
      Hold off until calibration shows it materially improves trace fit
      vs. experimental.
    - `_run_time_loop` signature: `a_sr, n_sr` scalars replaced with
      `tab_min_p, tab_max_p, tab_a, tab_n, n_tabs` arrays.
    - `update_cell_geometry` and `advance_endface_regression` signatures
      updated identically. `grain_geometry.py` got a private
      `_saint_robert_local` helper (duplicates burn_rate.py's lookup
      logic to keep grain_geometry.py a leaf module).
