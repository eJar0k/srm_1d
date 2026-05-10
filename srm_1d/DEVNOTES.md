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

### End-face injection uses a partition-of-unity hat function (v0.6.0+)
`update_cell_geometry` distributes each end-face's mass over the two
adjacent cells with weights `w = 1 − |x_face − x_center|/dx`, summing
to 1.0 by construction. This is the correct response to the v0.6.0
snapping discretization — `build_snapped_geometry` puts segment edges
on cell boundaries by design, which would cause the prior interval-
containment kernel (`x_lo ≤ x_face ≤ x_hi`) to silently double-count
because both adjacent cells satisfy the closed-interval test.

The new kernel is gated by `tests/test_endface_conservation.py` which
verifies <0.1% mass error against `2·ρ·r·A_face` for both snapped and
deliberately-unsnapped grids across resolutions 50, 100, 200, 500.

### Boundary clamp invariant in update_cell_geometry
The face-distribution kernel has a hard-edge clamp:
`if i == 0 and x_face < x: weight = 1.0` (and the symmetric N-1 clause).
This assumes the face is **inside or at the edge of the motor domain**.
If something pushes `x_fwd < 0` or `x_aft > L_motor` (e.g., a regression
overshoot or a malformed segment), the clamp will inject the face's
mass with weight=1 into the boundary cell, hiding the bug. Don't use
this kernel for grains whose faces could drift outside the domain
without a domain check upstream.

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

### Igniter (v0.7.0 pyrogen plenum model)
The v0.7.0 igniter is a forward hot-gas pyrogen plenum. A
`PyrogenChamber` burns a single-tab `Pyrogen`, vents through a
choked/subsonic orifice into cell 0, and feeds Goodman per-cell
solid heating. Grain cells ignite when `T_surf > T_ignition`.
- The v0.6.0 exponential igniter API was removed.
- Source coupling is mass plus temperature-weighted enthalpy; igniter
  momentum is intentionally deferred.
- Built-in pyrogen datasheets live under `srm_1d/motors/pyrogens/`.

### Zerox motor (LHS-calibrated, v0.6.0)
Forward-Finocyl + aft-BATES, ~1.45 kg "Risky Batman V3" propellant.
Static-fire data is the calibration ground truth — the openMotor
.ric values for `erosionCoeff` and propellant `a` were both
significantly off for this firing.

**Calibrated `.ric` values (`srm_1d/motors/zerox_LHS.ric`):**
- `erosionCoeff = 2.585e-10 m/(s·Pa)` — **2.34× the openMotor default**
  (1.105e-10). Risky Batman V3 / phenolic ablates much faster than
  the original .ric assumed. Real motor's throat opens to ~34.5 mm
  by burnout (Δr ≈ 9 mm) vs 29.65 mm with the original coefficient.
- `propellant a = 4.634e-6 m/s/Pa^n` — 0.917× the openMotor default
  (5.054e-6). The propellant is ~8% over-rated.

**v0.7.0 ignition calibration:**
- The v0.6.0 igniter knobs are gone. Zerox now uses pyrogen plenum
  parameters (`pyrogen_mass`, `pyrogen_throat_area`, `T_ignition`) and
  needs a fresh Phase 4 LHS.
- See `srm_1d/examples/zerox.py` for the current pyrogen-based run.

**Pinned-variant sensitivity** (LHS where one variable is held at the
Hasegawa-A inherited value; reveals which knobs are essential):

| Pinned at | Best MSE | Verdict |
|---|---|---|
| Main 6-var (no pin) | 0.071 | reference |
| `erosion_coeff_scale = 1.0` | 0.357 | **5× worse** — dominant lever |
| `a_scale = 1.0` | 0.149 | 2× worse — moderate |
| `pyrogen_throat_area` | TBD | v0.7.0 LHS pending |
| `T_ignition` | TBD | v0.7.0 LHS pending |
| `kappa = 0.45` | 0.104 | 1.5× — minor |
| `pyrogen_mass` | TBD | v0.7.0 LHS pending |

**Residual structural artifacts (NOT parametric — cannot be tuned away):**
- Spike behavior now depends on the pyrogen plenum and Goodman ignition
  coupling; Phase 4 validation will quantify the remaining residual.
- Sharp step at t≈1.9s — fin-burnout transition in the FMM finocyl
  model. Real motor has it smoothed/absent. The FMM table is
  axially uniform within each segment (per `fmm_grain.py:329-414`),
  which is geometrically reasonable for the user's fin layout but
  treats fin consumption as a single-radius event. Possible v0.7.0+
  fix: per-cell axial FMM table.

**Conservation diagnostic** (the trick that pinned the issue):
sweeping `a × {0.80..1.05}` showed `∫P dt ≈ 19.6 MPa·s` was nearly
invariant across the sweep (per 0D equilibrium ∫P dt ≈ m·c\*/Ā_t).
Experimental ∫P dt = 16.4 MPa·s → 20% gap could only come from
larger A_t in real motor (c\* drop unphysical at >15%). The
erosion sweep then settled it: erosion × 2.5–3.0 brackets the
experimental impulse-time integral. See `zerox_diagnostic.py`,
`zerox_a_sweep.py`, `zerox_erosion_sweep.py`, `zerox_lhs.py`.

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
- No inter-segment spacing in file; openMotor is 0-D and doesn't model gaps.
  srm_1d's adapter auto-applies `max(3mm, 5%·D_outer)` between segments.
- No gas transport properties in file; supplied via sibling
  `<motor>.transport.yaml` (mu, k, Cp) auto-discovered by `run_from_ric`,
  or via explicit `gas_props={...}`.
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
- v0.6.0: dynamic discretization, .ric data motors, kernel rewrites.
    - **Discretization API**: `convert_geometry(ric_grains, N_cells, spacing)`
      → `convert_geometry(ric_grains, target_propellant_cells)`.
      `cells_per_seg`/`cells_per_gap` defaults removed. `ric_to_sim_args`
      and `run_from_ric` lose `N_cells`/`spacing` kwargs the same way.
      Adapter routes through `build_snapped_geometry`.
    - **`build_snapped_geometry`** is the canonical builder: takes a
      `segments_spec` list (per-segment `D_bore_fwd`, `length`, optional
      `gap_after`/`fmm_table`) plus `target_propellant_cells`, returns a
      `MotorGeometry` with integer-snapped segment lengths and a
      Nyquist-CFD-clamped dx.
    - **End-face injection kernel**: replaced interval-containment
      (`x_lo ≤ x_face ≤ x_hi`) with a partition-of-unity hat function
      (`weight = 1 − |x_face − x|/dx`). See "Critical Gotchas". Gated
      by `tests/test_endface_conservation.py`.
    - **Volumetric overlap accumulator** in `update_cell_geometry`: the
      `break` that used to stop after the first matching segment is
      gone — cells straddling two segments now receive C_burn from
      both, which fixes a silent under-count at narrow gaps.
    - **Igniter model rewrite**: legacy Saint-Robert flame-front tracking
      was removed in v0.6.0, then the v0.6.0 decay placeholder was
      replaced in v0.7.0 by `PyrogenChamber` + Goodman ignition.
    - **Named motors moved to `srm_1d/motors/*.ric`**:
      `make_hasegawa_motor_A/B/C_geo/_nozzle`,
      `make_hasegawa_propellant_1`, `make_king_propellant_4525`,
      `make_bates_motor`, `make_single_cylinder`, `make_conical_grain`,
      `make_stepped_motor`, `make_example_bates` are all DELETED.
      Use `run_from_ric('srm_1d/motors/<motor>.ric')` or build geometry
      directly via `build_snapped_geometry`.
    - **Transport YAML sibling**: `<motor>.transport.yaml` next to a
      `.ric` file is auto-discovered by `run_from_ric` and supplies
      `mu`, `k`, `Cp` (which the openMotor schema doesn't carry).
      Override with `transport_path=...` or `gas_props=...`.
    - **`time_offset` removed** from `HASEGAWA_MOTOR_A_EXPERIMENTAL`
      global dict. Pass `time_offset=` to `plot_pressure`/`plot_summary`
      explicitly.
    - **Sensitivity tooling** lives at `srm_1d/tools/sensitivity.py`
      (`run_lhs`, `mse_fitness`, etc.). Example wrapper at
      `srm_1d/examples/hasegawa_a_lhs.py`.
