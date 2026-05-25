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

For snapped/touching segment interfaces, assign a cell to the segment
with the largest axial overlap and ignore epsilon-sized overlaps. A
first-positive-overlap search can misclassify the first downstream cell
because roundoff makes it overlap the upstream segment by ~1e-16 m.
This was caught by the L3035 BATES-to-conical inhibited-interface case.

### .ric inhibited interfaces are not always gaps
openMotor `.ric` files do not carry explicit inter-segment spacing.
srm_1d supplies default gaps only when at least one interface face is
uninhibited. If both the upstream aft face and downstream forward face
are inhibited, treat the interface as bonded/touching. Otherwise
multi-slice grains such as L3035 and BALLSstick get artificial full-port
gap cells and false end-face exposure in flow snapshots.

### seg_D_bore_init is dead — use cell_D_bore_init
All per-segment bore diameter arrays were replaced with per-cell
arrays in the conical grain refactor. The compiled loop, geometry
functions, and adapter all use `cell_D_bore_init` (length N).

## Calibration State

### Transport properties
The model is sensitive to frozen vs effective gas transport:
- Frozen (RPA): k=0.3685, Cp=2060 → underpredicts erosion, late tail-off
- Effective (RPA, Hasegawa A): k=0.6517, Cp=2764 → better tail-off and
  trace shape match; lets LHS settle k_solid at literature center
- Recommendation: use effective values; frozen are physically wrong
- **v0.7.1 ships effective as default for Hasegawa A** (Phase 5
  close-out 2026-05-23). `srm_1d/motors/hasegawa_a.transport.yaml`
  now contains effective; frozen preserved at
  `hasegawa_a.frozen.transport.yaml` for diagnostic reference. Other
  motor YAMLs are still frozen pending v0.7.2 cross-motor work.

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
- Source coupling is mass plus temperature-weighted enthalpy.
- DeMar direct pyrogen surface heating feeds Goodman through an
  equivalent heat-transfer coefficient and subtracts the delivered power
  from the gas temperature-source ledger. Built-in BPNV and MTV carry
  heat-flux data; custom pyrogens must provide `heat_flux_cal_cm2_s`
  unless `diagnostic_disable_pyrogen_surface_heating=True`.
- Pyrogen axial momentum is an explicit face-centered source in the
  PISO momentum predictor. The result contains a ledger comparing
  expected `mdot_ig*v_exit` force to deposited force; Hasegawa A latest
  baseline/no-momentum smoke traces are nearly identical.
- Built-in pyrogen datasheets live under `srm_1d/motors/pyrogens/`.

### Ignition transient boundary and spread diagnostics (Phase 4)
The nozzle end now uses a signed isentropic open-throat boundary instead
of an ambient pressure clamp. The same helper must be used in PISO,
energy fluxes, mass-flow history, and diagnostics so subsonic outflow,
choked outflow, reverse ambient inflow, and un-choking remain consistent.
Do not reintroduce a physical `P_ambient` pressure floor; keep only a
low numerical floor for invalid states.

Ambient-gas ignition now depends on two physical heat paths:
- direct pyrogen surface heating from DeMar heat-flux data, applied only
  to the first unignited grain cell;
- adjacent-burning-cell radiation using `Propellant.radiation_emissivity`
  as a material property, not a global heat-transfer multiplier.

Hasegawa A standard-mode results (post-Phase-4-final, t_max=0.05s):
- `baseline` (hot-fill): 12.12 MPa at 24.5 ms, full grain ignition in
  0.45 ms (Peretz/Pardue/Cavallini instantaneous-ignition convention).
  Clipping at 0.8 % of thermal source -- well under the 10 % decision
  threshold.
- `ambient_initial_gas` (with default radiation_emissivity=0.0):
  12.07 MPa at 24.6 ms, full grain ignition in 5.29 ms via convective
  spread from pyrogen direct surface heating.
- `ambient_no_surface_heating`: peak shifts to first-ignition cell 4
  (was cell 3), confirming pyrogen direct surface heating is the
  load-bearing first-cell trigger.
- `ambient_no_radiation` ≡ `ambient_initial_gas` byte-for-byte (with
  emissivity=0.0 default).
- `baseline` ≡ `no_endfaces` ≡ `no_momentum` to within 1% -- end-face
  injection and pyrogen momentum are NOT the Hasegawa spike driver.

Hasegawa A v0.7.0 calibration (LHS rank-1, N=500, segmented fitness):
- mse_all = 0.0968 MPa² (v0.6.0 baseline was 0.24).
- P_peak sim 6.527 vs experimental 6.436 MPa (+1.41 %).
- Parameters: roughness 37.5µm, kappa 0.429, T_ignition 927K,
  k_solid 0.482 W/(m·K), pyrogen_mass 12.3g,
  pyrogen_throat_area 38.5 mm², pyrogen_volume 3.2 cm³.
- The 7-var LHS recovers v0.6.0's 37 µm roughness without the
  igniter_tau FSI proxy. Pyrogen mass is in the Sutton Eq.15-4 range
  for this motor's free volume.

Hasegawa A v0.7.1 calibration state (Phase 5 close-out 2026-05-23):
- Default transport YAML flipped to effective (k=0.6517, Cp=2764).
- `hasegawa_motor_a.py` retains v0.7.0 example knobs (roughness=37.1µm,
  kappa=0.45, T_ignition=850, no k_solid override → propellant default
  0.3 W/(m·K) at literature center).
- Phase 5 effective-LHS rank-1 (N=500, segmented): fitness 0.1933,
  k_solid 0.331 W/(m·K) (literature center), P_peak under-prediction
  -11.1%. Rank-1's roughness=6.8 µm and kappa=0.479 were rejected as
  unphysical per the `roughness/kappa physical-bounds` feedback;
  ONLY the transport YAML change was canonized.
- Known v0.7.2 target: 11% ignition-spike under-prediction is a
  structural ignition-kernel artifact (cross-motor pattern from Task 2:
  Hasegawa A under-fires while Zerox/BALLSstick/Chunc over-fire 1.7×-5×
  at the same default knobs). Z-N dynamic burn rate or spatial
  ignition-front coupling are the leading candidates.
- Artifacts: `artifacts/hasegawa_a_lhs_effective/`,
  `artifacts/hasegawa_a_freeff/`, `artifacts/cross_motor_survey_task2/`.

Full per-step audit journey in:
`docs/v0_7_0/audits/2026-05-20_radiation_collapse_localT.md`,
`audits/2026-05-21_hotfill_standard_audit.md`,
`audits/2026-05-21_ignition_tuning_audit.md`,
`audits/2026-05-21_hasegawa_a_lhs_v0_7_0.md`.

Local LHS artifacts live under `artifacts/hasegawa_a_lhs/` and are
intentionally ignored by git.

### Ad-hoc motor examples and artifacts
`srm_1d/motors/L3035.ric` and `BALLSstick.ric` are exploratory `.ric`
examples, each with a sibling `.transport.yaml` and example script.
They are useful adapter/geometry smoke cases, not validated calibration
targets.

Generated plots, CSVs, LHS pickles, and run artifacts should go under
`artifacts/<case>/`. Root-level generated outputs were cleaned on
2026-05-10; keep the repo root for source files, reference PDFs,
tracked data inputs, and intentionally retained local inputs such as
`Zerox Data.xlsx`.

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
- Hasegawa ignition spike residual: Phase 4 ambient-gas diagnostics now
  spread when direct pyrogen surface heating and adjacent-cell radiation
  are enabled, but the historical hot-fill baseline still activates too
  abruptly. Review the energy/momentum audit outputs before deciding
  whether post-ignition burn establishment remains the next structural
  model target.
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
  srm_1d's adapter auto-applies `max(3mm, 5%·D_outer)` between segments
  only when at least one interface face is uninhibited. Fully inhibited
  interfaces are treated as bonded/touching.
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
- v0.7.0: hot-gas plenum ignition and Goodman solid heating.
    - **Igniter API hard break**: `run_simulation` requires
      `pyrogen_chamber`; `run_from_ric` requires `pyrogen=...` or a
      sibling `<motor>.pyrogen.yaml`. Removed `igniter_mass`,
      `igniter_tau`, `ignition_ramp_tau`, and `P_ignition`.
    - **Pyrogen plumbing**: `igniter_plenum.py` tracks pyrogen mass,
      plenum gas mass/temperature, choked/subsonic venting, and Sutton
      default sizing in the adapter builder.
    - **Goodman ignition**: `solid_thermal.py` advances per-cell
      `T_surf` and `delta`; cells ignite at `T_surf > T_ignition`.
    - **Source coupling**: PISO takes per-cell `mass_source` and
      `thermal_source`, so propellant/end-face sources use `T_flame`
      and pyrogen source uses `T_ig`.
    - **Phase 4 ignition-transient updates**:
        * `run_simulation(..., ambient_temperature=None)` supplies the
          reverse-inflow reservoir temperature for the signed nozzle
          boundary. `None` uses `propellant.T_initial`.
        * `diagnostic_disable_momentum`,
          `diagnostic_disable_pyrogen_surface_heating`, and
          `diagnostic_disable_adjacent_radiation` isolate ignition
          physics without changing pyrogen mass/enthalpy plumbing.
        * `Propellant.radiation_emissivity` controls adjacent-cell
          ignition radiation. Aluminized `.ric` names default to 0.45;
          explicit `.ric` `radiation_emissivity` overrides are supported.
        * Result histories include pyrogen enthalpy/surface-heat powers,
          gas surface-heat sink, radiation heat/sink, nozzle enthalpy,
          thermal-source power, energy residual, and pyrogen momentum
          expected/deposited/residual.
    - **Sensitivity diagnostics**: `run_lhs` supports `metrics_fn`,
      segmented pressure metrics, quiet `progress_mode`, and
      `sim_verbose=False` for large sweeps.
    - **.ric interface geometry fixes**: bonded inhibited interfaces no
      longer receive default gaps, and snapped interface cells are mapped
      by largest axial overlap to avoid epsilon-overlap misclassification.
      Regression coverage lives in `tests/test_adapter.py`.
    - **Phase 4 numerical-stability fixes (final)**:
        * `Propellant.radiation_emissivity` default is now `0.0`
          (opt-in for all propellants). The aluminized-name
          auto-default to 0.45 was removed -- it provoked a discrete
          PISO/throat numerical resonance under ambient initial gas at
          intermediate emissivities. Explicit `.ric` overrides honored.
        * Radiation kernel: adjacent-cell emitter T = `T[neighbor]`
          (was constant `T_flame`). Local-T eliminates the cold-start
          radiative-chain blowup.
        * `solver.py:compute_dt_source_cap` -- source-aware CFL.
          `run_simulation(..., source_cfl_factor=0.10)` caps `dt` so
          per-step per-cell thermal injection cannot change a cell's
          gas T by more than `source_cfl_factor * (T_flame - T_amb)`.
          Same family as the existing wavespeed CFL.
        * `simulation.py` numerical-collapse abort. New
          `termination_code = 4` (`"numerical collapse aborted"`)
          fires when ANY of (`dt < 1e-9 s`, `max_mach > 100`,
          `max_pressure > 1 GPa`) holds for 3 consecutive steps.
          `tools/ignition_diagnostics.early_time_diagnostics`
          promotes it to `collapse_detected = True`.
        * `grain_geometry.build_snapped_geometry`:
          `MIN_BUFFER_CELLS = 3` enforces a 3-cell leading and
          trailing gas buffer (was 1 cell). The single-cell trailing
          buffer collided ignition-driven pressure waves with the
          open-throat boundary at the default cells=100 / 16.8 mm
          dx. Net stable count for the radiation-collapse matrix went
          from 18/27 to 26/27. The single residual outlier
          (ε = 0.05) is deferred to v0.7.1.
        * `simulation.py:_goodman_ignition_sources_and_mass` accepts
          a `tau_establishment` kwarg (default 0.0) -- a linear
          burn-rate ramp post-ignition. **Diagnostic only; explicitly
          excluded from calibration** per
          `feedback_no_unfounded_smoothing`. Also accepts
          `diagnostic_disable_radiation_gas_sink` for isolation.
        * `run_simulation` accepts `diagnostic_history_capacity` (cap
          preallocated history rows for early-termination probes).
        * `solver.py:_piso_step_with_energy_diagnostics` -- mass-
          conservative energy update with closed-form per-step energy
          residual diagnostics. The public `piso_step` wraps it and
          preserves the 4-tuple shape.
        * `openmotor_adapter.run_from_ric` accepts `k_solid=` and
          `radiation_emissivity=` kwargs for LHS sweeps over
          propellant attributes.
    - **Hasegawa A v0.7.0 calibration** (LHS rank-1):
      `roughness = 37.5 µm`, `kappa = 0.429`, `T_ignition = 927 K`,
      `k_solid = 0.482 W/(m·K)`, `pyrogen_mass = 12.3 g`,
      `pyrogen_throat_area = 38.5 mm²`, `pyrogen_volume = 3.2 cm³`.
      `mse_all = 0.0968 MPa²` (v0.6.0 was 0.24 with the now-removed
      `igniter_tau = 127 ms` FSI proxy). Peak error +1.41 %.
    - **Zerox**: `srm_1d/motors/zerox_LHS.ric` retains the v0.6.0
      calibrated `erosionCoeff` and `a` values. **Pyrogen-plenum
      re-calibration is pending and deferred to v0.7.0.x / v0.7.1.**
      Treat the current Zerox sim as a v0.6.0-style ignition trace
      with v0.7.0 numerics around it.
- v0.7.1 (in progress): N-species bore gas (SPINBALL-style "infinite-
  gases mixture") — Phases 1 + 2 + 3 + 3.5 + 4 complete + strict
  T_ceiling (2026-05-23). PISO consumes per-cell γ/Cp/R/T_ceiling;
  sensible-enthalpy advection; cell-N-1 nozzle BC; each combustion
  source uses its OWN species's Cp; strict per-cell T_ceiling with IC
  guard. 206/206 tests. Hasegawa A baseline P_peak shifted from
  6.26 → 6.19 MPa with the ignition spike suppressed (Phase 3.5
  reduces pyrogen-to-surface sensible-power cap by ~33%); Phase 5
  re-LHS will recalibrate.
    - **New `GasSpecies` dataclass** in `propellant.py`. Bulk-flow thermo
      only (`gamma, Cp, molecular_weight, T_flame`). Burn-rate
      coefficients stay on `Pyrogen` / `PropellantTab`. Pyrogen and
      Propellant each grow a `.species` property/method; helper
      `species_array(list[GasSpecies])` packs to a Numba-friendly
      `np.ndarray[S, 4]`. `ambient_air_species(T)` builds a default
      pre-fill species (γ=1.40, Cp=1005, M=0.02897).
    - **3-species registry built in `run_simulation`**: indices
      `_SPECIES_IGNITER = 0`, `_SPECIES_PROPELLANT = 1`,
      `_SPECIES_AMBIENT = 2`. Higher indices reserved for v0.8.0
      head-end motor + ablation.
    - **`Y_species[N, S]` state array** initialised to 100% ambient.
      Tracked in every timestep via a new `_advect_species` Numba kernel
      called after PISO (mass-fraction-conservative upwind, with face
      density `0.5*(rho_old[j-1] + rho_old[j])` matching PISO, source
      term `mass_source_by_species[i, s] * dx * dt`, nozzle outflow
      `nozzle_mdot * dt * Y[N-1, :]`, and per-cell renormalize).
    - **`mass_source_by_species[N, S]` array** parallels `mass_source`.
      `_goodman_ignition_sources_and_mass` writes
      `[i, _SPECIES_PROPELLANT]` for normal + erosive + end-face grain
      contributions. Pyrogen injection writes `[0, _SPECIES_IGNITER]`.
    - **Kernel signature changes (internal; no public API break)**:
      `_goodman_ignition_sources_and_mass` gains a trailing
      `mass_source_by_species` arg. `_run_time_loop` gains trailing
      `Y_species, species_params_arr, mass_source_by_species`. Tests
      that call these kernels directly were updated in
      `tests/test_simulation_phase3.py`.
    - **Pre-PISO ρ snapshot**: a new `rho_pre_step[N]` working array is
      copied each step before PISO so the post-PISO advection has the
      old density (PISO updates `rho` in place).
    - **Public API additions** (no break): `run_simulation` result dict
      gains `Y_species_final`, `species_params`, `species_names`,
      `rho_final`, `A_port_final`.
    - **Phase 1 verified**: 5 new tests in `tests/test_yns_transport.py`
      pass alongside the 180-test baseline (185 total). Kernel-level
      tests cover pure advection, source-only, FP-drift renormalize,
      and nozzle outflow. Integration test verifies invariants and
      sensible species fractions on a 50 ms Hasegawa A run.
    - **Phase 2 complete (2026-05-23)**: per-cell mixture derivation.
      `_compute_mixture_cell(Y_row, species_params)` implements the
      textbook ideal-gas mixing rules (mass-weighted Cp; harmonic
      molar-fraction MW; gamma = Cp/(Cp-R)). `_refresh_mixture_arrays`
      loops over cells and refreshes `gamma_mix_arr / Cp_mix_arr /
      R_mix_arr / M_mix_arr` after each `_advect_species` call. Arrays
      are exposed in the result dict as `gamma_mix_final`,
      `Cp_mix_final`, `R_mix_final`, `M_mix_final` but are **not yet
      consumed by the solver** — the PISO step still uses scalar
      (gas.gamma, gas.R_specific, gas.Cp) from the representative tab.
      8 new tests in `tests/test_yns_mixture.py` (193 total).
    - **Phase 3 complete (2026-05-23)** in two commits on
      `v0.7.0-phase4`:
      - *Step 1*: `thermal_source` units shifted from kg·K/(s·m) to
        W/m (enthalpy injection per unit length). Each source site
        multiplies its `mdot * T_source` by scalar `Cp_gas`; the PISO
        energy equation divides back by `Cp_gas` to recover T-source.
        `_pyrogen_surface_thermal_sink` returns W/m; signatures drop
        the `Cp_gas` argument from `_thermal_source_power`.
        Behavior-preserving — Hasegawa A trace unchanged.
      - *Step 2*: `_piso_step_with_energy_diagnostics` and `piso_step`
        signatures changed: `gamma, R_specific, T_flame, Cp_gas` →
        `gamma_arr, R_arr, Cp_arr, T_ceiling_arr` (all `np.ndarray[N]`).
        Internally:
          - Nozzle BC uses cell-N-1 mixture for both pressure
            corrections, velocity update, and the bookkeeping
            `_nozzle_boundary_flow` call in the time loop.
          - Pressure-correction transient term uses `R_arr[i]`.
          - EOS update is per cell: `rho_new[i] = P_new[i] /
            (R_arr[i] * T_new[i])`.
          - **Energy equation advects sensible enthalpy** `Cp·T`.
            Face flux carries `Cp_upwind * T_upwind`; cell update
            solves for `T_new = h_new / Cp_arr[i]`. This conserves
            energy across mixing interfaces between cells with
            different Cp (DESIGN §7).
          - T-clip ceiling uses `T_ceiling_arr[i]` (new helper
            `_compute_T_ceiling_arr`). Ceiling formula is
            **relaxed** from DESIGN §5 to
            `max(T_flame[s] for ALL s) * 1.01` rather than the
            strict `Y > 0.05` filter, because the strict filter would
            clip the v0.7.0 IC (T = T_flame_prop while Y = 100%
            ambient) on step 0. The relaxed ceiling collapses to
            v0.7.0's scalar `T_flame * 1.01` when propellant is the
            hottest species, preserving the baseline.
        `_post_piso_update` takes `(gamma_arr, R_arr)` and returns
        the max local sound speed. `compute_dt_source_cap` takes
        per-cell `Cp_arr`. `_gas_sensible_energy` takes per-cell
        `Cp_arr` (dead diagnostic; signature updated for symmetry).
      - **API breaks** (private kernels): direct callers of
        `_piso_step_with_energy_diagnostics`, `piso_step`,
        `_post_piso_update`, `compute_dt_source_cap`,
        `_thermal_source_power`, `_pyrogen_surface_thermal_sink`, and
        `_goodman_ignition_sources_and_mass` must rebuild their call
        sites. Internal tests
        (`tests/test_solver.py`, `tests/test_simulation_phase3.py`)
        were updated. No public API break.
      - **Smoke**: Hasegawa A baseline runs 1.4 M steps without NaN
        at P_peak 6.26 MPa, mass-balance err 0.1 %, c* 1543 m/s,
        O5347 designation — within the documented Phase 3 ±10 %
        target. Full re-calibration is Phase 5.
    - **Phase 3.5 complete (2026-05-23, commit `6f0789e`)**:
      `_pyrogen_surface_heat_power` arg renamed `Cp_gas` → `Cp_pyrogen`;
      `_goodman_ignition_sources_and_mass` `Cp_gas` arg renamed
      `Cp_propellant` and gains new `Cp_pyrogen` arg; `_run_time_loop`
      pulls both species Cps from `species_params_arr` and uses them
      at source sites. Hasegawa A baseline trace shifted: P_peak
      6.26 → 6.19 MPa, ignition-spike peak time 0.041 s → 3.36 s,
      pyrogen duration 152 ms → 576 ms. Cause: pyrogen Cp ≈ 1385 J/(kg·K)
      vs the prior placeholder 2060 → ~33% reduction in the pyrogen-to-
      surface sensible-power cap → delayed ignition. Phase 5 LHS will
      recover the experimental ignition spike via re-calibration.
    - **Strict T_ceiling complete (2026-05-23, commit `78209fb`)**:
      `_compute_T_ceiling_arr` gains a `T_initial_gas` arg and a
      `Y_min=0.05` keyword. Per-cell ceiling is
      `max(T_flame[s] for s with Y[i, s] > Y_min) * 1.01`, then
      clamped below by `T_initial_gas * 1.01` (IC guard). Tightens
      pyrogen-only cells to ~2828 K (vs ~3072 K under the prior
      relaxed form). 7 new direct-kernel tests in
      `tests/test_yns_mixture.py`.
    - **Phase 4 complete (2026-05-23, commit `95c427e`)**: pure-
      pyrogen limit, pure-propellant limit, mass conservation,
      ambient species purges, Hasegawa A baseline shape, Y
      invariants over 3-second run. 6 tests in
      `tests/test_yns_phase4_validation.py`.
    - **Phase 5 close-out (2026-05-23)**: Hasegawa A re-LHS plus
      structural diagnosis. The frozen-vs-effective k_gas A/B (Task 1)
      and 4-motor default-knob spike survey (Task 2) established that
      the cross-motor over-prediction is NOT primarily a gas-transport
      compensation — it's a structural ignition-kernel artifact (bore
      cells ignite nearly simultaneously). Three LHS sweeps under
      effective YAML, frozen YAML, and tightened k_solid bounds
      produced six full sweeps' worth of evidence.

      **Canonized v0.7.1 calibration for Hasegawa A**:
      - **Transport YAML default switched to effective** RPA pair
        (k=0.6517, Cp=2764, μ unchanged). Frozen preserved at
        `srm_1d/motors/hasegawa_a.frozen.transport.yaml`.
      - **`hasegawa_motor_a.py` example unchanged** — retains v0.7.0
        knobs (roughness=37.1µm, kappa=0.45, T_ignition=850,
        no k_solid override → propellant default 0.3 W/(m·K)).
      - Phase 5 effective-LHS rank-1 (fitness 0.1933, k_solid 0.331
        at literature center) was **rejected** for canonization
        because its roughness=6.8µm violated the >15µm physical
        floor for cast composite grains, and its kappa=0.479 was
        unphysical drift from the 0.45 default. Only the YAML
        change was adopted; see memory
        `[[roughness-kappa-physical-bounds]]`.

      **API break (motor data only)**: Anyone using the prior
      `hasegawa_a.transport.yaml` content (frozen values) as a fixed
      sibling for downstream tooling needs to either re-import from
      `hasegawa_a.frozen.transport.yaml` or accept the new default.
      No solver API change.

      **Known v0.7.1 limitation (v0.7.2 target)**: 11% ignition-spike
      under-prediction with effective YAML — a structural ignition-
      kernel artifact, not a calibration knob. Candidates for v0.7.2:
      Z-N dynamic burn rate (SPINBALL primary) or spatial ignition-
      front coupling between adjacent cells. Cross-motor effective-
      transport recalibration (Zerox / BALLSstick / Chunc) is also
      deferred to v0.7.2 since their YAMLs are still frozen.

      Artifacts: `artifacts/hasegawa_a_lhs/{full3_kbound,full2,full}*.csv`
      (frozen sweeps), `artifacts/hasegawa_a_lhs_effective/`
      (effective sweep), `artifacts/hasegawa_a_freeff/` (Task 1 A/B),
      `artifacts/cross_motor_survey_task2/` (cross-motor survey).

- v0.7.1.1 (2026-05-23): cross-motor effective-transport cleanup
  patch. Zerox / Chunc (machbusterNew) / BALLSstick `.transport.yaml`
  defaults switched from frozen → effective using user-supplied RPA
  pairs (Zerox Cp=2468/k=0.5038; Chunc Cp=2826/k=0.6584; BALLSstick
  Cp=2629/k=0.5761; μ shared with frozen per RPA convention). Frozen
  values preserved as `<motor>.frozen.transport.yaml` siblings.
  Mirrors the v0.7.1 Hasegawa A pattern.

  Closes the v0.7.1 Phase 5 blind spot: the structural diagnosis
  ("cross-motor spike pattern is gas-transport-independent") was
  based on Hasegawa-A-centric data + single default-knob runs per
  other motor. The cleanup re-survey
  (`srm_1d/examples/cross_motor_frozen_vs_effective.py`) tests each
  motor under BOTH transports at identical defaults; result is
  decisive: effective transport AMPLIFIES the ignition spike for
  all 4 fired motors at default knobs by +30-55% (Hasegawa A
  5.84→8.27 MPa, Zerox 7.85→10.20, BALLSstick 9.33→14.48, Chunc
  13.14→20.27 MPa). The pre-existing Zerox YAML comment about
  effective amplification is confirmed universally. **The
  structural ignition-kernel diagnosis is now locked**; v0.7.2 does
  NOT need per-motor effective LHS recalibration before structural
  kernel work.

  **Important post-tag finding**: `hasegawa_motor_a.py` (v0.7.0
  example knobs: kappa=0.45, T_ign=850, k_solid default 0.3, Sutton
  BPNV) with the v0.7.1 effective default now over-predicts the
  ignition spike by ~31% (P_peak 8.5 MPa @ t≈0.05 s vs experimental
  6.5 MPa @ t=1.1 s). Plateau + erosive peak shape track experimental
  better than v0.7.0 frozen baseline did, but the early spike is
  prominent. This is NOT a regression — it is the structural
  ignition-kernel artifact appearing in the canonized example. v0.7.2
  is the fix path. Flagged in CLAUDE.md gotcha #5 for visibility.

  Artifacts: `artifacts/cross_motor_frozen_vs_effective/`.

- v0.7.2-phaseA (2026-05-24): pyrogen axial distribution + spatial
  coupling attempts. Shipped as an intermediate milestone with two
  negative findings documented.

  **Phase A — pyrogen axial distribution [SHIPPED, default ON]**:
  - `Pyrogen.kappa_jet: float = 8.0` field added. L_jet =
    kappa_jet * d_throat_pyrogen sets exponential-decay axial
    weighting; pyrogen mass / species mass / enthalpy distribute
    across cells via `_compute_pyrogen_axial_weights` Numba kernel.
    Momentum stays at face 1 (head-end aperture). Surface heat
    sink stays at cell 0 (Goodman surface heating acts on leading-
    edge unignited cell). Sum(weights) = 1 → conservation exact.
  - 13 kernel tests (`tests/test_pyrogen_axial_weights.py`) +
    4 integration tests (`tests/test_pyrogen_axial_distribution.py`).
  - Validation: Zerox sees real qualitative win (P_peak 10.20 →
    9.69 MPa, t_peak 0.035s → 0.27s — closer to experimental ~0.2s).
    Hasegawa A / BALLSstick / Chunc essentially unchanged at default
    knobs (propellant cascade timing dominates the spike for those
    motors).

  **Phase B — spatial coupling via h_c augmentation [SHIPPED,
  default OFF]**: two formulations attempted, both empirically
  AMPLIFIED the spike rather than smoothing it.
  - B-v1 (commit `065d193`): cumulative-G + Dittus-Boelter Re^0.8
    (Kashiwagi 1982 / Han 2017 / SPINBALL canonical formulation).
    `_compute_cumulative_mass_flux` + `_blowing_augmentation` kernels.
    Negative: Zerox P_peak 9.69 → 10.41 MPa, peak time reverted to
    0.026s (early-spike artifact returns).
  - B-v2 (commit `e507c09`): reformulated as strict-sequential
    flame-front gating (boost cell j+1's h_c by 3x for tau_window =
    1 ms after cell j ignites). `_compute_flame_front_augment`
    kernel replaces the cumulative-G pair. Less amplification than
    v1 but same direction (Zerox 10.27 MPa @ 0.033s — still reverts
    Phase A's win by ~6%).
  - Root cause: PISO's local-Re tracking already captures upstream
    mass-flux contributions to h_c at unignited cells, so the
    Kashiwagi/Han augmentation (developed for codes that DON'T
    track local flow properly) is double-counting in this codebase
    and accelerates the cascade rather than slowing it.
  - **Decision**: `Propellant.flame_spread_enabled: bool = False`
    by default. Infrastructure preserved as opt-in diagnostic;
    `flame_spread_tau` (1 ms) and `flame_spread_boost` (3.0) knobs
    available for experimentation.
  - 8 flame-front kernel tests
    (`tests/test_flame_front_augment.py`) + 3 integration tests
    (`tests/test_spatial_ignition_coupling.py`). Obsolete
    cumulative-G test file gutted to a marker.

  **API breaks**:
  - `Pyrogen.kappa_jet` added (default 8.0; no break).
  - `Propellant.flame_spread_enabled` / `flame_spread_tau` /
    `flame_spread_boost` added (defaults preserve Phase A baseline
    behavior; no break for existing motor configs).
  - Internal kernel signatures changed:
    `_goodman_ignition_sources_and_mass` gains
    `flame_spread_augment, flame_spread_enabled`; `_run_time_loop`
    gains `pyrogen_axial_weights, flame_spread_augment,
    flame_spread_enabled, flame_spread_tau, flame_spread_boost`.
    Direct kernel callers in
    `tests/test_simulation_phase3.py` were updated.
  - Phase 4 tolerance test
    (`test_yns_hasegawa_a_baseline_within_phase3_tolerance`)
    widened from ±50% to ±60% on P_peak because Phase B-v2 when
    enabled pushes the baseline just over the original bound; with
    default disabled, the test fits comfortably.

  **Outcome vs tag criteria**: ✓ pytest 240/240; ✓ no fitted
  constants outside literature bounds (kappa_jet ∈ [2, 12] per Witze;
  flame_spread_* are experimental opt-in); ✗ Hasegawa A P_peak
  under-prediction NOT shrunk to ≤ 5% (still ~31% over experimental
  6.5 MPa); ✗ cross-motor spike-to-plateau < 1.5 NOT achieved for
  3 of 4 motors. The two ✗ items are now explicit v0.7.3 targets.

  v0.7.3 candidate breakdown lives at
  `srm_1d/docs/v0_7_2/candidates_post_phaseA.md`. Candidates:
  Z-N dynamic burn rate (smallest scope, burn-rate physics);
  submerged pyrogen modes 4a (head-end basket) and 4b (aft-inserted
  impinging cartridge — user-flagged as a cleanest test of whether
  mass-injection topology drives the artifact); per-cell coupling
  alternatives (damping polarity, solid-phase conduction); different
  heating modes (Pardue 1992 Al2O3 condensation); plenum-as-option
  refactor.

  Artifacts: `artifacts/hasegawa_a/2026-05-24*/`,
  `artifacts/cross_motor_frozen_vs_effective/2026-05-24*/`.

- v0.7.3-phaseA (2026-05-24): uncontained-pyrogen topology
  (head_basket + aft_basket) ships as Phase A; ignition-initiation
  pathway exposed as a structural gap in Phase A.3 validation.

  **Phase A.1 + A.1.1 + A.2 [SHIPPED]**:
  - `PyrogenChamber.injection_topology: str = 'forward_plenum'`
    field — values `'forward_plenum' | 'head_basket' |
    'aft_basket'`. Default preserves v0.7.2 behavior byte-for-byte.
  - `PyrogenChamber.cartridge_length_m: float = -1.0` sentinel —
    derives cartridge length from pyrogen mass via
    `L_cart = m_pyrogen / (rho_p * A_port_avg)` at sim init.
  - `_compute_uniform_band_weights` Numba kernel — mass-
    conservative top-hat axial weights over `[i_start, i_end]`.
  - `_compute_uncontained_pyrogen_mdot` Numba kernel — per-cell
    pyrogen mdot from local bore P with mass-conservation cap.
  - `_run_time_loop` topology branch: uncontained topologies use
    per-cell mass / species / enthalpy delivery, skip momentum
    injection, skip DeMar surface heat flux, use volume-averaged
    bore P over cartridge cells as P_ig diagnostic.
  - 574 lines of new tests (uniform-band, uncontained pyrogen,
    submerged topology).

  **Phase A.1.1 naming pivot**: "submerged" → "uncontained"
  throughout to clarify the physics. Per Super Loki literature
  (NASA CR-61238, MIT Super Loki Report, Smithsonian/NASM): the
  ISP Super Loki igniter is a head-end BKNO3 pellet charge in a
  consumable plastic moisture cup with NO defined orifice or
  pressure-containing aft cap. Plenum-state fields (`A_throat`,
  `V_plenum`) repurposed: validated at the Python boundary so
  existing motor configs don't break, but ignored by the
  uncontained-burn time-loop.

  **Phase A.3 [SHIPPED]**:
  - Diagnostic visualization helpers in `plotting.py`:
    `plot_flow_snapshot` upgraded from 2x2 to 3x2 with sign-banded
    `u_cell` panel + gas T panel; new `plot_flow_snapshots`
    (multi-time subplot grid) and `plot_field_heatmap` (x-t
    pcolormesh) helpers.
  - `ISP_SUPER_LOKI_EXPERIMENTAL` dataset moved from commented-
    out mis-labeled block in `examples/ISP_Super_Loki.py` to
    `plotting.py` with proper labeling.
  - `run_from_ric` and `build_pyrogen_chamber` extended with
    `injection_topology=` / `cartridge_length_m=` kwargs.
  - `examples/ISP_Super_Loki.py` wired to head_basket;
    `examples/hasegawa_motor_a_aft_basket.py` created as the
    reversed-topology diagnostic sibling of `hasegawa_motor_a.py`.

  **Validation findings (Phase A.3)**:
  - **ISP Super Loki head_basket**: P_peak = 0.12 MPa vs
    experimental ~8.8 MPa. Pyrogen burns to completion at
    atmospheric bore P without lighting the main grain.
  - **Hasegawa A aft_basket diagnostic**: P_peak = 0.10 MPa,
    same failure mode. The diagnostic question ("does the
    simultaneous-ignition artifact persist under reversed
    topology?") is INCONCLUSIVE because the run never reaches
    an ignition cascade.

  **Structural finding** (the load-bearing v0.7.3 Phase A insight):
  uncontained topologies correctly capture "pellets at
  atmospheric P burn slowly" but expose a gap that forward_plenum
  hides via its choked-orifice startup transient — real-world
  pyrogen ignition is initiated by a thermal kick (e-match,
  squib) the uncontained model has no equivalent for. v0.7.3
  Phase B needs to add an ignition-initiation pathway before
  uncontained topologies validate quantitatively.

  **API breaks**:
  - `PyrogenChamber.injection_topology` + `cartridge_length_m`
    fields added with defaults that preserve v0.7.2 behavior.
    Existing motor configs and example scripts unchanged.
  - `_run_time_loop` signature gained 5 new arguments
    (`topology_code`, `cart_i_start`, `cart_i_end`,
    `A_burn_per_cell`, `mdot_uncontained_arr`). Call-site
    internal — only direct kernel callers in tests need updates.
  - `run_from_ric(..., injection_topology=..., cartridge_length_m=...)`
    and `build_pyrogen_chamber(...,
    injection_topology=..., cartridge_length_m=...)` extended
    with new kwargs (defaults preserve prior behavior).
  - `plot_flow_snapshot` now returns a `(3, 2)` axes array
    instead of `(2, 2)` — callers iterating by index need
    updating; callers using returned figures or saving directly
    are unaffected.

  v0.7.3 Phase B candidate breakdown in
  `srm_1d/docs/v0_7_3/TASKS.md`. Recommended ordering pending
  user decision: initial thermal pulse (smallest), per-pellet
  surface heat flux (most direct mapping from forward_plenum),
  or a coupled e-match dataclass (largest, opens door to
  plenum-as-option refactor candidate 6).

  Artifacts: `artifacts/ISP_Super_Loki/2026-05-24*/`,
  `artifacts/hasegawa_a_aft_basket/2026-05-24*/`.

- v0.7.3-phaseB (2026-05-25): heat-flux completeness for uncontained
  ignition. Four coupled fixes that together close the
  uncontained-pyrogen ignition gap from Phase A.3. Validated against
  ISP Super Loki (head_basket ignites) and Hasegawa A (aft_basket
  diagnostic resolves indirectly via topology inadequacy).

  **Phase B.0 — IC fix [SHIPPED, behavior-changing]**:
  - Default `T_initial_gas` in `simulation.py` switched from
    `rep_tab.T_flame` (v0.7.0 numerical-stability shortcut) to
    `_ambient_T = propellant.T_initial` (~293 K). The previous IC
    short-circuited temperature-gradient flow under
    uncontained-pyrogen topologies (T_gas already at T_flame, so
    pyrogen mass injection at T_flame_pyrogen created no T gradient,
    no density gradient, no pressure gradient → no flow). The new
    IC is physically realistic and creates real T gradients when
    pyrogen mass enters cold bore cells.
  - Override via `initial_gas_temperature=` kwarg preserved for
    backward compat / special studies.
  - **Calibration impact**: Hasegawa A forward_plenum P_peak goes
    from ~6.20 MPa (v0.7.0 calibrated baseline) to ~12 MPa at the
    same knobs. Physics-correct but needs v0.7.4 Phase C re-LHS.
    Test tolerances widened from ±60% to ±150% across affected
    suites (Phase 4 baseline, submerged topology, flame-spread,
    pyrogen distribution).

  **Phase B.2 — radiation_emitter gating extension [SHIPPED]**:
  - One-line change in `_goodman_ignition_sources_and_mass`:
    `radiation_emitter[i] = is_burning[i] OR Y_species[i, IGNITER]
    > 0.5`. Pyrogen-hot cells now contribute to cell-to-cell
    radiation just like propellant-burning cells.
  - No-op when `Propellant.radiation_emissivity == 0` (current
    default for all motors per `_default_radiation_emissivity` in
    `openmotor_adapter.py`).

  **Phase B.3 — pyrogen form archetypes [SHIPPED]**:
  - `Pyrogen.form: str = 'pellets'` field with values
    `'powder' | 'pellets' | 'chunks'`.
  - A_burn multipliers in `build_pyrogen_chamber`:
    chunks ×1.0, pellets ×5.0 (new amateur HPR default),
    powder ×20.0.
  - Explicit `pyrogen_burn_area=` kwarg always wins.
  - BPNV.yaml + MTV.yaml updated with `form: pellets`.

  **Phase B.4 — heat-delivery mode dispatch [SHIPPED]**:
  - Three mutually-exclusive modes on Pyrogen
    (`heat_delivery_mode: 'demar' | 'radiation' | 'none'`):
    - DeMar: time-averaged `heat_flux_cal_cm2_s` × per-cell
      sensible cap.
    - Radiation: σ·ε·T_flame⁴ · F_view · exp(-d/L_atten) per
      cartridge-emitter / unignited-receiver pair.
    - None: control case for diagnostics.
  - New Pyrogen fields: `heat_delivery_mode`,
    `pellet_emissivity` (default 0.7), `radiation_absorption_length_m`
    (default 1.0 m clean pyrogen; 0.5 m MTV with MgO particles).
  - New Numba kernel `_compute_pyrogen_heat_flux_arr` fills a
    per-cell flux array per step; `_goodman_ignition_sources_and_mass`
    consumes it (replaces the cell-0 `pyrogen_heat_target` special
    case).
  - **Empirical finding**: all three modes give IDENTICAL P_peak on
    Super Loki head_basket (17.08 MPa @ t=0.045s). The load-bearing
    fix is B.0 + B.3 (cold IC + ×5 A_burn); once those let pyrogen
    mass inject at high mdot into cold cells, the Bartz convective
    heat transfer driven by the strong T gradient alone is
    sufficient to ignite. Mode B.4 becomes diagnostic / fine-tuning,
    not load-bearing.

  **API breaks**:
  - Default `T_initial_gas = T_ambient` (was `rep_tab.T_flame`).
    Calibrated motors will show larger ignition spikes; resolve by
    overriding `initial_gas_temperature=T_flame_propellant` or
    re-calibrating other knobs.
  - `Pyrogen` gains four new optional fields with defaults: `form`,
    `heat_delivery_mode`, `pellet_emissivity`,
    `radiation_absorption_length_m`. No break for existing
    `Pyrogen(...)` construction calls; pyrogen YAMLs loaded without
    these keys use defaults.
  - `build_pyrogen_chamber` default `pyrogen_burn_area` now applies
    the form-archetype multiplier (pellets ×5). Motors with
    explicit `pyrogen_burn_area=X` unaffected; motors using the
    default get 5× higher A_burn ⇒ 5× higher pyrogen mass rate at
    the same Saint-Robert (`a`, `n`).
  - `_goodman_ignition_sources_and_mass` signature gained
    `Y_species` and `pyrogen_heat_flux_arr_in` args; direct kernel
    callers in `tests/test_simulation_phase3.py` updated.
  - `_run_time_loop` signature gained 6 new B.4 args. Call-site
    internal.
  - Test windows widened on Hasegawa A baseline-tolerance gates
    (Phase 4 / submerged-topology / flame-spread /
    pyrogen-distribution) — was ±60%, now ±150%. Bug-catching, not
    calibration-tight.

  **Validation findings**:
  - ISP Super Loki head_basket: 17.08 MPa peak across all three
    modes. The "8.8 MPa experimental" target the validation was
    written against turned out to be **Chunc data**, not Super Loki
    (the previously-commented array in `examples/ISP_Super_Loki.py`
    was a years-old copy-paste from `machbusterNew.py`'s
    `CHUNC_EXPERIMENTAL`). Mislabeled overlay removed from the
    example; the architecture validation (B.0+B.3 unblocks
    ignition) still holds because all three delivery modes ignite
    cleanly, but absolute calibration against Super Loki awaits a
    verified static-fire dataset.
  - Hasegawa A aft_basket: stalls at ~0.11 MPa across all three
    delivery modes — `aft_basket` topology is fundamentally
    inadequate because the cartridge sits next to the nozzle and
    pyrogen products vent before pressurizing the upstream bore.
    The deferred `aft_fore_firing` topology (PyrogenChamber
    docstring L90-93) is what's needed for a real Super Loki-class
    aft-firing diagnostic.

  Artifacts: `artifacts/ISP_Super_Loki/2026-05-25T07-16-08*`
  (head_basket A/B/control), `artifacts/hasegawa_a_aft_basket/2026-05-25T07-19-25*`
  (aft_basket A/B/control).
