# 1D SRM Internal Ballistics Simulator (srm_1d v0.5.0)

A transient 1D finite-volume solver for solid rocket motor internal
ballistics with the Ma et al. (2020) erosive burning model.

**Performance:** ~45–90k steps/s (compiled time loop with Numba JIT;
FMM grains run faster than cylindrical due to fewer per-step ops)
**Tests:** 98 via pytest | **Adapter:** reads openMotor .ric files
including all FMM grain types (Finocyl, Star, Moon, X, C, D, Custom)

## File Structure

```
srm_1d/
├── __init__.py              # Public API, v0.5.0
├── solver.py                # PISO: TDMA, pressure correction, CFL
├── burn_rate.py             # Ma et al.: Haaland → Gnielinski → bisection
├── grain_geometry.py        # BATES, conical, stepped, FMM; per-cell regress[i]
├── propellant.py            # Propellant (multi-tab) + gas properties + thermo utilities
├── nozzle.py                # openMotor-aligned Nozzle: thrust, Isp, CF, throat erosion
├── fmm_grain.py             # openMotor FmmGrain bridge — extracts FmmTable from any FMM grain
├── simulation.py            # Compiled time loop (_run_time_loop @njit)
├── plotting.py              # Pressure/thrust/snapshot/comparison plots
├── openmotor_adapter.py     # .ric reader, CSV export, grain metrics
├── tests/  (6 files, 98 tests)
└── examples/ (hasegawa_motor_a.py, bates_4seg.py)
```

## Architecture

### Compiled Time Loop
```
run_simulation(geo, propellant, nozzle, ...)   Python: setup, allocate
  → _run_time_loop(55 args)                    @njit: entire while loop
      ├─ advance_endface_regression   Numba→Numba direct calls
      ├─ advance_bore_regression      (mass-conserving burnout ramp)
      ├─ update_cell_geometry         (axial overlap fraction)
      ├─ compute_burn_rates           (every N steps)
      ├─ _ignition_source_and_mass    (fused: ignition+source+mass)
      ├─ piso_step                    (staggered grid, 2 corrections)
      └─ _post_piso_update            (fused: flow+friction+CFL)
  → wrap results: dict + summary    Python: snapshots, per-grain data
```

### Staggered Grid
Scalars (P,ρ,T) at N cell centers. Velocities at N+1 faces.
Head end: wall (u=0). Nozzle: choked flow BC (ṁ ∝ P·A_t).

### Geometry System
- `GrainSegment(D_bore_fwd, D_bore_aft, ..., fmm_table=None)` — per-segment
  config. Cylindrical/conical use bore diameters; FMM segments attach
  an `FmmTable` instead.
- Primary per-cell state: `regress[i]` — radial regression depth (m).
  Unifies cylindrical, conical, and FMM paths. `D_port` is derived per
  step (`D_bore_init + 2·regress` for analytic; `√(4·A_port/π)` for FMM).
- `cell_segment_type[i]` — `0` for analytic, `1` for FMM lookup.
- `cell_D_bore_init` array (length N) — per-cell initial bore diameter
  (analytic only; ignored for FMM cells).
- `cell_wall_web[i]` — radial regression at burnout (per-cell). For
  cylindrical/conical: `(D_outer − D_bore_init)/2`. For FMM: from the
  attached `FmmTable.wall_web`.
- `grain_frac` — axial overlap fraction for smooth end-face transitions.
- `f_active` — radial burnout ramp (applied to regression rate).
- Auto-inhibition: touching segments get shared faces inhibited.

### FMM Grain Support (v0.5.0)
- All 7 of openMotor's FMM grain types supported: Finocyl, Star Grain,
  Moon Burner, X Core, C Grain, D Grain, Custom Grain.
- `srm_1d.fmm_grain.from_openmotor(om_grain)` runs openMotor's
  regression-map pipeline on a populated FmmGrain instance, samples
  perimeter and port area at fine resolution, and returns an
  `FmmTable` ready for the @njit hot loop.
- Tables are packed CSR-style (`fmm_offset[]`, `fmm_reg_flat[]`,
  `fmm_perim_flat[]`, `fmm_port_flat[]`) per-grain; per-cell
  `cell_fmm_idx` pinpoints which segment each cell uses. This keeps
  the heavy data per-grain while preserving per-cell regression.
- Mixed-grain motors (BATES + Finocyl in the same MotorGeometry) work
  out of the box — `cell_segment_type[i]` switches paths per cell.
- Hydraulic diameter for FMM cells: `D_h = 4·A_port/perimeter`
  (correct for non-circular ports; reduces to `D` for circular).
- Optional dependencies: `scikit-fmm` (regression maps), `scikit-image`
  (openMotor Custom grain). The Cython `mathlib._find_perimeter_cy`
  is replaced with a Numba-JIT marching-squares perimeter, so MSVC
  build tools aren't required.

### Igniter Model
Pressure-dependent Saint-Robert: r_ign = a_ign × P^n_ign
Surface area tapers as (mass_remaining/mass_initial)^(2/3).
Parameters: igniter_mass, igniter_a, igniter_n, igniter_rho, igniter_A_burn.
Default coefficients approximate BKNO3.

### Throat Evolution (in-loop, not post-processing)
Erosion: rate = erosion_coeff [μm/(s·MPa)] × P [MPa] × 1e-6
Slag: rate = slag_coeff [(m·MPa)/s] / P [MPa]
dD/dt = 2 × (erosion - slag). Feeds back into nozzle BC each step.
Coefficients live on the `Nozzle` object passed to `run_simulation`.

### Nozzle (openMotor-aligned, v0.3.0)
- Field names mirror openMotor's `motorlib.Nozzle`: `D_throat`, `D_exit`,
  `efficiency`, `div_angle`, `conv_angle`, `throat_length`,
  `erosion_coeff`, `slag_coeff` (snake_case'd; units kept human-readable
  internally — adapter handles conversion).
- Thrust uses openMotor's adjusted-CF formula:
    CF_adj = divLoss × throatLoss × efficiency
           × (skinLoss × CF_ideal + (1 − skinLoss))
- Loss factors: `divergence_losses()` (cosine half-angle),
  `throat_losses(d_throat)` (RasAero throat-aspect ratio),
  `skin_losses()` (constant 0.99).

### Propellant (multi-tab, openMotor-aligned, v0.4.0)
- `Propellant` holds `tabs: list[PropellantTab]`. Each tab carries
  `(min_pressure, max_pressure, a, n, gamma, T_flame, molecular_weight)` —
  a 1:1 snake_case mirror of openMotor's `motorlib.PropellantTab`.
- `Propellant.select_tab(P)` — hard-switchover lookup matching openMotor's
  `getCombustionProperties` (strict containment, closest-boundary fallback).
- `Propellant.representative_tab(P_expected=None)` — pick one tab for
  sim-start gas thermo. Single-tab propellants always return their only tab.
- **v0.4.0 scope**: only `a` and `n` vary in the hot loop; gas thermo
  (γ, T_flame, MW) is frozen at sim start from the representative tab.
  See DEVNOTES "API Breaking Changes Log" for the full-per-step upgrade path.

### openMotor Adapter
- `load_ric()` → parse YAML with Python-tag handling
- `convert_propellant()` → Propellant with all .ric tabs preserved 1:1
  (MW g/mol→kg/mol per tab)
- `convert_geometry()` → MotorGeometry (auto-gap, auto-N_cells)
- `convert_nozzle()` → Nozzle (erosionCoeff m/(s·Pa) → μm/(s·MPa))
- `run_from_ric()` → result, perf, nozzle, geo, prop
- `save_csv()` → Time, Kn, P, Force, Mass Flow, per-grain regression/web
- Supported grains: BATES only (raises error for Finocyl/Star/etc.)

## Result Dict Keys

```python
result['time']       # ndarray, simulation times [s]
result['P_head']     # ndarray, head-end pressure [Pa]
result['P_exit']     # ndarray, nozzle-end pressure [Pa]
result['D_throat']   # ndarray, throat diameter [m]
result['Kn']         # ndarray, burning area / throat area
result['massflow']   # ndarray, nozzle mass flow [kg/s]
result['snapshots']  # list of dicts (P, u, Mach, T, r_total, r_erosive,
                     #   D_port, x, C_burn, endface_msource,
                     #   is_burning, is_grain)
result['grains']     # list of per-segment dicts (regression, web)
result['summary']    # dict (P_peak, t_peak, mass_balance_error, etc.)
```

## Known Issues

1. Igniter spike decay rate — 2/3 power taper approximate
2. Transport property sensitivity — frozen vs effective k_gas, Cp_gas
3. Micro-steps in BATES from spatial discretization (converge with N)
4. Multi-tab burn rate not yet interpolated (selects one tab)
5. Burnout ramp extends burn time ~30% (asymptotic f_active tail)

## Validated Parameters (Hasegawa Motor A)

Propellant: a=4.821e-5, n=0.3, ρ=1700, T_flame=3041K, γ=1.19, MW=0.0254
Gas (effective): μ=8.842e-5, k=0.55, Cp=1900 | (frozen): k=0.3685, Cp=2060
Igniter: mass=10g, a=1e-4, n=0.4, ρ=1800, A_burn=auto
Motor A: D_bore=40mm, D_outer=80mm, L=1680mm, D_throat=34mm, N=150

## Roadmap

1. FMM grain support (per-slice burning perimeters from regression maps)
2. Multi-tab burn rate interpolation
3. Sensitivity study tooling
4. Tapered finocyl geometry generation
5. openMotor front-end integration
