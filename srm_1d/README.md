# 1D SRM Internal Ballistics Simulator (srm_1d v0.7.0)

A transient 1D finite-volume solver for solid rocket motor internal
ballistics with the Ma et al. (2020) erosive burning model.

**Performance:** ~45–90k steps/s (compiled time loop with Numba JIT;
FMM grains run faster than cylindrical due to fewer per-step ops)
**Tests:** pytest | **Adapter:** reads openMotor .ric files
including all FMM grain types (Finocyl, Star, Moon, X, C, D, Custom)

## Quick start

```python
from srm_1d.openmotor_adapter import run_from_ric

result, perf, nozzle, geo, prop = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    roughness=37.1e-6,
    pyrogen='bpnv',
    T_ignition=850.0,
    P_cutoff=0.05e6,
)
```

`run_from_ric` auto-discovers a sibling `<motor>.transport.yaml`
alongside the .ric file and uses it for combustion gas transport
properties (mu, k, Cp). It also accepts `pyrogen='bpnv'`/`'mtv'` or
a sibling `<motor>.pyrogen.yaml` for v0.7.0 pyrogen ignition. For
parametric geometry construction without a .ric, use
`srm_1d.grain_geometry.build_snapped_geometry` directly.

## File Structure

```
srm_1d/
├── __init__.py              # Public API, v0.7.0
├── solver.py                # PISO: TDMA, pressure correction, CFL
├── burn_rate.py             # Ma et al.: Haaland → Gnielinski → bisection
├── grain_geometry.py        # GrainSegment, MotorGeometry, build_snapped_geometry
├── propellant.py            # Propellant (multi-tab) + gas properties + thermo utilities
├── nozzle.py                # openMotor-aligned Nozzle: thrust, Isp, CF, throat erosion
├── fmm_grain.py             # openMotor FmmGrain bridge — extracts FmmTable from any FMM grain
├── igniter_plenum.py        # Pyrogen chamber and choked/subsonic vent model
├── solid_thermal.py         # Goodman integral solid-heating ignition model
├── simulation.py            # Compiled time loop (_run_time_loop @njit)
├── plotting.py              # Pressure/thrust/snapshot/comparison plots
├── openmotor_adapter.py     # .ric reader, transport YAML loader, CSV export
├── motors/                  # Canonical motor + pyrogen data (v0.6.0+)
│   ├── hasegawa_a.ric, hasegawa_a.transport.yaml
│   ├── hasegawa_b.ric, hasegawa_b.transport.yaml
│   ├── hasegawa_c.ric, hasegawa_c.transport.yaml
│   ├── example_bates.ric, example_bates.transport.yaml
│   └── pyrogens/bpnv.yaml, mtv.yaml
├── tools/
│   └── sensitivity.py       # Latin Hypercube sweeps with parallel execution
├── tests/  (11 files, 133 tests)
└── examples/ (hasegawa_motor_a.py, bates_4seg.py, hasegawa_a_lhs.py)
```

## Architecture

### Compiled Time Loop
```
run_simulation(geo, propellant, nozzle, ...)   Python: setup, allocate
  → _run_time_loop(...)                        @njit: entire while loop
      ├─ advance_endface_regression   Numba→Numba direct calls
      ├─ advance_bore_regression      (mass-conserving burnout ramp)
      ├─ update_cell_geometry         (volumetric overlap + hat-function end-face kernel)
      ├─ compute_burn_rates           (every N steps)
      ├─ _goodman_ignition_sources_and_mass
      │   (Goodman ignition + mass/thermal source assembly)
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
- `build_snapped_geometry(segments_spec, D_outer, target_propellant_cells)`
  is the canonical builder. It computes a `dx` from
  `L_propellant / target_propellant_cells`, applies a Nyquist-CFD
  clamp so the smallest gap is ≥1 cell, then **integer-snaps** every
  segment length and inter-segment gap to multiples of `dx`. The
  result: cell boundaries align with segment edges by construction
  (eliminates fractional-cell artifacts at boundaries).
- Primary per-cell state: `regress[i]` — radial regression depth (m).
  Unifies cylindrical, conical, and FMM paths.
- `cell_segment_type[i]` — `0` for analytic, `1` for FMM lookup.
- End-face mass injection uses a **partition-of-unity hat function**:
  each face's mass splits over 2 adjacent cells with weights summing
  to 1.0. Coupled to snapping (which puts faces on cell edges); gated
  by `tests/test_endface_conservation.py` (<0.1% error vs analytic).
- Auto-inhibition: touching segments get shared faces inhibited.

### FMM Grain Support
- All 7 of openMotor's FMM grain types supported: Finocyl, Star Grain,
  Moon Burner, X Core, C Grain, D Grain, Custom Grain.
- `srm_1d.fmm_grain.from_openmotor(om_grain)` runs openMotor's
  regression-map pipeline on a populated FmmGrain instance, samples
  perimeter and port area at fine resolution, and returns an
  `FmmTable` ready for the @njit hot loop.
- Tables are packed CSR-style (`fmm_offset[]`, `fmm_reg_flat[]`,
  `fmm_perim_flat[]`, `fmm_port_flat[]`) per-grain; per-cell
  `cell_fmm_idx` pinpoints which segment each cell uses.
- Mixed-grain motors (BATES + Finocyl in the same MotorGeometry) work
  out of the box — `cell_segment_type[i]` switches paths per cell.
- Hydraulic diameter for FMM cells: `D_h = 4·A_port/perimeter`.
- Optional dependencies: `scikit-fmm` (regression maps), `scikit-image`
  (openMotor Custom grain). The Cython `mathlib._find_perimeter_cy`
  is replaced with a Numba-JIT marching-squares perimeter, so MSVC
  build tools aren't required.

### Igniter Model (v0.7.0)
Hot-gas pyrogen ignition uses a forward 0D plenum (`PyrogenChamber`)
venting through a choked/subsonic orifice into cell 0. The main solver
tracks pyrogen mass flow, plenum pressure/temperature, and per-cell
Goodman solid heating. Grain cells ignite when `T_surf > T_ignition`
(default 850 K). The old exponential igniter API is removed.

### Throat Evolution (in-loop, not post-processing)
Erosion: rate = erosion_coeff [μm/(s·MPa)] × P [MPa] × 1e-6
Slag: rate = slag_coeff [(m·MPa)/s] / P [MPa]
dD/dt = 2 × (erosion - slag). Feeds back into nozzle BC each step.
Coefficients live on the `Nozzle` object passed to `run_simulation`.

### Nozzle (openMotor-aligned, v0.3.0+)
- Field names mirror openMotor's `motorlib.Nozzle`: `D_throat`, `D_exit`,
  `efficiency`, `div_angle`, `conv_angle`, `throat_length`,
  `erosion_coeff`, `slag_coeff` (snake_case'd; units kept human-readable
  internally — adapter handles conversion).
- Thrust uses openMotor's adjusted-CF formula:
    CF_adj = divLoss × throatLoss × efficiency
           × (skinLoss × CF_ideal + (1 − skinLoss))

### Propellant (multi-tab, openMotor-aligned, v0.4.0+)
- `Propellant` holds `tabs: list[PropellantTab]`. Each tab carries
  `(min_pressure, max_pressure, a, n, gamma, T_flame, molecular_weight)` —
  a 1:1 snake_case mirror of openMotor's `motorlib.PropellantTab`.
- Gas transport (`mu_gas`, `k_gas`, `Cp_gas`) lives at the propellant
  level (not in tabs); supplied via sibling `<motor>.transport.yaml`
  or `gas_props={...}`.

### openMotor Adapter (v0.7.0)
- `load_ric()` → parse YAML with Python-tag handling
- `load_transport()` → parse sibling `<motor>.transport.yaml`
- `load_pyrogen()` → built-in (`'bpnv'`, `'mtv'`) or YAML pyrogen loader
- `convert_propellant()` → Propellant with all .ric tabs preserved 1:1
- `convert_geometry(grains, target_propellant_cells)` → MotorGeometry
  via `build_snapped_geometry`. Auto-applies inter-segment gap of
  `max(3mm, 5%·D_outer)` between segments.
- `convert_nozzle()` → Nozzle (erosionCoeff m/(s·Pa) → μm/(s·MPa))
- `run_from_ric()` → result, perf, nozzle, geo, prop. Auto-resolves
  sibling transport YAML and requires explicit `pyrogen=...` or a
  sibling `<motor>.pyrogen.yaml`. `verbose=False` suppresses setup and
  performance summary blocks for sweeps.
- `save_csv()` → Time, Kn, P, Force, Mass Flow, per-grain regression/web
- Supported grains: BATES, Conical (analytic) + all 7 FMM types.

### Sensitivity Tooling (v0.7.0)
`srm_1d.tools.sensitivity.run_lhs(motor_path, bounds, n_samples,
fitness_fn, **sim_kwargs)` runs an N-sample Latin Hypercube sweep with
`scipy.stats.qmc.LatinHypercube` and `concurrent.futures.ProcessPoolExecutor`.
Pluggable fitness factories: `mse_fitness`,
`impulse_error_fitness`, `peak_pressure_error_fitness`, and
`segmented_pressure_fitness`. Optional `metrics_fn` columns are persisted
to CSV. `progress_mode='brief'|'verbose'|'none'` controls terminal noise,
and `sim_verbose=False` is the sweep default. Example:
`srm_1d/examples/hasegawa_a_lhs.py`.

## Result Dict Keys

```python
result['time']       # ndarray, simulation times [s]
result['P_head']     # ndarray, head-end pressure [Pa]
result['P_exit']     # ndarray, nozzle-end pressure [Pa]
result['D_throat']   # ndarray, throat diameter [m]
result['Kn']         # ndarray, burning area / throat area
result['massflow']   # ndarray, nozzle mass flow [kg/s]
result['P_ig']       # ndarray, pyrogen plenum pressure [Pa]
result['T_ig']       # ndarray, pyrogen plenum temperature [K]
result['mdot_ig']    # ndarray, igniter mass flow into cell 0 [kg/s]
result['m_pyrogen']  # ndarray, remaining solid pyrogen mass [kg]
result['snapshots']  # list of dicts (P, u, Mach, T, r_total, r_erosive,
                     #   D_port, x, C_burn, endface_msource,
                     #   is_burning, is_grain, T_surf)
result['grains']     # list of per-segment dicts (regression, web)
result['summary']    # dict (P_peak, t_peak, mass_balance_error, etc.)
```

## Known Issues

1. Transport property sensitivity — frozen vs effective k_gas, Cp_gas
2. Multi-tab burn rate not yet interpolated (selects one tab)
3. Burnout ramp extends burn time ~30% (asymptotic f_active tail)
4. The `build_snapped_geometry` ±10mm warning fires on coarse grids
   (target_propellant_cells ≤ 100 with leading_gap=1mm); cosmetic only
5. Phase 3 Goodman ignition switches cells immediately to full burning.
   Hasegawa segmented diagnostics indicate the next structural fix should
   be a finite post-ignition burn-establishment/participation model.

## Calibration State (Hasegawa Motor A)

Motor: D_bore=40mm, D_outer=80mm, L=1680mm, D_throat=34mm
Propellant: a=4.821e-5, n=0.3, ρ=1700, T_flame=3041K, γ=1.19, MW=0.0254
Transport (RPA effective): μ=8.842e-5, k=0.3685, Cp=2060
Roughness: 37.1 μm remains the inherited v0.6.0 erosive baseline.
Igniter: v0.7.0 uses BPNV pyrogen plenum + Goodman `T_surf` ignition;
the old `igniter_tau`/`P_ignition` knobs are removed. Phase 4 segmented
LHS diagnostics are in `artifacts/hasegawa_a_lhs/` locally and point to
post-ignition burn establishment as the next modeling target.

## Roadmap

1. Per-step gas thermo (γ, T_flame, MW) for multi-tab propellants
2. Post-ignition burn participation/establishment after Goodman trigger
3. Squib stage before pyrogen ignition
4. RodTube grain support (PerforatedGrain extension to from_openmotor)
5. openMotor front-end integration (CFDSimulation subclass — deferred)
