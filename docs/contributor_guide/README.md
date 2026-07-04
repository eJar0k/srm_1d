# srm_1d — Contributor Guide

**Goal of this guide:** get a new contributor productive on the codebase
without having to reverse-engineer the sim core from scratch.

**Who it assumes you are:** you write good Python, and you know solid-rocket
internal ballistics (burn rate `r = a·Pⁿ`, Kn, c\*, nozzle expansion,
erosive burning, ignition transients). You do **not** need to be a CFD
person — you only need the high-level idea that the gas flow is solved
numerically. This guide fills the CFD gap and maps the code; it does **not**
re-teach SRM physics or Python.

> If you want the terse orientation the maintainers load every session, read
> the repo-root [`CLAUDE.md`](../../CLAUDE.md) first (what ships, quick-start
> commands, module one-liners). This guide is the *narrative* companion:
> CLAUDE.md tells you **what**, this tells you **how it works and why**.

---

## How this guide is organized

| # | Doc | What it covers | Status |
|---|-----|----------------|--------|
| — | **This README** | Orientation, the newcomer code map, the core-vs-opt-in feature catalog, and the reading order | ✅ |
| 01 | [`01_SIM_CORE.md`](01_SIM_CORE.md) | **The heart of the code, CFD-demystified**: finite volume, the staggered grid, the PISO pressure–velocity solve, sources, boundaries, time-stepping, and a guided tour of the main time loop | ✅ |
| 02 | [`02_BURN_AND_GEOMETRY.md`](02_BURN_AND_GEOMETRY.md) | The Ma-2020 erosive burn-rate chain (friction → heat transfer → bisection — a heat-transfer closure, **not** a flux power law) and how grain geometry regresses (analytic + FMM) each step | ✅ |
| 03 | [`03_IGNITION.md`](03_IGNITION.md) | The pyrogen igniter plenum + the Goodman solid-heating ignition kernel + the source-assembly step | ✅ |
| 04 | [`04_ADDINS.md`](04_ADDINS.md) | Deep-dive on each opt-in / experimental feature (Phase F/Z, radiation, topologies, diagnostic toggles) — the catalog below is the index | ✅ |
| 05 | [`05_IO_AND_OPENMOTOR.md`](05_IO_AND_OPENMOTOR.md) | `.ric` loading, the openMotor adapter, the plugin boundary, channels, plotting | ✅ |

**Recommended reading order for a newcomer:** this README → `01_SIM_CORE`
→ then whichever of 02/03 matches what you're touching. Keep
[`../../srm_1d/ARCHITECTURE.md`](../../srm_1d/ARCHITECTURE.md) (function-level
reference) and [`../../srm_1d/DEVNOTES.md`](../../srm_1d/DEVNOTES.md) (gotchas +
calibration state) open alongside.

---

## The 2-minute mental model

The simulator answers: *given a motor (grain geometry + propellant + nozzle +
igniter), what is the chamber-pressure and thrust history?*

It does that by treating the **bore/port as a 1-D duct** running head-end →
nozzle, chopping it into `N` axial cells, and marching a **transient
compressible gas solve** forward in time. Each timestep:

1. the **grain regresses** (burns back) a little, changing the port area;
2. the **burn rate** at each cell is computed (Saint-Robert + Ma erosive);
3. burning propellant and the igniter **inject mass, energy, and momentum**
   into the gas as *source terms*;
4. the **gas solver (PISO)** advances pressure / velocity / temperature /
   density one step, venting through the choked nozzle;
5. bookkeeping records the head pressure, mass flow, etc.

Loop until the propellant is gone (pressure drops below `P_cutoff`). The
whole loop is one Numba-JIT-compiled function for speed (~40–90 k steps/s).

**The single most important file to understand is
[`srm_1d/solver.py`](../../srm_1d/solver.py) (the gas solve) and the
`_run_time_loop` function in
[`srm_1d/simulation.py`](../../srm_1d/simulation.py) (the driver).**
`01_SIM_CORE` walks both.

---

## Newcomer code map

Modules, roughly in dependency order (leaves first). "Core" = you will read
it to understand a normal run; "peripheral" = I/O, tooling, or optional.

```
                       ┌─────────────────────────────────────────┐
 LEAVES (no project    │ solver.py       PISO gas solve (CFD)  ◄──┼── the core numerics
 dependencies — pure   │ burn_rate.py    Ma-2020 erosive rate     │
 numerics/data)        │ propellant.py   materials + gas thermo   │
                       │ grain_geometry.py  port regression       │
                       │ solid_thermal.py   Goodman ignition       │
                       └─────────────────────────────────────────┘
                       ┌─────────────────────────────────────────┐
 MID (import a leaf)   │ nozzle.py         thrust/Isp/CF, erosion │
                       │ igniter_plenum.py pyrogen chamber ODE     │
                       │ fmm_grain.py      openMotor FMM bridge    │
                       └─────────────────────────────────────────┘
                       ┌─────────────────────────────────────────┐
 DRIVER                │ simulation.py   run_simulation +          │  ◄── everything
                       │                 _run_time_loop (the loop) │      comes together
                       └─────────────────────────────────────────┘
                       ┌─────────────────────────────────────────┐
 EDGES (I/O, tooling)  │ openmotor_adapter.py  .ric ⇄ sim, CSV     │
                       │ plotting.py           matplotlib          │
                       │ srm1d_plugin.py       openMotor plugin     │
                       │ tools/                LHS sweeps, diags     │
                       └─────────────────────────────────────────┘
```

**Where a normal run flows** (the call graph you'll trace in `01_SIM_CORE`):

```
run_from_ric(<motor>.ric)                       # openmotor_adapter.py
  └─ load_ric / convert_{propellant,geometry,nozzle} / build_pyrogen_chamber
  └─ run_simulation(geo, propellant, nozzle, pyrogen_chamber, …)   # simulation.py
       ├─ (setup: build Numba arrays, initial conditions, history buffers)
       └─ _run_time_loop(… ~120 args …)          # the @njit while-loop
            └─ per step: geometry → burn → ignition/sources → PISO → post → record
                                             │
                                             └─ piso_step(…)        # solver.py
  └─ compute_motor_performance(result, nozzle, …)  # nozzle.py (thrust/Isp)
returns  (result, perf, nozzle, geo, prop)
```

The design decision to fold the entire loop into **one** `@njit` function
(rather than Python calling compiled kernels per step) is the main reason the
code reads as one very long function with many array arguments. It is a speed
tax on readability; `01_SIM_CORE` gives you the map so it's navigable.

---

## Core vs. opt-in: the feature catalog

A newcomer's biggest trap is mistaking an **experimental, default-off knob**
for load-bearing physics. The sim core (below) is always on. Everything in
the **opt-in** table is a research add-in, defaults to off/neutral, and can
be ignored while you learn the core.

### Always-on core physics
| Piece | Where | One-liner |
|---|---|---|
| PISO gas solve | `solver.py:piso_step` | Pressure–velocity–energy update per step |
| Saint-Robert + **Ma-2020 erosive** burn rate | `burn_rate.py:compute_burn_rates` | `r = a·Pⁿ` + erosive increment from local Re/roughness |
| Grain regression | `grain_geometry.py:update_cell_geometry` | Port area/perimeter grow as the grain burns back |
| **N-species mixture** | `simulation.py:_advect_species` / `_refresh_mixture_arrays` | Bore gas is a mix (igniter + grain products); γ/R/Cₚ are **per-cell** |
| Pyrogen igniter | `igniter_plenum.py` + source assembly | Hot-gas plenum vents into cell 0 to start the motor |
| Goodman ignition | `solid_thermal.py` | Per-cell solid heating; a cell "lights" when `T_surf > T_ignition` |
| Choked-nozzle BC | `solver.py:_nozzle_boundary_flow` | Signed isentropic throat outflow (with subsonic fallback) |
| Adaptive CFL | `solver.py:compute_dt_cfl` + `compute_dt_source_cap` | `dt` from wave speed + a source-rate cap |

### Opt-in / experimental add-ins (default OFF or neutral)
| Feature | Flag(s) & default | Gated in | What it does | Docs |
|---|---|---|---|---|
| **Phase F — flame-front** | `flame_front_enabled=False`, `flame_front_velocity` | `_advance_flame_front`, source assembly | Ignition spreads as a front at a finite speed instead of lighting the whole grain at once | `docs/v0_7_4/` |
| **Phase Z — Z-N dynamic burn** | `zn_enabled=False`, `kappa_zn` | `_advance_zn_burn_rate` (STEP 2) | Burn rate *lags* the quasi-steady Ma target with `τ = κ·α_s/r²` (adds burn-rate memory) | `docs/v0_7_4/` |
| Burn establishment lag | `tau_establishment=0.0` | source assembly | Post-ignition, ramps a cell's output 0→1 over `τ` | signature docstring |
| **`port_mach_cap`** | `port_mach_cap=0.0` (off) | `solver.py` PISO (Step 3c) | Clips interior face velocity to ~Mach·a — bounds the supersonic-fill artifact | `docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` |
| Adjacent-cell radiation | `propellant.radiation_emissivity=0.0` | Goodman source assembly | Burning cells radiate to unignited neighbors (ignition assist) | `docs/v0_7_0/` |
| Igniter topology | `injection_topology='forward_plenum'` (alt: `head_basket`/`aft_basket`) | topology_code branch | Where/how pyrogen enters (plenum-with-orifice vs uncontained pellets) | `docs/v0_7_3/` |
| Throat erosion | `nozzle.erosion_coeff` / `slag_coeff` | STEP 3b | Throat diameter evolves with the pressure history | `nozzle.py` |
| 6× `diagnostic_disable_*` | all `False` | throughout | Turn OFF a physics term (erosive, endfaces, momentum, pyrogen surface heat, adjacent radiation, radiation gas sink) to isolate its effect | signature docstring |

**Rule of thumb:** if a term is behind a `diagnostic_disable_*` /
`*_enabled` / a nonzero-to-activate coefficient, it is **not** part of the
default validated physics — it's a probe or a research feature. The
canonical validated run (Hasegawa A) uses the core table only, at the
v0.7.5 knob defaults (roughness 32 µm, kappa 0.44, T_ignition 756 K,
k_solid 0.271).

---

## Gotchas you will hit (read before editing the core)

Full list in [`../../srm_1d/DEVNOTES.md`](../../srm_1d/DEVNOTES.md). The three
that bite newcomers immediately:

1. **Numba cache staleness** — after editing any `@njit` function you MUST
   delete `srm_1d/__pycache__/` (and `.nbi`/`.nbc`) or you'll silently run the
   *old* compiled code. This is the #1 "my fix didn't work" cause.
2. **Mass conservation coupling** — the burnout ramp `f_active` must multiply
   **both** `C_burn` and the regression rate; either alone causes 3–40 % mass
   error. Don't "simplify" one without the other.
3. **The Python env** — use the pyenv 3.10.5 interpreter explicitly (it has
   numba/scipy/scikit-fmm); the system Python does not. See CLAUDE.md
   quick-start.

---

## Conventions this codebase follows (so your PR fits)

- **Defer to openMotor's data model** for new data structures (field names,
  semantics). Units are the documented exception — srm_1d keeps
  human-readable engineering units internally and converts at the adapter
  boundary.
- **`.ric` files are openMotor's format** — regenerate them from openMotor,
  never hand-edit.
- **Named motors are data**, not Python factories — add a `<motor>.ric` under
  `motors/` and load via `run_from_ric`.
- **Versioning is git tags**, and hard API breaks are acceptable (documented
  in DEVNOTES' breaking-change log) — no backward-compat shims.
- **No physics tuning without literature backing** — calibration knobs
  (roughness, kappa, k_solid) have published physical bounds; numerical knobs
  (CFL) are separate and freely adjustable.
