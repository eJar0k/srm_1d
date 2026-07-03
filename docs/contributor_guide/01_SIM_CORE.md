# 01 — The Simulation Core (CFD-demystified)

This is the document that makes the gas solver comprehensible. It assumes
you know SRM internal ballistics and Python, and that CFD is a black box to
you. By the end you should be able to open
[`solver.py`](../../srm_1d/solver.py) and the `_run_time_loop` in
[`simulation.py`](../../srm_1d/simulation.py) and know what every stage is
doing and why.

---

## 1. The physical model — a duct that makes its own gas

You already know the 0-D lumped-parameter picture: propellant burns, produces
gas at rate `ṁ_gen = ρ_p·A_burn·r`, chamber pressure rises until generation
balances nozzle outflow `ṁ_noz = P_c·A_t/c\*`, and Kn sets the equilibrium
pressure. That 0-D model has **no spatial dimension** — one pressure for the
whole chamber.

srm_1d keeps the axial dimension. The bore/port is modeled as a **quasi-1-D
duct** from the head end (`x=0`) to the nozzle (`x=L`). Gas can have different
pressure, temperature, density, and **velocity** at different axial stations.
That matters for exactly the phenomena this project cares about:

- **Erosive burning** — near the nozzle the fast axial gas flow enhances heat
  transfer to the burning surface and lifts the burn rate above the static
  `a·Pⁿ` value; the augmentation grows from ~zero at the head end toward the
  aft. srm_1d uses the **Ma (2020)** model, which captures this through the
  local *turbulent convective heat transfer* to the surface (a Reynolds-number
  / Gnielinski–Haaland closure, solved by bisection). Importantly it is **not**
  an explicit mass-flux (`G = ρu`) power law like the older, purely flux-based
  correlations (Lenoir–Robert and successors) — a distinction that matters,
  because `G`/Re appear only *inside* a heat-transfer closure, not as a direct
  `r ∝ Gⁿ` term. You cannot see any of this in 0-D. (Details in
  `02_BURN_AND_GEOMETRY.md`.)
- **The ignition transient** — the flame and pressure wave sweep down the
  bore over milliseconds; the head and aft ends light and pressurize at
  different times.

So the core solves the **1-D compressible Euler equations with area variation
and distributed source terms** (mass, momentum, and energy added along the
duct by combustion and the igniter). "Compressible" because Mach numbers reach
~0.3–0.4 in the plateau (and transiently higher), so density is not constant.

The three conservation laws, in words:

| Law | Plain statement | Sources in an SRM |
|---|---|---|
| **Mass** | gas in a cell changes by what flows across its faces **plus** what the propellant/igniter inject | burning surface adds `ṁ` |
| **Momentum** | gas accelerates due to pressure gradients, convection, and wall friction | igniter jet adds axial momentum |
| **Energy** | gas enthalpy changes by advection across faces **plus** injected combustion enthalpy | hot combustion products add `ṁ·Cp·T_flame` |

Plus the **ideal-gas equation of state** `P = ρRT` tying pressure, density,
and temperature together.

---

## 2. Discretization — finite volume on a staggered grid

### Finite volume, `N` cells
The duct is chopped into `N` equal cells of width `dx` (set from
`target_propellant_cells`; `N ≈ 50–200` typically). "Finite volume" means we
track the **cell-averaged** conserved quantities and update them by summing
**fluxes across the cell faces** — this is what makes the scheme
conservative (nothing is created or lost except by the explicit sources and
the nozzle). This is the natural framework for "how much mass/energy is in
this chunk of duct, and what crosses its boundaries."

### The staggered grid — the one CFD idea you must internalize
Two sets of locations:

```
      cell 0        cell 1        cell 2            cell N-1
   ┌──────────┬──────────┬──────────┬─  ...  ─┬──────────┐
   P,T,ρ,A    P,T,ρ,A    P,T,ρ,A              P,T,ρ,A       ← CELL CENTERS (scalars)
   u0        u1        u2        u3           u(N-1)  uN     ← FACES (velocity)
   ▲                                                  ▲
 head wall (u0=0)                              nozzle (uN)
```

- **Scalars** (pressure `P`, temperature `T`, density `rho`, port area
  `A_port`) live at **cell centers**, indices `0 … N-1`.
- **Velocity** `u` lives at **faces**, indices `0 … N`. Face `j` sits
  *between* cells `j-1` and `j`. Face `0` is the head-end wall
  (`u[0]=0`, no flow through a closed head); face `N` is the nozzle exit.

**Why stagger instead of putting everything at cell centers?** If pressure and
velocity shared the same location, the discrete pressure gradient at a point
would depend on its *non-adjacent* neighbors (`P[i+1]-P[i-1]`), skipping the
cell itself. That lets a checkerboard pressure pattern (high-low-high-low)
look "smooth" to the solver and produces spurious oscillations. Staggering
puts the face velocity `u[j]` directly between `P[j-1]` and `P[j]`, so the
pressure gradient driving it is `(P[j]-P[j-1])/dx` — adjacent cells only, no
skipping ([`solver.py:390`](../../srm_1d/solver.py#L390)). This tight coupling
is why the scheme is stable without artificial smoothing.

---

## 3. Why a "PISO" solver — the pressure–velocity coupling problem

Here is the crux that all the machinery exists to solve.

In compressible flow, pressure, velocity, and density are coupled: velocity is
driven by the pressure gradient (momentum), but the velocity field in turn
determines how mass piles up and therefore what the pressure/density *become*
(continuity + EOS). You cannot update velocity without knowing the new
pressure, and you cannot know the new pressure without knowing the velocity.
Acoustic waves make this coupling near-instantaneous (they carry pressure
information across a cell in `dx/a` seconds).

A naive fully-explicit scheme (use old pressure to push velocity, then update
density, then get new pressure from EOS) works but is **fragile and requires
tiny time steps**, and doesn't enforce the continuity constraint cleanly.

**PISO** (Pressure-Implicit with Splitting of Operators, Issa 1986) resolves
the coupling by *operator splitting*:

1. **Predict** the velocity using the *old* pressure gradient — call it `u*`.
   It won't satisfy continuity yet.
2. **Correct**: solve a pressure-correction equation `P'` such that when you
   nudge the velocities by `u' ∝ -∂P'/∂x`, mass is conserved in every cell.
   This is an implicit (whole-domain) solve — a tridiagonal system, because
   each cell's `P'` couples only to its two neighbors.
3. **Correct again** (2nd corrector — the "I" in PISO) to tighten the
   split-operator error.

The pressure-correction equation is the discrete analog of a Poisson/Helmholtz
equation. On our 1-D grid it's **tridiagonal**, so it's solved directly with
the **Thomas algorithm** (`thomas_solve`, O(N)) rather than iteratively — fast
and exact.

> **This is the seam where the acoustic-CFL performance work lives.** The
> pressure solve is implicit, but because it's *non-iterated* (2 correctors,
> old temperature in the coefficients), the scheme is only stable up to an
> **acoustic** Courant number ≈ 0.5. See
> [`../core_loop_opt/DESIGN_LEVER_B.md`](../core_loop_opt/DESIGN_LEVER_B.md) if
> you're touching the time-step limit.

---

## 4. The PISO step, annotated

All in [`solver.py:_piso_step_with_energy_diagnostics`](../../srm_1d/solver.py#L276)
(the public wrapper `piso_step` just drops the energy-diagnostic returns).
Inputs: current `rho, u, P, T`, geometry (`A_port, D_hyd`), the three source
arrays, friction `f_darcy`, `dt`, and **per-cell** gas thermo arrays
(`gamma_arr, R_arr, Cp_arr` — see §7). Stages:

**Step 1 — Momentum predictor** ([L368](../../srm_1d/solver.py#L368)), faces
`j=1…N-1`:
```
u*[j] = u[j] + (dt/ρ_face) · ( −∂P/∂x        # pressure gradient, OLD P
                               + convection    # upwind ρu·u flux
                               + friction       # −f/(2D)·ρ|u|u  (Darcy)
                               + momentum_source[j] )  # igniter jet
```
`u*[0]=0` is the closed head wall. The nozzle face is extrapolated for now and
fixed up by the pressure step.

**Step 2 — Pressure correction #1** ([L430](../../srm_1d/solver.py#L430)).
Assemble a tridiagonal system for `P'` from discrete continuity in each cell:
```
(ρ_new − ρ_old)·A·dx/dt  +  (ṁ*_east − ṁ*_west)  =  mass_source·dx
```
with `ρ_new = (P+P')/(R·T)` (the transient-density term — the diagonal
coefficient `a_t = A·dx/(R·T·dt)`, [L454](../../srm_1d/solver.py#L454)) and the
velocity correction `u'[j] = −d_face·(P'[j]−P'[j-1])/dx`. The **nozzle
boundary** contributes its flow sensitivity `∂ṁ/∂P` to the last cell's diagonal
([L473](../../srm_1d/solver.py#L473)). Solve with `thomas_solve`, then update
`P += P'` and correct the face velocities.

**Step 3 — Pressure correction #2** ([L507](../../srm_1d/solver.py#L507)).
Recompute density with the updated pressure and repeat the correction — this
is the second PISO corrector that reduces the operator-splitting error. Then
set the nozzle face velocity from the final boundary mass flow
([L564](../../srm_1d/solver.py#L564)).

**Step 3c — `port_mach_cap`** ([L574](../../srm_1d/solver.py#L574), *opt-in*,
default off). Clips interior face velocities to `≤ Mach·a` to bound the
grid-divergent supersonic-fill artifact. No-op in the plateau. See the
ignition-spike docs.

**Step 3b — Energy** ([L602](../../srm_1d/solver.py#L602)). Advect **sensible
enthalpy** `h = Cp·T` (not `T` directly — across a face between cells with
different `Cp`, `h` is the conserved scalar) using upwind face fluxes, in a
mass-conservative form:
```
m_new·h_new = m_old·h_old + dt·(h_flux_west − h_flux_east + thermal_source·dx)
T_new = h_new / Cp,  then clipped to [T_floor, T_ceiling]
```
The `T_ceiling` is a per-cell physical cap (a cell can't exceed the flame
temperature of the gas actually in it). The clip's energy correction is
tracked in the diagnostics for conservation auditing.

**Step 4 — Density / EOS** ([L684](../../srm_1d/solver.py#L684)):
`ρ_new = P_new/(R·T_new)`, with a pressure floor.

Returns updated `(rho, u, P, T)` (+ a bundle of energy-balance diagnostics the
driver records).

---

## 5. Source terms — how the SRM physics enters the gas solve

The PISO step is pure gas dynamics; it knows nothing about propellant. All the
SRM physics reaches it through **three per-cell source arrays** the driver
assembles each step:

| Array | Units | Filled by | Meaning |
|---|---|---|---|
| `mass_source[i]` | kg/(m·s) | burning grain + igniter | gas mass injected per unit length |
| `thermal_source[i]` | W/m | same | enthalpy injected per unit length (`ṁ·Cp·T_source`) |
| `momentum_source[j]` | N/m³ | igniter orifice jet | axial momentum injected at faces |

The grain contribution is `mass_source[i] = ρ_p · C_burn[i] · r_total[i]`
(propellant density × effective burning perimeter × burn rate → kg/(m·s)),
assembled in the ignition/source step; the igniter contribution comes from the
pyrogen plenum (`igniter_plenum.py`) venting into cell 0 (or, for the
uncontained topologies, into the cartridge cells directly). A cell only
contributes once it has **ignited** (`T_surf > T_ignition`, tracked by the
Goodman kernel). This is covered in `03_IGNITION.md`; for the core, the key
idea is: **combustion is a boundary condition on the gas solve, delivered as
source terms.**

---

## 6. Boundary conditions

- **Head end** (face 0): a solid wall, `u[0] = 0`. No flow through the closed
  forward closure.
- **Nozzle** (face N):
  [`_nozzle_boundary_flow`](../../srm_1d/solver.py#L231) computes a **signed
  isentropic throat mass flow** from the last cell's `P, T` and the throat
  area — choked when the pressure ratio is above critical (the normal firing
  case), with a subsonic-outflow and even reverse-inflow branch for the
  low-pressure fill/tail. Its pressure derivative `∂ṁ/∂P` is what couples the
  nozzle into the pressure-correction matrix, so the chamber "feels" the
  nozzle implicitly rather than as a lagged outflow.

This is the spatially-resolved version of the `ṁ_noz = P_c·A_t/c\*` you know
from 0-D.

---

## 7. The N-species gas mixture (per-cell thermo)

A subtlety that surprises newcomers: the bore gas is **not** a single gas. It's
a mixture of **igniter products** (species 0) and **main grain combustion
products** (species 1), and their proportions vary by cell and time (the
igniter dominates cell 0 at startup; grain products dominate later). Because
`γ`, `R`, and `Cp` differ between the two, the solver carries **per-cell
arrays** `gamma_arr / R_arr / Cp_arr`, not scalars.

- [`_advect_species`](../../srm_1d/simulation.py#L939) transports the species
  mass fractions `Y` with the flow each step (after PISO updates `rho`).
- [`_refresh_mixture_arrays`](../../srm_1d/simulation.py#L326) recomputes the
  per-cell `γ/R/Cp` from the updated `Y` for the next step.

The PISO step and the nozzle BC consume these per-cell arrays (the nozzle uses
cell `N-1`'s mixture). If you see `_arr` suffixes on gas properties, this is
why. (Full design: `docs/v0_7_1/`.)

---

## 8. Time stepping — the adaptive `dt`

The step is chosen each iteration ([`simulation.py:1617`](../../srm_1d/simulation.py#L1617))
as the **minimum** of:

1. **Wave-speed CFL** ([`compute_dt_cfl`](../../srm_1d/solver.py#L741)):
   `dt = cfl_target · dx / (|u|_max + a_max)`, where `a=√(γRT)` is the sound
   speed. This keeps information from crossing more than a fraction
   `cfl_target` (default 0.3) of a cell per step. Because it includes `a`, this
   is an **acoustic** CFL — the dominant cost, and the target of the Lever-B
   perf work.
2. **Source-rate cap** ([`compute_dt_source_cap`](../../srm_1d/solver.py#L785)):
   limits `dt` so a violent per-cell energy injection can't change a cell's
   temperature by more than `source_cfl_factor` of the flame–ambient range in
   one step. Needed to survive the ignition cascade.
3. **`dt_max`** — a hard ceiling for quiet phases.

Understanding this is background for the profiling story in
[`../core_loop_opt/README.md`](../core_loop_opt/README.md): ~99.9 % of steps
are the acoustic-CFL-limited plateau.

---

## 9. Guided tour of the main loop

[`_run_time_loop`](../../srm_1d/simulation.py#L1316) is one `@njit` function.
Ignore its ~120 arguments (they're the price of keeping the whole loop
compiled); the **body is a clean 6-step sequence** you can read top to bottom:

| Step | Code | What happens | Deeper doc |
|---|---|---|---|
| **1. Geometry** | [L1621](../../srm_1d/simulation.py#L1621) | Advance regression (bore + end-faces); every `burn_update_interval` steps recompute per-cell `A_port`, perimeter `C_burn`, `D_hyd` via `update_cell_geometry` | 02 |
| **2. Burn rates** | [L1657](../../srm_1d/simulation.py#L1657) | `compute_burn_rates` (Saint-Robert + Ma erosive) per ignited cell; **Phase Z** relaxation if `zn_enabled` | 02 |
| **3. Ignition + sources** | [L1708](../../srm_1d/simulation.py#L1708) | Step the pyrogen plenum; advance Goodman solid heating; light cells past `T_ignition`; assemble `mass/thermal/momentum` source arrays (branching on `topology_code`); **Phase F** front gate if enabled | 03 |
| **3b. Throat** | [L1965](../../srm_1d/simulation.py#L1965) | Evolve `D_throat` from erosion/slag if `throat_is_evolving` | — |
| **4. PISO** | [L1985](../../srm_1d/simulation.py#L1985) | `_piso_step_with_energy_diagnostics` — the §4 gas update | this doc |
| **5. Post-PISO** | [L2003](../../srm_1d/simulation.py#L2003) | `_post_piso_update`: recompute face→cell velocities, Re, Mach, friction, max sound speed (feeds next step's burn rate + `dt`) | this doc |
| **6. Bookkeeping** | [L2008](../../srm_1d/simulation.py#L2008) | `_advect_species` + `_refresh_mixture_arrays`; nozzle mass-flow tally; record histories; snapshots at `snapshot_interval` | this doc |

Loop exits on `t ≥ t_max`, `P_head < P_cutoff` after ignition, burnout, or a
numerical-health abort.

---

## 10. How to trace and experiment

- **Entry point:** run a motor with
  `python -m examples.hasegawa_motor_a` (the canonical validated case). It
  calls `run_from_ric` → `run_simulation` → `_run_time_loop`.
- **The result dict** (returned by `run_simulation`) is plain data — keys like
  `time`, `P_head`, `dt`, `n_ignited`, `max_pressure`, `snapshots`,
  `summary`, plus the energy/momentum audit histories. Everything downstream
  (plots, CSV, performance) reads this dict, so it's the clean inspection
  surface.
- **You can't breakpoint inside `@njit`.** To debug the loop, either (a) read
  the recorded per-step histories, (b) temporarily comment `@njit` off a
  function to run it in pure Python (slow but debuggable — remember the cache
  gotcha when you put it back), or (c) add a diagnostic history array.
- **Isolate a physics term** with the `diagnostic_disable_*` flags (§ feature
  catalog in the README) — e.g. `diagnostic_disable_erosive=True` gives you
  the Saint-Robert-only baseline to compare against.
- **After any `@njit` edit:** delete `srm_1d/__pycache__/` and rerun
  `python -m pytest tests/` (≈406 tests).

---

## 11. Where to go next

- **The burn rate & geometry** (Step 1–2): `02_BURN_AND_GEOMETRY.md` (planned)
  — the Ma-2020 chain and FMM grain regression.
- **Ignition & sources** (Step 3): `03_IGNITION.md` (planned) — plenum +
  Goodman + source assembly.
- **Function-level reference:** [`../../srm_1d/ARCHITECTURE.md`](../../srm_1d/ARCHITECTURE.md).
- **Gotchas, calibration state, API-break log:**
  [`../../srm_1d/DEVNOTES.md`](../../srm_1d/DEVNOTES.md).
- **The design docs** under `docs/v0_7_x/`, `docs/v0_8_0/`, and
  `docs/core_loop_opt/` record *why* each subsystem is the way it is — read the
  one matching the area you're changing.
