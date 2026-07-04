# 02 — Burn Rate & Geometry (loop Steps 1–2)

Covers how the code turns "the grain is burning" into the `mass_source` the
gas solver consumes: the **Ma-2020 erosive burn rate**
([`burn_rate.py`](../../srm_1d/burn_rate.py)) and the **port regression**
([`grain_geometry.py`](../../srm_1d/grain_geometry.py)). You know SRM burn-rate
physics; this doc is about *how srm_1d represents it spatially* and — most
importantly — *what the Ma erosive model actually is* (a common point of
confusion).

---

## Part A — Burn rate

Each cell's total burn rate is
```
r_total = r₀ + r_e
```
`r₀` is the normal (Saint-Robert) rate; `r_e` is the erosive increment. These
are computed per ignited cell in
[`compute_burn_rates`](../../srm_1d/burn_rate.py#L512) → `burn_rate_cell`,
every `burn_update_interval` steps.

### A.1 — Normal rate `r₀`: multi-tab Saint-Robert
`r₀ = a(P)·P^n(P)` — the burn law you know, but with **multiple tabs** (a/n
pairs valid over different pressure ranges), matching openMotor's
`getCombustionProperties`. Tab selection
([`select_tab_idx`](../../srm_1d/burn_rate.py#L308)) is a **hard switchover**:
strict containment (`min < P < max`) first, else the closest-boundary tab.
No interpolation between tabs — this is deliberate, to match openMotor exactly.

### A.2 — Erosive rate `r_e`: the Ma-2020 heat-transfer model

> **Read this if you take nothing else from the doc.** The Ma (2020) model is
> **not** an explicit mass-flux correlation. It is a **convective-heat-transfer
> closure**. The older erosive models (Lenoir–Robert 1954 and successors)
> write the erosive rate as an explicit power law in the mass flux,
> `r_e ∝ G^0.8·…` with fitted constants. **srm_1d does not use that form.** In
> Ma-2020, the crossflow enters *only* through the Reynolds number that drives
> a named heat-transfer correlation — there is no `r_e ∝ Gⁿ` term and no fitted
> erosive constant. Keep this straight when reading or explaining the code.

`r_e` comes from a **surface energy balance** (Ma Eqs. 2–7): the extra heat the
crossflow delivers to the surface drives extra pyrolysis.
```
r_e = (T_flame − T_surface) / (ρ_p · C_ps · (T_surface − T_initial)) · h
```
where `h` is the **convective heat-transfer coefficient**. Everything hard is
in computing `h`. The chain ([`burn_rate_cell`](../../srm_1d/burn_rate.py#L342)):

1. **Reynolds number** `Re = ρuD_hyd/μ` — computed *post-PISO* from the current
   flow (`_post_piso_update`) and passed in. This is the *only* place the
   crossflow (`G = ρu`) enters, and it enters as `Re`, buried inside the
   correlations below. Cells with `Re < 100` get **no erosion** (`r_e = 0`).
2. **Haaland friction factor** `f`
   ([`haaland_friction`](../../srm_1d/burn_rate.py#L82), Ma Eq. 15) — Darcy
   factor across laminar (`64/Re`), turbulent (Haaland's explicit
   Colebrook–White approximation using roughness `ε/D`), and a linear
   transition blend. Feeds the turbulent Nusselt formula.
3. **Gnielinski Nusselt number** `Nu`
   ([`gnielinski_nusselt`](../../srm_1d/burn_rate.py#L153), Ma Eqs. 8–10) —
   convective heat transfer for laminar / transition / turbulent flow, with an
   **entrance correction** (`D/L`, heat transfer is enhanced near the head end)
   and a **temperature-ratio correction** `(T_gas/T_surf)^κ`. Then
   `h₀ = Nu·k/D_hyd` — the "bare" HTC before blowing.
4. **Transpiration (blowing) correction**
   ([`transpiration_correction`](../../srm_1d/burn_rate.py#L257), Ma Eq. 16) —
   the burning surface injects mass into the boundary layer, thickening it and
   *reducing* heat transfer: `h = h₀ · β/(exp(β)−1)`, with the blowing
   parameter `β = ρ_p·r·Cp/h₀`. This is a **self-limiting feedback**: more burn
   → more blowing → less heat → less erosion, so erosive burning saturates
   instead of running away.

Because `β` depends on `r`, the equation is **implicit in `r`**. The code
solves it by **bisection** on `F(r) = r − r₀ − r_e(r)`
([L459](../../srm_1d/burn_rate.py#L459)), 30 iterations. Bisection (not
fixed-point) because the transpiration feedback makes the fixed-point Jacobian
exceed 1 at moderate erosion (it would diverge); `F` is monotincreasing so
bisection always converges. *(This per-cell root-find is the hot spot the
deferred "A2 bisection warm-start" perf idea targets — see
[`../core_loop_opt/`](../core_loop_opt/).)*

**Zero tuning constants.** Every term traces to a named correlation (Haaland,
Gnielinski) or a measured physical property (from CEA/RPA gas transport or
propellant characterization). The only calibration knobs are `roughness`
(surface roughness `ε`) and `kappa` (the Gnielinski temperature-ratio
exponent) — both with published physical bounds (see the
`feedback_roughness_kappa_physical_bounds` calibration notes; don't push them
outside physical ranges to chase a fit). This "no arbitrary constants" property
is *why* Ma-2020 is the model of record here and the flux-power-law models are
not — see `feedback_keep_ma_erosive_model`.

**Transport-property caveat:** `Pr`, `k_thermal`, `Cp_gas`, `μ` for the gas
strongly affect `Re` and `h`, hence the erosive rate. srm_1d uses **effective**
(not frozen) RPA/CEA transport for the fired motors; if these aren't in the
`.ric`, the code refuses to fabricate them (it hard-faults) rather than
guessing — see `05`.

---

## Part B — Grain geometry / regression

The geometry side answers: *as the grain burns back, how do the port area,
burning perimeter, and hydraulic diameter change per cell?* All in
[`grain_geometry.py`](../../srm_1d/grain_geometry.py).

### B.1 — Primary state is regression depth
The authoritative per-cell state is `regress[i]` — how far the surface has
burned back (meters). Everything else (`A_port`, perimeter, `D_hyd`) is
*derived* from it each update. Advanced each step by
`advance_bore_regression`: `regress[i] += r_total[i] · f_active · dt`.

### B.2 — Two cell kinds: analytic vs FMM
`update_cell_geometry` branches per cell on `cell_segment_type`:

- **Analytic (type 0)** — cylindrical/conical BATES ports:
  `D_port = D_bore_init + 2·regress`, `A_port = π/4·D²`, `perimeter = π·D`,
  `D_hyd = D_port`. Cheap, closed-form.
- **FMM (type 1)** — complex cross-sections (Finocyl, Star, Moonburner, C/D,
  X-core, Custom). srm_1d can't evaluate those analytically, so it **borrows
  openMotor's Fast Marching Method**: at setup, [`fmm_grain.py`](../../srm_1d/fmm_grain.py)
  runs openMotor's regression-map generation once and samples
  `(reg_depth → perimeter, port_area)` into flat CSR-packed tables. In the hot
  loop, each FMM cell looks up its perimeter/area from its table by regression
  depth. `D_hyd = 4·A_port/perimeter` (correct for non-circular ports).

This is why `fmm_grain.py` exists and lazily imports the sibling openMotor
checkout — srm_1d reuses openMotor's grain geometry rather than reimplementing
it.

### B.3 — `C_burn`: the burning perimeter (the link to `mass_source`)
The quantity that couples geometry back into the gas solve is
```
C_burn[i] = base_perimeter · grain_frac · f_active
```
the **effective burning perimeter** in cell `i`, from which
`mass_source[i] = ρ_p · C_burn[i] · r_total[i]` (Step 3). Two modulating
factors:

- **`grain_frac`** — axial overlap: how much of the cell actually contains
  grain (handles cells straddling a segment edge or gap; cells spanning two
  segments accumulate `C_burn` from both — the v0.6.0 volumetric-overlap
  accumulator).
- **`f_active`** — the **radial burnout ramp**: as a cell's web is consumed,
  its contribution ramps to zero. **Mass-conservation gotcha:** `f_active` must
  multiply **both** `C_burn` *and* the regression rate (`advance_bore_regression`
  applies the same ramp). Applying it to one but not the other causes 3–40 %
  mass error. Do not "simplify" this.

### B.4 — End-face burning
Inhibited/uninhibited grain *ends* regress axially too
(`advance_endface_regression`, Saint-Robert only — no erosive term on an end
face). Their mass injection uses a **partition-of-unity end-face kernel**: each
end face's mass is split across its 2 adjacent cells with weights summing to
1.0 (coupled to the geometry snapping that puts segment edges on cell
boundaries). Gated by `tests/test_endface_conservation.py`. Can be isolated
off with `diagnostic_disable_endfaces`.

### B.5 — Update cadence & taper (extensions)
- `update_cell_geometry` and `compute_burn_rates` run every
  `burn_update_interval` flow steps (geometry changes slowly vs. the acoustic
  timestep, so recomputing every step is wasteful). Between updates the rates
  and areas are held.
- **Axial taper** (v0.8.0): grains can taper along their length (bore and/or
  outer diameter), giving per-cell `cell_D_outer[i]` and per-station FMM
  tables. This is a solver-agnostic extension layered on the above; see
  `project_taper_feature` notes and `motorlib/taper.py`. You can ignore it
  while learning the core.

---

## How Steps 1–2 chain into the rest

```
r_total (Ma, §A) ──► advance_bore_regression ──► regress[i]
                                                    │
                          update_cell_geometry ◄────┘
                                    │
             A_port, D_hyd, C_burn ─┼─► Step 3: mass_source = ρ_p·C_burn·r_total
                                    └─► PISO consumes A_port/D_hyd (Step 4)
                          Re (post-PISO) ──► feeds next step's Ma burn rate
```

The loop is a feedback: flow sets `Re` → `Re` sets the erosive burn rate →
burn rate regresses the grain and injects mass → mass drives the flow. Next:
`03_IGNITION.md` (how a cell *starts* burning, and the source assembly).
