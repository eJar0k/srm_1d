# 03 — Ignition & Source Assembly (loop Step 3)

Step 3 is where a motor *starts*. It answers two questions each timestep:
**(1)** what does the igniter inject, and **(2)** which cold propellant cells
have heated enough to start burning — and then it packs everything into the
`mass_source / thermal_source / momentum_source` arrays the PISO step consumes.

Two modules plus the in-loop assembly:
[`igniter_plenum.py`](../../srm_1d/igniter_plenum.py) (the pyrogen igniter),
[`solid_thermal.py`](../../srm_1d/solid_thermal.py) (the Goodman ignition
kernel), and the `_goodman_ignition_*` / source-assembly block inside
`_run_time_loop`.

---

## Part A — The pyrogen igniter

A cold motor won't self-sustain: the bore is at ambient pressure and the
surface is cold, so `r₀ = a·Pⁿ` is tiny and no heat reaches the grain. The
**igniter** provides the initial hot-gas kick that pressurizes the bore and
heats the surface until the grain takes over. srm_1d models real pyrogen
igniters (a small pyrotechnic charge — `bpnv`, `mtv`, …), not an abstract
"pressure ramp."

### A.1 — The plenum model (`forward_plenum`, the default)
[`PyrogenChamber`](../../srm_1d/igniter_plenum.py) is a **0-D plenum**: a small
chamber holding a pyrogen charge that burns, pressurizes, and vents hot gas
through an orifice into **cell 0** of the main bore. Its state is
`[m_pyrogen, m_gas, T_gas]`, advanced each timestep by
[`_step_plenum_ode`](../../srm_1d/igniter_plenum.py) (RK4):

1. the solid pyrogen burns (its own Saint-Robert law) → gas generation;
2. the plenum gas gains mass + energy and loses it out the orifice;
3. venting mass flow is [`_choked_orifice_mdot`](../../srm_1d/igniter_plenum.py)
   — choked when the plenum/bore pressure ratio is above critical, subsonic
   otherwise.

It returns `mdot_igniter`, plenum pressure `P_ig`, and temperature `T_ig` for
the step. Those become cell-0 (or head-distributed) source terms:

- **mass:** `mass_source[0] += mdot_igniter / dx` (as igniter species, `s=0`);
- **enthalpy:** `thermal_source[0] += mdot_igniter · Cp_pyrogen · T_ig`
  (note: the pyrogen's **own** `Cp`, not the grain gas's — species-specific);
- **axial momentum:** the orifice exit velocity injected at face 1, scaled by
  `igniter_axial_momentum_fraction` (default 1.0 = a pure head-end axial jet).

The igniter also **heats the grain surface directly** (DeMar convention) — see
Part B.

Sizing defaults (charge mass, plenum volume, orifice area) come from Sutton
correlations in `build_pyrogen_chamber`; you normally just pass
`pyrogen='bpnv'`.

### A.2 — Uncontained topologies (opt-in)
`injection_topology` can instead be `head_basket` or `aft_basket`: loose
pyrogen pellets sitting in the bore that burn at the **local** bore pressure
with **no plenum, no orifice, no momentum jet**, injecting mass/enthalpy
directly into their host cells. This is a research feature for modeling
basket igniters; details in `04_ADDINS.md`. The default `forward_plenum` is
what the validated motors use.

---

## Part B — Igniting the grain: the Goodman kernel

A propellant cell starts burning when its **surface temperature** reaches the
ignition threshold `T_ignition`. To know the surface temperature we must model
the transient heat conduction into the cold solid as hot gas washes over it —
that's [`solid_thermal.py`](../../srm_1d/solid_thermal.py).

### B.1 — Heat-balance integral (why δ, not a full PDE)
Solving the full 1-D heat-conduction PDF into each cell's solid every timestep
would be expensive. The **Goodman heat-balance integral method** instead
assumes a **cubic temperature profile** in the solid over a **thermal
penetration depth `δ`**, and tracks just that one state variable per cell.
Given `δ`, the surface temperature is *algebraic*
([`_compute_T_surf`](../../srm_1d/solid_thermal.py#L27)):
```
T_surf = (3·k_solid·T_initial + h_c·δ·T_gas) / (3·k_solid + h_c·δ)
```
and `δ` advances by an ODE ([`_goodman_rhs`](../../srm_1d/solid_thermal.py#L46))
integrated with RK4 (with adaptive sub-stepping when the step is stiff). As the
gas keeps heating the surface, `δ` grows and `T_surf` climbs from `T_initial`
toward `T_gas`.

Here `h_c` is the **convective heat-transfer coefficient from the hot bore gas
to the cold surface**, `k_solid` is the propellant thermal conductivity, and
`alpha = k_solid/(ρ_p·C_ps)` is its thermal diffusivity.

### B.2 — The ignition criterion
When [`_surface_has_ignited`](../../srm_1d/solid_thermal.py#L115) fires
(`T_surf > T_ignition`), the cell flips to `is_burning = True`, its
`ignition_time` is recorded, and from then on it contributes a grain
`mass_source` (Part C). This **surface-temperature threshold** criterion is the
one the project validated against (Keller 1966); see the
`project_ignition_fidelity_candidates` notes.

The ignition calibration knobs are **`k_solid`** and **`T_ignition`** (v0.7.5
defaults 0.271 W/(m·K) and 756 K). Extra heat paths into the surface, all
optional, are covered in `04`: **DeMar** direct pyrogen-plume heating (the
igniter radiates/convects onto nearby grain), and **adjacent-cell radiation**
(a burning cell radiates to unignited neighbors).

> **Note for readers:** the time loop calls `_step_goodman_ode` +
> `_surface_has_ignited` **directly** from `solid_thermal.py` each step for
> unignited cells (imported at the top of `simulation.py`), so this module is
> the live ignition kernel — not just a reference. (Its older docstring
> claiming it's "not wired in yet" was stale and has been corrected.)

---

## Part C — Source assembly (the output of Step 3)

With the igniter stepped and grain ignition updated, Step 3 fills the three
source arrays that Step 4 (PISO) consumes:

```
for each cell i:
    if igniter injects here:   add pyrogen mass + enthalpy (+ momentum at face)
    if grain cell is burning:  add ρ_p · C_burn[i] · r_total[i]  (mass)
                               add that ṁ · Cp_grain · T_flame     (enthalpy)
```
branching on `topology_code` for where the igniter mass lands. The block also
handles some **energy-balance bookkeeping** that trips people up:

- the pyrogen's enthalpy is injected as gas **and** it radiates to the walls —
  the radiated part is **debited** from the gas so it isn't double-counted;
- the bore gas loses heat to still-cold walls (a convective wall-loss sink).

These are small but correct; they're recorded in the energy-audit histories so
conservation can be checked. Each contribution can be isolated with the
`diagnostic_disable_*` flags (`04`).

**Result:** by the end of Step 3, `mass_source`, `thermal_source`, and
`momentum_source` fully describe what combustion + ignition add to the gas this
step. Step 3b evolves the throat, Step 4 runs PISO, and the cycle repeats. Next:
`04_ADDINS.md` — the opt-in features (Phase F/Z, radiation, topologies, the
diagnostic toggles) in mechanism-level detail.
