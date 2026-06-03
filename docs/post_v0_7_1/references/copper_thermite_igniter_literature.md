# Copper Thermite as Alternative Igniter — Literature Review

**Research date**: 2026-05-22 (during v0.7.1 Phase 5)
**Authoring agent**: general-purpose subagent (haiku) delegated by the
v0.7.1 Phase 5 session. Distilled report from the agent's full
transcript; raw transcript was transient and is not preserved.

## Why this exists

The amateur rocketry and industrial development communities are
increasingly using copper thermite (CuO/Al) as an alternative to
gas-rich pyrogens like BPNV. The user raised whether srm_1d's
current pyrogen plenum model would accommodate it, and whether
copper thermite should be on the v0.7.x or v0.8 roadmap.

This document captures the open-literature answer.

## Chemistry and combustion

| Property                | Value                                     |
|-------------------------|-------------------------------------------|
| Reaction                | 2 Al(s) + 3 CuO(s) → Al₂O₃(l) + 3 Cu(l→v) |
| Adiabatic flame temp.   | ~2843 K (Reese et al. 2015)               |
| Action time             | ~75 ms (mass-independent)                 |
| Combustion regime       | Deflagration                              |
| Cu vapor yield          | ~77.6 g/kg thermite                       |
| Gas mass fraction       | ~10% of products (truly low-gas)          |
| Solid/liquid products   | Al₂O₃ particles (5-30 µm), molten Cu      |

For comparison, BPNV's T_flame is ~2800 K (similar) but its product
gas fraction is ~80%, with N₂ and CO₂ dominating — completely
different thermochemistry from thermite's condensed-phase
combustion.

## The "no gas" claim — quantified

The often-cited "thermite produces no gas" characterization is
**partially true but misleading** for modeling purposes:

- Copper thermite produces ~10% gaseous products by mass (Cu vapor +
  O₂ from CuO decomposition) — small but non-zero.
- Al/Fe₂O₃ thermite variants ARE explicitly gasless (rated as
  vacuum-compatible).
- The salient point for srm_1d modeling: **~90% of the chemical
  energy is delivered as condensed-phase mass** (molten metal
  droplets + Al₂O₃) and radiation, NOT as hot gas advected through a
  plenum.

## Heat-flux delivery

| Quantity                                            | Value / source              |
|-----------------------------------------------------|-----------------------------|
| Minimum effective flux for SRM ignition             | 0.47 MW/m² (~20 cal/cm²/s)  |
| Reese 2015 dual-criteria ignition model             | Critical T + critical E∫    |
| Droplet size (high-speed imaging)                   | 10-50 µm typical            |
| Combustion time                                     | ~50-75 ms                   |
| Heat-transfer mode                                  | Radiative + droplet contact |
| Typical standoff distance (amateur designs)         | 5-20 mm                     |

Quantitative thermite-specific heat-flux numbers are sparse in open
literature. Reese et al. 2015 use a **dual-criteria ignition model**
requiring both:

1. **Critical temperature**: grain surface temperature raised above
   ignition threshold.
2. **Critical energy**: subsurface energy integral sufficient for
   sustained burn after the thermite extinguishes.

The "no bore pressurization" property of thermite is operationally
significant — unlike BPNV's axial gas jet, thermite does not
pressurize the chamber before the propellant ignites. This is a
safety advantage but means pressure-fed systems may require a
secondary igniter.

## srm_1d modeling implications

The current v0.7.0+ pyrogen plenum is structurally a **0D chamber +
choked nozzle that injects (mass + enthalpy + axial momentum) into
bore cell 0**. Mass conservation through that nozzle is the integrity
assumption.

Thermite breaks the assumption. Its primary energy delivery is
condensed-phase droplets with negligible bore-gas contribution. The
bore sees **heat flux without mass injection** until the propellant
itself ignites.

### Minimal viable model (v0.8 Phase 1, ~1 week effort)

- Add a `thermite_surface_heat_flux(t)` boundary condition applied
  directly to grain surface in cell 0.
- Parameterize by: thermite mass, burn rate, droplet size
  distribution, emissivity.
- **No mass injection** into bore — treat as a pure surface flux
  source (similar to the current adjacent-cell radiation pathway).
- Apply Reese's dual-criteria ignition check post-hoc.

### Full treatment (v0.8 Phase 2+, 2+ months effort)

- 0D droplet trajectory module (ballistic + cooling).
- Radiosity solver coupling droplet-cloud temperature to grain
  surface.
- Probabilistic impingement enforcing droplet mass conservation.
- Likely requires 2D axisymmetric geometry to capture droplet
  scatter — significant departure from srm_1d's 1D PISO core.

## Roadmap placement

| Phase             | Action                                              |
|-------------------|-----------------------------------------------------|
| v0.7.x            | Not in scope. Document as deferred igniter type.    |
| v0.8 Phase 1      | Surface-flux-only thermite (minimal viable model)   |
| v0.8 Phase 2+     | Droplet physics / radiosity if needed for validation |

The architectural separation (different boundary-condition path)
means thermite doesn't risk disturbing the v0.7.x pyrogen
calibration. Adding it post-v0.7.1 is a natural expansion that
serves the amateur rocketry / CubeSat communities (where thermite is
common for cost/safety reasons).

## Primary literature

1. **Reese et al. 2015**, "CuO/Al Thermites for Solid Rocket Motor
   Ignition," AIAA Journal of Propulsion and Power.
   **Gold standard.** Experimental characterization, dual-criteria
   ignition model, hot-fire validation. Safety testing showing
   CuO/Al safer than legacy igniters. This is the canonical
   reference for any v0.8 implementation.
   https://arc.aiaa.org/doi/10.2514/1.B34771

2. **Nishii et al. 2024**, "Application of an Al/Fe₂O₃ Thermite
   Reaction to an Igniter of a Hybrid Rocket." Only published 1D
   transient modeling of a thermite igniter coupled to a fuel
   grain. Simplified one-step kinetics + FDM. Confirms gasless
   combustion of the Al/Fe₂O₃ variant. The closest existing
   reference for what a v0.8 srm_1d implementation would look like.
   https://www.researchgate.net/publication/368302735_Application_of_an_AlFe2O3_Thermite_Reaction_to_an_Igniter_of_a_Hybrid_Rocket

3. **Richard Nakka's Experimental Rocketry — Thermite Experiments**.
   Authoritative amateur reference. Documents copper vs iron
   thermite, burn rates, safety considerations. Aligned with
   srm_1d's open-access philosophy. Useful for validation against
   amateur-scale data once the model lands.
   https://www.nakka-rocketry.net/thermites.html

4. **U.S. Patent 4,464,989** (1984), "Integral Low-Energy Thermite
   Igniter." Historical formulation reference; consolidates
   thermite-formulation space (CuO, MoO₃, Fe₂O₃ variants).
   Establishes "low-gas-output" criterion as a deliberate design
   feature.

5. **INL Digital Library** — "Quantification of Heat Flux from a
   Reacting Thermite Spray." Quantitative thermite flux
   measurements (one of the few published).
   https://inldigitallibrary.inl.gov/sites/sti/sti/4310589.pdf

6. **Frontiers in Chemistry 2018** — "Combustion Characteristics of
   Physically Mixed 40 nm Aluminum/Copper Oxide Nanothermites."
   Nano-scale CuO/Al variant; relevant only if v0.8+ explores
   nano-thermites which have different combustion kinetics.
   https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2018.00465/full

## Conclusion

> Copper thermite is a viable v0.8+ feature. It requires a separate
> boundary-condition path (surface-flux-only, no mass injection) and
> cannot be cleanly bolted onto the v0.7.x 0D plenum model.
> Reese 2015's dual-criteria ignition framework is the canonical
> reference. Out of scope for v0.7.1.

## Related memory

- `[[srm-1d-copper-thermite-igniter-research]]` — distilled memory
  entry used by future srm_1d sessions for quick lookup.
- `[[srm-1d-pyrogen-heat-flux-literature-bounds]]` — sibling
  research on conventional pyrogen heat flux; thermite extends
  those ranges and shifts the architecture.
- `[[feedback_igniter_conventions]]` — user's amateur-rocketry
  igniter intuition; thermite is a natural extension they'll want
  to model eventually.
