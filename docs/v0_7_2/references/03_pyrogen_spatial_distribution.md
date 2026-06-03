# Pyrogen-Injection Spatial Distribution — Extended Literature Digest

*Subagent literature dive, 2026-05-23. Extends what srm_1d already had
documented (`docs/post_v0_7_1/references/pyrogen_heat_flux_literature.md`,
`copper_thermite_igniter_literature.md`, and the v0.7.0 plenum/Goodman
design docs).*

## What srm_1d already has documented

Existing references:
- **DeMar 1995 / Sandia LDRD 2022 / AIAA JSR BPN pyrogen** — heat-flux
  *magnitude* bounds (0.47 MW/m² functional minimum to 23.2 MW/m²
  sustained, 1 GW/m² peak near-field). Time profile is open work;
  **axial profile is not addressed.**
- **Reese 2015 + Nishii 2024 (thermite)** — dual-criteria ignition;
  surface-flux boundary condition concept. Both still treat the
  igniter as a head-end region, not a distributed axial source.
- **v0.7.0 DESIGN.md** explicitly defers "multi-cell impingement-
  region distribution (Cavallini SPINBALL)" pending physical
  jet/basket geometry inputs. Cell-0 injection is acknowledged as a
  simplification.

## New findings from this dive

**1. Peretz/Kuo/Caveny/Summerfield (1973) is the foundational SRM
ignition-transient framework.** It treats the head-end igniter as a
*boundary condition at x=0* (mass + enthalpy injection through the
head-end face) but then computes axial energy deposition implicitly
via **a flame-spreading sub-model**: the convective heat flux to the
unignited propellant surface is computed cell-by-cell using a
**Dittus-Boelter Nu ~ Re^0.8 Pr^0.4** correlation against the *local*
bore gas state. Surface ignition occurs cell-by-cell when the
integrated surface energy reaches a threshold. This is the canonical
"lumped at head-end, axially distributed via convective transport"
architecture that nearly every academic 1D ignition code follows,
including Han (NASA, 1990-1993), Cai & Yang 2001, Salita's SPP, and
Cavallini SPINBALL.

**2. Cavallini SPINBALL ("impingement region") is the most explicit
published axial-distribution model.** The PhD thesis (~2009) and the
SPINBALL papers (Favini/Cavallini/Di Giacinto, AIAA papers and JANNAF
reports) define an **impingement region** of finite axial length
L_imp downstream of the head end. Inside L_imp, igniter mass and
enthalpy are deposited per-cell weighted by local geometry; outside
L_imp, only the convective flame-spreading mechanism operates. L_imp
is parameterized by *igniter physical hardware* (basket geometry,
axial-vs-radial nozzle split, jet angle) rather than by a
first-principles jet correlation.

**3. Jet-in-crossflow (JICF) is the right physical analog for
radial-vent pyrogens, but axial-vent pyrogens are coaxial jets.**
Soo-Young No 2015 (review in *Int. J. Spray & Combustion Dynamics*)
catalogs JICF penetration correlations of the form
**y/d_j = A · q^a · (x/d_j)^b** with q the momentum flux ratio,
A ~ 1.1-2.5, a ~ 0.4-0.5, b ~ 0.25-0.34. For Hasegawa A-class
pyrogens venting into a quiescent bore at ignition, the head-end
pyrogen plume is effectively a **coaxial axisymmetric jet** (not
JICF) until the bore flow develops. Coaxial-jet theory (Tollmien,
Reichardt, Witze) gives a **potential-core length
L_pc ≈ 6-7 · d_throat** beyond which centerline velocity decays as
U_c/U_0 ≈ 6.3 · d_throat/x. Lateral half-width grows linearly:
r_1/2 ≈ 0.086 · x. Both give a physically motivated decay length
for energy deposition.

**4. NTRS 19710018794 (Hersch/Rieser, NASA TM, 1971)** —
"Correlation of secondary sonic and supersonic gaseous jet
penetration" provides directly applicable empirical jet-penetration
data for sonic injection into a confined flow at SRM-representative
Mach numbers (0.5-4) and pressure ratios. Their correlations are
L_pen/d_j as a function of P_j/P_inf and M_inf and are the cleanest
open-source dataset for sizing an axial-deposition length.

**5. Submerged/aft-pyrogen variants (US Patent 4503773, Aerojet)**
confirm that when industry needs the energy delivered anywhere other
than the head-end, they redesign the igniter hardware — they don't
model the head-end plume as automatically distributing. This
reinforces that the cell-0 lumping is a *modeling* artifact, not a
known physical fidelity ceiling.

## Implementation guidance for srm_1d

**Simplest defensible profile: exponential decay weighted by a
characteristic jet length L_jet.**

For per-cell *i* with axial center x_i, distribute the pyrogen
mdot / enthalpy / momentum as:
```
w_i = exp(-x_i / L_jet) · dx_i / sum_j[exp(-x_j / L_jet) · dx_j]
mdot_i = w_i · mdot_plenum
```
This conserves total mass exactly and reduces to current cell-0-only
behavior in the limit L_jet → 0. **L_jet** is the one new YAML knob,
parameterized as `L_jet = κ_jet · d_throat` where κ_jet is the
dimensionless decay length. Defensible defaults from the literature:
- **κ_jet ≈ 6-10** for choked sonic axial-vent pyrogens (Witze +
  Hersch/Rieser, potential-core + early-decay band).
- **κ_jet ≈ 2-4** for predominantly radial-vent pyrogens (impingement
  is local; JICF penetration y/d_j is modest at q < 5).

For ROC-A class hardware where d_throat ≈ 7 mm and bore length
≈ 0.3 m, L_jet ≈ 0.04-0.07 m — meaning energy is realistically spread
over the **first 10-25% of the bore**, not concentrated in cell 0.
This is the qualitative change you're after.

A second-tier upgrade adds a **lateral-impingement weight**: weight
w_i by the local available burning surface area (already in
`_pyrogen_surface_thermal_sink`), so the same axial profile properly
distributes heat-flux among segments of different perimeter.

## Known limitations

- Open literature is much weaker on heat-flux **vs axial position**
  measurements than on flux magnitude — quantitative validation of
  κ_jet will require either CFD calibration or Hasegawa-A-class
  flame-spreading delay matching.
- Exponential decay is an empirical shape, not derived from JICF or
  coaxial-jet theory exactly. Witze gives ~1/x decay beyond the
  potential core, which is heavier-tailed; exponential is the safe
  interpolant that doesn't put unphysical energy at the aft end.
- Effect on the v0.7.1 calibration (k_solid, pyrogen heat flux) is
  not free — a positive L_jet will reduce cell-0 hot spot and may
  move the rank-1 LHS basin. Expect re-calibration scope.
- Coupling to the v0.7.1 N-species mixture is straightforward: each
  cell receives w_i fraction of the pyrogen species mass; Phase 3.5's
  per-species Cp applies cell-by-cell automatically.

## Primary citations (with URLs)

1. **Peretz, Caveny, Kuo, Summerfield (1973)** — "Starting Transient
   of Solid Propellant Rocket Motors with High Internal Gas
   Velocities," NASA TN. Foundational head-end-boundary + axial-flame-
   spreading framework. https://ntrs.nasa.gov/citations/19740005393
2. **Hersch & Rieser (1971)** — "Correlation of Secondary Sonic and
   Supersonic Gaseous Jet Penetration," NASA TM. Empirical L_pen/d_j
   vs pressure ratio and Mach number for sonic injection.
   https://ntrs.nasa.gov/api/citations/19710018794/downloads/19710018794.pdf
3. **Han (1990 / 1991)** — "Ignition Transient Analysis of Solid
   Rocket Motor," NASA NTRS. 1D and 2D axisymmetric SRM ignition
   modeling with head-end boundary igniter and Dittus-Boelter flame
   spreading. https://ntrs.nasa.gov/citations/19910009672 and
   https://ntrs.nasa.gov/citations/19920006646
4. **Bai, Han, Pardue (1993)** — "2D Axisymmetric Analysis of SRM
   Ignition Transient," AIAA 29th JPC. Higher-dimensional confirmation
   that head-end igniter + convective spreading captures the physics.
   https://ntrs.nasa.gov/citations/19930066098
5. **Soo-Young No (2015)** — "A Review on Empirical Correlations for
   Jet/Spray Trajectory of Liquid Jet in Uniform Cross Flow,"
   *Int. J. Spray & Combustion Dynamics*. Power-law penetration
   correlations in momentum flux ratio.
   https://journals.sagepub.com/doi/10.1260/1756-8277.7.4.283
6. **Witze (1974) — generic centerline velocity decay** for turbulent
   axisymmetric jets; standard reference for coaxial-jet potential-
   core length and 1/x decay. Summary at NIST:
   https://www.nist.gov/publications/generic-centerline-velocity-decay-curve-initially-turbulent-axisymmetric-jets
7. **Cavallini, Favini, Di Giacinto** — SPINBALL SRM internal
   ballistics papers, including impingement-region multi-cell
   mass/enthalpy distribution. PhD thesis:
   https://www.studocu.com/row/document/yildirim-beyazit-universitesi/grain-simulation/phdthesis-cavallini-padis/67377396
   SPINBALL summary:
   https://www.researchgate.net/publication/268483069_SRM_internal_ballistic_numerical_simulation_by_SPINBALL_model
8. **Parametric igniter study (2025)** — recent ScienceDirect study
   with one axial + twelve radial nozzles, examining mass-flow /
   temperature / pressure / angle effects on induction interval and
   flame spreading.
   https://www.sciencedirect.com/science/article/abs/pii/S1290072925006453
