# Favini, Cavallini, Di Giacinto, Serraglia (2008) — AIAA 2008-5141

**Full title:** "An Ignition-to-Burn Out Analysis of SRM Internal Ballistic and Performances"

**Venue:** 44th AIAA/ASME/SAE/ASEE Joint Propulsion Conference, 20-23 July 2008, Hartford, Connecticut

**PDF:** `srm_1d/docs/post_v0_7_0/digiacinto2008.pdf` (19 pp)

**Authorship note:** First author Favini; "digiacinto2008" filename matches the project's naming convention but Favini is principal author. Cavallini is the PhD student at this point.

---

## 1. Scope and significance

This is the **first public SPINBALL paper** — the SPIT → SPINBALL transition announcement. Its body, however, is **80% devoted to GREG** (the 3D level-set grain burnback module). The flowfield-physics changes from SPIT to SPINBALL are introduced briefly in §II.A; the underlying IT sub-models (igniter, ignition criterion, conduction, heat transfer) are explicitly **deferred to the SPIT references** (Refs 1, 3, 6, 7):

> "For the other sub-models present in SPINBALL and heritage of the SPIT code, we address the interested reader to the following works (Ref. 1,3,6,7)." (p. 2)

So this paper is the **bridge document**, not the model document. For spike-taildown research:

- It is decisive about what changed (the gasdynamics formulation, GREG, TP ablation models)
- It is silent about IT sub-model equations
- It explicitly validates Z9 head-end pressure for the entire burn (Fig. 7) — see §4 below

## 2. What changed from SPIT to SPINBALL (§II.A)

The single architectural change the paper foregrounds: the **"variable properties mixtures of gases" formulation**. Each cell carries its own thermophysical state (γ, Cp, R, MW), computed as a **weighted average** of:

- mixture inflows from the two adjacent cells (L, R) at the current timestep, and
- source contributions in the current cell from grain combustion, igniter, cavities, and TP ablation.

Weights are species concentrations evaluated at that cell-timestep. This is fundamentally different from:

- **SPIT (≤2007):** 3 explicit species mass equations (igniter / pressurizing gas / propellant gas) with each species's thermophysical properties held **constant in space and time**.
- **SPINBALL (2008+):** "infinite-gases" mixture — properties vary in space and time per cell, but no species-specific transport equation is solved.

> "This chamber flowfield model defines a set of differential equations, composed of: one mass conservation equation, one momentum conservation equation and one energy conservation equation for the gases mixture present in the bore; closed by the perfect gases equation." (p. 2)

The paper does **not** give the explicit equations here — those appear in the 2009 conference paper (Eq. 1) and the 2009 thesis. The 2008 paper only describes the philosophy.

## 3. Gasdynamics numerics

> "The numerical discretization is represented by a finite volume Godunov scheme first/second order accurate in space and time, coupled with an exact or approximate Riemann solver. The time discretization is characterized by a CFL-like condition that guarantees the stability of the numerical scheme adopted." (p. 3)

First-order space/time as default; second-order with ENO + MinMod limiter + Heun stepping available. Exact Riemann solver. Standard Godunov mixture-of-gases construction.

**Critical contrast with srm_1d:** srm_1d v0.7.0 uses **PISO** (pressure-implicit splitting of operators) with TDMA + adaptive CFL. PISO is a pressure-based, incompressible-leaning method; Godunov is a density-based, fully-compressible shock-capturing method. SPINBALL's choice is dictated by needing to capture igniter-seal-rupture shocks and strong IT-region discontinuities. For Hasegawa A's smooth IT (no nozzle seal, gentle pyrogen ramp), this difference is not first-order; PISO is fine.

## 4. Ablation models (§II.B)

Simple semi-empirical TP/nozzle ablation models tied to convective + radiative heat flux from the bore. Bartz correlation is the example. Mass added back to the bore.

> "These models are based on modeling the mass flow rate from these materials, and consequently their regression velocity, as dependent directly on the heat fluxes due to radiation and convection from chamber hot gases" (p. 3)

**srm_1d v0.7.0:** Has throat erosion in the openMotor-aligned Nozzle (Bartz-flavored). Does not have liner/case TP ablation — for Hasegawa A this is irrelevant; the motor is short, the run is short, the case is non-ablative.

## 5. GREG (Grain REGression, §II.C)

The 3D level-set burnback module. ~9 of the paper's 19 pages. Key features:

- **Level Set on Eulerian rectangular or cylindrical structured grids** (Eq. 1). Hamilton-Jacobi PDE: `φ_t + r_b · |∇φ| = 0`.
- **Banded SDF initialization from STL files** (Eq. 5–8). Cuts initialization cost by 1-2 orders of magnitude vs. full SDF.
- **First-order Godunov-type numerical Hamiltonian with exact Riemann solver** (Eq. 2). CFL: `Δt ≤ min[CFL/r_b · (|n_x|/Δx + |n_y|/Δy)]` (Eq. 4).
- **Grain geometry extraction**: integrate Heaviside (volume, Eq. 9), Dirac delta (area, Eq. 11) over the level set field — smeared regularizations per Peskin 1977 (Eq. 13–14). Burning surface from volume budget (Eq. 15): `S_b(t^{n+1/2}) = [V_b(t^{n+1}) − V_b(t^n)] / [r_b · Δt]`. Burn perimeter analogously (Eq. 16).
- **Off-line vs. on-line coupling** with the flowfield solver. Off-line: GREG runs once with constant r_b, produces lookup tables of (P_b, P_w, A_p) vs. web; flowfield interpolates. On-line: full coupling. 2008 paper uses off-line.

**Relevance to srm_1d:** srm_1d v0.7.0 uses FmmGrain from openMotor — scikit-fmm-driven 2D fast-marching with table lookup. Mathematically related (level set = front propagation), but srm_1d ports openMotor's structure rather than re-deriving it. GREG's 3D level-set is more general (handles full 3D finocyl + submergence); srm_1d's 2D FMM + table approach is sufficient for cylindrical bates / Hasegawa A.

## 6. Z9 ignition-to-burnout validation (§III.B)

The lone full-run validation result. SPINBALL HEP trace vs experimental for Z9, the third Vega solid stage (HTPB, 11 t propellant, finocyl + aft star, submerged nozzle).

**Result quote** (Fig. 7):

> "Figure 7 shows a good agreement of the numerical predicted head-end pressure respect to the experimental one, until almost 0.7 adimensional time, longer that point the disagreement could be ascribed to the following facts."

The disagreement after τ ≈ 0.7 is attributed to:

1. **Inadequate burning-surface evolution coupling** — at this point GREG is not yet coupled on-line to SPINBALL, so they use SPP-generated tables. Late-time geometry effects are missed.
2. **TP and ablative material sub-models not yet refined** — relevant because grain burning surface is shrinking, so TP gas addition becomes a larger fraction of the total mass flow rate.

**Critical for spike-taildown question:** No mention of any IT-specific disagreement. SPINBALL/SPIT captures the Z9 IT to the resolution of the plot — peak P, rise rate, and equilibration look correct. The 30% time-window of disagreement is at the **tail-end** of the run, not the IT spike taildown.

## 7. Authoritative SPIT pointers (refs 1, 3, 6, 7)

The four SPIT references this paper points at for "what we actually do for ignition transient":

| Ref | Citation | Likely content |
|---|---|---|
| 1 | Di Giacinto et al. 2007, AIAA, 43rd JPC, Cincinnati | VEGA SRM IT firing test prediction + post-firing analysis. Likely the most complete SPIT description in the modern era. |
| 3 | Favini et al. 2006, AIAA, 42nd JPC, Sacramento (Paper 5141 in 2006 series) | Z9 SFT pre-test prediction with SPIT. |
| 6 | Favini et al. 2005, EUCASS, 1st European Conf, Moscow | "Pressuring Gas Effects on Pressure Oscillations during IT" — multi-species effects, this is likely where the 3-gas SPIT model is documented mathematically. |
| 7 | Serraglia et al. 2004, ESA/DLR, Cologne | "Gas Dynamic Features in SRM with Finocyl Grain during Ignition" — finocyl-specific IT phenomenology. |

Additionally from the bibliography (not pointed at by the deferral, but relevant):

- Ref 10: Di Giacinto & Serraglia 2001, AIAA 2001-3448, "Modeling of Solid Motor Start-up" — likely the earliest published SPIT model description.
- Ref 13: Serraglia 2003, PhD thesis, "Modeling and Numerical Simulation of Ignition Transient of Large SRMs" — the doctoral parent of SPIT.

**Gap status:** If we ever need to port SPINBALL's exact igniter sub-model or ignition criterion equations, **Ref 10 (2001) and Ref 13 (2003 thesis)** are the targets. Neither is in repo. **Not blocking** for the current spike-taildown analysis — see [[spinball_walkthrough]] for why.

## 8. What this paper contributes uniquely (not in 2009 paper or thesis)

1. The **explicit SPIT-to-SPINBALL transition narrative** — establishes that variable-properties mixtures is *the* architectural innovation.
2. **GREG-specific equation set** (Level Set numerics + STL-banded SDF). The 2009 thesis covers this but in less compact form.
3. **Z9 entire-burn validation figure** — useful as a cross-comparison datum if we ever want to compare our adapter-loaded motors to a published SPINBALL run.

## 9. What this paper does NOT give us

- The Q1D flowfield equations in explicit form (these appear in 2009, see [[extraction_spinball_2009]]).
- Igniter sub-model equations (deferred to SPIT references — see §7).
- Ignition criterion equations (deferred to SPIT references).
- Solid-phase conduction equations (deferred).
- Heat transfer model details (Bartz mentioned, no further detail).
