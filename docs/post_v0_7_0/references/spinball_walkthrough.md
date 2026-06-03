# SPINBALL vs srm_1d v0.7.0 — Mechanism-by-Mechanism Walkthrough

**Purpose.** This document synthesizes [[extraction_spit_to_spinball_2008]], [[extraction_spinball_2009]], and the older [[extraction_peretz_pardue_cavallini]] (with the correction in §10 of the 2009 extraction) into a decision instrument for the next physical-model addition to srm_1d. The motivating problem is the Hasegawa A spike-taildown residual (the dominant remaining MSE contributor in v0.7.0; weight 0.35 in the segmented fitness).

**Bottom line up front (BLUF).** SPINBALL has **no specific mechanism aimed at spike-taildown improvement** that we don't already have or already exceed in v0.7.0. Its principal innovation over SPIT — per-cell variable thermophysical properties from mixture-weighted source history — addresses *QSS/tail-off* accuracy, not IT spike-taildown. For Hasegawa A's residual, the candidates remain [[project_next_session_focus|the three from session-end discussion]]:

1. **Zeldovich-Novozhilov dynamic burn rate** (primary recommendation),
2. **Al2O3 two-phase thermal lag** (secondary; physically motivated for Hasegawa Propellant 1's 14% Al loading),
3. **Peretz-aligned participation fraction** (lowest-leverage; dovetails with Z-N).

A fourth candidate — SPINBALL's **per-cell variable γ/Cp/R** — is in scope but is unlikely to dominate the spike-taildown signal. See §4 "Verdict" below.

---

## 1. Mechanism comparison table

| Mechanism | SPINBALL (2008/2009) | srm_1d v0.7.0 | Verdict |
|---|---|---|---|
| **Gasdynamics: governing equations** | Q1D Euler, area-varying duct, conservative form (Cavallini 2009 Eq. 1): 1 mass + 1 momentum + 1 energy. Igniter momentum injected via ṁ_ig·v̄_inj. | Q1D compressible Navier-Stokes (1 mass + 1 momentum + 1 energy), PISO/TDMA on staggered grid. Igniter mass + enthalpy + axial momentum into cell 0 / face 1 (`igniter_axial_momentum_fraction` knob; `_orifice_exit_velocity` computes `v̄_inj` from plenum state). | Tie. Both inject mass + enthalpy + momentum. srm_1d adds a geometric axial-fraction knob for non-purely-axial igniters. |
| **Numerical method** | Finite-volume Godunov, exact Riemann solver, first-order in space/time for Z23. Density-based, shock-capturing. CFL ~0.5. | PISO (pressure-implicit, splitting of operators) + TDMA + adaptive CFL. Pressure-based. ~45-90k steps/s with Numba. | Different families. PISO is appropriate for smooth IT; Godunov is appropriate for nozzle-seal-rupture shocks. **Hasegawa A has no such shocks — PISO is fine.** |
| **Per-cell thermophysical state** | γ, Cp, R, MW **vary per cell, per step** — weighted average of L/R adjacent cell mixtures + source terms ("infinite-gases" formulation). | Single set of (γ, Cp, R) per propellant tab (PropellantTab), with **frozen vs effective transport** as a tunable knob. Igniter gas treated with same γ/Cp/R as propellant gas. | **SPINBALL is more comprehensive.** srm_1d's frozen-vs-effective toggle is a single-cell-global approximation of the same physics. See §3.1. |
| **Igniter mass-flow source** | ṁ_ig(t) **prescribed** from experiment / igniter design tabulated trace. T_ig and H_ig prescribed. No igniter chamber thermodynamics. Multi-cell impingement region for VEGA-class motors with radial igniter nozzles. | **0D pyrogen chamber** (Saint-Robert burn + choked/subsonic vent) computes ṁ_ig(t) from pyrogen burning. Injected into cell 0 with hat-function partitioning. Sutton sizing defaults. | **srm_1d is more physically grounded** here — we model the igniter chamber, SPINBALL does not. SPINBALL has axial-momentum injection (we don't); we have plenum thermodynamics (they don't). |
| **Ignition criterion** | T_surf-based critical-temperature criterion with **pressure-dependent T_ign(p)**. Equations referred to SPIT lineage (2001-2007 papers). | T_surf ≥ T_ign (constant), via Goodman cubic-polynomial integral method on the solid side. Peretz 1973 form. | **Tie on form; SPINBALL has T_ign(p) variation.** Pressure-dependent T_ign is a small refinement that has not been characterized for Hasegawa Propellant 1. Low priority. |
| **Solid-phase conduction** | 1D PDE in y normal to surface. Likely uses Crank-Nicolson or similar finite difference. Explicit equations referred to SPIT. | **Goodman cubic-polynomial heat-balance integral** (DESIGN.md Eq. 8). ODE form, ~30 μs penetration depth, RK4-style. | srm_1d's form is **identical to Peretz 1973** and Pardue 1992. SPINBALL likely uses a full PDE, which is more accurate but ~10× more expensive. **Goodman is sufficient.** See [[equations_goodman_integral]]. |
| **Erosive burning** | **Lenoir-Robillard 1957** with Lawrence 1968 hydraulic-diameter modification. Semi-empirical, β fitted. | **Ma 2020** — Haaland friction factor + Gnielinski Nusselt + transpiration correction + bisection solver. Zero arbitrary constants. | **srm_1d is materially ahead of SPINBALL** on erosive burning physics. Ma 2020 (2020) postdates the SPINBALL formulation (Lenoir-Robillard 1957). |
| **Dynamic burning (Zeldovich-Novozhilov)** | **Not implemented.** Discussed in Cavallini 2009 thesis §2.2.3, dismissed for VEGA scale ("minor effect at IT and tail-off"). | **Not implemented.** | **Both lack it.** Candidate #1 for Hasegawa A spike-taildown. See §4. |
| **Two-phase / Al2O3 particles** | **Not modeled.** Explicit assumption: "non-reactive mixture of perfect gases." | **Not modeled.** | **Both lack it.** Pardue 1992 has it (two-fluid IPSA). Candidate #2 for Hasegawa A spike-taildown. See §4. |
| **Multi-species transport** | Variable mixture properties — see "per-cell thermophysical state" row. No explicit species transport equation in SPINBALL (was in SPIT). | Not modeled. Single-gas thermo per tab. | SPINBALL's approach is per-cell variable properties, not species advection. See §3.1. |
| **Grain burnback** | **GREG**: full-matrix 3D level-set on rectangular or cylindrical grids. Initialized as banded SDF from STL CAD files. Off-line or on-line coupling to flowfield. | **FmmGrain** ported from openMotor — 2D scikit-fmm fast-marching on a per-slice basis. Numba lookup. Effectively off-line: web-indexed tables. | GREG is more general (full 3D, submergence, finocyl). srm_1d's port is sufficient for cylindrical / bates / Hasegawa A. Not relevant to spike-taildown. |
| **Cavity / submergence** | Cavity model — 0D ODE bolt-on per cavity (slots, submergence). Source ṁ_s into the flowfield. | Not implemented. | Not relevant to Hasegawa A. |
| **Nozzle throat erosion / TP ablation** | Bartz-style convective heat flux drives prescribed throat regression. TP ablation as mass source back into bore (semi-empirical). | openMotor-aligned Nozzle has throat erosion (Bartz-flavored). No TP ablation. | srm_1d is comparable for what Hasegawa A needs. |
| **Calibration architecture** | Separate 0DQSS pre-extracts HUMP(web), η_c*, η_cF, A_t(t) from SFT data; SPINBALL consumes as inputs. | Coupled LHS sweep of model parameters against experimental trace (`tools/sensitivity.py`). | Different philosophies; both legitimate. srm_1d's tight coupling is appropriate for a research code; SPINBALL's separation is appropriate for industrial use. |

## 2. Comparison summary: what SPINBALL does that we don't

### 2.1. Igniter axial-momentum injection

SPINBALL adds an explicit `ṁ_ig·A_p·v̄_inj/V` term to the momentum equation (Cavallini 2009 Eq. 1, momentum). **srm_1d v0.7.0 already does this** ([simulation.py:645-656](../../../simulation.py#L645)) — `_orifice_exit_velocity()` computes `v̄_inj` from the plenum state (choked or subsonic), and `pyrogen_momentum_expected = ṁ_ig · v̄_inj · igniter_axial_momentum_fraction` is deposited at face 1 (between cell 0 and cell 1). The `fraction` parameter (default 1.0 = purely axial head-end jet) is a refinement over SPINBALL for non-purely-axial igniter geometries. Test coverage at [test_ignition_diagnostics.py::test_disabling_momentum_removes_igniter_momentum_source](../../../tests/test_ignition_diagnostics.py).

**Status:** Not a gap. Documented for completeness.

### 2.2. Per-cell variable γ, Cp, R

This is the headline SPINBALL innovation over SPIT. Each cell carries its own thermophysical state, evaluated as a weighted average of the mixtures flowing in from adjacent cells + the source terms in that cell. The weights are species concentrations at that cell-step.

**Impact on Hasegawa A:** The spike-taildown phase is precisely when the igniter gas (frozen-like, k≈0.37) is being displaced by propellant gas (effective-like, k≈0.65). The current v0.7.0 model uses a single k for both — typically the propellant tab value, which over-predicts heat transfer when igniter gas is still resident. A per-cell variable γ/Cp/R model would handle this transition.

**Estimated leverage on spike-taildown residual:** Moderate. The frozen-vs-effective sensitivity in v0.7.0 calibration is documented (`feedback` memory `feedback_no_unfounded_smoothing` and Hasegawa A calibration state). The peak is at frozen-end (k~0.37) but the QSS approaches effective-end (k~0.65). A spatial+temporal blend implicitly does what the LHS chose with a global toggle.

**Action:** Add to v0.7.1 roadmap as a SECONDARY candidate. Implementation cost is moderate (2 extra state variables per cell: Y_ig and Y_prop, plus a mixing rule for γ and Cp). Predicted dominant gain is in chamber-fill realism, not specifically spike taildown.

### 2.3. Pressure-dependent ignition temperature T_ign(p)

SPINBALL uses `T_ign(p)`, srm_1d uses a constant T_ign. Calibrated v0.7.0 value is 927 K.

**Impact on Hasegawa A:** Marginal. The ignition front passes any given cell at p ~ 0.5-5 MPa over a time interval shorter than the spike-taildown window. T_ign(p) variation is at most a few tens of K over this range. Localized in time to the ignition transition, not to the spike taildown.

**Action:** Deprioritized; not a spike-taildown solution.

## 3. What srm_1d does that SPINBALL doesn't

### 3.1. Pyrogen plenum thermodynamics (igniter chamber model)

srm_1d v0.7.0's `igniter_plenum.py` is a 0D ODE for a Saint-Robert pyrogen burning into a choked/subsonic vent. SPINBALL prescribes ṁ_ig(t) from experimental measurement and does not model the igniter chamber.

**Implication:** srm_1d can predict igniter behavior from BPNV / KNO3 / KNSU formulation data without needing an experimental ṁ_ig(t) trace. This is the *amateur-rocket convention* the project has adopted (see `feedback_igniter_conventions` memory).

### 3.2. Ma 2020 erosive burning model

srm_1d uses Haaland → Gnielinski → bisection. SPINBALL uses Lenoir-Robillard with Lawrence modification. Ma 2020 is more physically grounded (no fitted β) and was published 11 years after the SPINBALL formulation was frozen.

**Implication:** When the impingement-region erosive burn rate in Z23 reaches "same order of magnitude or larger than APN" (Cavallini 2009, p. 18), this is computed with Lenoir-Robillard. The corresponding number in srm_1d would be Ma 2020's prediction, which has lower epistemic uncertainty.

### 3.3. Goodman integral solid heating

Identical to Peretz 1973; both srm_1d v0.7.0 and Peretz 1973 use the cubic-polynomial integral method. SPINBALL likely uses a full PDE (not explicitly stated, but the SPIT lineage suggests this).

**Implication:** srm_1d is computationally cheaper here (~30× cheaper per ignition event) at ~5% accuracy cost — see [[equations_goodman_integral]]. This is a deliberate trade-off and not a deficiency.

## 4. Verdict — which post-v0.7.0 candidate to pursue first

The spike-taildown residual is the problem. Candidates from session-end discussion (see [[project_next_session_focus]]):

### Candidate A: Zeldovich-Novozhilov dynamic burn rate ★ PRIMARY

**Physical mechanism.** During the spike, P ramps faster than the propellant surface thermal layer can equilibrate. The classical Saint-Robert form `r = a·P^n` assumes quasi-equilibrium — the surface temperature has fully adjusted to the current pressure. In reality, during a fast pressure rise, the surface is colder than the quasi-equilibrium temperature for that pressure, so the instantaneous burning rate is **lower than `a·P^n`**. As P plateaus, the surface catches up; burning rate transient-overshoots and then settles. The result is a **sharper spike with a faster taildown** — exactly the asymmetry we need.

**SPINBALL position.** Discussed in Cavallini 2009 thesis §2.2.3 (Eq. 2.42, the Z-N formulation). **NOT implemented in SPINBALL.** Dismissed for VEGA scale as "minor effect at IT and tail-off." But:

- VEGA-class motors are L/D ~ 6-12 with long fill times (~100 ms) and slow pressure rise rates.
- Hasegawa A is L/D ~ 3-5 with a fast pressure rise (~5-20 ms to peak).
- The Z-N relaxation timescale `τ_ZN ~ α / r²` is **fixed by propellant properties**, not motor scale.
- For HTPB-AP at peak burn rate, `α ~ 1.5e-7 m²/s`, `r ~ 5-10 mm/s`, giving `τ_ZN ~ 1.5-6 ms` — **comparable to Hasegawa A's spike duration but ≪ VEGA's**.

So Cavallini's "minor effect" judgment is correct for VEGA and likely wrong for Hasegawa A.

**Implementation cost.** Moderate. One extra state per cell (the dynamic burning rate `r_dyn[i]`) plus a relaxation ODE coupling `r_dyn` to the steady-state `r_ss = a·P^n + r_erosive` over timescale `τ_ZN`. Numba-friendly. No fitted constants if `τ_ZN = κ·α/r²` with κ a documented O(1) constant from Novozhilov 1973 / Strand-Razdan-Strand 1986.

**Sources.**
- Novozhilov, B. V., *Nonstationary Combustion of Solid Propellants*, USSR Acad. Sci. Press, 1973 (translated to English).
- Strand, L. D., Razdan, M. K., and Strand, J. C., "Computer modeling of solid propellant ignition transients," AIAA J. 24(3), 1986.
- Zarko, V. E., "Nonstationary combustion of energetic materials," AIAA Comb. Conf., 1989.

**Predicted Hasegawa A impact.** Sharpens spike peak (may overshoot — re-calibration of erosion_coeff likely needed) and accelerates taildown. **Single highest-leverage physical addition for our specific residual.**

### Candidate B: Al2O3 two-phase thermal lag ★ SECONDARY (physically motivated)

**Physical mechanism.** Hasegawa Propellant 1 is 14% Al by mass. Al combustion at the surface produces ~3-10 μm Al2O3 particles that are entrained in the gas. These particles have a thermal mass and a heat-transfer relaxation time `τ_p ~ d²/α_p ~ 1-3 ms` — they leave the surface at gas temperature, then act as a transient heat sink as the gas above the surface heats up. During the spike, this caps the local gas temperature → flattens the spike → slower buildup but lower peak.

**SPINBALL position.** **NOT modeled** (single-phase explicit). Pardue 1992 has it; their two-fluid Shuttle SRM simulation reduces single-fluid error by ~50% on thrust but is mostly insensitive on head-end pressure (Pardue 1992 Fig. 9 — see [[extraction_peretz_pardue_cavallini]] §2.6).

**Implementation cost.** Moderate-high. One extra state per cell (T_p[i], the particle temperature) with relaxation `dT_p/dt = (T_gas − T_p)/τ_p`. Two-phase momentum coupling can be skipped if `u_p ≈ u_gas` (Stokes number small — for our particle size and bore velocity, this holds in QSS but may not during the spike when bore velocity is highest).

**Predicted Hasegawa A impact.** Reduces spike peak (currently +1.41% — could go to -1% or worse, helping or hurting depending on sign), changes taildown shape. Pardue 1992 explicitly notes that **two-phase mostly affects thrust, not pressure** — but Pardue's test was a large, low-Mach SRM. For Hasegawa A's high-velocity bore, the particle lag may matter more.

**Caution.** This is the candidate most likely to require *re-doing* the v0.7.0 LHS calibration from scratch. Z-N (Candidate A) is additive to existing physics; Al2O3 changes the underlying energy balance.

### Candidate C: Peretz-aligned participation fraction — LOW PRIORITY

**Physical mechanism.** Same Z-N physics viewed from the solid side. Naive Goodman gives a too-fast surface temperature ramp (~30 μs saturation); a surface-kinetics timescale from Z-N matches the right ~1-10 ms order. In a multi-cell ignition flame-spread context, this means a fraction of cells that have nominally ignited (T_surf > T_ign per Goodman) but haven't yet reached steady-state burn rate.

**SPINBALL position.** Not addressed explicitly. SPIT-lineage refs may have it.

**Implementation cost.** Low if Z-N is implemented (Candidate A); essentially free as a side effect. **Best treated as part of Candidate A, not a separate candidate.**

### Candidate D: Per-cell variable γ/Cp/R (SPINBALL's "infinite-gases") — TERTIARY

**Physical mechanism.** Igniter gas vs propellant gas have different (γ, Cp, R, k). During spike taildown the cell composition transitions from igniter-dominated to propellant-dominated. Currently v0.7.0 uses a single set of thermophysical properties, calibrated effectively (k ≈ 0.5 by LHS) as a compromise.

**SPINBALL position.** **Implemented as the headline change from SPIT.**

**Implementation cost.** Moderate. Two extra state arrays per cell (Y_ig[i], Y_prop[i] — passive scalars advected with the flow, with source terms from igniter and grain combustion). Cell γ and Cp computed as mixture-weighted averages.

**Predicted Hasegawa A impact.** Improves chamber-fill realism (igniter dominance transition) and may help spike-taildown indirectly. But the calibrated frozen-effective compromise in v0.7.0 already represents this physics approximately. Estimated leverage on spike-taildown residual: **lower than Candidates A and B**.

## 5. Recommendation

**Implement Z-N dynamic burn rate (Candidate A) first**, as a v0.7.0.x or v0.7.1 addition.

Reasoning:
- It is the candidate most specifically targeted at the spike-taildown problem (P-rise asymmetry → r_b asymmetry).
- It is additive: existing r = a·P^n + r_erosive remains as the steady-state target; Z-N adds a relaxation toward it.
- It has no fitted constants (τ_ZN derives from α/r²).
- The implementation is Numba-friendly.
- It is exactly the physics SPINBALL DOES NOT have, so we are not duplicating their work.

**Sequence after Z-N:**

1. **Z-N implementation** (v0.7.1): one cell-level state, one ODE, one new term in burn_rate.py.
2. **Re-calibrate Hasegawa A** LHS with the new model (3-7 var depending on whether κ is treated as free).
3. **If residual still present**: implement Al2O3 two-phase (Candidate B) and re-calibrate. This is the more invasive step.
4. **Variable γ/Cp/R per cell** (Candidate D): only if Z-N + (optionally) Al2O3 don't close the residual. Lower priority because it is more about chamber-fill than spike-taildown.

## 6. SPIT-lineage gap — when (and whether) to chase it

The SPIT igniter sub-model and ignition criterion equations are not in repo. They live in:

- Di Giacinto & Serraglia 2001 (AIAA 2001-3448) — "Modeling of Solid Motor Start-up"
- Serraglia 2003 PhD thesis — "Modeling and Numerical Simulation of IT of Large SRMs"
- Favini et al. 2005 (EUCASS) — "Pressuring Gas Effects on Pressure Oscillations during IT"

**Do we need these for spike-taildown?** No — see the verdict above. The dominant candidates (Z-N, Al2O3) come from Novozhilov/Pardue, not SPIT.

**Do we need these ever?** Only if we decide to port SPINBALL's variable-properties mixture formulation (Candidate D) or its multi-cell impingement region (irrelevant for Hasegawa A). Both are tertiary in our priority order.

**Action:** Defer. If/when Candidate D becomes priority, add SPIT acquisition to the agenda.

## 7. Links

- [[extraction_spinball_2009]] — full equation extraction from the 2009 SPINBALL conference paper
- [[extraction_spit_to_spinball_2008]] — SPIT-to-SPINBALL transition narrative + GREG details
- [[extraction_peretz_pardue_cavallini]] — earlier extraction with Peretz 1973 + Pardue 1992 + Cavallini 2009 thesis. Note the §3.1 6-species claim correction in [[extraction_spinball_2009]] §10.
- [[equations_goodman_integral]] — Goodman cubic-polynomial integral method (Peretz 1973 / srm_1d v0.7.0 solid-heating kernel)
- [[project_next_session_focus]] — session-end context that motivated this research
- [[project_hasegawa_calibration_state]] — v0.7.0 LHS rank-1 parameter set
