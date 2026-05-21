# Cavallini, Favini, Di Giacinto, Serraglia (2009) — AIAA 2009-5512

**Full title:** "SRM Internal Ballistic Numerical Simulation by SPINBALL Model"

**Venue:** 45th AIAA/ASME/SAE/ASEE Joint Propulsion Conference, 2-5 August 2009, Denver, Colorado

**PDF:** `srm_1d/docs/post_v0_7_0/cavallini2009.pdf` (22 pp)

**Companion to:** the 2009 Cavallini PhD thesis (74323997.pdf, ~330 pp, already in repo). Where the 2008 paper announced the SPIT→SPINBALL transition, this is the **canonical equation-set publication for SPINBALL**.

---

## 1. The SPINBALL Q1D flowfield equations (Eq. 1, p. 4) — VERBATIM

```
∂(ρ·A_p)/∂t + ∂(ρ·u·A_p)/∂x = r_b·P_b·ρ_p + (ṁ_s·A_p)/V + (ṁ_ig·A_p)/V          (mass)

∂(ρ·u·A_p)/∂t + ∂[(ρ·u² + p)·A_p]/∂x − p·(∂A_p/∂x) =
                            (ṁ_ig·A_p·v̄_inj)/V + (1/2)·ρ·u²·c_f                  (momentum)

∂(ρ·E·A_p)/∂t + ∂[(ρ·E + p)·u·A_p]/∂x =
                r_b·P_b·ρ_p·H_f + (ṁ_ig·A_p·H_ig)/V + (ṁ_s·A_p·H_s)/V            (energy)
```

Closed by perfect gas EOS. ṁ_s = ablation/cavity source mass; ṁ_ig = igniter mass flow; v̄_inj = igniter combustion-product axial injection velocity (modern improvement over Peretz's mass-only injection). H_f, H_ig, H_s = propellant, igniter, and source total enthalpies.

**Critical observation:** this is **ONE mass equation**, not the six-species system that the 2009 thesis description suggested. The earlier extraction at `srm_1d/docs/v0_7_0/references/extraction_peretz_pardue_cavallini.md` §3.1 has this wrong and should be corrected. The "infinite-gases" formulation is **not** about transporting species — it is about per-cell variable thermophysical state evaluated from a weighted mixture of inflows, as documented below.

## 2. SPIT vs. SPINBALL formulation contrast (Eq. 2, p. 4)

**SPIT — 3-gas model (Eq. 2):**

```
∂(ρ_pr·A_p)/∂t + ∂(ρ_pr·u·A_p)/∂x = r_b·P_b·ρ_p + ṁ_{s,pr}·A_p/V
∂(ρ_ig·A_p)/∂t + ∂(ρ_ig·u·A_p)/∂x = ṁ_{s,ig}·A_p/V + ṁ_ig·A_p/V
∂(ρ_in·A_p)/∂t + ∂(ρ_in·u·A_p)/∂x = ṁ_{s,in}·A_p/V

∂(ρu·A_p)/∂t + ... = ṁ_ig·A_p·v̄_inj/V + (1/2)ρu²·c_f
∂(ρE·A_p)/∂t + ... = r_b·P_b·ρ_p·H_f + ṁ_ig·A_p·H_ig/V + ṁ_s·A_p·H_s/V
```

Three species — propellant gas (pr), igniter gas (ig), pressurizing gas (in) — each with its own mass conservation equation. Each species has **constant** thermophysical properties.

**SPINBALL — infinite-gases mixture (Eq. 1, repeated):**

Single total-mass equation. No species transport. Thermophysical properties (γ, Cp, R, MW) are **per-cell, per-timestep**, evaluated as a weighted average of:

> "the mass fluxes of the mixture coming from adjacent cells, located at the left and right of the considered one, and from the sources terms from the grain combustion reactions, the igniter and the cavity model" (p. 4)

Weights = concentrations in the cell at the current step. The chemical equilibrium tables of the propellant supply pressure-dependent properties for the source-term gas; combustion efficiency η_c* (a global scaling factor) modifies these.

**Why the switch:** SPIT-style 3-species transport with constant per-species properties **cannot represent**:
- pressure variation of propellant-gas properties (real-equilibrium chemistry → γ(p), Cp(p), MW(p))
- TP-ablation-product mixing
- combustion-inefficiency-driven offset from ideal equilibrium

SPINBALL achieves the *physics result* of multi-species transport (locally-correct thermodynamics) without the *cost* of advecting many species. Per-cell mixture state is updated from the weighted-average operation each timestep.

## 3. Burning rate model (Eq. 3, p. 5) — VERBATIM

```
r_b = a(T_i)·(p/p_ref)^n   +   k_eb·h_c · exp[-β(D_h, P)·r_b·ρ_p / (ρu)]
       └─────APN─────┘         └────────────erosive (Lenoir-Robillard)────────────┘
```

**Identical to SPIT, identical to Peretz 1973.** Lenoir-Robillard 1957 erosive form with Lawrence 1968 hydraulic-diameter modification.

**Comparison to srm_1d:** srm_1d v0.7.0 uses **Ma 2020** (Haaland → Gnielinski → bisection with a friction-factor closure on a roughness Reynolds number). Lenoir-Robillard is a 1957-vintage semi-empirical correlation; Ma 2020 has more physics. **srm_1d is ahead of SPINBALL on erosive burning physics.** Both share the additive `r = r_APN + r_erosive` split.

> "Note that, certainly, some difficulties for the quasi steady state 0D model is related ... while instead, even a rough calibration of the Q1D unsteady model gives, **as strong heritage of the SPIT model**, a good accordance with the experimental data during the ignition transient." (p. 15)

The IT fidelity is credited to **SPIT lineage**, not SPINBALL-era changes. So whatever IT physics SPINBALL has, it inherited from the 2001-2003 SPIT generation. (See [[extraction_spit_to_spinball_2008]] §7 for the SPIT-lineage paper trail we'd need to chase if we ever want those equations.)

## 4. Numerical method (p. 4)

> "The numerical discretization method is represented by a finite volume Godunov scheme **first order accurate in space and time**, coupled with an exact Riemann solver. The time discretization is characterized by a CFL-like condition that guarantees the stability of the numerical scheme adopted." (p. 4)

First-order in space and time **for the 2009 SPINBALL paper specifically**. The 2008 paper claimed first/second-order; 2009 paper reports first-order in actual runs. ENO + MinMod + Heun second-order is in GREG (Eq. 10 area) but not necessarily activated for the flowfield Z23 results shown.

## 5. 0DQSS sub-model for SFT reconstruction (§II.C)

A separate 0D quasi-steady model used to **pre-extract calibration parameters** from static-firing-test (SFT) data:
- η_c* (combustion efficiency)
- η_cF (nozzle efficiency, when thrust is available)
- HUMP law h(web) — local r_b enhancement vs. nominal
- nozzle throat area evolution A_t(t)

The 0DQSS algorithm (Eqs. 4-9) is a closed-form inverse problem given experimental ˜p(t), ˜F(t), and a monotonic throat-erosion model. SPINBALL uses these calibrated parameters as **inputs**, not free parameters.

**Relevance to srm_1d:** srm_1d does not have an SFT reconstruction tool. Calibration in srm_1d v0.7.0 is via LHS sweep against the experimental trace directly (see `srm_1d/tools/sensitivity.py`, Hasegawa A LHS). Cavallini's approach decouples calibration from the dynamic simulation; srm_1d's approach couples them. Both are legitimate.

## 6. Coupling GREG↔SPINBALL (Eq. 11-12, p. 7-8)

Off-line coupling in the 2009 results. GREG runs once with constant r_b → tables of {A_p, P_b, P_w} vs. web. SPINBALL interpolates these tables at runtime.

```
P_b,j = (1/Δx) · Σ_{i(j)} S_b,i(j)        (Eq. 11)
ṁ_b,j = ρ_p · r_b,j · P_b,j · Δx          (Eq. 12)
```

On-line coupling (full bidirectional) is discussed as future work — the 2009 paper acknowledges that off-line is only valid when r_b is "slowly varying" in space, i.e. low-erosive QSS.

**Relevance to srm_1d:** srm_1d's FmmGrain port from openMotor uses lookup-table coupling (constant-velocity geometry tables indexed by web). Functionally analogous to GREG off-line. **srm_1d does not currently do on-line bidirectional grain↔flowfield coupling.** For Hasegawa A this is not needed (small motor, simple geometry).

## 7. Z23 results (§III) — what the paper actually demonstrates

**Z23 motor:** Vega second stage, HTPB1912, 23 t propellant, finocyl, submerged nozzle, ~6 m long, t_b ≈ 77 s, MEOP ≈ 106 bar.

**Run 1 — no calibration** (Fig. 10): linear throat-area evolution, no HUMP, no η_c*. SPINBALL HEP and 0DQSS HEP both significantly below experimental. **Critical detail:** "the interesting fact consists, however, in the very small displacement for the prevision of the two different models" — that is, Q1D and 0DQSS give nearly the same answer at QSS without calibration.

**Run 1 IT detail** (Fig. 11): even uncalibrated, SPINBALL matches experimental HEP during IT to plotting accuracy — credited to SPIT heritage.

**Run 2 — full calibration** (Fig. 12): uses HUMP(web), η_c*, A_t(t) from the 0DQSS SFT reconstruction. Now SPINBALL HEP matches experimental closely over the entire combustion time.

**Critical findings for spike-taildown:**

a. Z23 has very low ρu in the bore → **erosive burning contribution is "totally negligible" during QSS** (p. 16). The Z23 spike taildown is therefore *not driven by erosive burning evolution* in their case. For Hasegawa A which has high ρu, this is the opposite limit.

b. **IT spike is captured well; the remaining HEP disagreement is in the QSS** (small total pressure drops not in 0DQSS — Fig. 14, 15). The QSS-trace residual is ~0.5% of instantaneous ṁ — much smaller than srm_1d's spike-taildown problem.

c. The paper does **not report** any Z23 IT post-spike taildown residual. The IT they show looks good and they move on.

d. **IT geometry observation** (Fig. 16, p. 18-19):

> "the flame spreading and the igniter jets impinging directly on the motor burning surface have some effects in terms of non-uniform regression of the motor grain geometry, underlined in the web field curves in time during the IT. This phenomenon is located in particular, in the plume of the igniter radial nozzles jets, where velocities are high and in the stagnation point of the jet on the propellant grain surface, so that the burning rate is high, **with an erosive contribution with the same order of magnitude, or even predominant, with respect to the APN term**."

For Z23 IT, the impingement region exhibits **erosive-dominant** burn rate locally where igniter jets hit. This is a multi-cell spatial phenomenon driven by igniter jet geometry. For Hasegawa A's head-end pyrogen-into-cell-0 setup, no analogous mechanism exists.

## 8. What this paper confirms is NOT in SPINBALL

- **Dynamic burning (Zeldovich-Novozhilov):** not mentioned anywhere; absent.
- **Two-phase / particles:** explicit assumption "non-reactive mixture of perfect gases" (p. 4). Single phase.
- **Solid-phase explicit equations:** referred to SPIT lineage (refs 6-13).
- **Igniter sub-model details:** referred to SPIT lineage. SPINBALL adds the v̄_inj momentum term to the momentum equation but the igniter chamber thermodynamics (e.g. pyrogen plenum) are NOT modeled — ṁ_ig(t) is prescribed.
- **Ignition criterion equations:** referred to SPIT lineage.

## 9. What's authoritative in this paper that nothing else gives us

1. **The actual Q1D conservation system used in SPINBALL** (Eq. 1) — this is the SPINBALL equation set, in concise form, vetted against Z23. The thesis description was either compressed or differs from the published equations.
2. **The SPIT-to-SPINBALL formulation reasoning** (Eq. 1 vs. Eq. 2 contrast). Establishes that the gain is per-cell variable thermophysical state, not species transport.
3. **Z23 entire-burn quantitative results.** No prior reference has these in this form.
4. **The 0DQSS-as-precalibration-tool architecture.** Conceptually different from srm_1d's coupled-LHS calibration.

## 10. Cross-reference correction

The earlier extraction at [extraction_peretz_pardue_cavallini.md](../../v0_7_0/references/extraction_peretz_pardue_cavallini.md) §3.1 states:

> "Six mass-conservation equations track six gas species independently."

This is **incorrect** per Eq. 1 of the 2009 conference paper. The correct statement is:

> SPIT has 3 species mass equations (propellant, igniter, pressurizing gas). SPINBALL has 1 total mass equation plus per-cell variable thermophysical properties evaluated from a mixture-weighted average of inflows + sources.

The thesis text the original extraction was based on likely described variable-properties tables indexed by something like 6 species but the *transport equations* are not 6 — they are 3 (SPIT) or 1 (SPINBALL). This correction is captured in [[spinball_walkthrough]] and should be noted in the v0.7.0 extraction file when convenient (low priority; that file is frozen as a v0.7.0 artifact).
