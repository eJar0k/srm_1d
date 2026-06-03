# Spatial Ignition-Front Coupling — Extended Literature Digest

*Subagent literature dive, 2026-05-23. Extends what srm_1d already had
documented (`docs/v0_7_0/references/equations_goodman_integral.md` +
the 18-variant radiation-collapse audit in
`docs/v0_7_0/audits/2026-05-21_ignition_tuning_audit.md`).*

## What srm_1d already has documented

- **Per-cell Goodman cubic-polynomial integral method** (Peretz 1973
  §III.G, Eqs. III-44 to III-50) — each unignited bore cell holds
  its own scalar `(δ, T_surf)`, advanced by a 4th-order RK ODE. **No
  spatial coupling between adjacent cells.**
- **Convective h_c** is computed locally from a Bartz/Dittus-Boelter-
  style closure on the cell's own (ρ, u, T) state; "blowing"
  suppression goes h → 0 after ignition.
- **Radiation** (per-cell ε·σ·T⁴ emitter + 3-cell buffer) is wired
  and tested, but explicitly noted as *non-load-bearing for spike
  shape* in the 18-variant Cartesian audit (Δ peak < 0.05 MPa across
  ε ∈ {0, 0.45}).
- **Pyrogen plenum** (v0.7.0+) injects mass + sensible enthalpy +
  axial momentum into cell 0 only; no axial distribution of the
  impingement footprint. The Phase 5 cross-motor signature (Hasegawa
  A under-fires while Zerox/BALLSStick/Chunc over-fire 1.7-5×) is
  fully consistent with a simultaneous-bore-ignition artifact rather
  than a transport-property issue.

## New findings from this dive

**1. Salita "Modern SRM ignition transient modeling" (AIAA
2001-3443 / -3447, parts I and V)** is the most consequential
reference srm_1d's docs don't already cite. Salita explicitly
catalogues the post-Peretz lineage: the major industrial 1D codes
(SPP/SAIB at SEA, Tetra Research's solid-motor ignition tool) all
converged on the same architecture — Peretz Q1D gas dynamics +
Goodman per-cell solid + an **additional axial impingement-region
sub-model** that *redistributes* the prescribed igniter mdot/heat
flux over multiple head-end cells rather than dumping it all in cell
0. The flame spread is still emergent from per-cell T_surf > T_ign
crossings, but the head-end heat-flux is no longer a delta function
at x=0.

**2. The Sapienza SPIT/SPINBALL impingement-region sub-model**
(Cavallini 2009, Di Giacinto / Favini / Serraglia heritage, also
EUCASS proceedings) is the canonical formulation: a **separate
control volume** spanning the geometrically calculated jet footprint
distributes igniter mass + enthalpy across a contiguous range of
bore cells with an axial weighting derived from jet cone angle /
radial impingement geometry, plus an **enhanced convective h_c** in
that band (the "erosive contribution" Figure 7.9b in Cavallini 2009).
Cells outside the band see only standard mass-flux-driven h_c. This
is the only published spatial mechanism that demonstrably reproduces
sequential bore ignition without an empirical flame-speed law.

**3. Convective flame-spread theory (Glick, Most, Kashiwagi
1971-1982; Combust. Flame & Combust. Sci. Tech.)** treats the
propagating front as a **turbulent-boundary-layer-coupled ignition
problem**: the local h_c at an unignited cell is *augmented* by the
mass addition from all *upstream burning* cells, because mass flux
G(x) = ∫₀ˣ ρ_p·r_b·P_b dx' grows along the bore. Kashiwagi 1982
(CST 28, "Flame-Spreading over the Surface of a Solid Propellant,
Part II: Simplified Model") gives a closed-form spread velocity
V_fs = f(G, P, Y_ox) using a critical-T criterion identical to
Peretz's — but with h_c(x,t) that scales with cumulative upstream
G. **This is the structural fix for the simultaneous-ignition
artifact**: replace local-state h_c with one that integrates
upstream mass addition.

**4. Han & Cai 2017, "Numerical Modeling and Studies of Ignition
Transients in End-Burning-Grain Solid Rocket Motors" (JPP, DOI
10.2514/1.B36024)** is the cleanest published demonstration of
cell-coupled flame spreading in a 1D code. They show that *without*
enhanced-h_c coupling in the impingement region, the ignition is
essentially instantaneous and the spike unphysical — the same
artifact srm_1d exhibits. Their fix: a piecewise axial h_c
distribution from a prescribed igniter-jet profile.

**5. Peretz follow-ups — post-1973 lineage.** The Frazer-Hicks
heat-flux-integral criterion (an ∫q²dt threshold) and the
Vilyunov-Zarko gas-phase-Arrhenius criterion are the two main
alternatives that emerged but are *not* dominantly better than
Peretz's critical-T criterion for engineering predictions per
Cavallini's own assessment. Several authors (Johnston AIAA-1995,
"Solid rocket motor internal flow during ignition"; Yıldız et al.
recent ScienceDirect 2025) all retain Peretz's critical-T criterion
and instead invest in better gas-side h_c distributions. **The
literature consensus is that the per-cell criterion is fine; the
load-bearing physics gap is the axial heat-flux distribution.**

## Implementation guidance for srm_1d

The lowest-footprint, most-defensible coupling mechanism for srm_1d's
existing Goodman kernel is **mass-addition-augmented convective h_c**,
in the Kashiwagi 1982 / Han 2017 style. Concretely: where srm_1d
currently computes `h_c[i]` from local cell-i state, replace with
`h_c[i] = h_c_local(state[i]) · f_blow(G_cum[i] / G_ref)` where
`G_cum[i] = G_igniter + Σ_{j<i, burning} ρ_p·r_b[j]·P_b[j]·Δx[j] / A_p[i]`
is the cumulative upstream mass flux.

This requires **one extra cumulative-sum pass per timestep** before
the Goodman ODE step — Numba-friendly, no new fitted constants if
`f_blow` is taken as Dittus-Boelter (`∝ Re^0.8 ∝ G^0.8`). Cells far
from the igniter naturally see weaker h_c at t=0 (only the pyrogen
contributes to G_cum) and ramp up as upstream cells light, producing
sequential propagation.

The complementary higher-effort enhancement is a **multi-cell pyrogen
footprint** (Sapienza impingement-region style): change cell-0-only
pyrogen injection to a Gaussian or top-hat axial weighting over the
first N cells set by jet cone geometry; this fixes the head-end
concentration without touching the Goodman solver. Doing both
together is the published architecture in SPIT/SPINBALL and Han 2017.

## Known limitations

- Kashiwagi-style empirical V_fs correlations are calibrated on
  flat-strand experiments, not motor bores; the constants will need
  re-fitting against Hasegawa A.
- Cumulative-G coupling is not strictly physical when the boundary
  layer separates / recirculates (sudden expansions, fin slots —
  Unnikrishnan 2001 / Mukunda-Paul JPP 2007 case).
- Mass-addition h_c augmentation will interact with srm_1d's
  existing Ma 2020 erosive burning term (which also uses Re_local
  from G) — care needed not to double-count enhancement.
- Han 2017 explicitly notes order-of-magnitude sensitivity of
  ignition-transient shape to the prescribed h_c axial profile —
  i.e. this approach trades one calibration knob (k_solid) for
  another (axial profile shape). The benefit is the new knob is
  *physically anchored* to igniter geometry rather than a fudge.
- None of the cited sources have published *open* implementations;
  all are proprietary (SPP, SAIB, Tetra) or research codes
  (SPIT/SPINBALL).

## Primary citations (with URLs)

1. Salita, M. (2001). "Modern SRM ignition transient modeling. I —
   Introduction and physical models." AIAA 2001-3443.
   https://arc.aiaa.org/doi/10.2514/6.2001-3443
2. Salita, M. (2001). "Modern SRM ignition transient modeling. V —
   Prospective developments in CFD simulation." AIAA 2001-3447.
   https://www.researchgate.net/publication/269216223_Modem_SRM_ignition_transient_modeling_V_-_Prospective_developments_in_CFD_simulation
3. Han, S. & Cai, W. (2017). "Numerical Modeling and Studies of
   Ignition Transients in End-Burning-Grain Solid Rocket Motors."
   J. Propulsion & Power. https://doi.org/10.2514/1.B36024
4. Kashiwagi, T. (1982). "Flame-Spreading over the Surface of a
   Solid Propellant, Part II: Simplified Model." Combustion Science
   & Technology 28(1-2).
   https://www.tandfonline.com/doi/abs/10.1080/00102208208952536
5. Unnikrishnan, C. et al. (2001). "Effect of flame spread mechanism
   on starting transients of solid rocket motors." AIAA 2001-3854.
   https://arc.aiaa.org/doi/abs/10.2514/6.2001-3854
6. Johnston, W. A. (1995). "Solid rocket motor internal flow during
   ignition." J. Propulsion & Power.
   https://arc.aiaa.org/doi/abs/10.2514/3.23869
7. Mukunda, H. S. & Paul, P. J. — "Flame-Spreading Process in Solid
   Rocket Motor with Fin Slots." J. Propulsion & Power.
   https://arc.aiaa.org/doi/10.2514/1.B40101
8. Yıldız, B. et al. (2025). "Parametric study of igniter design on
   ignition transient performance in solid rocket motors."
   https://www.sciencedirect.com/science/article/abs/pii/S1290072925006453
