# AP/HTPB+Al Thermal Conductivity Literature Bounds

**Research date**: 2026-05-22 (during v0.7.1 Phase 5 LHS calibration)
**Authoring agent**: general-purpose subagent (haiku) delegated by the
v0.7.1 Phase 5 session. Distilled report from the agent's full
transcript; raw transcript was transient and is not preserved.

## Why this exists

The v0.7.1 Phase 5 N=1200 LHS sweep with literature-bounded heat
flux produced a top-5 ensemble whose fitted `k_solid` values
spanned **0.175 → 0.776 W/(m·K)** — a 4.4× range for a single
propellant composition (69 % AP / 17 % HTPB / 14 % Al, the
Hasegawa Motor A propellant). The user (correctly) objected to
the spread on physical-fidelity grounds and asked for literature
to pin a defensible range.

This document captures the open-literature answer.

## Published k_solid for AP/HTPB(+Al) composites

| Component / formulation                            | k_solid [W/(m·K)]   |
|----------------------------------------------------|---------------------|
| Neat HTPB binder (Rajoriya et al.)                 | ~0.199              |
| AP/HTPB composites (60-70 % AP / 17 % HTPB)        | 0.25 - 0.35         |
| Hasegawa A class (69 AP / 17 HTPB / 14 Al)         | 0.25 - 0.35 (est.)  |
| Thakre & Yang baseline (calibrated fit)            | 0.297               |
| Graphene-silver nano-enhanced HTPB (ceiling)       | ~0.30               |
| **Conservative LHS upper for AP/HTPB**             | **0.40**            |

### Aluminum contribution

14 % aluminum loading in Hasegawa A does **NOT** measurably raise
k_solid above the binder + AP baseline. Bulk aluminum k is
~237 W/(m·K) but dispersed Al particles in a polymer matrix
follow near-rule-of-mixtures damping. One nanocomposite study
saw only ~52 % k improvement even with full surface treatment
of the additives, bringing the matrix from ~0.20 to ~0.30
W/(m·K). 14 % mass-loaded Al particles contribute ≤ 0.05 W/(m·K)
above the unaluminized composite baseline.

## Temperature dependence: physically plausible, not quantified

Open literature **does not provide** quantitative k_solid(T)
functions or tabulated k(T) data for AP/HTPB composites between
300 K and ~1000 K. What is documented:

- AP undergoes orthorhombic → cubic phase transition at ~513 K
  (~240 °C). Lattice rearrangement should produce measurable
  thermal property changes, but they have not been systematically
  measured for the composite.
- AP decomposition onset effects above ~650 K documented for
  RP-1 liquid propellant, NOT for AP/HTPB composites.

**Implication for srm_1d**: implementing a temperature-dependent
k_solid(T) function would require extrapolation from generic
polymer behavior, not direct measurement. Treating k_solid as a
scalar is defensible until/unless T-dependent data exists. Adding
T-dependence as a "fit a polynomial to make traces match" exercise
is the kind of unfounded smoothing the project explicitly avoids
(see `feedback_no_unfounded_smoothing` memory).

## "Effective" vs "measured" k_solid

The ballistic-modeling literature routinely **uses k_solid as a
fit parameter** to match experimental ignition transients
(Broyden-Fletcher-Goldfarb-Shanno parameter fitting in
Thakre & Yang; sensitivity demos in EUCASS 2019). This is
documented standard practice — srm_1d is not unique in this
respect.

**Critical observation**: every fitted value found in the surveyed
literature clusters within the **0.20-0.40 W/(m·K) physical band**.
No source supports k_solid > 0.6 W/(m·K) for AP/HTPB-class
composites except as a fit parameter compensating for absent
physics. The v0.7.1 Phase 5 LHS rank-4 (k_solid = 0.776) and
similar high-k entries are therefore *physically pathological*
and should be excluded from a calibration that aims to be
honestly bounded.

## What the literature implies about srm_1d's spike-vs-tail tension

User-observed tension: rank-1 (low-k 0.175) fits tail well,
spike acceptably; rank-4 (high-k 0.776) fits spike well, tail
poorly. Re-ranking under priority-aligned weights produces
basin shifts but doesn't resolve the conflict.

The k_solid literature establishes that **this tension cannot be
resolved by k_solid alone** when k_solid is bounded to the
physical 0.20-0.40 range. Predicted outcomes from a re-sweep
with `k_solid ∈ (0.26, 0.32)`:

- **Fit holds tightly**: calibration is honest; the segment
  tradeoff is intrinsic to the current model — pick the entry
  that best matches user's stated peak-P + erosive priorities.
- **Fit degrades sharply**: model gap confirmed. The right next
  step is NOT to relax k_solid back. The right step is to
  evaluate v0.7.2/v0.7.3 candidates:
  1. **Al2O3 two-phase thermal lag** (Pardue 1992) — molten
     oxide condensation extends tail-off pressure decay.
     Already on the v0.7.3 candidate list.
  2. **Z-N dynamic burn rate** — relaxation ODE on steady r_b;
     SPINBALL walkthrough's leading spike-taildown candidate.
     Slotted as v0.7.2 candidate.
  3. **Temperature-dependent Cps** — AP phase changes around
     510 K should produce a Cps step we currently treat as
     scalar.
  4. **Better throat erosion physics** — late-stage P decay
     rate depends sensitively on throat-area evolution; we
     treat erosion as a calibrated scalar (`erosion_coeff`).

## Recommended LHS bounds

- **Defensible (broad)**: `k_solid ∈ (0.20, 0.40)`. Captures the
  full literature range. Sensible default for a generic AP/HTPB
  calibration sweep.
- **Literature-center (tight)**: `k_solid ∈ (0.26, 0.32)`. The
  central AP/HTPB+Al estimate band. Use this for Hasegawa A
  specifically when the goal is "calibrate honestly within
  physical truth."

The next Phase 5 re-sweep should adopt the **literature-center
range** to honestly test whether the spike-vs-tail tension is
fundamentally a k_solid issue (it isn't) or a missing-physics
issue (likely yes).

## Primary literature

1. **Thakre & Yang** — "A Model of AP/HTPB Composite Propellant
   Combustion in Rocket-Motor Environments." Uses
   k_solid = 0.297 W/(m·K) as the baseline ignition-modeling
   value for AP/HTPB. References parameter calibration as
   standard practice.
   https://www.researchgate.net/profile/Piyush-Thakre-2/publication/237293623_A_Model_of_APHTPB_Composite_Propellant_Combustion_in_Rocket_Motor_Environments

2. **AIAA classic** — "Thermal Conductivity — A Parameter in
   Solid Propellant Burning." Establishes k_solid sensitivity
   in burn-rate transients. Canonical older reference.
   https://arc.aiaa.org/doi/abs/10.2514/3.3743

3. **EUCASS 2019** — "Ignition Study at Small-Scale Solid Rocket
   Motor." Recent transient-ignition model; explicit k_solid
   parametric-sensitivity demo.
   https://www.eucass.eu/doi/EUCASS2019-0757.pdf

4. **AIAA 1972** — "Thermal Diffusivity of Ammonium Perchlorate."
   Direct measurement of AP-only thermal properties.
   https://arc.aiaa.org/doi/abs/10.2514/3.3505

5. **ResearchGate** — "Enhanced thermal conductivity of HTPB
   composites by graphene-silver hybrid." Establishes neat HTPB
   baseline ~0.199 and nano-enhancement ceiling ~0.30.
   https://www.researchgate.net/publication/348924610_Enhanced_the_thermal_conductivity_of_hydroxyl_terminated_polybutadiene_HTPB_composites_by_graphene-silver_hybrid

6. **ScienceDirect** — "Thermal conductivity estimation of high
   solid loading particulate composites: A numerical approach."
   Numerical framework; references Rajoriya et al. baselines.
   https://www.sciencedirect.com/science/article/abs/pii/S1290072917310244

7. **ACS / Energy & Fuels** — "Effect of RP-1 Compositional
   Variability on Thermal Conductivity at High Temperatures
   and High Pressures." Bounds compositional-variability effects
   to ~10 %, not 3× ranges.
   https://pubs.acs.org/doi/abs/10.1021/ef900435b

8. **AP phase transition reference**:
   https://www.sciencedirect.com/science/article/abs/pii/0167577X95001298

## Conclusion

> The v0.7.1 Phase 5 LHS k_solid spread of 0.175-0.776 W/(m·K)
> is **3× wider than the published physical bounds** for AP/HTPB+Al
> composite propellants (0.20-0.40 W/(m·K), centered 0.25-0.30).
> The high-k LHS basin entries are not physically defensible
> calibrations — they are k_solid acting as a free parameter
> absorbing errors from absent transient mechanisms.
>
> Recommended next step: re-sweep with `k_solid ∈ (0.26, 0.32)`,
> per the literature center. If the fit degrades sharply, the gap
> is missing physics (Al2O3 lag, Z-N dynamic burn rate,
> temperature-dependent Cps) and the right response is to slate
> those candidates for v0.7.2/v0.7.3 rather than re-widening
> k_solid.

## Related memory

- `[[srm-1d-ap-htpb-k-solid-literature-bounds]]` — quick-lookup
  memory equivalent.
- `[[srm-1d-pyrogen-heat-flux-literature-bounds]]` — sibling
  heat-flux lit dive.
- `[[srm-1d-v0-7-1-progress-state]]` — Phase 5 calibration arc.
- `[[project_v0_7_1_post_phase35_trace_assessment]]` — user's
  validation priorities that motivate the segment-weight
  discussion.
- `[[project_spinball_research_state]]` — Z-N dynamic burn rate
  candidacy flagged here as the natural follow-up if a tight
  k_solid sweep degrades the fit.
- `[[feedback_no_unfounded_smoothing]]` — relevant if a temporary
  k_solid(T) polynomial is tempting; literature doesn't support
  fitting one.
