# Z-N Dynamic Burn Rate — Extended Literature Digest

*Subagent literature dive, 2026-05-23. Extends what srm_1d already had
documented (`docs/post_v0_7_0/references/spinball_walkthrough.md` and the
Cavallini 2009 + Novozhilov 1973 + Strand 1986 + Zarko 1989 chain
extracted in `extraction_spinball_2009.md`).*

## What srm_1d already has documented

The `spinball_walkthrough.md` recommended Z-N as the primary spike-
taildown candidate with the form `τ_ZN = κ·α/r²` and an order-of-
magnitude estimate (`τ_ZN ≈ 1.5–6 ms` for HTPB/AP at peak burn,
comparable to Hasegawa A's spike width). The SPINBALL extractions
(Cavallini 2009 thesis §2.2.3) confirm Z-N was discussed but
*dismissed* for VEGA-scale motors, citing Novozhilov 1973,
Strand-Razdan-Strand 1986, and Zarko 1989. The post-v0.7.1 README
lists Z-N as queued v0.7.2 work.

## New findings from this dive

**1. Greatrix 2008 is the canonical 1D-IB implementation reference.**
"Transient Burning Rate Model for Solid Rocket Motor Internal Ballistic
Simulations" (Int. J. Aero. Eng., 2008, paper 826070) is the single
most directly applicable paper to srm_1d's situation. Greatrix builds
a general Z-N model where the *integrated solid-phase temperature
distribution* (not the surface gradient) drives instantaneous r_b.
Critically, he validates on **ignition pressure spikes in low-L\***
motors and finds dynamic burning is the *dominant* mechanism producing
the spike — exactly the regime where the QSS Saint-Robert prediction
"misses the pressure spike completely." This matches the structural
artifact in srm_1d's Hasegawa A / Zerox / BALLSstick cross-motor
survey (over-firing 1.7-5× at default knobs).

**2. The model is robust enough for production internal-ballistics
codes.** Multiple groups have wired Z-N into CFD-coupled IB
simulators with reported transient-pressure accuracy within 5% of
experimental measurement, including recent 2025 ScienceDirect work
on overload-condition instability. SEA's SPP'12 (the JANNAF reference
code) includes an "IGT — 1-D time dependent true ignition transient
module" with Z-N-class dynamic burn rate, confirming this is
industry-validated, not just academic.

**3. A 2023 MDPI Aerospace paper revisits the formulation** ("A
Phenomenological Model for the Unsteady Combustion of Solid
Propellants from a Zel'dovich-Novzhilov Approach," 10(9), 767). The
authors derive explicit linear and lowest-order quadratic nonlinear
closures for the propellant response, both rooted in the Z-N
quasi-steady approximation — useful if srm_1d wants a fully
closed-form alternative to the integral approach.

**4. Alternative formulations have been surveyed and are mostly
inferior for this use case.** Zarko & Gusachenko's 2010 "Critical
Review of Phenomenological Models for Studying Transient Combustion
of Solid Propellants" (Int. J. Spray and Combust. Dynamics, 2(2),
151-167) concludes Zeldovich's 1942 phenomenological framework
remains the most defensible foundation; many proposed
generalizations rely on assumptions that narrow applicability.
**Levine-Culick (QSHOD response-function form), Beckstead-Derr-Price
(BDP, kinetics-based), and Krier-Tien generalizations** all exist,
but Zarko/Gusachenko flag them as either limited to narrow
operating windows or built on questionable assumptions. There is no
strong contender other than Z-N for the time-domain transient case
(BDP is more appropriate for frequency-domain combustion-instability
work).

**5. Greatrix explicitly reports a "burning rate limiting function"
for numerical stability** — this is the key engineering detail. The
relaxation ODE can become stiff under fast P-rise; Greatrix imposes
a cap/limiter to keep the integration stable. This is reported as
both a numerical-stability device AND a calibration knob aligned
with observed combustion-response data.

## Implementation guidance for srm_1d

**Functional form**: treat the steady-state burn rate
`r₀ = a·P^n + r_erosive(Ma 2020)` as a *target* and relax `r_b`
toward it via an ODE per cell. The Greatrix integrated-temperature
form is most defensible; the simpler first-order relaxation
`dr_b/dt = (r₀ − r_b)/τ_ZN` with `τ_ZN = κ·α/r_b²` and κ = O(1) is
the cheap entry point. Use `α ≈ 1.0–2.0 × 10⁻⁷ m²/s` for AP/HTPB+Al,
consistent with the existing `ap_htpb_k_solid_literature.md` bounds
(k_solid 0.20-0.40, ρ ≈ 1700 kg/m³, Cp ≈ 1500 J/(kg·K)). For
Hasegawa A at r_b ≈ 7 mm/s, this gives τ_ZN ≈ 2-4 ms — exactly the
spike duration.

**Numerics**: implement Greatrix's limiter from the start; pair it
with srm_1d's source-CFL cap. If stiffness shows up (likely during
the first 1-2 ms when r_b is increasing fastest), sub-step the Z-N
ODE inside one PISO step rather than implicit-coupling — the
per-cell ODE is decoupled, so split-operator is clean. Numba-friendly.

## Known limitations

- **Heterogeneity blind spot**: Z-N treats the propellant as
  homogeneous; AP/HTPB's bimodal/trimodal AP distribution produces
  local hot-spot effects Z-N can't represent. Acceptable for global
  P(t), but don't expect it to capture small-scale local kinetics.
- **Over-correction risk on slow-rise / long-L\* motors**: Greatrix
  and SPINBALL both note Z-N's effect is "minor" for VEGA-scale
  motors with τ_rise ≫ τ_ZN. Won't hurt these motors, but the LHS
  may converge to nearly the QSS solution — no spike-taildown gain
  to harvest there.
- **Initial-condition sensitivity**: choice of `r_b(t=0)` matters.
  Igniter-side initialization (r_b ≈ 0 at unignited cells) is
  correct but means the relaxation timescale during first ignition
  can be longer than τ_ZN suggests.
- **Pressure-coupled response only, not strain-rate**. Erosive-burning
  strain-rate coupling exists (some 2010s work extends Z-N with a
  strain-rate term in solid-phase energy) but is not standard.
- **Calibration burden if κ is freed**. No-fitted-constant claim
  depends on κ ≈ 1; treating κ as fit-free is a strong assumption —
  Greatrix's limiter is effectively a hidden fit parameter.
- **Re-calibration required**. Existing v0.7.1 Hasegawa A LHS rank-1
  will not transfer — `erosion_coeff` in particular will likely need
  re-fitting downward because Z-N sharpens the spike.

## Primary citations (with URLs)

1. Greatrix, D.R., "Transient Burning Rate Model for Solid Rocket
   Motor Internal Ballistic Simulations," *Int. J. Aerospace Eng.*,
   2008. https://onlinelibrary.wiley.com/doi/10.1155/2008/826070 —
   **canonical IB-coupled Z-N reference; explicit low-L\* ignition-
   spike validation**
2. Zarko, V.E. and Gusachenko, L.K., "Critical Review of
   Phenomenological Models for Studying Transient Combustion of
   Solid Propellants," *Int. J. Spray and Combustion Dynamics*,
   2(2), 151-167, 2010.
   https://journals.sagepub.com/doi/10.1260/1756-8277.2.2.151 —
   **alternative-models survey; Z-N comes out best for time-domain**
3. "A Phenomenological Model for the Unsteady Combustion of Solid
   Propellants from a Zel'dovich-Novzhilov Approach," *Aerospace*
   (MDPI), 10(9), 767, 2023.
   https://www.mdpi.com/2226-4310/10/9/767 — **modern explicit-
   closure Z-N derivation**
4. Greatrix, D.R., "Scale Effects on Solid Rocket Combustion
   Instability Behaviour," *Energies*, 4(1), 90, 2011.
   https://www.mdpi.com/1996-1073/4/1/90/htm — **Z-N + scale-effect
   validation; same Greatrix family**
5. "A Review of Calculations for Unsteady Burning of Solid
   Propellant," AIAA Journal review.
   https://arc.aiaa.org/doi/abs/10.2514/3.4980 — **foundational
   survey covering Novozhilov, Culick, Levine, BDP**
6. "Review of dynamic burning of solid propellants in gun and rocket
   propulsion systems," *Symposium (Int.) on Combustion*, 1977.
   https://www.sciencedirect.com/science/article/pii/S0082078477804062
   — **historical anchor; defines τ₁ (solid) and τ₂ (gas) relaxation
   times explicitly**
7. "Numerical investigation of combustion instability in solid rocket
   motors under overload conditions," *Int. J. Thermal Sciences*,
   2025.
   https://www.sciencedirect.com/science/article/abs/pii/S1359431125019829
   — **recent CFD-coupled Z-N application, reports <5% pressure-trace
   error**
8. Coats, D.E. et al., "IHPRPT Improvements to the Solid Performance
   Program (SPP)," AFRL-PR-ED-TR-2003-0061.
   https://apps.dtic.mil/sti/tr/pdf/ADA467974.pdf — **SPP IGT module
   (1-D true ignition transient with Z-N-class dynamic burn rate);
   JANNAF reference**
