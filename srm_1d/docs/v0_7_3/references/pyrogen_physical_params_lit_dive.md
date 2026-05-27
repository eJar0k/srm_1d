# Pyrogen Physical Parameters — Literature Dive

**Research date**: 2026-05-25 (during v0.7.3 Phase B follow-on)
**Scope**: defensible Saint-Robert burn-rate constants, density, flame
temperature, and pellet geometry for three pyrogen materials whose
YAMLs currently carry seed values flagged "pending measured data."

**Materials covered**:

1. **BPNV** — Boron / KNO₃ / Viton (binder), 25:60:15 nominal mass.
2. **MTV** — Magnesium / Teflon / Viton (binder), Mg-rich variants.
3. **Cu/Al thermite** — CuO + Al stoichiometric; no current YAML.

**Confidence legend**:

| Tier   | Meaning                                                        |
|--------|----------------------------------------------------------------|
| HIGH   | Primary-source experimental measurement, peer-reviewed journal |
| MEDIUM | Derived from related composition / single-source value         |
| LOW    | Industry-standard estimate / amateur-practice consensus        |

**Units convention**: srm_1d's `PropellantTab` stores Saint-Robert as
`r_b [m/s] = a · P^n  with P in Pa`. The published BKNO₃ result
below is reported as `r [mm/s] = 71.1 · P^0.589 with P in MPa`. The
conversion is:

```
a_SI [m/(s·Pa^n)] = (a_pub_mm_s_MPa) · 1e-3 · (1e-6)^n
                  = 71.1 · 1e-3 · (1e-6)^0.589
                  = 71.1 · 1e-3 · 4.92e-4
                  ≈ 3.50e-5
```

This factor governs how the YAML's `a` should be reported in SI; the
conversion is shown for every entry below.

---

## 1. BPNV (Boron / Potassium Nitrate / Viton, 25:60:15)

| Property                       | Value                              | Units             | Citation                                                     | Confidence |
|--------------------------------|------------------------------------|-------------------|--------------------------------------------------------------|------------|
| Saint-Robert form              | `r = 71.1 · P^0.589`               | mm/s, MPa         | Mizushima et al. 2016 (J. Japan Soc. Aeron. Sp. Sci. 64-2)   | HIGH       |
| → SI form                      | `a ≈ 3.50e-5, n = 0.589`           | m/(s·Pa^n), —     | derived from Mizushima                                       | HIGH       |
| Validity range                 | 5 kPa – 7 MPa                      | Pa                | Mizushima et al. 2016                                        | HIGH       |
| Low-pressure floor             | `r = 3.0 ± 0.5` (P < 5 kPa)        | mm/s              | Mizushima et al. 2016                                        | HIGH       |
| Test composition               | B : KNO₃ : binder = 28 : 70 : 2    | mass %            | Mizushima et al. 2016                                        | HIGH       |
| Mass-ratio difference vs BPNV  | binder 2% vs 15% Viton             | —                 | YAML composition                                             | MEDIUM     |
| Pressed-pellet density (Viton) | 1700                               | kg/m³             | DeMar ProCast/QuickBurst formulation; current YAML value     | LOW        |
| Pressed-pellet density (NC)    | 1850–1900                          | kg/m³             | Island Pyrochemical Industries pellet spec (MIL-P-46994B)    | MEDIUM     |
| Adiabatic flame temperature    | 2800                               | K                 | NASA-CR-205274; current YAML; ICT pyrotechnic handbook       | MEDIUM     |
| Gas-product mole weight        | 0.030 (estimated, K-rich exhaust)  | kg/mol            | current YAML; ICT/Applied Sciences 2020 chemistry of B+KNO₃  | LOW        |
| γ (effective)                  | 1.25                               | —                 | current YAML estimate (K-rich, ~25% condensed)               | LOW        |
| Impetus W                      | 5000                               | psi·in³/g         | DeMar 2021 (commercial BPNV spec)                            | MEDIUM     |
| Time-averaged heat flux        | 69.4                               | cal/(cm²·s)       | DeMar 1995/2021 (canonical calorimeter measurement)          | HIGH       |
| Typical pellet diameter (industry) | 3.2                            | mm                | Mizushima et al. 2016 (this is the standard NASA test size)  | HIGH       |
| Typical pellet L/D (industry)  | 0.5–1.0                            | —                 | Mizushima geometry; NASA basket-igniter convention           | MEDIUM     |
| Typical pellet diameter (amateur, ProCast) | 1.5–3.0               | mm                | QuickBurst/Apogee ProCast straw-cast practice (5/32" molds)  | LOW        |
| Typical pellet L/D (amateur)   | 1–3                                | —                 | amateur HPR practice (cut-to-length straws)                  | LOW        |

### Key BPNV findings

- **The Mizushima et al. 2016 value `r = 71.1 P^0.589` (mm/s, MPa)
  is the highest-confidence Saint-Robert fit for any BKNO₃-class
  pyrogen.** It was measured on 3.2 mm-diameter pressed pellets in
  vacuum/low-pressure chambers — directly relevant to the
  amateur-rocketry pyrogen regime where bore P starts near
  atmospheric.
- The Mizushima composition (B:KNO₃:NC = 28:70:2) is binder-light
  vs BPNV's 25:60:15. Viton binder absorbs ~5–10% of the heat-of-
  reaction and modestly suppresses `a`; expect BPNV's effective `a`
  to be ~10–20% lower than the 3.50e-5 SI value at the same `n`.
  **Recommended YAML update**: a = **2.8e-5 to 3.5e-5**, n = **0.59**
  (replacing seed 2.0e-5 / 0.50).
- The current YAML `a = 2.0e-5` is ~43% below the literature
  central value when `n = 0.59` is adopted. At 1 MPa, current YAML
  gives `r = 0.63 mm/s`; literature gives `r = 1.10 mm/s` for
  Mizushima composition (~1.0 mm/s for BPNV-binder-corrected).
- Density 1700 kg/m³ in the YAML is consistent with Viton-bonded
  amateur pyrogen; industry NC-bonded pellets (e.g., IPI MIL-P-46994B)
  press to 1850–1900 kg/m³. **A `binder_kind` switch in the YAML
  would be cleaner than a single ρ value.**
- Particle dimensions: **amateur HPR practice** is 5/32" (~4 mm) cast
  cylinders cut to ~6 mm length (L/D ≈ 1.5). **Industry** standard
  is the NASA 3.2 mm × 3.2 mm right cylinder (L/D = 1).

### BPNV primary literature

- **Mizushima et al. 2016**, "Burning Rate Measurements of 3.2-mm-dia.
  BKNO₃ Pellets in Low Pressures," *J. Japan Soc. Aeron. Sp. Sci.*
  Vol 64 No 2. Source for `r = 71.1 P^0.589`.
  https://www.jstage.jst.go.jp/article/jjsass/64/2/64_64_139/_article/-char/en
- **Mizushima et al. 2024**, "Measurements of Burning Rates of BKNO₃
  Pellets with a Diameter of 3.2 mm under Vacuum Conditions,"
  *Trans. JSASS Aerospace Tech. Japan* Vol 22.
  https://www.jstage.jst.go.jp/article/tastj/22/0/22_22.25/_pdf
- **Berdoyes 1991** (Steinz-era ref), "Pressure Transients for Boron-
  Potassium Nitrate Igniters in Inert, Vented Chambers," NASA CR
  20150020422. Validates pressure-rise dynamics that bound `a, n`.
  https://ntrs.nasa.gov/citations/20150020422
- **Applied Sciences 2020** (MDPI), "Thermal Analysis and Stability
  of B/KNO₃." T_flame chemistry, two-stage combustion.
  https://www.mdpi.com/2076-3417/9/17/3630
- **AIAA JSR (Sloan 1979)**, "Performance Prediction of BPN Pyrogen-
  Type Igniters for Rocket Motors."
  https://arc.aiaa.org/doi/pdf/10.2514/3.57183
- **Island Pyrochemical Industries pellet spec** (MIL-P-46994B / AS-
  4362C compliance). Industrial pellet dimensions + density.
  https://www.islandpyrochemical.com/boron-potassium-nitrate-pellets/

---

## 2. MTV (Magnesium / Teflon / Viton, Mg-rich)

| Property                       | Value                              | Units             | Citation                                                     | Confidence |
|--------------------------------|------------------------------------|-------------------|--------------------------------------------------------------|------------|
| Saint-Robert form              | `r ≈ 2.0 – 6.3 · P^n`, n = 0.30 ± 0.10 | mm/s, MPa     | Koch 2002 (MTV-III/IV); Kuwahara/Koch IPB ref               | MEDIUM     |
| → SI form (central)            | `a ≈ 2.0e-5`, n ≈ 0.30             | m/(s·Pa^n), —     | derived; Koch class data                                     | MEDIUM     |
| Burn rate, Mg=58%, ρ=1.60      | 2.02                               | mm/s              | Kuwahara/Koch PEP 24:65                                      | HIGH       |
| Burn rate, Mg=60%, ρ=1.80      | 6.27 (peak)                        | mm/s              | Koch PEP 27 (MTV-IV, 2002)                                   | HIGH       |
| Validity range                 | 0.02 – 0.1 MPa (Vieille law fit)   | MPa               | Koch 2002                                                    | HIGH       |
| Extrapolation note             | Pressure exponent stays low up to ~5 MPa; weak P-dependence is the defining MTV trait | — | Koch & Clément 2007 review                | MEDIUM     |
| Pressed-pellet density (peak)  | 1800                               | kg/m³             | Koch 2002; current YAML value                                | HIGH       |
| Pressed-pellet density (loose) | 1600                               | kg/m³             | Kuwahara/Koch (low-consolidation)                            | MEDIUM     |
| Adiabatic flame temperature    | 2700–3100                          | K                 | Koch PEP 27 (varies with Mg%); current YAML 3000 is mid-range| HIGH       |
| Gas-product mole weight        | 0.030–0.040 (Mg-rich → much condensed MgF₂) | kg/mol  | Kosanke pyrotechnic handbook; current YAML 0.032             | LOW        |
| γ (effective)                  | 1.22                               | —                 | current YAML estimate (MgF₂ + C particles)                   | LOW        |
| Impetus W                      | 4700 (estimated)                   | psi·in³/g         | current YAML; consistent with Koch's heat-of-explosion       | LOW        |
| Time-averaged heat flux        | 110 (current YAML)                 | cal/(cm²·s)       | "DeMar 2021" reference in YAML; Williams 2013 implies ~30% radiative loss → consistent | MEDIUM     |
| Williams 2013 radiative fraction | ~30% of Q_chem                   | —                 | Williams 2013 PEP                                            | HIGH       |
| Fireball surface T (1 kg scale) | ~1800                             | K                 | Williams 2013                                                | HIGH       |
| Typical pellet diameter (industry) | 6.0 (decoy flare); 3.2 (igniter) | mm              | Koch & Clément 2007                                          | MEDIUM     |
| Typical pellet L/D (igniter)   | 1                                  | —                 | Brazilian/Egyptian MTV igniter qualification papers          | MEDIUM     |

### Key MTV findings

- **MTV's pressure exponent is genuinely low** (0.2–0.4 depending on
  Mg% and consolidation), well below the current YAML's `n = 0.5`
  seed. This is the *defining* property of MTV that makes it popular
  for flares (insensitive burn rate). For srm_1d simulation:
  **lowering n to ~0.30 will substantially reduce predicted P-peak
  during the ignition transient** vs the current n=0.5 seed,
  because chamber P rises into the multi-MPa regime while MTV's burn
  rate barely accelerates.
- **No primary source gives a Saint-Robert fit with the Pa/m units
  srm_1d uses directly.** The Koch / Kuwahara data are reported as
  point burn-rate values at fixed atmospheric pressure, or as
  Vieille-law fits over the 0.02–0.1 MPa range. **Extrapolating to
  the ~1–10 MPa SRM bore-pressure regime is a real uncertainty** —
  most MTV literature is flare-grade (1 atm).
- The 110 cal/(cm²·s) heat-flux value in the YAML is consistent
  with Williams 2013 (~30% radiative fraction of a high-flame-T
  pyrotechnic) but is not a direct calorimeter measurement from
  that paper. **Recommended YAML update**: a = **2.0e-5**, n = **0.30**
  (replacing seed 3.0e-5 / 0.50). Keep ρ, T_flame, heat_flux as-is.
- Particle dimensions for igniters are typically 3.2 mm × 3.2 mm
  pellets (similar to BKNO₃ basket-igniter convention). MTV is not
  commonly used in amateur HPR (electrostatic-sensitivity concern).

### MTV primary literature

- **Koch 2002**, "Metal-Fluorocarbon-Pyrolants IV: Thermochemical
  and Combustion Behaviour of MTV," *Propellants, Explosives,
  Pyrotechnics* Vol 27. Burn-rate-vs-composition curve.
  https://onlinelibrary.wiley.com/doi/abs/10.1002/prep.200290004
- **Kuwahara & Koch**, *PEP* Vol 24 No 1 (1999), "Pyrotechnic
  Igniter Compositions." Includes 58% Mg / 4% Viton burn rate of
  2.02 mm/s @ ρ = 1.60.
  https://www.ibb.ch/publication/Igniters/PEP_24_65.pdf
- **Williams 2013**, "Heat Flux Measurement from Bulk MTV Flare
  Composition Combustion," *PEP*. Radiative fraction, fireball T.
  https://onlinelibrary.wiley.com/doi/10.1002/prep.201200111
- **Koch & Clément 2007**, "Special Materials in Pyrotechnics: V.
  Military Applications of MTV." Review.
- **Aulia et al. 2017** (Brazilian/Egyptian rocket motor qualification),
  "Qualification of MTV Pyrotechnic Composition Used in Rocket Motors
  Ignition System." Mass-fraction sensitivity for ignition delay.
  https://www.academia.edu/34566946/
- **DTIC ADA243244**, "A Theoretical Study of the Combustion of
  Magnesium/Teflon/Viton." Heat-of-reaction tables, condensed-phase
  fraction.

---

## 3. Cu/Al Thermite (CuO + Al, stoichiometric)

**Architectural disclaimer up front**: Cu/Al thermite is **NOT a
Saint-Robert pyrogen.** Its primary energy delivery is radiative
flux + molten droplet impingement; chemistry is condensed-phase
(~90% condensed products by mass). The current
`PropellantTab(a, n)` rate law does not describe the underlying
physics. See architectural recommendation at end.

| Property                          | Value                          | Units            | Citation                                          | Confidence |
|-----------------------------------|--------------------------------|------------------|---------------------------------------------------|------------|
| Reaction                          | 2 Al + 3 CuO → Al₂O₃ + 3 Cu    | —                | Reese et al. 2013                                 | HIGH       |
| Adiabatic flame temperature       | 2843                           | K                | Reese et al. 2013                                 | HIGH       |
| Heat of reaction                  | ~4.0                           | MJ/kg            | Reese et al. 2013; Fischer & Grubelich 1996       | HIGH       |
| Pressed-pellet density (TMD ~5.1) | 3500–4500 (typical 70–90% TMD) | kg/m³            | Reese et al. 2013; Frontiers Chem. 2018           | HIGH       |
| Gas-product mass fraction         | ~10% (Cu vapor + residual O₂)  | —                | Reese et al. 2013                                 | HIGH       |
| Gas-product mole weight           | ~64 (Cu vapor dominant)        | g/mol            | derived from Reese product analysis               | MEDIUM     |
| γ (effective, gas-only)           | 1.30 (monatomic-rich, Cu)      | —                | estimate                                          | LOW        |
| Linear flame propagation rate     | 5–25 (loose-pressed); 1–10 (high-TMD) | mm/s      | Reese 2011/2013; Frontiers Chem. 2018 (Al/CuO)    | HIGH       |
| Pressure exponent n               | ~0 (deflagration is approximately P-independent in tested range) | — | Reese 2013; consensus thermite literature        | HIGH       |
| Action time (effective)           | 75                             | ms (mass-indep.) | Reese et al. 2013                                 | HIGH       |
| CuO mesh                          | -325 mesh (~44 μm and finer)   | —                | Reese et al. 2013; user's amateur convention      | HIGH       |
| Al mesh                           | -325 mesh (~44 μm and finer)   | —                | Reese et al. 2013                                 | HIGH       |
| Pre-burn particle size            | 10–50                          | μm               | Reese 2013 SEM characterization                   | HIGH       |
| Post-burn agglomerate (droplet)   | 10–50                          | μm               | Reese 2013 high-speed imaging                     | HIGH       |
| Heat-of-reaction radiative fraction | ~50–70% (depends on standoff) | —                | INL "Quantification of Heat Flux from a Reacting Thermite Spray" | MEDIUM     |
| Peak surface heat flux at 5 mm standoff | ~10–30                  | MW/m² (≈ 250-750 cal/cm²/s) | INL spray-flux measurement              | MEDIUM     |
| Critical ignition T (Reese)       | ~625                           | K (AP/HTPB surface) | Reese 2013                                     | HIGH       |
| Critical ignition energy (Reese)  | ~50                            | J/cm² (∫q dt over heat-up phase) | Reese 2013 dual-criterion       | HIGH       |

### Is Saint-Robert even the right rate law for Cu/Al thermite?

**No.** Three independent reasons:

1. **Pressure independence.** Thermite deflagration in the 0.1–10 MPa
   range is governed by Al particle melting + CuO decomposition
   kinetics, not pressure-coupled mass diffusion. Reese 2013, the
   nanothermite literature, and Fischer/Grubelich 1996 all report
   `n ≈ 0` (or n undefined because no pressure dependence was
   detectable across the measured range).
2. **No surface-regression geometry.** A Cu/Al pellet does not have
   a propagating planar flame front the way a BKNO₃ pellet does;
   the reaction front is irregular, with copper droplets ejected
   from the bulk before full conversion. "Linear burn rate" is a
   convenient average, not a regression-rate physics quantity.
3. **Energy delivery is condensed-phase.** ~90% of chemical energy
   reaches the grain as molten Cu/Al₂O₃ droplet impingement +
   thermal radiation. srm_1d's current pyrogen plenum injects
   *gaseous* mass + enthalpy + momentum. Thermite breaks both the
   mass-balance (no gas to advect) and the heat-balance (no enthalpy
   to deliver via gas).

### Architectural recommendation

The Phase A.2 uncontained-pyrogen kernel **can be coerced to model
Cu/Al thermite** with these compromises:

- **Burn rate**: set `n = 0` and `a = 1.5e-5` (gives r ≈ 15 mm/s at
  any P, the Reese loose-pressed central value). Acceptable
  approximation for "how fast does the thermite mass disappear."
- **Mass-injection mode**: the existing `mdot_pyrogen` kernel will
  inject 100% of the consumed thermite mass into the bore as gas.
  This is wrong by ~9× — actual gas yield is ~10%. **Recommended
  workaround**: scale the YAML `rho_pyrogen` down by 10× from the
  TMD value (e.g., effective ρ ≈ 400 kg/m³) so the swept-volume
  *mass-flow* is right even though the depletion timing is wrong by
  the same factor. This is a hack; document it.
- **Heat delivery**: the v0.7.3-phaseB `heat_delivery_mode = 'demar'`
  with `heat_flux_cal_cm2_s = 250` (= 10 MW/m², INL central) gives
  a defensible peak-flux value, but the time profile is wrong (real
  thermite is a 75 ms pulse, not the steady-state-averaged
  approximation `demar` implements).
- **Flame temperature**: T_flame = 2843 K (Reese) is fine.
- **What CANNOT be modeled with the current architecture**: the
  ~90% mass delivered as molten droplets that don't advect with
  the bore gas; secondary ignition lag if the propellant grain
  doesn't directly intercept the droplet spray.

**Cleaner path (v0.8+)**: add a `Thermite` dataclass parallel to
`Pyrogen` with `mass_delivery_mode='gas'|'droplet'` and
`heat_delivery_mode='surface_flux_only'`. This was already
recommended in `srm_1d/docs/post_v0_7_1/references/copper_thermite_igniter_literature.md`.
The Phase B uncontained-pyrogen scaffolding has built much of what
Thermite needs already — adding the gas-yield decoupling is the
remaining work.

### Cu/Al thermite primary literature

- **Reese et al. 2013**, "CuO/Al Thermites for Solid Rocket Motor
  Ignition," *J. Propulsion & Power* Vol 29 No 5 (often cited as
  "Reese 2015" because the print issue postdates the online
  appearance). The canonical reference. Dual-criteria ignition,
  flame T = 2843 K, action time = 75 ms, hot-fire validation.
  https://arc.aiaa.org/doi/abs/10.2514/1.B34771
  Open mirror: https://nakka-rocketry.net/articles/Journal-of-Propulsion-and-Power_-CuO-Al-Thermites-for-Solid-Rocket-Motor-Ignition.pdf
- **Reese 2011 JPC paper** (precursor): "CuO/Al Igniters for Solid
  Rocket Motor Ignition." Earlier sizing data + safety testing.
  https://www.davidree.se/content/2.talks-papers/5.JPC-2011/1.%20CuO-Al-Igniters-for-Solid-Rocket-Motor-Ignition.pdf
- **Fischer & Grubelich 1996**, "Theoretical Energy Release of
  Thermites, Intermetallics, and Combustible Metals," Sandia
  SAND96-2826C. Heat-of-reaction tables.
- **INL Digital Library**, "Quantification of Heat Flux from a
  Reacting Thermite Spray." One of the few published spray-flux
  measurements.
  https://inldigitallibrary.inl.gov/sites/sti/sti/4310589.pdf
- **Frontiers in Chemistry 2018**, "Combustion Characteristics of
  Physically Mixed 40 nm Aluminum/Copper Oxide Nanothermites."
  Nano-scale variant; bounds the size-dependence trend.
  https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2018.00465/full
- **Reese — Purdue dissertation 2013** (Wright advisor). Complete
  characterization including SEM, high-speed video, calorimeter.
  Referenced in JPP paper.

---

## Summary table: YAML changes most worth making

| Material      | Field                       | Current YAML        | Recommended         | Rationale                                                    |
|---------------|-----------------------------|---------------------|---------------------|--------------------------------------------------------------|
| BPNV          | `n`                         | 0.50 (seed)         | **0.59**            | Mizushima 2016 measurement on 3.2 mm BKNO₃ pellets           |
| BPNV          | `a` (m/(s·Pa^n))            | 2.0e-5 (seed)       | **3.0e-5** (range 2.8e-5 to 3.5e-5) | Mizushima 71.1·P^0.589 minus ~10-15% binder correction |
| BPNV          | `rho` (kg/m³)               | 1700                | 1700 (Viton) or 1850 (NC) — add `binder_kind` field | NC-bonded industrial pellets press denser     |
| MTV           | `n`                         | 0.50 (seed)         | **0.30**            | Koch's defining low-pressure-exponent property               |
| MTV           | `a` (m/(s·Pa^n))            | 3.0e-5 (seed)       | **2.0e-5**          | Calibrated against 2.02 mm/s @ 0.1 MPa, ρ=1.60 (Kuwahara)    |
| MTV           | (no other changes)          | —                   | —                   | ρ, T_flame, M, γ, heat_flux all reasonable                   |
| Cu/Al thermite | (no YAML — new file needed) | absent              | **opt-in opt-out**  | See architectural caveat: Saint-Robert is wrong rate law     |

**Highest-impact single change**: BPNV `n: 0.50 → 0.59` + `a: 2.0e-5 →
3.0e-5`. This roughly doubles the burn rate at 1–5 MPa where srm_1d
spends most of the ignition transient and brings the pyrogen mass-
flow into the regime that NASA/JAXA igniter design correlations were
calibrated against. The MTV update has lower impact because MTV is
not the default and the calibration data are sparser.

---

## Open questions deferred

1. **BPNV binder correction factor.** Mizushima used B:KNO₃:NC =
   28:70:2; BPNV is 25:60:15. The 13 percentage-point binder
   difference is large enough that the burn rate should be
   re-measured for the Viton-rich formulation. Until then, the
   recommended `a` carries a 10–15% binder-effect uncertainty band.

2. **MTV high-pressure extrapolation.** All MTV burn-rate data is
   at flare-grade 1 atm. Extrapolation to SRM bore P (1–10 MPa) is
   physically defensible (the low-`n` property is fundamental, not
   a low-P artifact) but not directly measured.

3. **Cu/Al thermite gas yield as a function of confinement.** Reese
   measured ~10% in open burns; confined burns (inside an SRM
   port) may shift toward higher gas fraction due to Cu re-
   condensation suppression. No literature found.

4. **Particle dimensions for amateur vs industry.** Industry pellets
   are 3.2 mm right cylinders (L/D=1, NASA standard). Amateur ProCast/
   QuickBurst practice is wider variance (1.5–5 mm, L/D 1–3). The
   v0.7.3-phaseB `form` field captures this categorically (pellets vs
   chunks vs powder) — recommend documenting the diameter choice in
   each motor's `.ric` instead of the pyrogen YAML.

---

## Related repo references

- `srm_1d/docs/post_v0_7_1/references/pyrogen_heat_flux_literature.md`
  — heat-flux bounds (DeMar / Sandia LDRD / functional thresholds).
- `srm_1d/docs/post_v0_7_1/references/copper_thermite_igniter_literature.md`
  — prior thermite literature review; this document supersedes the
  rate-law discussion and confirms the architectural recommendation.
- `srm_1d/motors/pyrogens/bpnv.yaml` — receives BPNV YAML updates.
- `srm_1d/motors/pyrogens/mtv.yaml` — receives MTV YAML updates.
- `srm_1d/propellant.py:Pyrogen` — Saint-Robert dataclass. No
  schema changes needed for BPNV/MTV updates. A new dataclass
  (`Thermite` or `Pyrogen(rate_law='thermite')`) would be needed
  to cleanly support Cu/Al; out of scope for v0.7.3.
