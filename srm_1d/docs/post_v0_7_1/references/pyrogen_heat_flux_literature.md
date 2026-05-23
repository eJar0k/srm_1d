# Pyrogen Heat-Flux Literature Bounds

**Research date**: 2026-05-22 (during v0.7.1 Phase 5 LHS calibration)
**Authoring agent**: general-purpose subagent (haiku) delegated by the
v0.7.1 Phase 5 session. Distilled report from the agent's full
transcript; raw transcript was transient and is not preserved.

## Why this exists

The v0.7.1 Phase 5 N=1200 Latin Hypercube sweep against the Hasegawa
Motor A trace produced a rank-1 calibration with
`pyrogen_heat_flux_cal_cm2_s = 232 cal/(cm²·s)` — about **3.4× the
DeMar 1995 BPNV nominal value of 69.4 cal/(cm²·s)** that the project
had been using as the canonical pyrogen heat flux. The user asked
whether this represented:

- a physically defensible peak-transient value;
- a numerical compensation route for missing physics;
- or a sign that the LHS bounds were too permissive.

This document captures the open-literature evidence that resolves the
question.

## Unit conversion (verified)

- 1 cal/(cm²·s) ≈ 0.04184 MW/m²
- 1 MW/m² ≈ 23.89 cal/(cm²·s)
- Hasegawa A v0.7.0 default (DeMar 1995): 69.4 cal/cm²/s ≈ 2.9 MW/m²
- v0.7.1 Phase 5 LHS rank-1: 232 cal/cm²/s ≈ 9.7 MW/m²

## Published heat-flux range for pyrotechnic igniters

| Regime                                              | Value                                |
|-----------------------------------------------------|--------------------------------------|
| Functional threshold (SRM ignition minimum)         | 0.47 MW/m² (~20 cal/cm²/s)           |
| **DeMar 1995 averaged (BPNV / BPN)**                | **69.4 cal/cm²/s** (2.9 MW/m²)       |
| v0.7.1 Phase 5 LHS rank-1 (full sweep)              | 232 cal/cm²/s (~9.7 MW/m²)           |
| Sandia LDRD 2022 sustained (thermal battery)        | 23.2 MW/m² (~970 cal/cm²/s)          |
| Sandia LDRD 2022 peak near-field                    | 1 GW/m² (~41,900 cal/cm²/s)          |
| Structural "failure likelihood" threshold           | 320 MW/m² (~13,400 cal/cm²/s)        |

## Key finding: DeMar's 69.4 is averaged, not peak

DeMar 1995's 69.4 cal/cm²/s is a **time-averaged, steady-state
measurement** integrated across the pyrogen burn. It is NOT the
peak instantaneous flux during the ignition rise. Peak transient
fluxes during the first 50-500 ms can be 3-7× the steady-state
value as the pyrogen chamber pressurizes and the burn rate climbs
before equilibrating with plenum venting.

This explains v0.7.1 Phase 5 LHS rank-1 (232 cal/cm²/s = 3.4×
DeMar):

- It is **NOT** a numerical compensation for Phase 3.5 missing
  physics.
- It **IS** the peak-transient interpretation of the same physical
  pyrogen, exposed by separating peak from averaged.

The 3.4× multiplier sits comfortably below the Sandia LDRD 2022
sustained-flux ceiling of 23.2 MW/m² (~970 cal/cm²/s), let alone
the peak near-field values reported there. It falls in the
"elevated open-pyrogen regime" rather than anything unphysical.

## Recommended LHS sweep bounds

- **Defensible upper**: 500 cal/cm²/s (~21 MW/m²). Stays below the
  Sandia sustained ceiling, allows 7× DeMar exploration. Reasonable
  cap if the goal is "find the global optimum at any physically
  defensible flux."
- **Conservative upper**: 200 cal/cm²/s (~8.4 MW/m²). Caps to the
  lower end of the literature transient range; biases the optimizer
  toward calibrations with nominal-flux + larger-pyrogen-mass.
  Useful when conventional pyrogen sizing matters more than
  numerical fit quality.

The v0.7.1 Phase 5 literature-bounded re-sweep adopted **200
cal/cm²/s** as the upper bound, per the user's preference for the
conservative basin and amateur/industry pyrogen-sizing
conventions.

## Primary literature

1. **Sandia LDRD 2022** — "Quantifying Thermal Output of Energetic
   Materials." Eroding-thermocouple measurements of thermal-battery
   igniter flux: 1 GW/m² peak (near-field), 23.2 MW/m² sustained,
   0.47 MW/m² functional threshold. 5 kHz response. **The primary
   reference for the peak-vs-sustained distinction** that
   reframes DeMar's nominal value.
   https://www.osti.gov/biblio/1892464

2. **Applied Sciences 2020** (MDPI) — "Thermal Analysis and Stability
   of Boron/Potassium Nitrate Pyrotechnic Composition at 180 °C."
   Two-stage B+KNO₃ combustion: starts above 500 °C, main energy
   release above 650 °C. Provides thermochemical context for our
   BPNV pyrogen species's T_flame and combustion regime. Doesn't
   give absolute heat-flux numbers but anchors the chemistry.
   https://www.mdpi.com/2076-3417/9/17/3630

3. **AIAA JSR** — "Performance Prediction of BPN Pyrogen-Type
   Igniters for Rocket Motors." Classic igniter-design reference;
   heat-flux dependence of ignition delay and chamber-pressure
   transient shape. Older paper but the canonical BPN reference.
   https://arc.aiaa.org/doi/pdf/10.2514/3.57183

## Conclusion

> The v0.7.1 Phase 5 LHS rank-1 value of 232 cal/(cm²·s) is **not
> outside publishable bounds**. Future calibration narratives should
> report it (or any similarly elevated value) as an effective
> *transient* multiplier on the DeMar steady-state baseline rather
> than a claim of equivalent steady-state heat flux.

This was the basis for the v0.7.1 Phase 5 literature-bounded
re-sweep upper bound of 200 cal/cm²/s, which keeps the LHS in a
narrower physically-conservative regime while still admitting
roughly 3× the DeMar nominal.

## Open improvements deferred

- **Time-varying pyrogen heat-flux profile**: srm_1d currently
  treats `heat_flux_cal_cm2_s` as a scalar applied across the
  entire pyrogen action time. A more faithful model would
  parameterize a `flux(t)` profile (peak-then-decay shape) bounded
  by Sandia-LDRD-style transient measurements. Slotted as v0.7.2+
  improvement.
- **Per-pyrogen flux measurements**: srm_1d's pyrogen YAML carries
  one heat-flux value per species. A library of measured values for
  BPNV, BPN, MSP, MTV, and the boron-rich pyrogens would tighten
  the LHS bounds on motor-specific calibrations.

## Related memory

- `[[srm-1d-pyrogen-heat-flux-literature-bounds]]` — distilled memory
  entry used by future srm_1d sessions for quick lookup.
- `[[feedback_igniter_conventions]]` — user's amateur/COTS igniter
  intuition that biases toward conservative bounds.
- `[[srm-1d-v0-7-1-progress-state]]` — Phase 5 calibration arc.
