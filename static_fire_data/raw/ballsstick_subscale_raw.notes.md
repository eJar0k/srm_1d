# ballsstick_subscale_raw.csv — provenance & cleaning notes

**Status: RAW. Not usable for validation as-is. Do not load directly.**

## Provenance
Raw DAQ export from a **2″ subscale** firing of the BALLSStick design.
**No 3″-scale static fire exists** — so validating `motors/BALLSstick.ric` (3″)
against this requires both **cleaning** and a **justified subscale→3″ scaling**.

Columns: `Relative Time, Date, Time Stamp UTC, Volt, Volt, PSI, Pounds`
(2 voltage channels; PSI = pressure transducer; Pounds = load cell).

## Measured state (verified 2026-06-19)
- **417 rows**; raw time **[1388.16, 1392.32] s** — this is the **DAQ power-on
  clock**, not zeroed to ignition.
- **Firing window** (PSI > baseline + 10 %): t = [1388.38, 1390.69] →
  **duration ≈ 2.31 s**.
- **Pre-fire baseline PSI = 60.7**, where ambient should read ≈ **14.7** →
  **≈ +46 psi discrepancy** on the pressure channel.
- **Peak = 1329.7 psi = 9.168 MPa** raw; **≈ 8.85 MPa** if the 46 psi is treated
  as a pure zero-offset.
- **Load cell:** peak ≈ **2800 N**; baseline-subtracted integrated impulse
  ≈ **3971 N·s** vs the expected **≈ 4670 N·s** → **0.85×** (≈15 % low).

## Known issues (user-reported + confirmed above)
1. **Time not zeroed** — starts at the DAQ power-on timestamp (~1388 s).
2. **Pressure channel calibration suspect.** Ambient reads 60.7 psi, not ~14.7.
   **Unresolved:** is this a constant **zero-offset** (subtract ~46 psi) or a
   **gain/scale error** (or both)? The two give different peak pressures and
   different trace *shapes* — do not assume offset-only without a check
   (e.g. does the corrected plateau match the expected Kn·c\* for the 2″ grain?).
3. **Load cell suspect.** Motor mass is **not** subtracted (and the static weight
   *drifts* as propellant burns, so a single tare is not exact); integrated
   impulse comes out 0.85× the expected ~4670 N·s.
4. **Subscale.** 2″ → 3″ scaling is required and must be **physically justified**
   (Kn / c\* / L\* reasoning), not eyeballed — see `feedback_no_unfounded_smoothing`.

## Cleaning TODO (in order)
1. Zero the time axis to ignition (trim to the firing window).
2. Resolve the pressure calibration (offset vs gain — verify against expected
   chamber pressure for the 2″ motor before committing to a correction).
3. Resolve the load cell (tare/motor-mass handling; reconcile the impulse gap).
4. Decide + justify the subscale→3″ scaling.
5. Emit a cleaned `../ballsstick_subscale.csv` and add its load recipe to
   `../README.md`; record every step and its justification here.

**Until steps 1–4 are resolved, this trace cannot support a claim about the 3″
BALLSStick.** It is committed for provenance and so the cleaning is reviewable.
