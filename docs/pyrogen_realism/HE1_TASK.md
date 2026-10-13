# HE1 — head-end igniter-jet heat transfer (opt-in)

**User decision 2026-10-12** ("merge now, head-end next"). With the DeMar flux off
(IG3), the auto-sized igniters leave head-end dead zones that never light:
- BALLSstick 4.83 g: cells 3-11, x 36-117 mm;
- Hasegawa A 14.98 g: cells 3-10, x 59-176 mm. Its late peak is 5.36 MPa,
  0.83x the 6.44 MPa experiment.

SM1 showed this is gas-side, not a solid-model artifact. Once the downstream
grain lights, the 1-D flow at the closed head end stalls (u ~ 0.3-6 m/s), and
Peretz's bore h, built from the LOCAL velocity, collapses 5-15x. The user's
point: a 1-D model cannot represent the jet impingement and recirculation of a
real head end, which feed hot gas to the surface far faster. Several dead cells
end only 1-12 K short of 756 K.

Read first:
- `HEAD_END_NOTES.md`: Peretz Eqs. III-28..III-30 read from the page image, and
  how Peretz used them;
- `R9_SOURCES.md`: the Carlson & Seader 1967 summary (`docs/references/carlson1967.pdf`);
- `SM1_REPORT.md`; `IG2_REPORT.md` section 2 (the bore closure).

## The model (opt-in, default off)

New kwarg, e.g. `head_end_jet_h=False`. It only acts in the `forward_plenum`
topology, because there is a jet. The head basket has no jet: state that and
leave it alone.

While the igniter flows (plenum vent mdot > 0), grain cells within the duct-data
range (x/d <= ~10; report the exact cut you use and why) get a jet-driven h from
the Cohen-Derr-Price correlations:
- III-28 at impingement;
- III-29 upstream of impingement;
- III-30 downstream of impingement.

Use them exactly as published, with no new constants. Input definitions:
- **Re_d** from the **igniter mass flow** and the local port diameter,
  Re = 4 mdot_ig / (pi D mu). This is Carlson & Seader's statement. Peretz does
  not define Re_d for III-28..30; say so.
- **mu and properties:** the same film-temperature treatment the Peretz bore
  closure uses (`_peretz_heat_transfer_coeff`, Peretz III-24..26).
- **x, x_im, x_r:** x from the igniter exit (the forward closure); x_im = one port
  diameter (Peretz, citing Love et al.); x_r = x - x_im.
- **Ap/At:** the local reading, consistent with the default `peretz_ap_reference`.
- **Pr:** III-28..30 as quoted carry no Pr term, while Carlson & Seader's reference
  form has Pr^0.4. Use the published form, report the size of the difference, and
  flag it.
- **Driving temperature:** the local gas T, as now. Peretz prescribed an
  effective T_ig (adiabatic minus 5 % and 25-35 %), which the sim does not need.
  Note it.

**Combination with the local bore h.** This is a structural choice the user
should see. Primary: h = max(h_bore_local, h_jet) in the jet region while the
igniter flows. The heat transfer is at least what the jet drives, and at least
what the local flow drives; Peretz built III-31 to equal III-30 when the only
flow is the igniter's. Also run the variant "h_jet replaces h_bore in the jet
region" as a sensitivity. After igniter burn-out the jet term is off; nothing in
the data covers that period.

## Optional check against the original data

Render Carlson & Seader Figs. 8-10: single pages, <= 110 dpi, cropped. Compare
the III-29 / III-30 multiplier on Dittus-Boelter at matching Ap/At and X/Dp.
Label it as read by eye: it is not a digitization. The user may digitize the
figures properly later.

## Measure (solid models: `conduction` and `peretz`; `goodman` as reference)

1. **Dead zones** at fixed charges: BALLSstick 4.83 g and Hasegawa A 14.98 g,
   forward plenum, DeMar off. Report:
   - cells never lit and the per-cell ignition time;
   - startup peak and plateau;
   - BALLSstick peak/trough and peak/mid vs ~1.50 / ~1.35;
   - Hasegawa A late peak vs 6.44 MPa.
2. **Auto-sizing (IG3 default path)** on the three canonical motors: final
   charge, attempts, cells never lit, startup ratios vs experiment.
3. **Chunc:** canonical forward plenum, and confirm the head basket is unchanged.
4. **Forward-plenum ignition thresholds**, as in IG2 section 4.1.
5. **Gates:** B0 PASS with defaults; full suite green.

## Deliverables

- `docs/pyrogen_realism/HE1_REPORT.md`, with a recommendation. The user decides
  defaults.
- Opt-in implementation, tests (hand-checked III-28..30 values, continuity at
  x_im, default identity) and probes on your branch.
- Paywalled-reference rule applies. The Cohen-Derr-Price paper (CPIA 220) is
  not online; mark anything that depends on its definitions "pending".
