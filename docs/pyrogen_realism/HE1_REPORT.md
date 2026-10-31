# HE1 -- head-end igniter-jet heat transfer (opt-in) -- report

Branch `he1/head-end` (from `en1/total-energy` bd0f898, i.e. with the EN1 and RT1 options available), 2026-10-31.
Spec: `HE1_TASK.md`, `HEAD_END_NOTES.md`, `R9_SOURCES.md`. Probes: `docs/solver_verification/probes/{he1_summary,he1_threshold,en1_motor_tables}.py`;
raw rows `docs/solver_verification/worklog/HE1_results/`. Tests: `tests/test_he1_jet.py`. Figures (read by eye): `figures/carlson1967_fig8.png`, `carlson1967_fig9_10.png`.
All runs through `guarded_run`; measured with `solid_heating_model='conduction'` (passed explicitly) and `goodman` as the reference.

## 1. The model (opt-in, default off, bit-identical)

`run_simulation(..., head_end_jet_h=False|True|'max'|'replace')`. Only in the `forward_plenum` topology (the head baskets have no jet and are untouched: test
`test_head_basket_is_untouched` asserts the pressure trace and per-cell ignition times are identical with the option on). While the igniter flows
(`mdot_igniter > 0`), UNLIT grain cells with `x <= 10 D_hyd` get the Cohen-Derr-Price coefficients as published in Peretz et al. 1973 (Eqs. III-28..III-30):

    III-29  x < x_im :  Nu_d = 0.7 Re_d^0.8 (Ap/At)^0.4 (0.022 + 0.066 x/x_im)
    III-28  x = x_im :  Nu_d = 0.0615 Re_d^0.8 (Ap/At)^0.4     (equals III-29 at x_im to 0.2 %; the code uses III-29 below x_im and III-30 from x_im)
    III-30  x >= x_im:  Nu_d = 0.7 Re_d^0.8 (Ap/At)^0.4 [0.031 + 0.057 exp(-0.5 x_r/d)],   x_r = x - x_im

and `h = Nu_d k / D`. No new constants. Choices and caveats:

* **Re_d = 4 mdot_igniter / (pi D mu)**: the IGNITER mass flow (plenum vent flow at that step) and the local port diameter `D_hyd`, as Carlson and Seader
  state. Peretz does not define Re_d for III-28..30; Cohen-Derr-Price (CPIA 220 vol. II) is not available, so the definition is **pending** that source.
* mu, Pr, k at the film temperature exactly as the Peretz bore closure (III-24..III-26): `mu = 0.87e-7 sqrt(W) T_af^0.65`, `Pr = gamma/(1.77 gamma - 0.45)`, `k = Cp mu/Pr`.
* x is the distance from the igniter exit (head end), `x_im = D` (one port diameter, Peretz citing Love et al.), `x_r = x - x_im`.
* `Ap/At` is read as the bore closure reads it (`peretz_ap_reference`, default local port area).
* Range cut: **x/D <= 10**, because Peretz states the duct data extend only to x/d ~ 10; beyond it the jet term is zero (the bore closure alone). Nothing covers the head end after igniter burn-out, so the term is off then.
* **Pr factor**: III-28..30 as quoted carry no Pr term; the Carlson and Seader Dittus-Boelter reference form has Pr^0.4. The published form is used. Including Pr^0.4 would lower h by
  a factor 0.72^0.4 = 0.88 (-12 %) at the sim's Pr ~ 0.72: small against the 2-5x effects below, flagged.
* **Driving temperature**: the cell's local gas T (or the RT1 recovery temperature if selected); Peretz's prescribed effective T_ig (adiabatic minus 5 % and 25-35 %) is not needed and not used.
* **Combination with the bore h**: primary `'max'`: `h = max(h_bore, h_jet)` in the jet region while the igniter flows; sensitivity `'replace'`: `h = h_jet` there.
* The jet h at ignition is 5.6-6.7 kW/(m^2 K) (max over cells; `ignition_diag` column 14).

## 2. Carlson and Seader Figs. 8-10 vs III-29/III-30 (READ BY EYE, not a digitisation)

Nu ratio to Dittus-Boelter vs X/Dp, t = 0 s curves. The formulas give, relative to DB (0.023): 0.7*mult/0.023 * (Ap/At)^0.4 with mult 0.022+0.066 x/D (x<D) or 0.031+0.057 exp(-0.5 (x-D)/D):

| X/Dp | formula, Ap/At = 2 | figure, Ap/At = 2 (Dp/Di 5-15) | formula, Ap/At = 3 | figure, Ap/At = 3 (Dp/Di 5-10) |
|---|---|---|---|---|
| 1 | 3.5 (the peak) | 3.5-5.2 (the peak, 1.1-2.5 for Dp/Di 5-15) | 4.2 (peak) | 2.5-3.5, still rising |
| 2 | 2.6 | 4.0-5.0 | 3.1 | 3.4-4.5 (peak at 2-2.7) |
| 3 | 2.1 | 3.3-3.6 | 2.4 | 3.0-4.0 |
| 10 | 1.27 | 1.3 | 1.49 | 1.3 |

Reading: same order of magnitude and the same far-field (X/Dp = 10) value; the formulas peak at X/D = 1 by construction, whereas the figures peak later (1.1-2.7, later for larger Dp/Di and Ap/At)
and are higher at 2-5 D (the formulas are 25-40 % below the figures there) and lower at the nozzle exit (X/Dp < 0.7). Fig. 8 (Ap/At = 1.1) has no resolved peak and values of about 3-4 near X/Dp = 0.3-0.5, far above the
formula's 1.3 there. So III-28..30 look like a reduced collapse of the figure family (as `HEAD_END_NOTES.md` guessed) rather than a transcription; for the dead cells (x/D 3-10) the formulas, if anything, UNDER-state the
measured ratio. A proper digitisation is left to the user.

## 3. Results (conduction solid; sensible/legacy/static = "A", total energy + conservative momentum + recovery driver = "B")

### 3.1 Dead zones at the fixed charges (forward plenum, DeMar off)

| motor / charge | form | jet | first / t90 [ms] | never lit | peak [MPa] @ ms | plateau / trough, mid | ratio |
|---|---|---|---|---|---|---|---|
| Chunc 2.62 g | A | off | 5.31 / 7.01 | 10 (cells 3-12) | 8.211 @ 26.3 | plateau 6.766 | 1.214 |
| Chunc 2.62 g | A | **max** | 2.26 / 4.73 | **0** | 11.514 @ 21.4 | plateau 8.148 | 1.413 |
| Chunc 2.62 g | A | replace | 2.26 / 4.84 | 0 | 11.480 @ 21.6 | 8.148 | 1.409 |
| Chunc 2.62 g | B | off | 4.91 / 6.81 | 9 | 8.529 @ 25.4 | 6.894 | 1.237 |
| Chunc 2.62 g | B | max | 1.87 / 4.49 | 0 | 11.746 @ 20.2 | 8.144 | 1.442 |
| BALLS 4.83 g | A | off | 6.02 / 8.08 | 9 (cells 3-11) | 6.812 @ 24.4 | 3.523 / 3.884 | pk/tr 1.933, pk/mid 1.754 |
| BALLS 4.83 g | A | **max** | 1.22 / 5.08 | **0** | 9.841 @ 20.2 | 4.345 / 4.887 | 2.265 / 2.014 |
| BALLS 4.83 g | B | off | 5.15 / 7.22 | 8 | 7.442 @ 22.2 | 3.604 / 3.986 | 2.065 / 1.867 |
| BALLS 4.83 g | B | max | 1.17 / 4.75 | 0 | 10.336 @ 18.8 | 4.341 / 4.887 | 2.381 / 2.115 |
| Hasegawa A 14.98 g | A | off | 5.88 / 10.05 | 8 (cells 3-10) | late 5.363 @ 2.69 s (0.83x exp) | | |
| Hasegawa A 14.98 g | A | **max** | 2.32 / 6.51 | **0** | late **6.107** @ 2.40 s (**0.95x**) | | |
| Hasegawa A 14.98 g | A | replace | 2.32 / 6.83 | 0 | late 6.107 @ 2.40 s | | |
| Hasegawa A 14.98 g | B | off | 5.31 / 9.55 | 8 | late 5.359 @ 2.69 s | | |
| Hasegawa A 14.98 g | B | max | 1.96 / 6.32 | 0 | late 6.103 @ 2.40 s | | |

Goodman (reference) rows behave the same: every dead zone is lit (0 never lit in all 12 goodman rows with the jet), first light 0.8-2.1 ms, Chunc ratio 1.46-1.52, BALLS pk/tr 2.37-2.50,
Hasegawa A late peak 6.10-6.11 MPa (goodman/off already lights all cells in form B with the RT1 driver: 6.07).
Per-cell ignition times (Chunc, A, max): cells 3..13 light at 2.3, 2.5, 2.7, 2.9, 3.1, 3.4, 3.5, ... ms: an ignition front propagating from the head end outward; without the jet the same cells never light.
Experiment: Chunc ratio 1.027, plateau about 8.43; BALLS pk/tr about 1.50, pk/mid about 1.35; Hasegawa A late peak 6.44 MPa.

### 3.2 `max` vs `replace`
Practically identical (Chunc first light identical, peak 11.51 vs 11.48 MPa; t90 differs by 0.1-0.3 ms; Hasegawa A late peak identical): in the jet region while the igniter flows h_jet exceeds the stagnating bore h by 5-15x,
so `max` is `replace`. The structural choice does not matter for the result.

### 3.3 Auto-sizing path (IG3 default path; conduction; the charge it finds, cells never lit, startup ratios)

| motor | form | jet | final charge [g] (attempts) | first / t90 [ms] | never lit | peak [MPa] | plateau / trough, mid | ratio |
|---|---|---|---|---|---|---|---|---|
| Chunc | A | off | 2.10 (4) | 8.34 / 11.48 | 1 | 9.745 @ 35 ms | 7.991 | 1.219 |
| Chunc | A | max | **1.34** (2) | 7.02 / 454.8 | 0 | 9.942 @ 489 ms | 8.112 | 1.226 |
| Chunc | B | off | 2.10 (4) | 8.48 / 11.75 | 1 | 9.805 | 7.988 | 1.227 |
| Chunc | B | max | 1.34 (2) | 6.28 / 454.8 | 0 | 9.955 @ 488 ms | 8.109 | 1.228 |
| BALLS | A | off | 3.86 (4) | 9.59 / 12.93 | 8 | 5.815 @ 33 ms | 3.607 / 3.971 | pk/tr 1.612, pk/mid 1.464 |
| BALLS | A | max | **3.09** (3) | 5.17 / 12.62 | 0 | 7.291 @ 34 ms | 4.345 / 4.855 | 1.678 / 1.502 |
| BALLS | B | off | 3.86 (4) | 8.44 / 11.69 | 7 | 6.213 | 3.690 / 4.075 | 1.684 / 1.525 |
| BALLS | B | max | 3.09 (3) | 3.57 / 10.41 | 0 | 7.562 | 4.342 / 4.859 | 1.742 / 1.556 |
| Hasegawa A | A | off | 14.98 (3) | 5.88 / 10.05 | 8 | late 5.363 @ 2.69 s | | 0.83x exp |
| Hasegawa A | A | max | **9.58** (1) | 5.60 / 238.9 | 0 | late **6.113** @ 2.64 s | | **0.95x** |
| Hasegawa A | B | off | 14.98 (3) | 5.31 / 9.55 | 8 | late 5.359 | | |
| Hasegawa A | B | max | 9.58 (1) | 4.99 / 238.8 | 0 | late 6.109 | | |

With the jet the fallback finds a 21-36 % SMALLER charge (Chunc 2.10 -> 1.34 g, BALLS 3.86 -> 3.09 g, Hasegawa A 14.98 -> 9.58 g) that lights the whole grain. The ignition is then progressive: the fore part lights
in 2-6 ms and the rest cascades (t90 240-455 ms for Chunc and Hasegawa A). Startup ratios at these charges: Chunc 1.23 (exp 1.027; unchanged from 1.22 off), BALLS pk/tr 1.68 and pk/mid 1.50 (exp 1.50 / 1.35; 1.61 / 1.46 off),
Hasegawa A late peak 6.11 MPa (0.95x of 6.44) against 5.36 off.

### 3.4 Forward-plenum ignition thresholds (IG2 section 4.1 criterion: minimum charge with >= 99 % of the grain lit within 0.25 s; conduction; A; bisected to 3 %)

| motor | jet off [g] | jet max [g] | ratio |
|---|---|---|---|
| Chunc | 2.00 | 1.38 | 0.69 |
| BALLSstick | 8.35 | 2.83 | 0.34 |
| Hasegawa A | 21.2 | 9.31 | 0.44 |

The transition is sharp (fraction lit jumps from 0.02 to 1.0 between consecutive bisection steps with the jet: the head-end front either runs through the grain or stalls). Real charges: Chunc 4.5-6 g.

## 4. What this says

1. **The head-end jet h removes the dead zones everywhere**: 24 of 24 fixed-charge rows with the jet have 0 cells never lit (conduction and goodman, A and B), against 7-12 without it. The Hasegawa A late peak returns to
   6.10-6.11 MPa (0.95x experiment) for both forms and both solid models; the RT1 recovery driver only did this under goodman.
2. **At the dead-zone charges the startup is far too violent**: first light 1.2-2.3 ms and the whole grain lit in ~5 ms; Chunc ratio 1.41-1.44 (exp 1.03), BALLS pk/tr 2.3-2.4 (exp 1.5). The charges were sized under the old, weak head-end h,
   so they are now over-sized for this physics; the auto-size path, run with the jet, shrinks the charge by 21-36 % and then the startup ratios are about as good as, or slightly worse than, without it (1.23, 1.68/1.50) while ALL cells light and the Hasegawa late peak is 0.95x.
3. The ignition thresholds drop by 31-66 %, i.e. with the jet the model needs 1.4 / 2.8 / 9.3 g (Chunc / BALLS / Hasegawa A) where the real charges are 4.5-6 g for Chunc: the over-sized-charge spike is therefore a consequence of any realistic charge, as in IG2/IG3, not removed by this term.
4. The recovery driver (RT1) and total energy (EN1) matter little once the jet is on (A and B give nearly the same rows): the jet term dominates the head end.

## 5. Recommendation (the user decides the defaults)

* Ship `head_end_jet_h` opt-in (`'max'`). It is a sourced, constant-free term, it removes the last dead-zone artifact, and it is what makes the Hasegawa A late peak right under the current default solid model.
* If it is adopted as the default, the canonical charges must be re-derived through the auto-size path (it does this by itself; expect 1.3 / 3.1 / 9.6 g), B0 re-recorded, and the knobs (roughness/kappa/T_ignition/k_solid) re-checked, since they were fit with the dead zones.
* Unresolved: the Re_d definition (pending CPIA 220), the post-burn-out head-end h (no data), and the x/D <= 10 cut (data limit). A digitisation of Figs. 8-10 would settle whether the figures (higher than the formulas at 2-5 D) should be used instead.

## 6. Gates

Defaults are unchanged by construction (`head_end_jet_h=False`: `h_jet = 0`, no branch taken). Tests added/updated: `tests/test_he1_jet.py` (59 passed together with the RT1 and phase-3 kernel-call tests); three kernel-call tests updated for the extra argument and the (N, 15) `ignition_diag`.
B0 and the full suite are run once after the rebase onto the conduction-default base (see the final report).
