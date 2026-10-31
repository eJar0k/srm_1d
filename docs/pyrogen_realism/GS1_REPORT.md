# GS1 report -- grid and contact-front convergence of the startup

Branch `gs1/grid` (base `4d68573`, conduction default). Spec: `GS1_TASK.md`. No defaults were changed.
Probes: `probes/gs1_row.py` (one row), `probes/gs1_run.py` (guarded driver), `probes/gs1_analyze.py` (tables).
Raw results: `results/gs1/{rows,profiles,full}.jsonl` (+ `*.guard.jsonl`: guard wall / peak GB per job). Full tables:
`results/gs1/tables.md`.

Configuration: the defaults of `4d68573` (conduction solid model, Peretz ignition h, DeMar off, MIL-P-46994B I-A pellets,
`ma2020`), roughness 32 um, kappa 0.44, T_ignition 756 K, explicit charges (the ignition fallback never resizes them). The
EN1 configuration (total energy / conservative momentum / recovery driver) is not on main and was skipped.
Grids: N = 100, 200, 400 and 800 target propellant cells (dx 8.1 / 4.1 / 2.0 / 1.0 mm for Chunc; 10.2 / 5.1 / 2.5 / 1.3 mm BALLSstick;
16.8 / 8.4 / 4.2 / 2.1 mm Hasegawa A). The 800-cell grid fit trivially (60-180 s per short row), so it is included.

Rows. **Clean** (charge well above threshold): Chunc forward plenum 10 g, BALLSstick 10 g, Hasegawa A 30 g, Chunc head basket
4.5 g (the real igniter), `t_max` 0.15 s. **Threshold** (the fixed IG3 charges): Chunc 2.62 g, BALLSstick 4.83 g, Hasegawa A
14.98 g, `t_max` 0.3 s (0.2 s at N = 800). **Full burns** at the threshold charges for N = 100 and 200: Chunc 1 s, BALLSstick 2 s,
Hasegawa A 10 s. **Profile rows**: clean rows stopped at 0.5 ms and 0.9 ms, final T / species profiles (contact front).
All numbers below come from these runs.

## 0. Summary

1. **The startup is not grid-converged, even at 800 cells.** Every startup metric moves monotonically with refinement at an
   observed order p of about 0.4-1.4 (first-order-like). The 100-cell grid under-resolves the startup far more than T8's ~10 %:
   Chunc forward-plenum 10 g peak 23.6 / 27.8 / 30.9 / 33.0 MPa for N = 100 / 200 / 400 / 800 (Richardson 36.6); first light 1.55 / 1.33 /
   1.20 / 1.13 ms (Richardson 1.06); peak/window 2.79 / 3.28 / 3.67 / 3.97. **The over-spike gets WORSE with refinement**
   (Chunc experiment pk/plateau 1.027): at the fixed Chunc charge 1.26 / 1.37 / 1.45 / 1.50 (Richardson 1.59). The 100-cell
   calibrations were fitted to an under-resolved startup.
2. **Plateau and late-peak values are grid-robust; startup values are not.** Startup-window medians of the clean rows move < 2.5 % over 100 -> 800 cells; at the fixed threshold charges the 0.2-0.3 s window keeps rising at
   400/800 cells as the dark zone shrinks (Chunc 6.50 / 6.60 / 7.22 / 7.58 MPa), so it is only robust from 100 to 200 there (+1.6 %). Late peaks: Chunc 7.249 -> 7.331
   (+1.1 %, N 100 -> 200), Hasegawa A 5.363 -> 5.274 (-1.7 %); BALLSstick 5.556 -> 5.894 (+6 %) is not.
3. **Contact smearing is NOT the grid-limited effect.** The gas-T / hot-gas fronts of the clean plenum rows are about 350-390 mm wide (10-90 %) at
   0.9 ms on all four grids (width ratio 0.96-0.99 per halving of dx; first-order numerical diffusion would give 1.41 and, at 3.62 sqrt(D t), a
   width of only 43-126 mm). The front is a wide physical ramp of the jet-fed fill. What does change with the grid is where the front is: its speed rises
   by 17-21 % from N = 100 to 800 (Chunc 360 -> 435 m/s).
4. **The ignition front runs ahead of the gas front** (600-1000 m/s vs 360-435 m/s on the clean plenum rows) and the head end lights last; the order
   in which grain cells light is not a travelling front (the last quartile has negative speed on most grids: the head end lights after the aft end).
5. **Dark-zone position** (not the cell count) is the meaningful knife-edge metric: Hasegawa A's stays at 172-177 mm on all four grids
   (converged); Chunc's shrinks 102 -> 87 -> 50 -> 38 mm and BALLSstick's 117 -> 79 -> 50 -> 39 mm (Richardson ~33 mm for both, p 1.4-1.7).
6. **No grid up to 800 cells meets a 5 % tolerance on the startup metrics; 800 cells costs ~27 min per Chunc full burn.** See the LHS
   recommendation (Section 5).

## 1. Startup metrics vs grid (clean rows, `t_max` 0.15 s)

Window median = median over 0.1-0.15 s (a startup window, not the plateau). p = observed order from the 200/400/800 triple, R = Richardson value.

| row | metric | N=100 | N=200 | N=400 | N=800 | p | R |
|---|---|---|---|---|---|---|---|
| Chunc FP 10 g | first light [ms] | 1.55 | 1.33 | 1.20 | 1.13 | 0.96 | 1.06 |
| | t90 [ms] | 2.36 | 2.13 | 2.01 | 1.92 | 0.61 | 1.77 |
| | peak [MPa] | 23.58 | 27.76 | 30.92 | 32.95 | 0.64 | 36.63 |
| | t_peak [ms] | 13.44 | 11.46 | 10.37 | 9.77 | 0.87 | 9.05 |
| | peak/window | 2.79 | 3.28 | 3.67 | 3.97 | 0.38 | 4.95 |
| | dark cells / end [mm] | 0 | 0 | 1 / 7 | 7 / 10 | | |
| BALLSstick FP 10 g | first light [ms] | 2.80 | 2.28 | 1.96 | 1.78 | 0.82 | 1.54 |
| | peak [MPa] | 13.54 | 17.47 | 20.72 | 23.12 | 0.44 | 29.91 |
| | peak/window | 2.49 | 3.15 | 3.77 | 4.32 | 0.17 | 8.8 |
| Hasegawa A FP 30 g | first light [ms] | 2.09 | 1.84 | 1.74 | 1.70 | 1.76 | 1.69 |
| | peak [MPa] | 8.13 | 8.95 | 9.27 | 9.45 | 0.84 | 9.68 |
| | peak/window | 1.56 | 1.68 | 1.76 | 1.80 | 1.05 | 1.84 |
| Chunc head basket 4.5 g (the real igniter) | first light [ms] | 4.40 | 3.60 | 3.09 | 2.83 | 1.00 | 2.58 |
| | t90 [ms] | 6.05 | 5.22 | 4.66 | 4.36 | 0.89 | 4.00 |
| | peak [MPa] | 15.83 | 18.83 | 21.07 | 22.39 | 0.77 | 24.24 |
| | t_peak [ms] | 20.42 | 16.82 | 14.79 | 13.62 | 0.80 | 12.04 |
| | peak/window | 1.83 | 2.21 | 2.45 | 2.63 | 0.47 | 3.09 |

Head basket row: error vs Richardson N = 100 / 200 / 400 / 800: first light 70 / 40 / 20 / 10 %; peak 35 / 22 / 13 / 8 %; t_peak 70 / 40 / 23 / 13 %.

## 2. Threshold rows (fixed IG3 charges; the knife edge under refinement)

| row | metric | N=100 | N=200 | N=400 | N=800 | p | R |
|---|---|---|---|---|---|---|---|
| Chunc 2.62 g | first light [ms] | 5.32 | 4.39 | 3.94 | 3.72 | 1.06 | 3.52 |
| | dark cells | 10 | 19 | 22 | 35 | | |
| | dark-zone end [mm] | 101.6 | 87.4 | 49.8 | 38.1 | 1.69 | 32.8 |
| | peak [MPa] / peak over window | 8.21 / 1.26 | 9.07 / 1.37 | 10.48 / 1.45 | 11.38 / 1.50 | 0.65 / 0.63 | 12.96 / 1.59 |
| BALLSstick 4.83 g | first light [ms] | 6.01 | 4.39 | 3.80 | 3.56 | 1.33 | 3.41 |
| | dark cells | 9 | 13 | 17 | 28 | | |
| | dark-zone end [mm] | 117.1 | 78.9 | 49.7 | 38.8 | 1.44 | 32.5 |
| | peak [MPa] | 6.81 | 8.49 | 10.14 | 11.39 | 0.41 | 15.2 |
| Hasegawa A 14.98 g | first light [ms] | 5.87 | 5.33 | 4.96 | 4.74 | 0.72 | 4.40 |
| | dark cells | 8 | 18 | 39 | 82 | | |
| | dark-zone end [mm] | 176.4 | 172.2 | 174.3 | 177.5 | | converged |
| | peak [MPa] | 4.47 | 4.45 | 4.49 | 4.51 | 1.37 | 4.51 |

The dark-cell count is a grid artefact (it scales with N). The dark-zone POSITION converges for Hasegawa A and keeps shrinking for Chunc and BALLSstick
(about 38 mm at 800 cells, 33 mm extrapolated), so the head-end dead zone of those two motors, 100-117 mm on the 100-cell grid, is roughly 3x too long. These
rows show that the knife edge is partly discretisation: the auto-charge fallback sizes a charge to a threshold that itself moves with the grid.

## 3. Full burns (N = 100 and 200; plateau and late peak)

| row | N | steps | wall [s] | peak GB | first light [ms] | dark cells | startup peak [MPa] | pk/window | late peak [MPa] | t late [s] |
|---|---|---|---|---|---|---|---|---|---|---|
| Chunc 2.62 g, 1 s | 100 | 537172 | 24 | 1.53 | 5.32 | 10 | 8.206 | 1.21 | 7.249 | 1.00 |
| | 200 | 1073603 | 75 | 2.00 | 4.39 | 19 | 9.072 | 1.32 | 7.331 | 1.00 |
| BALLSstick 4.83 g, 2 s | 100 | 851541 | 48 | 2.01 | 6.01 | 9 | 6.813 | 1.55 | 5.556 | 2.00 |
| | 200 | 1695749 | 117 | 2.82 | 4.39 | 13 | 8.489 | 1.77 | 5.894 | 1.93 |
| Hasegawa A 14.98 g, 10 s | 100 | 2347930 | 92 | 2.42 | 5.87 | 8 | 5.363 (late; t 2.69 s) | - | 5.363 | 2.69 |
| | 200 | 4696087 | 488 | 3.79 | 5.33 | 18 | 5.274 (late; t 2.71 s) | - | 5.274 | 2.71 |

Late-peak changes 100 -> 200: Chunc +1.1 %, Hasegawa A -1.7 %, BALLSstick +6.1 % (the unlit head-end block of BALLSstick is large and moves). The Hasegawa A
late peak (6.44 MPa experiment) stays 0.82-0.83x on both grids; the early (startup) peak of that run is below the late one on every grid.

## 4. Contact front

Profile rows (clean rows stopped at 0.5 and 0.9 ms, before first light): 10-90 % width of the gas-T rise `w`, position of the 50 % level `x50`, and the
leading edge (1 % level) `lead`. `u_f` is the speed of x50 between 0.5 and 0.9 ms. Numerical prediction: first-order upwind
D = u dx (1 - C)/2 (C = 0.3 = the CFL target), width 3.62 sqrt(D t) at t = 0.9 ms. Full table in `tables.md`.

| row | N | dx [mm] | u_f [m/s] | measured w at 0.9 ms [mm] | predicted numerical w [mm] | measured / predicted | w(coarser)/w |
|---|---|---|---|---|---|---|---|
| Chunc FP 10 g | 100 | 8.13 | 360 | 353 | 110 | 3.2 | - |
| | 200 | 4.06 | 392 | 367 | 81 | 4.5 | 0.96 |
| | 400 | 2.03 | 419 | 378 | 59 | 6.4 | 0.97 |
| | 800 | 1.02 | 435 | 391 | 43 | 9.2 | 0.97 |
| Hasegawa A FP 30 g | 100 | 16.8 | 347 | 266 | 155 | 1.7 | - |
| | 800 | 2.10 | 349 | 341 | 55 | 6.2 | 0.96 |
| Chunc head basket 4.5 g | 100 / 200 / 400 / 800 | 8.1 / 4.1 / 2.0 / 1.0 | 284 / 291 / 319 / 342 | 94 / 67 / 104 / 134 | 98 / 70 / 52 / 38 | 1.0 / 1.0 / 2.0 / 3.5 | - / 1.40 / 0.65 / 0.78 |

- **Is the smearing real?** In the numerical-diffusion sense, no: the width does not shrink as sqrt(dx) (ratio ~0.97, not 1.41), and the
  measured width is 3-9x the first-order prediction. It is a wide (about 0.35-0.39 m, nearly the bore length) physical ramp set by the jet-fed fill
  (exponentially weighted plenum injection and the expansion), the same in mm on all grids. The hot-gas species front matches (w 267-338 mm).
  The leading edge (1 % of the rise = +30 K, harmless to ignition) runs 190-245 mm ahead of x50 on every grid.
- **Head basket (uncontained, the realistic igniter).** The front is much sharper (50-130 mm) and is not converged: x50 at 0.5 ms moves 93 -> 114 -> 143 -> 164 mm
  and u_f 284 -> 342 m/s from N = 100 to 800 (+20 %), and the width grows with refinement (a steepening front, not a smearing one).
  This is the case that matters for the LHS, and its speed has about 20 % grid error at N = 100 and 8 % at 400-800.
- **How much of the ignition-front speed does it explain?** Little. The grain lights far ahead of the 50 % gas front (Chunc FP: Q1-Q3 ignition speeds
  665-1210 m/s vs 360-435 m/s for x50) and the gas front speed itself is only mildly grid-dependent (+17-21 %). The first-light delay, whose order is ~1
  (Chunc FP 1.55 -> 1.06 ms R), is set at the head end near the injection cell, which suggests dilution of the injected energy in a coarse first cell, not
  smearing of a contact front far downstream. That is a hypothesis, not tested here (it would need the cell-0 energy budget at two grids).
  Q1-Q3 ignition speeds move 5-30 % across grids; the last quartile is negative on most grids (the head end lights after the aft end).

## 5. Recommended LHS grid and its cost

Error of each grid vs the Richardson-extrapolated value (worst case over the rows with a clean 200-400-800 triple, observed p >= 0.3; indicative,
since the observed p of 0.4-1.4 is not asymptotic): first light 82 / 48 / 27 / 15 % (N = 100 / 200 / 400 / 800); t90 72 / 40 / 24 / 15 %; startup
peak 55 / 44 / 33 / 25 %; t_peak 115 / 78 / 57 / 41 %; peak/window 44 / 34 / 26 / 20 %. For the realistic head-basket row alone (Section 1): first light
70 / 40 / 20 / 10 %, peak 35 / 22 / 13 / 8 %.

Cost per run (explicit charge, no fallback; guard figures). Measured: short rows at all four grids; full burns at N = 100 and 200. N = 400 and 800
full-burn figures are extrapolated by N^2 (steps ~ N, cost per step ~ N) and are NOT measured.

| N | Chunc 0.15 s (measured) | BALLSstick 0.15 s | Chunc 1 s full | BALLSstick 2 s full | Hasegawa A 10 s full | peak memory |
|---|---|---|---|---|---|---|
| 100 | 4 s | 18 s | 24 s | 48 s | 92 s | 1.5-2.4 GB |
| 200 | 12 s | 23 s | 75 s | 117 s | 488 s | 2.0-3.8 GB |
| 400 | 34 s | 42 s | ~300 s (est.) | ~470 s (est.) | ~2000 s (est.) | est. 3-6 GB (6M-row history regrow for Hasegawa A) |
| 800 | 117 s | 114 s | ~1600 s (est.) | ~1900 s (est.) | ~8000 s (est.) | est. 6 GB+ |

Each LHS point adds a cold fallback charge search (4-5 attempts for auto charges, 3.5-12 s each on the 100-cell grid, x(N/100)^2 on finer grids); explicit charges avoid it.

Recommendation (the user picks the tolerance):

- **Plateau / late-peak targets (T1: within ~3 %)**: N = 100 already meets it for Chunc and Hasegawa A (late peak +1.1 % / -1.7 % from 100 to 200). BALLSstick (+6 %) does
  not. If the fit uses only plateau / late-peak values, N = 100-200 is enough.
- **Startup targets (T2: within 10 % of the converged value)**: needs N = 800 (first light and peak 10 % / 8 % for the realistic igniter, worst rows 15-25 %), i.e. ~27 min
  per Chunc full burn: not an overnight LHS (1000 samples x 3 motors would be ~50 days). **N = 400 gives 13-33 % and ~5 min per Chunc run.**
- **Practical choice: N = 200** as the LHS grid, with the startup metrics (first light, peak, peak/plateau, dark-zone end) given a stated grid-error band of
  about -22 / -40 % (peak / first light on the realistic igniter) or, better, dropped from the fit objective in favour of the converged ones (plateau, late peak, burn
  time), until a higher-order advection (Section 6) cuts the startup error. At N = 200 a Chunc full burn is 75 s and a BALLSstick 117 s (Hasegawa A 8 min, 3.8 GB: above
  a 4 GB guard, so it needs the larger budget). If the user insists on fitting the startup, **N = 400** is the smallest grid whose realistic-igniter errors (first light 20 %,
  peak 13 %) are within the same order as the experimental uncertainty, at ~4x the cost.

## 6. Recommendation: higher-order advection

The startup metrics converge at ~first order, which is what a first-order upwind scheme gives for the quantities set by the head-end injection, the hot-gas arrival
and the front speed (head-basket x50 +20 % from 100 to 800 cells). A second-order MUSCL/TVD scheme (van Leer or minmod limiter) for T and the species mass
fractions (and, in a second step, rho and u) is the natural follow-up. Expected benefit: the observed order rises to ~2 where the front is smooth, so the N = 200 values
should land near the first-order N = 800-1600 values (a ~4-8x cost saving at equal startup error). Caveats from this study: (a) the wide plenum-row front is already grid-independent, so MUSCL
will not sharpen it; the gain is in the head-end injection cell and the head-basket front. (b) The test to run first: MUSCL on T and Y at N = 200 against the first-order
N = 800 row values (Chunc FP: first light 1.13 ms, peak 32.95 MPa; head basket: 2.83 ms, 22.39 MPa). GS1 measured only; nothing was implemented.

## 7. Process notes

- 70 guarded jobs in total, ~48 min of job wall (serial equivalent; the actual wall was shorter with 2-3 workers). Peak per job 1.1-3.8 GB, all at or under the cap; at most 3
  concurrent jobs.
- Warm-up: one cold-compile job (Chunc FP 10 g, N = 100, 61 s including the compile) was run alone before the parallel jobs; no LLVM / compile contention was seen.
- **Stalls:** no job stalled or timed out (the unexplained B0 timeout of 3878 s was not reproduced). Two Hasegawa A full-burn jobs hit the 4 GB guard because the probe
  stored a snapshot every 0.2 ms for 10 s (~1 GB of frames); they were re-run with a 5 ms snapshot interval and completed (N = 100: 93 s / 2.4 GB; N = 200: 488 s / 3.8 GB, which needed
  the 6M-row history regrow). N = 200 Hasegawa A is the only job above ~3 GB.
- Cut / not done: the species-front check uses the final profile at two times only (not a time series); no 800-cell full burns; no 400-cell full burns (the extrapolated costs above are
  estimates); the observed orders (0.4-1.4) mean the Richardson values carry an uncertainty of order 10-20 %; the EN1 configuration was skipped (not on main).
