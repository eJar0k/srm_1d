# SM2 report -- controlled comparison: Peretz surface-temperature ODE vs exact conduction

Branch `sm2/compare` (from `openmotor-frontend` 0eada1e). Goodman is out of the running (user, 2026-10-12); it appears below only as a reference
and as one of the three forcing sources. "Peretz ODE" = `solid_heating_model='peretz'`; "Peretz bore h" = `ignition_heat_transfer='peretz'` (the default
ignition h); "conduction" = `solid_heating_model='conduction'`. All tables are generated from raw results by `probes/sm2_tables.py`
(full output with a few extra tables: `results/sm2/sm2_tables.md`; raw per-run JSON lines `results/sm2/sm2_*.jsonl`). **No defaults changed; the user decides.**

## 0. What was built

1. **Open-loop replay.** `run_simulation(record_forcing_steps=n)` records, per step and per grain cell, the surface forcing the ignition gate saw
   (`h_total`, `T_driver`; flux = `h_total (T_driver - T_s)`, the same combined Robin forcing for all three models). Cells that have already lit record the bare
   convective h of their local gas (the forcing they would see if still unlit), so a replay can run past the reference run's own ignition time.
   `srm_1d/tools/solid_replay.py` pushes a recorded forcing through any solid model cell by cell. **Validation: replaying a run's own forcing through its own
   model reproduces its per-cell ignition times with max |dt| = 0.0 s (every case in table 2, and `tests/test_sm2_replay.py`).** The recorder does not change the trajectory (tested bit-for-bit).
2. **Ignition-event diagnostics (always on, bookkeeping only).** `result['ignition_diag']`: per cell at its lighting step, h_total, T_driver, q_pre, T_s, bore-gas T, u, local Mach, P, t.
3. The recorded forcing (~3 GB, 100 cells x 60-120 k steps per run) lives outside git in `C:/tmp/sv/sm2_data/*_forcing.npz`; regenerate with `probes/sm2_sweep.py A|X|H|D|C <out.jsonl>` then `probes/sm2_tables.py`.
4. Forcing was recorded from a Goodman, a Peretz-ODE and a conduction run of every case, so the forcing itself does not favour a model; the 3x3 matrix (forcing source x replayed model) is table 2.

Cases (fixed explicit charges, forward plenum unless noted, library BPNV, DeMar off, Peretz bore h, k_solid 0.271, T_ignition 756 K, 100 propellant cells):
Chunc 2.62 g, BALLSstick 4.83 g, Hasegawa A 14.98 g (IG3's auto charges under Goodman), Chunc 4.5 g in the head basket. **Added beyond the spec (high-charge rows):** Chunc 10 g,
head basket 8 g, BALLSstick 10 g, Hasegawa A 30 g -- see the first finding for why.

## 1. Findings

1. **The fixed IG3 charges sit at the ignition threshold, where the solid models cannot be compared cleanly.** 2.62 g / 4.83 g / 14.98 g were sized under Goodman as the first charge that lights;
   each solid model has its own threshold, so at the same charge Goodman, Peretz ODE and conduction are at different distances from it. The result is threshold-hypersensitive: e.g. Chunc 2.62 g,
   conduction leaves 10 cells dark (plateau 6.77 MPa, peak 8.21), the Peretz ODE leaves 1 (plateau 8.00, peak 10.37), Goodman 0 (peak 8.70 at plateau 8.07); a +/- one-cell change moves t90 by tens of ms.
   First light and t50 are the robust metrics there. Hence the high-charge rows, where every model lights (almost) every cell and the answer is clean.
2. **Clean regime (high charge): the Peretz ODE and conduction differ very little.** First light: the Peretz ODE is later by 0.13-0.23 ms (3.5-8 % of the 1.5-3.0 ms lighting time); t50 by 0.19-0.42 ms,
   t90 by 0.19-0.39 ms, the 90 %-lighting duration by +0.06 .. +0.20 ms. Startup peaks agree to 0.1 % (Chunc), -1.2 % (head basket), +1.8 % (BALLSstick), and Hasegawa A's late peak to 1.6 % (6.106 vs 6.009 MPa; conduction leaves 1 cell dark).
   Peak/plateau, pk/trough, pk/mid agree to the same level (tables 1b).
3. **Feedback.** For first light the open-loop gap equals the closed-loop gap (feedback ~0.00-0.03 ms): the difference in first light is purely the solid model. For t50/t90 the closed-loop gap is 2-4x the open-loop
   gap at high charge (e.g. head basket t90: closed 0.39 ms, open 0.08, i.e. ~80 % of it is flow feedback: the later-lighting model changes the flow the remaining cells see). At the near-threshold charges feedback is larger and of either sign
   (head basket t50: closed 1.14 ms vs open 0.29 ms), because a one-cell change in what is lit changes the downstream flow.
4. **Calibration moves the lighting time far more than the model choice does** (high charge, table 11 and 10c). T_ignition 756 -> 950 K (the v0.7.5 LHS band is 750-950 K) shifts first light by 0.8-2.2 ms for either model, i.e. 4-10x the model gap;
   k_solid 0.26 -> 0.32 shifts it by 0.12-0.31 ms, about the size of the gap; the spread in peak from T_ignition is 1.1-2.7 MPa, vs a model difference in peak of <= 0.31 MPa. The two models also respond to T_ignition in different ways near threshold:
   at T_ign = 850-950 K with the fixed charges, **conduction fails to light the grain at all (Chunc, BALLSstick: 6 of 9 (k, T) combinations dark; Hasegawa A 7 of 9) while the Peretz ODE still lights it (1, 1 and 3 of 9 dark), with t90 of 14-110 ms** (table 10b).
   At high charge conduction leaves more cells dark at high T_ignition (BALLSstick at 950 K, k 0.32: 11 cells vs 3).
5. **The ignition criterion differs in kind.** (Table 4b.) Conduction cells light by a sustained flux: the time spent within 100 K of threshold before lighting is median 0.11-0.67 ms, typically at most 2-7 ms (one 11.7 ms cell, BALLSstick 10 g).
   The Peretz ODE has the long creep tails SM1 predicted: median 0.17-1.08 ms but max 8.2 ms (Chunc 2.62 g), 194 ms (BALLSstick 4.83 g), 10 ms (Hasegawa A 14.98 g), 29 ms (BALLSstick 10 g). Cells that never light under conduction reach a peak T_s of 680-748 K (median) and 755 K (max),
   i.e. near misses of the 756 K threshold; the Peretz ODE lights such cells by holding the temperature (SM1 section 2.3: it overshoots T_s by ~+275 K after an h drop).
6. **Other confounders are 5-100x larger than the solid-model gap** (table 11, first light): ignition-h closure (Gnielinski vs Peretz bore h) 3.4-4.9 ms; charge (high vs fixed) 2.3-3.8 ms and the auto-sizing path itself (auto charges
   differ per model: Chunc 1.68 g Peretz ODE vs 2.10 g conduction vs 2.62 g Goodman; first light 11.6 vs 8.3 ms); DeMar flux on 0-1.8 ms; grid 0.5-1.6 ms. None of these changes which model lights first by more than the model gap's own sign flips near threshold.
7. **Grid (200 cells).** The model difference survives refinement in magnitude, not always in sign: first-light gap (Peretz ODE - conduction) at 200 cells: Chunc -0.23, head basket +0.21, BALLSstick -0.04, Hasegawa A -0.41 ms (100 cells: -0.05, +0.39, -0.14, -0.27); peak gap at 200 cells (head basket 18.49 vs 18.83 MPa, -1.8 %) is as at 100 cells (-2.2 %). The 200-cell startup peaks are higher for both models (head basket 18.5-18.8 vs 15.5-15.8: the known ~10-20 % under-resolution, T8).
8. **Fill Mach.** (Table 5; thread-F confounder, not fixed here.) The bore Mach reaches 1.04-1.09 at 1 ms and 1.4-1.8 at 2 ms at the fixed charges (Chunc 1.69, BALLSstick 1.42, Hasegawa A 1.76) and 0.7-1.7 at 5 ms in the Peretz-ODE and conduction runs (Goodman, which has finished lighting: 0.4-1.4), and 1.3-1.5 at 1 ms and up to 2.4 at 2 ms (Hasegawa A 30 g) at the high charges. At the first-lit cell the local Mach is 0.43-0.93 and the gas 2290-2560 K at 0.7-1.2 MPa at the fixed charges (0.66-0.83, 2425-2640 K, 1.3-2.1 MPa at the high charges; table 5 gives each model). These are the numbers to compare after a total-energy change.
   Because h uses the local velocity, the supersonic fill inflates the early h; it acts equally on all three models (the open-loop forcing is identical), so it cannot create a Peretz-vs-conduction difference, but it does shift the absolute lighting times and the threshold charges.

## 2. Closed-loop tables (fixed charges, then high charge)

### 1. Closed loop, fixed charges (the base matrix: Peretz bore h, DeMar off, k_solid 0.271, T_ign 756 K, 100 cells)

| case | model | charge [g] | first [ms] | t50 | t90 | t90-first | last | never | peak [MPa] / plateau etc. |
|---|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | goodman | 2.62 | 4.45 | 4.82 | 16.32 | 11.87 | 96.66 | 0 | 8.696 @ 24.5 ms | pl 8.066, pk/pl 1.08 |
| Chunc fwd plenum | peretz | 2.62 | 5.27 | 7.21 | 7.87 | 2.61 | 13.79 | 1 | 10.375 @ 26.5 ms | pl 7.997, pk/pl 1.30 |
| Chunc fwd plenum | conduction | 2.62 | 5.32 | 6.54 | 7.72 | 2.40 | 7.72 | 10 | 8.206 @ 26.3 ms | pl 6.766, pk/pl 1.21 |
| Chunc head basket | goodman | 4.50 | 3.41 | 3.75 | 4.15 | 0.75 | 8.78 | 0 | 16.516 @ 18.8 ms | pl 8.444, pk/pl 1.96 |
| Chunc head basket | peretz | 4.50 | 4.79 | 6.79 | 7.23 | 2.44 | 10.04 | 0 | 15.483 @ 21.4 ms | pl 8.476, pk/pl 1.83 |
| Chunc head basket | conduction | 4.50 | 4.40 | 5.65 | 6.05 | 1.65 | 11.26 | 0 | 15.833 @ 20.4 ms | pl 8.464, pk/pl 1.87 |
| BALLSstick | goodman | 4.83 | 5.02 | 5.42 | 5.78 | 0.76 | 5.80 | 9 | 7.194 @ 22 ms | pk/tr 2.04, pk/mid 1.85 |
| BALLSstick | peretz | 4.83 | 5.88 | 7.79 | 9.20 | 3.32 | 1580.07 | 1 | 7.574 @ 26 ms | pk/tr 1.83, pk/mid 1.64 |
| BALLSstick | conduction | 4.83 | 6.01 | 7.15 | 8.35 | 2.34 | 8.95 | 9 | 6.813 @ 24 ms | pk/tr 1.93, pk/mid 1.75 |
| Hasegawa A | goodman | 14.98 | 4.20 | 4.93 | 6.04 | 1.85 | 6.25 | 8 | late 5.362 (0.83x 6.44) @ 2.68 s | early 4.65 |
| Hasegawa A | peretz | 14.98 | 5.60 | 8.63 | 10.46 | 4.85 | 14.53 | 2 | late 5.915 (0.92x 6.44) @ 2.48 s | early 5.22 |
| Hasegawa A | conduction | 14.98 | 5.87 | 8.72 | 10.32 | 4.45 | 10.97 | 8 | late 5.363 (0.83x 6.44) @ 2.69 s | early 4.47 |


### 1b. Closed loop, HIGH charge (added; far above every model's ignition threshold: Chunc 10 g, head basket 8 g, BALLSstick 10 g, Hasegawa A 30 g)

| case | model | charge [g] | first [ms] | t50 | t90 | t90-first | last | never | peak [MPa] / plateau etc. |
|---|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | goodman | 10.00 | 1.14 | 1.43 | 1.76 | 0.63 | 6.41 | 0 | 25.128 @ 12.5 ms | pl 8.163, pk/pl 3.08 |
| Chunc fwd plenum | peretz | 10.00 | 1.67 | 2.21 | 2.55 | 0.88 | 17.88 | 0 | 23.598 @ 13.5 ms | pl 8.157, pk/pl 2.89 |
| Chunc fwd plenum | conduction | 10.00 | 1.55 | 2.02 | 2.36 | 0.81 | 16.92 | 0 | 23.577 @ 13.4 ms | pl 8.156, pk/pl 2.89 |
| Chunc head basket | goodman | 8.00 | 1.62 | 1.97 | 2.30 | 0.68 | 7.18 | 0 | 26.196 @ 13.5 ms | pl 8.348, pk/pl 3.14 |
| Chunc head basket | peretz | 8.00 | 2.26 | 3.17 | 3.47 | 1.21 | 10.00 | 0 | 25.314 @ 14.2 ms | pl 8.360, pk/pl 3.03 |
| Chunc head basket | conduction | 8.00 | 2.07 | 2.75 | 3.08 | 1.01 | 9.27 | 0 | 25.624 @ 13.9 ms | pl 8.356, pk/pl 3.07 |
| BALLSstick | goodman | 10.00 | 2.06 | 2.43 | 7.10 | 5.04 | 1408.92 | 0 | 14.290 @ 17 ms | pk/tr 3.37, pk/mid 2.98 |
| BALLSstick | peretz | 10.00 | 3.04 | 3.84 | 4.70 | 1.67 | 51.05 | 0 | 13.776 @ 17 ms | pk/tr 3.18, pk/mid 2.81 |
| BALLSstick | conduction | 10.00 | 2.80 | 3.43 | 4.40 | 1.60 | 20.07 | 1 | 13.536 @ 17 ms | pk/tr 3.19, pk/mid 2.83 |
| Hasegawa A | goodman | 30.00 | 1.49 | 2.17 | 2.84 | 1.35 | 3.44 | 0 | late 6.105 (0.95x 6.44) @ 2.38 s | early 8.86 |
| Hasegawa A | peretz | 30.00 | 2.23 | 3.13 | 3.91 | 1.68 | 4.59 | 0 | late 6.106 (0.95x 6.44) @ 2.38 s | early 8.04 |
| Hasegawa A | conduction | 30.00 | 2.09 | 2.90 | 3.68 | 1.59 | 3.88 | 1 | late 6.009 (0.93x 6.44) @ 2.42 s | early 8.13 |


## 3. Open-loop replay (3 forcing sources x 3 models)

### 2. Open-loop replay on the base-matrix forcing

Rows = forcing source (the closed-loop run whose recorded surface forcing is replayed); columns = solid model replayed on it. Cells = first / t50 / t90 / t90-first [ms], never-lit in the replay window. The diagonal equals the closed-loop run (self-replay check: max |dt| shown).

**Chunc fwd plenum**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 207 | 4.45 / 4.82 / 16.32 / 11.87, never 0 | 5.31 / 5.86 / 9.59 / 4.28, never 1 | 5.07 / 5.46 / 17.60 / 12.54, never 1 | 0.0e+00 |
| peretz | 206 | 4.39 / 5.76 / 6.64 / 2.25, never 0 | 5.27 / 7.21 / 7.87 / 2.61, never 1 | 5.29 / 6.83 / 7.62 / 2.33, never 2 | 0.0e+00 |
| conduction | 205 | 4.39 / 5.79 / 6.52 / 2.13, never 2 | 5.27 / 6.87 / 7.43 / 2.15, never 2 | 5.32 / 6.54 / 7.72 / 2.40, never 10 | 0.0e+00 |

**Chunc head basket**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 200 | 3.41 / 3.75 / 4.15 / 0.75, never 0 | 4.47 / 4.64 / 4.87 / 0.40, never 0 | 4.03 / 4.35 / 4.75 / 0.71, never 0 | 0.0e+00 |
| peretz | 200 | 3.32 / 5.30 / 6.17 / 2.85, never 0 | 4.79 / 6.79 / 7.23 / 2.44, never 0 | 4.38 / 6.47 / 7.01 / 2.64, never 0 | 0.0e+00 |
| conduction | 200 | 3.34 / 4.87 / 5.57 / 2.23, never 0 | 4.87 / 5.92 / 6.22 / 1.35, never 0 | 4.40 / 5.65 / 6.05 / 1.65, never 0 | 0.0e+00 |

**BALLSstick**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 251 | 5.02 / 5.42 / 5.78 / 0.76, never 9 | 5.73 / 6.26 / 9.30 / 3.56, never 2 | 5.24 / 5.88 / 20.31 / 15.07, never 10 | 0.0e+00 |
| peretz | 247 | 4.93 / 6.39 / 7.09 / 2.17, never 3 | 5.88 / 7.79 / 9.20 / 3.32, never 2 | 5.92 / 7.39 / 8.94 / 3.02, never 3 | 0.0e+00 |
| conduction | 249 | 4.93 / 6.50 / 6.98 / 2.04, never 4 | 5.90 / 7.50 / 8.59 / 2.69, never 2 | 6.01 / 7.15 / 8.35 / 2.34, never 9 | 0.0e+00 |

**Hasegawa A**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 406 | 4.20 / 4.93 / 6.04 / 1.85, never 8 | 5.89 / 6.91 / 8.35 / 2.46, never 4 | 5.93 / 7.33 / 82.70 / 76.76, never 9 | 0.0e+00 |
| peretz | 407 | 4.02 / 7.24 / 9.14 / 5.12, never 3 | 5.60 / 8.63 / 10.46 / 4.85, never 2 | 5.73 / 8.29 / 10.10 / 4.37, never 3 | 0.0e+00 |
| conduction | 405 | 4.05 / 7.83 / 9.46 / 5.41, never 6 | 5.62 / 9.04 / 10.64 / 5.02, never 2 | 5.87 / 8.72 / 10.32 / 4.45, never 8 | 0.0e+00 |


### 2b. Open-loop replay on the high-charge forcing

Rows = forcing source (the closed-loop run whose recorded surface forcing is replayed); columns = solid model replayed on it. Cells = first / t50 / t90 / t90-first [ms], never-lit in the replay window. The diagonal equals the closed-loop run (self-replay check: max |dt| shown).

**Chunc fwd plenum**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 211 | 1.14 / 1.43 / 1.76 / 0.63, never 0 | 1.53 / 1.71 / 2.42 / 0.89, never 0 | 1.41 / 1.65 / 2.40 / 0.98, never 0 | 0.0e+00 |
| peretz | 210 | 1.11 / 1.61 / 2.14 / 1.03, never 0 | 1.67 / 2.21 / 2.55 / 0.88, never 0 | 1.54 / 2.14 / 2.51 / 0.97, never 0 | 0.0e+00 |
| conduction | 211 | 1.11 / 1.61 / 2.07 / 0.96, never 0 | 1.66 / 2.09 / 2.40 / 0.74, never 0 | 1.55 / 2.02 / 2.36 / 0.81, never 0 | 0.0e+00 |

**Chunc head basket**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 200 | 1.62 / 1.97 / 2.30 / 0.68, never 0 | 2.14 / 2.38 / 2.62 / 0.48, never 0 | 1.95 / 2.29 / 2.58 / 0.62, never 0 | 0.0e+00 |
| peretz | 200 | 1.58 / 2.63 / 3.09 / 1.51, never 0 | 2.26 / 3.17 / 3.47 / 1.21, never 0 | 2.06 / 3.03 / 3.38 / 1.32, never 0 | 0.0e+00 |
| conduction | 200 | 1.59 / 2.44 / 2.82 / 1.23, never 0 | 2.25 / 2.88 / 3.16 / 0.91, never 0 | 2.07 / 2.75 / 3.08 / 1.01, never 0 | 0.0e+00 |

**BALLSstick**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 252 | 2.06 / 2.43 / 7.10 / 5.04, never 1 | 2.71 / 2.96 / 6.26 / 3.55, never 0 | 2.46 / 2.81 / 7.41 / 4.95, never 1 | 0.0e+00 |
| peretz | 250 | 2.01 / 3.20 / 3.79 / 1.78, never 1 | 3.04 / 3.84 / 4.70 / 1.67, never 0 | 2.79 / 3.64 / 4.51 / 1.72, never 1 | 0.0e+00 |
| conduction | 251 | 2.02 / 3.02 / 3.52 / 1.50, never 1 | 2.97 / 3.58 / 4.27 / 1.30, never 1 | 2.80 / 3.43 / 4.40 / 1.60, never 1 | 0.0e+00 |

**Hasegawa A**

| forcing from | window [ms] | replay goodman | replay peretz | replay conduction | self-replay max abs dt [s] |
|---|---|---|---|---|---|
| goodman | 408 | 1.49 / 2.17 / 2.84 / 1.35, never 0 | 2.05 / 2.56 / 3.13 / 1.08, never 0 | 1.94 / 2.49 / 3.09 / 1.15, never 1 | 0.0e+00 |
| peretz | 409 | 1.44 / 2.38 / 3.40 / 1.96, never 0 | 2.23 / 3.13 / 3.91 / 1.68, never 0 | 2.08 / 3.03 / 3.84 / 1.76, never 0 | 0.0e+00 |
| conduction | 409 | 1.45 / 2.39 / 3.36 / 1.91, never 0 | 2.22 / 3.00 / 3.74 / 1.52, never 0 | 2.09 / 2.90 / 3.68 / 1.59, never 1 | 0.0e+00 |


## 4. Feedback share

### 3. Feedback share: closed-loop vs open-loop Peretz-ODE - conduction gap (base matrix)

Gap = (Peretz ODE) - (conduction), ms. closed = full runs; open = replay of the SAME recorded forcing (mean and range over the Peretz-ODE and conduction forcing sources; the range including the Goodman forcing is added); feedback = closed - mean(open).

| case | metric | closed gap | open gap (mean) [min..max over sources] | feedback share (closed - open) |
|---|---|---|---|---|
| Chunc fwd plenum | first | -0.05 | -0.03 [-0.04 .. -0.02]; incl. Goodman forcing [-0.04 .. 0.25] | -0.02 |
| Chunc fwd plenum | t50 | 0.68 | 0.35 [0.33 .. 0.38]; incl. Goodman forcing [0.33 .. 0.40] | 0.32 |
| Chunc fwd plenum | t90 | 0.16 | -0.02 [-0.29 .. 0.25]; incl. Goodman forcing [-8.01 .. 0.25] | 0.17 |
| Chunc fwd plenum | dur90 | 0.20 | 0.01 [-0.25 .. 0.28]; incl. Goodman forcing [-8.26 .. 0.28] | 0.19 |
| Chunc head basket | first | 0.39 | 0.44 [0.42 .. 0.46]; incl. Goodman forcing [0.42 .. 0.46] | -0.05 |
| Chunc head basket | t50 | 1.14 | 0.29 [0.26 .. 0.33]; incl. Goodman forcing [0.26 .. 0.33] | 0.84 |
| Chunc head basket | t90 | 1.18 | 0.19 [0.17 .. 0.22]; incl. Goodman forcing [0.13 .. 0.22] | 0.99 |
| Chunc head basket | dur90 | 0.79 | -0.25 [-0.29 .. -0.20]; incl. Goodman forcing [-0.31 .. -0.20] | 1.04 |
| BALLSstick | first | -0.14 | -0.08 [-0.12 .. -0.05]; incl. Goodman forcing [-0.12 .. 0.50] | -0.05 |
| BALLSstick | t50 | 0.64 | 0.37 [0.35 .. 0.39]; incl. Goodman forcing [0.35 .. 0.39] | 0.26 |
| BALLSstick | t90 | 0.85 | 0.25 [0.23 .. 0.26]; incl. Goodman forcing [-11.01 .. 0.26] | 0.60 |
| BALLSstick | dur90 | 0.99 | 0.33 [0.30 .. 0.35]; incl. Goodman forcing [-11.50 .. 0.35] | 0.66 |
| Hasegawa A | first | -0.27 | -0.19 [-0.25 .. -0.13]; incl. Goodman forcing [-0.25 .. -0.04] | -0.08 |
| Hasegawa A | t50 | -0.09 | 0.33 [0.32 .. 0.33]; incl. Goodman forcing [-0.42 .. 0.33] | -0.42 |
| Hasegawa A | t90 | 0.14 | 0.34 [0.32 .. 0.36]; incl. Goodman forcing [-74.34 .. 0.36] | -0.21 |
| Hasegawa A | dur90 | 0.40 | 0.53 [0.49 .. 0.57]; incl. Goodman forcing [-74.30 .. 0.57] | -0.12 |


### 3b. Same, high charge

Gap = (Peretz ODE) - (conduction), ms. closed = full runs; open = replay of the SAME recorded forcing (mean and range over the Peretz-ODE and conduction forcing sources; the range including the Goodman forcing is added); feedback = closed - mean(open).

| case | metric | closed gap | open gap (mean) [min..max over sources] | feedback share (closed - open) |
|---|---|---|---|---|
| Chunc fwd plenum | first | 0.13 | 0.12 [0.11 .. 0.13]; incl. Goodman forcing [0.11 .. 0.13] | 0.00 |
| Chunc fwd plenum | t50 | 0.19 | 0.07 [0.07 .. 0.07]; incl. Goodman forcing [0.06 .. 0.07] | 0.12 |
| Chunc fwd plenum | t90 | 0.19 | 0.04 [0.04 .. 0.04]; incl. Goodman forcing [0.02 .. 0.04] | 0.15 |
| Chunc fwd plenum | dur90 | 0.06 | -0.08 [-0.10 .. -0.07]; incl. Goodman forcing [-0.10 .. -0.07] | 0.15 |
| Chunc head basket | first | 0.19 | 0.19 [0.18 .. 0.20]; incl. Goodman forcing [0.18 .. 0.20] | 0.00 |
| Chunc head basket | t50 | 0.42 | 0.14 [0.12 .. 0.15]; incl. Goodman forcing [0.09 .. 0.15] | 0.28 |
| Chunc head basket | t90 | 0.39 | 0.08 [0.08 .. 0.09]; incl. Goodman forcing [0.04 .. 0.09] | 0.31 |
| Chunc head basket | dur90 | 0.20 | -0.11 [-0.11 .. -0.10]; incl. Goodman forcing [-0.15 .. -0.10] | 0.31 |
| BALLSstick | first | 0.23 | 0.21 [0.17 .. 0.25]; incl. Goodman forcing [0.17 .. 0.26] | 0.03 |
| BALLSstick | t50 | 0.41 | 0.18 [0.16 .. 0.19]; incl. Goodman forcing [0.14 .. 0.19] | 0.23 |
| BALLSstick | t90 | 0.30 | 0.03 [-0.14 .. 0.19]; incl. Goodman forcing [-1.15 .. 0.19] | 0.27 |
| BALLSstick | dur90 | 0.07 | -0.18 [-0.30 .. -0.06]; incl. Goodman forcing [-1.41 .. -0.06] | 0.25 |
| Hasegawa A | first | 0.14 | 0.14 [0.13 .. 0.15]; incl. Goodman forcing [0.12 .. 0.15] | 0.00 |
| Hasegawa A | t50 | 0.23 | 0.10 [0.10 .. 0.10]; incl. Goodman forcing [0.07 .. 0.10] | 0.13 |
| Hasegawa A | t90 | 0.24 | 0.07 [0.06 .. 0.08]; incl. Goodman forcing [0.05 .. 0.08] | 0.17 |
| Hasegawa A | dur90 | 0.09 | -0.07 [-0.08 .. -0.07]; incl. Goodman forcing [-0.08 .. -0.07] | 0.17 |


## 5. The ignition criterion

### 4. The ignition criterion: what lit each cell (base matrix; `ignition_diag` at the ignition step)

q_pre = h_total (T_driver - T_s) at the start of the lighting step [MW/m2]; "h" = h_total at that step [W/m2K]. "low-flux lit" = lit cells with q_pre < 0.5 MW/m2 (temperature held/crept to threshold rather than driven by a sustained flux; the Peretz-ODE failure mode of SM1).

| case | model | lit | q_pre median [p10 .. p90] | h median [p10 .. p90] | low-flux lit | T_s at ign median [K] |
|---|---|---|---|---|---|---|
| Chunc fwd plenum | goodman | 100 | 4.71 [2.52 .. 4.85] | 2828 [1161 .. 3876] | 0 | 759.7 |
| Chunc fwd plenum | peretz | 99 | 10.35 [4.66 .. 13.32] | 4406 [2893 .. 5606] | 0 | 756.1 |
| Chunc fwd plenum | conduction | 90 | 9.93 [6.47 .. 11.90] | 4171 [3159 .. 4999] | 0 | 756.2 |
| Chunc head basket | goodman | 100 | 5.40 [5.14 .. 5.53] | 2898 [2512 .. 3711] | 0 | 760.3 |
| Chunc head basket | peretz | 100 | 12.99 [8.18 .. 16.38] | 5655 [3751 .. 6902] | 0 | 756.1 |
| Chunc head basket | conduction | 100 | 11.23 [8.35 .. 14.91] | 4862 [3621 .. 6262] | 0 | 756.2 |
| BALLSstick | goodman | 91 | 4.70 [4.59 .. 4.74] | 2707 [2534 .. 3076] | 0 | 759.2 |
| BALLSstick | peretz | 99 | 9.68 [4.56 .. 11.22] | 4347 [2934 .. 5064] | 2 | 756.1 |
| BALLSstick | conduction | 91 | 8.88 [6.84 .. 9.95] | 4068 [3056 .. 4483] | 0 | 756.1 |
| Hasegawa A | goodman | 92 | 4.64 [4.30 .. 4.88] | 3430 [2965 .. 3952] | 0 | 758.1 |
| Hasegawa A | peretz | 98 | 7.73 [4.68 .. 9.87] | 3338 [2516 .. 4261] | 0 | 756.2 |
| Hasegawa A | conduction | 92 | 6.86 [5.10 .. 9.08] | 2961 [2421 .. 3918] | 0 | 756.3 |


### 4b. Lighting by sustained flux or by creeping temperature: time spent within 100 K of T_ignition before lighting (own-forcing replay)

creep = t_ignite - (time T_s first exceeded T_ign - 100 K), over lit cells, ms: median / p90 / max. A flux-driven rise crosses the last 100 K in a fraction of a millisecond; a long creep means the surface hovers near threshold (held temperature). "dark cells peak T_s" = highest T_s reached by cells that never lit (median / max) in the replay window.

| case | charge set | model | lit | creep median / p90 / max [ms] | dark cells | dark-cell peak T_s median / max [K] |
|---|---|---|---|---|---|---|
| Chunc fwd plenum | fixed | goodman | 100 | 0.084 / 13.795 / 94.045 | 0 | - |
| Chunc fwd plenum | fixed | peretz | 99 | 0.595 / 1.979 / 8.237 | 1 | 637 / 637 |
| Chunc fwd plenum | fixed | conduction | 90 | 0.300 / 1.430 / 2.105 | 10 | 733 / 755 |
| Chunc head basket | fixed | goodman | 100 | 0.019 / 0.602 / 5.430 | 0 | - |
| Chunc head basket | fixed | peretz | 100 | 0.431 / 1.229 / 5.296 | 0 | - |
| Chunc head basket | fixed | conduction | 100 | 0.285 / 0.823 / 6.873 | 0 | - |
| BALLSstick | fixed | goodman | 91 | 0.032 / 2.374 / 2.564 | 9 | 735 / 756 |
| BALLSstick | fixed | peretz | 98 | 0.795 / 2.193 / 193.950 | 2 | 665 / 737 |
| BALLSstick | fixed | conduction | 91 | 0.413 / 2.257 / 2.485 | 9 | 728 / 755 |
| Hasegawa A | fixed | goodman | 92 | 0.280 / 1.593 / 1.717 | 8 | 733 / 754 |
| Hasegawa A | fixed | peretz | 98 | 1.079 / 2.081 / 10.112 | 2 | 676 / 729 |
| Hasegawa A | fixed | conduction | 92 | 0.666 / 2.372 / 2.431 | 8 | 733 / 755 |
| Chunc fwd plenum | high | goodman | 100 | 0.041 / 0.204 / 3.502 | 0 | - |
| Chunc fwd plenum | high | peretz | 100 | 0.165 / 0.338 / 7.006 | 0 | - |
| Chunc fwd plenum | high | conduction | 100 | 0.110 / 0.328 / 4.394 | 0 | - |
| Chunc head basket | high | goodman | 100 | 0.032 / 0.314 / 2.167 | 0 | - |
| Chunc head basket | high | peretz | 100 | 0.200 / 0.497 / 2.764 | 0 | - |
| Chunc head basket | high | conduction | 100 | 0.146 / 0.369 / 2.032 | 0 | - |
| BALLSstick | high | goodman | 99 | 0.137 / 1.745 / 10.503 | 1 | 738 / 738 |
| BALLSstick | high | peretz | 100 | 0.321 / 0.723 / 29.343 | 0 | - |
| BALLSstick | high | conduction | 99 | 0.217 / 0.654 / 11.671 | 1 | 680 / 680 |
| Hasegawa A | high | goodman | 100 | 0.024 / 0.259 / 0.910 | 0 | - |
| Hasegawa A | high | peretz | 100 | 0.176 / 0.400 / 1.401 | 0 | - |
| Hasegawa A | high | conduction | 99 | 0.113 / 0.360 / 0.756 | 1 | 748 / 748 |


## 6. Fill Mach and bore-gas state at first light

### 5. Fill Mach and bore-gas state at the first-lit cell (thread-F confounder: the solver transports sensible enthalpy only)

Max |Mach| over the bore at fixed times, and at the first-lit cell the local Mach, bore-gas temperature, velocity, pressure and q_pre. For a later total-energy check these are the numbers to compare.

| case | model | max Mach @1 ms | @2 | @5 | @10 | @20 ms | first-lit: t [ms] | T_gas [K] | Mach | u [m/s] | P [MPa] | h | q_pre [MW/m2] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | goodman | 1.04 | 1.70 | 0.77 | 0.35 | 0.33 | 4.45 | 2288 | 0.93 | 653 | 0.74 | 3132 | 4.80 |
| Chunc fwd plenum | peretz | 1.04 | 1.69 | 1.45 | 0.36 | 0.33 | 5.27 | 2337 | 0.93 | 654 | 0.68 | 2957 | 4.67 |
| Chunc fwd plenum | conduction | 1.04 | 1.69 | 1.45 | 0.37 | 0.33 | 5.32 | 2340 | 0.93 | 654 | 0.68 | 2946 | 4.67 |
| Chunc head basket | goodman | 1.00 | 1.50 | 0.39 | 0.35 | 0.32 | 3.41 | 2430 | 0.60 | 435 | 1.03 | 3298 | 5.52 |
| Chunc head basket | peretz | 1.00 | 1.49 | 1.70 | 0.39 | 0.33 | 4.79 | 2559 | 0.60 | 434 | 1.04 | 3330 | 6.01 |
| Chunc head basket | conduction | 1.00 | 1.49 | 1.56 | 0.36 | 0.32 | 4.40 | 2531 | 0.60 | 434 | 1.04 | 3342 | 5.93 |
| BALLSstick | goodman | 1.01 | 1.42 | 1.36 | 0.43 | 0.40 | 5.02 | 2326 | 0.72 | 518 | 1.23 | 3021 | 4.74 |
| BALLSstick | peretz | 1.01 | 1.42 | 1.37 | 0.79 | 0.41 | 5.88 | 2311 | 0.72 | 517 | 1.18 | 2936 | 4.56 |
| BALLSstick | conduction | 1.01 | 1.42 | 1.37 | 0.42 | 0.40 | 6.01 | 2310 | 0.72 | 517 | 1.17 | 2920 | 4.54 |
| Hasegawa A | goodman | 1.09 | 1.77 | 0.91 | 0.54 | 0.45 | 4.20 | 2442 | 0.50 | 354 | 0.96 | 2936 | 4.95 |
| Hasegawa A | peretz | 1.09 | 1.76 | 0.73 | 0.62 | 0.46 | 5.60 | 2538 | 0.44 | 317 | 0.87 | 2468 | 4.40 |
| Hasegawa A | conduction | 1.09 | 1.76 | 0.73 | 0.59 | 0.46 | 5.87 | 2548 | 0.43 | 309 | 0.85 | 2383 | 4.27 |
| Chunc fwd plenum (high charge) | goodman | 1.45 | 0.90 | 0.37 | 0.33 | 0.30 | 1.14 | 2642 | 0.72 | 523 | 1.75 | 5130 | 9.68 |
| Chunc fwd plenum (high charge) | peretz | 1.45 | 1.28 | 0.36 | 0.34 | 0.30 | 1.67 | 2620 | 0.83 | 607 | 2.11 | 6781 | 12.64 |
| Chunc fwd plenum (high charge) | conduction | 1.45 | 1.43 | 0.36 | 0.33 | 0.30 | 1.55 | 2626 | 0.82 | 594 | 2.05 | 6487 | 12.13 |
| Chunc head basket (high charge) | goodman | 1.32 | 1.87 | 0.40 | 0.34 | 0.30 | 1.62 | 2485 | 0.74 | 541 | 1.35 | 4663 | 8.06 |
| Chunc head basket (high charge) | peretz | 1.32 | 1.86 | 0.38 | 0.34 | 0.30 | 2.26 | 2474 | 0.73 | 529 | 1.73 | 5581 | 9.59 |
| Chunc head basket (high charge) | conduction | 1.32 | 1.86 | 0.38 | 0.34 | 0.30 | 2.07 | 2474 | 0.73 | 532 | 1.63 | 5346 | 9.19 |
| BALLSstick (high charge) | goodman | 1.34 | 1.59 | 0.49 | 0.46 | 0.41 | 2.06 | 2617 | 0.67 | 488 | 1.82 | 4030 | 7.50 |
| BALLSstick (high charge) | peretz | 1.34 | 1.59 | 1.11 | 0.46 | 0.41 | 3.04 | 2588 | 0.72 | 518 | 2.10 | 4742 | 8.69 |
| BALLSstick (high charge) | conduction | 1.34 | 1.58 | 1.00 | 0.46 | 0.41 | 2.80 | 2600 | 0.71 | 514 | 2.06 | 4633 | 8.54 |
| Hasegawa A (high charge) | goodman | 1.48 | 2.38 | 0.46 | 0.44 | 0.43 | 1.49 | 2475 | 0.73 | 543 | 1.34 | 4905 | 8.43 |
| Hasegawa A (high charge) | peretz | 1.48 | 2.37 | 0.56 | 0.44 | 0.44 | 2.23 | 2443 | 0.66 | 477 | 2.09 | 6728 | 11.35 |
| Hasegawa A (high charge) | conduction | 1.48 | 2.37 | 0.46 | 0.44 | 0.44 | 2.09 | 2425 | 0.66 | 481 | 1.98 | 6472 | 10.81 |


## 7. Other factors (closed loop, per model)

### 6. Ignition-h closure: Gnielinski instead of Peretz bore h (fixed charges)

| case | model | first | t50 | t90 | t90-first | never | peak / late | (base: t90, peak) |
|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | peretz | 1.03 | 1.75 | 2.28 | 1.25 | 0 | 12.257 (charge 2.62 g) | 7.87, 10.375 |
| Chunc fwd plenum | conduction | 0.98 | 1.69 | 2.23 | 1.25 | 0 | 12.258 (charge 2.62 g) | 7.72, 8.206 |
| Chunc head basket | peretz | 1.06 | 1.95 | 2.47 | 1.41 | 0 | 17.500 (charge 4.50 g) | 7.23, 15.483 |
| Chunc head basket | conduction | 0.98 | 1.83 | 2.36 | 1.38 | 0 | 17.552 (charge 4.50 g) | 6.05, 15.833 |
| BALLSstick | peretz | 1.22 | 1.99 | 2.68 | 1.46 | 0 | 10.475 (charge 4.83 g) | 9.20, 7.574 |
| BALLSstick | conduction | 1.14 | 1.93 | 2.67 | 1.53 | 0 | 10.444 (charge 4.83 g) | 8.35, 6.813 |
| Hasegawa A | peretz | 1.56 | 2.84 | 3.84 | 2.28 | 0 | 6.106 (charge 14.98 g) | 10.46, 5.915 |
| Hasegawa A | conduction | 1.47 | 2.70 | 3.69 | 2.22 | 0 | 6.106 (charge 14.98 g) | 10.32, 5.363 |


### 7. DeMar flux ON (library demar mode; default is off)

| case | model | first | t50 | t90 | t90-first | never | peak / late | (base: t90, peak) |
|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | peretz | 5.20 | 7.09 | 7.73 | 2.53 | 0 | 10.830 (charge 2.62 g) | 7.87, 10.375 |
| Chunc fwd plenum | conduction | 4.80 | 6.06 | 6.65 | 1.86 | 0 | 10.995 (charge 2.62 g) | 7.72, 8.206 |
| Chunc head basket | peretz | 3.47 | 6.46 | 6.88 | 3.42 | 0 | 15.659 (charge 4.50 g) | 7.23, 15.483 |
| Chunc head basket | conduction | 3.14 | 5.37 | 5.75 | 2.62 | 0 | 16.049 (charge 4.50 g) | 6.05, 15.833 |
| BALLSstick | peretz | 5.88 | 7.70 | 9.14 | 3.26 | 0 | 8.640 (charge 4.83 g) | 9.20, 7.574 |
| BALLSstick | conduction | 5.45 | 6.89 | 8.07 | 2.62 | 0 | 8.599 (charge 4.83 g) | 8.35, 6.813 |
| Hasegawa A | peretz | 4.49 | 7.90 | 9.84 | 5.35 | 0 | 6.108 (charge 14.98 g) | 10.46, 5.915 |
| Hasegawa A | conduction | 4.09 | 6.69 | 8.30 | 4.21 | 0 | 6.107 (charge 14.98 g) | 10.32, 5.363 |


### 8. Auto-sized charge (IG3: DeMar 300 psi sizing + x1.25 fallback), per model

| case | model | first | t50 | t90 | t90-first | never | peak / late | (base: t90, peak) |
|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | peretz | 11.61 | 18.59 | 20.13 | 8.52 | 1 | 9.587 (charge 1.68 g) | 7.87, 10.375 |
| Chunc fwd plenum | conduction | 8.34 | 10.72 | 11.52 | 3.17 | 1 | 9.745 (charge 2.10 g) | 7.72, 8.206 |
| BALLSstick | peretz | 13.10 | 20.75 | 23.80 | 10.70 | 1 | 6.722 (charge 3.09 g) | 9.20, 7.574 |
| BALLSstick | conduction | 9.59 | 11.57 | 13.28 | 3.69 | 8 | 5.815 (charge 3.86 g) | 8.35, 6.813 |
| Hasegawa A | peretz | 14.51 | 21.82 | 23.90 | 9.40 | 2 | 5.916 (charge 11.98 g) | 10.46, 5.915 |
| Hasegawa A | conduction | 5.88 | 8.72 | 10.32 | 4.44 | 8 | 5.363 (charge 14.98 g) | 10.32, 5.363 |


### 9. Grid: 200 propellant cells instead of 100 (fixed charges)

| case | model | first | t50 | t90 | t90-first | never | peak / late | (base: t90, peak) |
|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | peretz | 4.16 | 5.90 | 6.53 | 2.37 | 5 | 10.691 (charge 2.62 g) | 7.87, 10.375 |
| Chunc fwd plenum | conduction | 4.39 | 5.69 | 6.39 | 1.99 | 19 | 9.072 (charge 2.62 g) | 7.72, 8.206 |
| Chunc head basket | peretz | 3.81 | 5.69 | 6.05 | 2.24 | 0 | 18.493 (charge 4.50 g) | 7.23, 15.483 |
| Chunc head basket | conduction | 3.60 | 4.85 | 5.22 | 1.61 | 1 | 18.828 (charge 4.50 g) | 6.05, 15.833 |
| BALLSstick | peretz | 4.35 | 6.01 | 7.20 | 2.84 | 4 | 9.488 (charge 4.83 g) | 9.20, 7.574 |
| BALLSstick | conduction | 4.39 | 5.51 | 6.64 | 2.25 | 13 | 8.489 (charge 4.83 g) | 8.35, 6.813 |
| Hasegawa A | peretz | 4.92 | 7.49 | 9.27 | 4.35 | 6 | 5.822 (charge 14.98 g) | 10.46, 5.915 |
| Hasegawa A | conduction | 5.33 | 7.84 | 9.43 | 4.11 | 18 | 5.274 (charge 14.98 g) | 10.32, 5.363 |


## 8. Calibration

### 10b. Calibration spread vs the model gap, fixed (near-threshold) charges

Spreads are max - min over the 3 k levels (T_ign = 756) and over the 3 T_ign levels (k = 0.271), counting only combinations in which the grain lights at all within the window.

| case | model | k spread at T=756 (0.26..0.32): first / t50 / t90 [ms], peak [MPa] | T_ign spread at k=0.271 (756..950): first / t50 / t90 [ms], peak | combos that never light | gap P-C at nominal first / t50 / t90 |
|---|---|---|---|---|---|
| Chunc fwd plenum | peretz | 1.01 / 1.70 / 1.80, 0.25 | 25.04 / 99.73 / 102.25, 0.78 | 1 of 9 | |
| Chunc fwd plenum | conduction | 1.27 / 1.78 / 2.03, 0.31 | - / - / - | 6 of 9 | |
| Chunc fwd plenum | gap P-C | | | | -0.05 / 0.68 / 0.16 |
| Chunc head basket | peretz | 0.59 / 0.92 / 0.99, 0.35 | 4.42 / 7.04 / 7.48, 2.45 | 0 of 9 | |
| Chunc head basket | conduction | 0.56 / 0.74 / 0.74, 0.30 | 4.79 / 6.07 / 6.34, 2.34 | 0 of 9 | |
| Chunc head basket | gap P-C | | | | 0.39 / 1.14 / 1.18 |
| BALLSstick | peretz | 1.18 / 1.76 / 2.13, 0.43 | 26.72 / 86.96 / 92.17, 1.02 | 1 of 9 | |
| BALLSstick | conduction | 1.58 / 1.89 / 2.19, 0.41 | - / - / - | 6 of 9 | |
| BALLSstick | gap P-C | | | | -0.14 / 0.64 / 0.85 |
| Hasegawa A | peretz | 1.20 / 2.49 / 2.14, 0.14 | 6.29 / 10.17 / 10.57, 0.15 | 3 of 9 | |
| Hasegawa A | conduction | 0.34 / 0.72 / 0.61, 0.13 | - / - / - | 7 of 9 | |
| Hasegawa A | gap P-C | | | | -0.27 / -0.09 / 0.14 |


### 10c. Calibration spread vs the model gap, HIGH charge (Chunc 10 g, head basket 8 g, BALLSstick 10 g, Hasegawa A 30 g)

Same definitions as 10b.

| case | model | k spread at T=756 (0.26..0.32): first / t50 / t90 [ms], peak [MPa] | T_ign spread at k=0.271 (756..950): first / t50 / t90 [ms], peak | combos that never light | gap P-C at nominal first / t50 / t90 |
|---|---|---|---|---|---|
| Chunc fwd plenum | peretz | 0.13 / 0.17 / 0.20, 0.34 | 0.83 / 1.06 / 1.17, 1.83 | 0 of 9 | |
| Chunc fwd plenum | conduction | 0.12 / 0.14 / 0.16, 0.32 | 0.76 / 0.86 / 1.14, 1.96 | 0 of 9 | |
| Chunc fwd plenum | gap P-C | | | | 0.13 / 0.19 / 0.19 |
| Chunc head basket | peretz | 0.23 / 0.30 / 0.31, 0.28 | 1.36 / 1.82 / 1.87, 1.98 | 0 of 9 | |
| Chunc head basket | conduction | 0.21 / 0.26 / 0.24, 0.24 | 1.22 / 1.43 / 1.45, 1.69 | 0 of 9 | |
| Chunc head basket | gap P-C | | | | 0.19 / 0.42 / 0.39 |
| BALLSstick | peretz | 0.31 / 0.42 / 0.49, 0.32 | 2.22 / 2.87 / 3.44, 2.11 | 0 of 9 | |
| BALLSstick | conduction | 0.30 / 0.33 / 0.43, 0.34 | 2.21 / 2.48 / 11.71, 2.70 | 0 of 9 | |
| BALLSstick | gap P-C | | | | 0.23 / 0.41 / 0.30 |
| Hasegawa A | peretz | 0.16 / 0.19 / 0.20, 0.19 | 0.98 / 1.25 / 1.34, 1.20 | 0 of 9 | |
| Hasegawa A | conduction | 0.15 / 0.17 / 0.15, 0.15 | 0.87 / 1.03 / 1.05, 1.09 | 0 of 9 | |
| Hasegawa A | gap P-C | | | | 0.14 / 0.23 / 0.24 |


(Per-combination detail: `results/sm2/sm2_tables.md` tables 10 and 10d.)

## 9. Factor-by-factor attribution

### 11. Factor-by-factor attribution (first light and t50, ms; Peretz ODE / conduction)

Each factor entry is the absolute change of that metric caused by the factor, for the Peretz ODE and for conduction ("P / C"). "model gap" is the Peretz-ODE minus conduction difference itself (closed loop, then the open-loop part, then feedback = closed - open). Calibration = k_solid 0.26..0.32 at T=756 and T_ign 756..950 K at k=0.271 (high-charge runs, KX); h closure = Gnielinski bore h instead of Peretz (H); DeMar on (D); charge = high charge (X) and auto-sized (C); grid = 200 cells (G).

**first light**

| case | model gap closed | model gap open (mean) | feedback | k 0.26..0.32 | T_ign 756..950 | h closure | DeMar on | charge: high | charge: auto | grid 200 |
|---|---|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | -0.05 | -0.03 | -0.02 | 0.13 / 0.12 | 0.83 / 0.76 | 4.24 / 4.34 | 0.07 / 0.52 | 3.59 / 3.77 | 6.34 / 3.03 | 1.11 / 0.92 |
| Chunc head basket | 0.39 | 0.44 | -0.05 | 0.23 / 0.21 | 1.36 / 1.22 | 3.73 / 3.42 | 1.33 / 1.27 | 2.53 / 2.33 | - / - | 0.98 / 0.80 |
| BALLSstick | -0.14 | -0.08 | -0.05 | 0.31 / 0.30 | 2.22 / 2.21 | 4.66 / 4.87 | 0.00 / 0.56 | 2.84 / 3.21 | 7.23 / 3.58 | 1.52 / 1.63 |
| Hasegawa A | -0.27 | -0.19 | -0.08 | 0.16 / 0.15 | 0.98 / 0.87 | 4.05 / 4.40 | 1.11 / 1.78 | 3.37 / 3.79 | 8.91 / 0.01 | 0.68 / 0.55 |

**t50**

| case | model gap closed | model gap open (mean) | feedback | k 0.26..0.32 | T_ign 756..950 | h closure | DeMar on | charge: high | charge: auto | grid 200 |
|---|---|---|---|---|---|---|---|---|---|---|
| Chunc fwd plenum | 0.68 | 0.35 | 0.32 | 0.17 / 0.14 | 1.06 / 0.86 | 5.47 / 4.85 | 0.12 / 0.48 | 5.01 / 4.52 | 11.37 / 4.18 | 1.31 / 0.85 |
| Chunc head basket | 1.14 | 0.29 | 0.84 | 0.30 / 0.26 | 1.82 / 1.43 | 4.85 / 3.83 | 0.33 / 0.29 | 3.62 / 2.90 | - / - | 1.10 / 0.80 |
| BALLSstick | 0.64 | 0.37 | 0.26 | 0.42 / 0.33 | 2.87 / 2.48 | 5.80 / 5.22 | 0.08 / 0.26 | 3.95 / 3.72 | 12.96 / 4.42 | 1.78 / 1.64 |
| Hasegawa A | -0.09 | 0.33 | -0.42 | 0.19 / 0.17 | 1.25 / 1.03 | 5.79 / 6.02 | 0.73 / 2.03 | 5.49 / 5.81 | 13.19 / 0.00 | 1.13 / 0.87 |


Reading: in the clean regime the solid-model gap in first light is 0.13-0.23 ms (open loop; feedback adds nothing) and in t50/t90 0.2-0.4 ms (of which feedback is most). Every other factor is larger: calibration T_ignition 0.8-2.2 ms (4-10x),
h closure 3.4-4.9 ms (>20x), charge 2.3-3.8 ms, grid 0.5-1.6 ms, k_solid within the LHS range ~ the gap. Near the threshold the models differ in threshold position and in the number of cells lit, which dominates the metrics there.

## 10. Recommendation (the user decides)

- **Peak and lighting time barely depend on the choice** once the igniter is above threshold (first light 0.13-0.23 ms, t50/t90 0.2-0.4 ms, peak <= 1.8 %, Hasegawa A late peak 1.6 %). The decision should therefore rest on correctness and on the near-threshold behaviour, not on those metrics.
- **Conduction is the exact model of the stated problem** (verified to 0.3-0.8 K, SM1 section 2), is flux-driven (no long creep tails: typically at most 2-7 ms within 100 K of threshold, one 11.7 ms cell, vs up to 194 ms for the Peretz ODE), has no new parameters, and costs 7-20 % wall time (SM1 section 3.4).
  The Peretz ODE is free but is only exact for constant h, and in the marginal cells it lights by holding temperature: it lights grain at T_ign 850-950 K in cases where the exact solution does not light it at all. If the question is "which one will I be able to defend", it is conduction.
- **If the Peretz ODE is chosen** (e.g. for cost), the evidence here says it is acceptable above the ignition threshold and unreliable near it; its marginal lighting events should not be trusted, and the dead-zone cells (3-11 BALLSstick, 3-10 Hasegawa A) are exactly such events (SM1).
- **Either way T_ignition (and k_solid) must be re-fitted for the chosen model, and the dominant lever is T_ignition**, not the solid model: its band (750-950 K) moves first light 4-10x more than the model choice. Re-fit targets: first-light timing / ignition pressure rise, not the plateau.
- **The IG3 auto-sizing is not model-neutral:** it picks a different charge per solid model (Chunc 1.68 / 2.10 / 2.62 g for Peretz ODE / conduction / Goodman) and sits at the threshold, where the startup is hypersensitive (IG3 DEVNOTES). Any model switch changes the auto charges and the startup peaks with them. A fixed charge, or a margin above threshold, is needed for fair validation comparisons.

## 11. Gates

- **B0** (`baseline_gate.py`, default `solid_heating_model='goodman'`, against IG3's recorded baseline, no re-record): **PASS**
  (`results/sm2_B0_output.txt`): Chunc peak 8.698446 (rel 4e-15), plateau 8.066193 (9e-16), n_steps 533513 exact; BALLSstick peak/trough/mid rel 4e-16 .. 3e-15, n_steps 850245 exact;
  Hasegawa A late peak 5.362271 (1.6e-11), t_peak 2.681822 s, n_steps 2347427 exact. The recorder and the diagnostics do not alter the default trajectory.
- **Full suite** (`results/sm2_suite.txt`): **827 passed** (the post-SM1 base was 821; +6 new in `tests/test_sm2_replay.py`: recorder leaves the trajectory untouched, ignition-diagnostic contents,
  self-replay reproduces the closed-loop ignition times for all three models, lighting metrics). The direct callers of `_goodman_ignition_sources_and_mass` in two test files were extended with the 7 new trailing arguments (internal API change).

## 12. Not done / limits

- Open-loop replay cannot give the peak (it needs the flow); the peak differences are closed-loop only, and near threshold they are dominated by the number of cells lit.
- The calibration grid is 3 x 3 (k, T_ign) in short windows (0.2-0.3 s), at the fixed and high charges; no finer grid, no re-fit.
- The Hasegawa A late peak in the high-charge row is the 30 g charge's, not the validated 14.98 g one; the 0.83x / 0.92x numbers (table 1) are the fixed-charge ones.
- Forcing recorded only for the first 120000 steps (~0.2-0.4 s) of the base and high-charge runs and 60000 steps for the H / D / C runs; cells lit after that count as never in the open-loop tables (the windows are listed in table 2).
- Fill Mach / total-energy (thread F) was recorded, not fixed.
- The Goodman rows are references only.
