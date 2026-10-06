# R1 - Review of the friend's Codex commits (15) on srm_1d

Reviewer: R1 agent, 2026-10-05. Their code: `C:/tmp/sv/cx` (detached at `codexfork/main`, base `db0a980`,
which IS an ancestor of our `openmotor-frontend`: `git merge-base db0a980 openmotor-frontend` = `db0a980`, 33 commits behind).
Their openMotor: `C:/tmp/sv/cx-om`. Neither repo was changed. Nothing committed or pushed.
Their pytest suite was NOT run (spec allows only the Q1 grid and an optional Q5 probe). Their tests were read, not executed.

Runs performed (all in `C:/tmp/sv/cx`, env as in spec; preflight printed
`srm_1d C:\tmp/sv/cx\srm_1d\__init__.py`, `motorlib C:\tmp\sv\cx-om\motorlib\__init__.py`).
Scripts are in `C:/tmp/sv_cx/` (`q1_run.py`, `q_extract.py`, `q5_run.py`). Every run is Chunc (`motors/machbusterNew.ric` of THEIR tree),
50 ms startup, `cfl 0.3` unless stated. Total wall about 6 min. Windows/py3.10.5: no portability problem hit
(their lockfile is macOS-only and was not used; README commands use `../.venv-srm/bin/python`, which does not exist here).

## 0. Headline findings

1. **Their code is numerically identical to ours at default, canonical settings.** Their tree + canonical knobs
   (roughness 32 um, kappa 0.44, T_ign 756, k_solid 0.271) + library BPNV, 100 cells gives peak `11.472831280605046` MPa,
   t_peak `0.016454109339050758` s, max pre-peak Mach `3.944021091805623` - the exact B0 / T7 values (`T6.md` gates, `T7.md` base row).
   So their diagnostics and bracket change did not perturb the default trajectory.
2. **Q1 is fully explained and it is NOT a physics disagreement.** Their fill Mach is `max(max_mach[t <= 3 ms])` on a STALE, non-canonical
   Chunc config (MTV, 35 um, kappa 0.45, T_ign 900, k_solid 0.4, P_cutoff 50 kPa). Ours is `max(max_mach[t < t_peak])` on canonical knobs.
   Reproduced 2.2002 (their "current BPNV 2.200") and 3.944 (ours) in the same tree, differing only in knobs + window.
3. **Their "legacy BPNV is the cause of the supersonic fill" conclusion is overstated.** Over the full pre-peak window the legacy/library ratio is
   4.82/3.94 = 1.22 (canonical) not 2.2x, and library BPNV and MTV are both still supersonic (3.69 / 3.21) later in the fill than their 3 ms window.
4. **Real bug in OUR tree, found by their commit `e74a937`: slag coefficient unit conversion.** openMotor's `slagCoeff` is `(m*Pa)/s` and
   `slagRate = slagCoeff/P[Pa]`; ours divides by P in MPa but takes the .ric value unconverted (`srm_1d/openmotor_adapter.py:1151` on `openmotor-frontend`),
   so any nonzero slag is 1e6x too large. Default results do not change (all 28 committed `.ric` have `slagCoeff: 0.0`).
5. **Q5: the 806-cell pressure-floor events are TEMPORAL (timestep) in both their config and our canonical config** - reproduced in canonical:
   193 floor cell-steps at CFL 0.3 (raw P -734 Pa, matches T3) -> 0 at CFL 0.15 and 0.075. Details in Q5.
6. **Codex's "29.84 us" refresh delay IS reproducible** (we said it was not in T5/T7): it is a per-ignition-event delay to the next burn-rate refresh,
   not the max refresh gap. Reproduced `2.9842e-5 s` on their 100-cell MTV startup. Their energy numbers also reproduce (Q3).

## 1. Per-commit review (in order)

Legend: default change = does the default simulation result change.

### 1. `3a6ff95` Add reproducible baseline registry and workspace locks - **COMPARE** (adopt nothing wholesale)
- Does: `cases/` (JSON registry of 8 validation cases with sha256 of inputs + measurement schemas + quality flags, a loader that refuses raw/missing traces),
  `cases/baseline_runner.py` (runs a frozen config, writes `run_manifest.json` with git state, dependency versions, input hashes, metrics, CSV, plots),
  `requirements/` (macOS-x86_64 pip lock, `workspace-lock.json` pinning srm_1d `db0a980` and openMotor `92e425f`), `scripts/validate_case_registry.py`, `scripts/verify_workspace.py`.
- Correctness: default change: none (no solver edits). `baseline_configs.json` freezes the Chunc case from the stale template `examples/run_chunc.py`:
  `pyrogen mtv, roughness 35e-6, kappa 0.45, T_ignition 900, k_solid 0.4, P_cutoff 50 kPa` (`cases/baseline_configs.json:39-46`) - NOT the v0.7.5 canonical
  defaults (32 um / 0.44 / 756 / 0.271). Every downstream script (spatial extension, limiting cases, floor probes, fill comparison "current") inherits that.
  The runner reports metrics only; there is no pass/fail comparison to a recorded value. Registry sha256 pins are brittle: against our committed tree
  8 of 35 checksums differ (7 distinct `.ric`: hasegawa_a, zerox, machbusterNew, BALLSstick, ST2.0, hlbltest, FatBoy; our igniter refresh `5b7db7d` re-saved them; line-ending effects not ruled out).
  `workspace-lock.json` pins an openMotor commit that predates the igniter provider.
- Tests: `tests/test_baseline_runner.py` (config matches example constants; dry-run manifest has no reproduction claim), `tests/test_case_registry.py`. Metadata-level only; no solver run.
- Overlap: our B0 gate (`probes/baseline_gate.py`). See Q6.
- Verdict: COMPARE. Optional small adopt (see port plan P7): manifest provenance fields. SKIP lockfile, workspace verifier, checksum registry.

### 2. `e74a937` Add numerical convergence verification tooling - **PORT** (slag fix) + COMPARE (rest)
- Does: (a) `convert_nozzle`: `slag_ours = slagCoeff / 1e6` + test (`srm_1d/openmotor_adapter.py:974`, `tests/test_adapter.py::test_slag_unit_conversion`);
  (b) `run_simulation(history_capacity=)` exact preallocation override; (c) `scripts/run_numerical_verification.py` (Chunc grid/CFL study, startup + full profiles,
  `_extract_metrics`); README/lock text.
- Correctness: (a) CORRECT. openMotor `motorlib/motor.py:441-447`: `slagRate = (1/pressure) * slagCoeff` with pressure in Pa and units `(m*Pa)/s` (`nozzle.py:25`);
  srm_1d `simulation.py` step 3b: `s_rate = slag_coeff / P_throat_MPa`; hence `slag_ours = slagCoeff/1e6`. Our tree is unfixed (line 1151, docstring lines 32/1145 say "same convention").
  Default change in theirs: none for committed motors (all slag 0). (b) harmless; ours already has `diagnostic_history_capacity`. (c) fine as a probe; the Chunc config is the stale one (above).
- Tests: slag test is exact (2e-3 -> 2e-9). `history_capacity` tests fine.
- Overlap: none for (a). (b) equivalent to our `diagnostic_history_capacity`. (c) overlaps our T7 probe.
- Verdict: **PORT the slag conversion + test (bug fix, XS)**. SKIP (b). (c) COMPARE (feeds P5).

### 3. `d08bcc5` Add per-cell numerical limit diagnostics - **PORT (adapted)**
- Does: passes `limit_activation_counts/duration/abs_energy` arrays (3 x N) into `_piso_step_with_energy_diagnostics`; counts every activation of the **temperature floor, temperature ceiling and pressure floor** per cell,
  with duration and |correction energy| (T limits); returns `result['numerical_limits']`.
- Correctness: read-only tallies; the clip logic is unchanged (`solver.py` diff only adds a `limit_index` and counters). Default change: none. In my canonical 100-cell/50 ms run:
  `limit_temperature_ceiling_total_activations = 254483` (8.3% of the 3.05e6 cell-steps), `temperature_floor 1975`, `clipping_active_steps 10140`, `integrated |clip energy| 532 J`
  (their MTV config: 62215 / 13944 / 7381 / 307 J). Our tree has NO per-cell or per-limit T-clip accounting (only the net signed `clipping_correction_power` and the T3 scalar pressure-floor counters),
  so this is new information relevant to thread F. Risk: None-default optional arrays under numba (first/last-time arrays are indexed whenever `limit_activation_counts is not None`; passing counts without the time arrays would fail to compile). Internal API only.
- Tests: `test_temperature_ceiling_clip_is_localized` (single-cell PISO step with forced ceiling clip: count/duration/energy) - good, direct. Run-level test checks shapes only.
- Overlap: our T3 (`pressure_floor_*` scalars + mass/energy added + `pressure_floor_min_raw_P`) covers only the pressure floor. Theirs adds spatial/temporal resolution and the T limits; ours adds the floor's mass/energy bookkeeping (theirs does not account the mass the clamp adds). See Q4.
- Verdict: PORT (adapted to our 21-tuple return / T6 histories), default-neutral.

### 4. `d773b55` Track numerical limit activation windows - **PORT (fold into #3)**
- Does: first/last activation time per cell for the three limits; per-**face** activation count/duration/first/last for the optional `port_mach_cap`; passes `step_start_time_s`.
- Correctness: read-only. Face indexing `j-1` for interior faces matches the `N-1` sized arrays. Default change: none. `port_mach_cap` defaults off, so the Mach-cap part is dead weight at default
  (our REOPENED thread E/A used it as a no-op probe).
- Tests: first/last time asserted in the single-cell test; `test_enabled_port_mach_cap_reports_face_activations` (cap 1e-6 forces activations). Good.
- Overlap: T3 has first-activation time (scalar) only.
- Verdict: PORT first/last times; Mach-cap face counters optional (cheap, same code).

### 5. `e46562b` Add source and geometry cadence diagnostics - **COMPARE (one piece to PORT)**
- Does: `geometry_update_interval` kwarg (None = follow burn interval), `ValueError` if burn or geometry interval <= 0, summary keys, ARCHITECTURE text on step ordering;
  in the script `_ignition_refresh_delays` (per-ignition delay to the next burn refresh).
- Correctness: identical mechanism to our T5 (`simulation.py:1713` `if step % geometry_update_interval == 0`). Default change: none. They validate `burn_update_interval <= 0`; ours does not
  (`burn_update_interval=0` would be a `step % 0` inside the njit loop) - tiny hardening to port.
  `_ignition_refresh_delays` (`scripts/run_numerical_verification.py:139`): for each ignited cell, `time[next refresh step] - ignition_time`, `next = (ignition_step // interval + 1) * interval`.
  This is the statistic behind Codex's "29.84 us" and it reproduces: `ignition_rate_refresh_delay_max_s = 2.9842e-5` (MTV, 100 cells, 100 ignitions, median 1.33e-5 s); canonical BPNV: 2.875e-5 s, median 9.9e-6 s.
  Our T5/T7 "max refresh gap" (143.8 us Chunc) is a different, larger statistic (max time between refreshes), so the "not reproduced" note in T5/T7 should be amended: not a discrepancy.
- Tests: validation tests plus summary-key assertions. Fine.
- Overlap: our T5 + T7. Verdict: COMPARE; PORT the validation + the delay statistic as a post-hoc function (P4).

### 6. `e749f71` Verify and expand erosive closure brackets - **COMPARE** (details in Q2)
- Does: `_burn_rate_closure_residual` helper; keep `r_hi = 20*r0` when it brackets, otherwise double up to 60 times, else analytic `r0 + ebf*h0`; `D_hyd <= 0` returns `(r0, 0)`; tab-boundary tests.
- Correctness: sound. Differences from our T2 in Q2. Default change: none for non-clamped states (see Q2 for the one caveat).
- Tests: `test_bisection_expands_upper_bracket_and_closes_residual` (2 states, asserts `r_total > 20 r0` and residual < 1e-7) is a good direct regression; the `TestPressureTabs` tests are 2 cases vs our 81 (T4).
- Verdict: COMPARE; port only the `D_hyd <= 0` guard (hardening bundle P1).

### 7. `2b2c776` Add manifested limiting-case verification - **SKIP** (covered by T4/T7)
- Does: Chunc 100-cell 3 s baseline vs `diagnostic_disable_erosive`, a 0.6-1.8 s stationarity screen (CV <= 5%, |slope| <= 5%/s), closure checks at zero crossflow / Re 5e5.
- Correctness: the closure check "finite crossflow: normal rate preserved" is tautological (`r_erosive = r_total - r0` inside `burn_rate_cell`, so `r_total - r_erosive == r0` by construction except for the `max(.,0)` clamp);
  our T4 does the same identity over 56 grid states with a non-vacuity guard. The stationarity screen duplicates T7 finding 7 (Chunc [0.4,0.9] s: +11.5 %/s; theirs fails on 0.6-1.8 s, consistent). Default change: none.
- Tests: helper-only unit tests (stationarity metric, relative change, resume checksum refusal); no simulation test.
- Overlap: T4 (b)(c)(d) and T7. Our T4 worklog already compared 12.1574 -> 8.1357 MPa. Verdict: SKIP.

### 8. `79867c0` Audit sensible enthalpy balance - **COMPARE** (details in Q3)
- Does: new `convective_wall_heat_sink_power` channel (debit previously unrecorded), `energy_convention` summary dict, `sensible_enthalpy_audit()` (source ledger, nozzle telescoping check, undecomposed between-step storage jump).
- Correctness: same ledger as our T6 (we added the same channel as `convective_wall_sink_power`). Default change: none (bookkeeping). Reproduced numbers in Q3.
- Gap: ledger is blind to the left-neighbour radiation sink leak that T6 found (they only ran emissivity 0).
- Tests: synthetic-array test of the audit function (reconstruction arithmetic), summary key assertions. No physics test (ours has 14 incl. PISO identity).
- Verdict: COMPARE (ours is a superset). Nothing to port.

### 9. `cab26d0` Reconcile Chunc startup fill configurations - **COMPARE** (Q1; informs T8-B)
- Does: 2x2x3 factorial: pyrogen {mtv, library bpnv, legacy bpnv (`replace(load_pyrogen('bpnv'), T_flame=2800, M=0.030, gamma=1.25, gas_mass_fraction=1.0)`)} x knobs {current stale, "historical" 32 um/0.44/756/P_cutoff 10 kPa} at 50/100/200 cells.
- Correctness: reproducible and cleanly separates igniter from knobs - the legacy-marker construction is useful. Flaws: (i) the "historical knobs" variant still leaves `k_solid = 0.4` and the old `.ric`; (ii) the Mach metric window is 3 ms (see Q1); (iii) "current" is non-canonical. Default change: none.
- Tests: variants/legacy marker unit tests only. Verdict: COMPARE; reuse variant construction for T8-B (P5).

### 10. `8b0619c` Extend Chunc startup spatial verification - **PORT (adapt)** for T8-A (Q7)
- Does: ladder 100/200/400/800/1600 target cells, 50 ms, metrics peak/t_peak/fill Mach/min P, successive-change screen, Richardson order for monotone triplets (`_approximate_observed_order`), floor-censoring flag, resume from JSON with checksum guard.
- Correctness: order formula is the standard `p = log2((f2-f1)/(f3-f2))`, extrapolation `f + d/(2^p-1)` - correct for ratio 2 (they note target-cell vs actual-cell caveat). Withholds order for non-monotone series (good). Default change: none.
- Tests: pure-function tests of the order/changes. Good enough for the helper.
- Overlap: our queued T8 (`TASKS.md:401`). Verdict: PORT/adapt.

### 11. `86c587f` Instrument pressure-floor sensitivity - **PORT (adapt, small)**
- Does: `run_simulation(pressure_floor_pa=1e3)` (validated finite > 0) threaded into the PISO kernel as a runtime value; per-cell `maximum_raw_deficit_pa`; `scripts/run_pressure_floor_sensitivity.py` (800/1600 cells x floors 100/1000/10000 Pa).
- Correctness: at the default 1000 Pa the clamp is identical. Overriding the floor changes the physics of the run (it is a limiter-sensitivity knob) - correctly labelled diagnostic. Default change: none.
- Tests: validation parametrized (0, -1, NaN, inf), PISO test with floor 750 and deficit 250. Good.
- Overlap: T3 counters (we log min raw P globally; they log max deficit per cell). Verdict: PORT the kwarg + validation (it is the cleanest experiment for T3's -734 Pa and thread F). SKIP the sweep script as written (it runs the stale config).

### 12. `55b19cb` Trace pressure-floor correction events - **SKIP in-kernel; keep as a possible out-of-kernel replay**
- Does: 225 lines inside `_piso_step_with_energy_diagnostics` STEP 4 recording a 30-field ledger (west/east mass flows, continuity RHS, local diagonal, nozzle state, predictor velocities, densities) at each cell's worst-deficit event.
- Correctness: read-only, but most fields are RECONSTRUCTED from `u_star`, `P_prime`, `d_face`, `rho_new_1` rather than copied from the solver's own intermediates (e.g. `first_west_mdot = rho_w*u_star*A_face`), so they can drift from what PISO actually used; I did not verify each against the solver algebra. Their own doc says the row cannot attribute a global tridiagonal correction to one cell. Default change: none. Maintenance cost: high (solver-internal coupling, many positional args).
- Tests: shape and one single-cell value test.
- Overlap: none (we have no event trace). Verdict: SKIP for now; if thread F test #1 needs a local ledger, write an out-of-kernel replay probe instead.

### 13. `2d664ff` Record Chunc pressure-floor event data - **COMPARE (data point, not portable)**
- Docs + `docs/v0_8_0/data/chunc_pressure_floor_event_probe_806.json`. MTV stale config. Event at gap cell 803, x = 0.8164 m (aft bore -> 47.63 mm gap expansion, 2 cells upstream of the nozzle), raw pressure -3882 Pa, first correction -1161.6 Pa, nozzle in the subsonic-outflow state. Verdict: cite, do not port the file.

### 14. `32d9230` Record pressure-floor CFL sensitivity - **COMPARE (data point)**
- `chunc_pressure_floor_cfl_probe_806.json`: 806 cells, CFL 0.3/0.15/0.075 -> floor activations 110 / 0 / 0, min P 1.000 / 10.097 / 12.081 kPa. Config stale (MTV, 35 um, 0.45, 900 K, k_solid 0.4). See Q5. Verdict: cite.

### 15. `9ee1a91` Record floor-free finer CFL probe - **COMPARE (data point)**
- Adds the CFL 0.075 row and the interpretation block (`floor_free_solution_established: false`, next test = cadence + grid from the floor-free 0.075 baseline). Verdict: cite.

Summary table

| # | commit | verdict | changes default result? |
|---|---|---|---|
| 1 | 3a6ff95 registry / locks / runner | COMPARE | no |
| 2 | e74a937 slag fix + history_capacity + verification script | PORT (slag), COMPARE (rest) | no (committed motors have slag 0) |
| 3 | d08bcc5 per-cell limit diagnostics | PORT (adapted) | no |
| 4 | d773b55 activation windows + Mach-cap faces | PORT (with #3) | no |
| 5 | e46562b cadence diagnostics | COMPARE (+ port validation, delay statistic) | no |
| 6 | e749f71 bracket expansion | COMPARE (+ D_hyd guard) | no |
| 7 | 2b2c776 limiting cases | SKIP | no |
| 8 | 79867c0 enthalpy audit | COMPARE (T6 superset) | no |
| 9 | cab26d0 fill configurations | COMPARE (informs T8-B) | no |
| 10 | 8b0619c spatial extension | PORT (adapt for T8-A) | no |
| 11 | 86c587f floor override + deficit | PORT (kwarg, small) | no at default |
| 12 | 55b19cb floor event trace | SKIP | no |
| 13 | 2d664ff event data | COMPARE (cite) | n/a |
| 14 | 32d9230 CFL data | COMPARE (cite) | n/a |
| 15 | 9ee1a91 CFL data | COMPARE (cite) | n/a |

## 2. Answers

### Q1 Fill Mach: configuration and metric behind their numbers
Source: `scripts/run_fill_configuration_comparison.py` (`cab26d0`) + `_extract_metrics` in `scripts/run_numerical_verification.py`.
- **Metric:** `max_fill_mach = max(result['max_mach'][time <= FILL_WINDOW_S])`, `FILL_WINDOW_S = 0.003` (`run_numerical_verification.py:35,159`).
  `result['max_mach']` is the max over **cells** of `abs(Mach[i])` (cell-centred), `simulation.py:2088-2108`. Same history array as ours.
- **Window:** first 3 ms of the run. Ours (T7) is `max_mach` over `t < t_peak` (whole pre-peak window; Chunc t_peak = 16.5 ms).
- **Igniter source:** `pyrogen='bpnv'` -> library `pyrogens/bpnv.yaml` (post-`bfc2f3f`, has `gas_mass_fraction`) for "current"; `LEGACY_BPNV_MARKER` -> in-script `replace(load_pyrogen('bpnv'), T_flame=2800, M=0.030, gamma=1.25, gas_mass_fraction=1.0)`.
  The embedded `.ric` igniter block is never used (a `pyrogen` argument is always passed). The legacy fixture matches the stale embedded BPNV we found in T1.
- **Other config:** `target_propellant_cells` 50/100/200, `cfl 0.3`, `t_max 0.05`, and the baseline_configs Chunc knobs: roughness 35 um, kappa 0.45, T_ign 900 K, k_solid 0.4, P_cutoff 50 kPa (`cases/baseline_configs.json:39-46`) - not canonical.
  Ours: canonical knobs 32 um / 0.44 / 756 K / 0.271.
- **Reproduction (100 cells, their tree, `C:/tmp/sv_cx/q1_run.py`):**

| variant (100 cells) | their metric (t<=3 ms) | pre-peak max (our metric) | time of max | P_peak MPa |
|---|---|---|---|---|
| their `current_bpnv` (stale knobs) | **2.2002** (theirs: 2.200) | 3.694 | 3.69 ms | 11.981 |
| their `current_legacy_bpnv` | **4.8207** (theirs 4.821) | 4.8207 | 1.55 ms | 13.313 |
| their `current_mtv` | **0.4976** (theirs 0.498) | 3.2105 | 6.94 ms | 12.157 |
| canonical knobs + library BPNV | 3.9440 | 3.9440 | 2.48 ms | 11.4728 |

  Canonical row equals our T7/B0 values to all printed digits (peak 11.472831280605046, t_peak 0.016454109339050758, Mach 3.944021091805623).
- **Why 2.2 vs 3.94:** (a) with the stale knobs the Mach excursion happens at 3.69 ms, just outside their 3 ms window, so they report the earlier, smaller value (3.694 would be the max); (b) knobs (T_ign 900 vs 756) delay ignition by about 1.2 ms. Neither is a cell/face or igniter-source difference.
  Same effect for MTV: 0.4976 inside the window but 3.21 at 6.9 ms.
- **Consequence:** their statement "the Section-6 supersonic fill is mainly the pre-`bfc2f3f` BPNV gas" needs qualification. Legacy vs library ratio on the full pre-peak window: 4.82/3.69 = 1.31 (their knobs), 1.22 vs canonical.
  The fill remains supersonic at 100 cells with every igniter (3.2-4.8). Grid divergence: canonical BPNV at 806 cells reaches 17.5 at CFL 0.3 (Q5).
  Ours (embedded legacy 4.81 at 1.38 ms and library 3.94, `TASKS.md` T8 text) is consistent with these.
- Ambiguity: I did not rerun their 50/200-cell points; the 100-cell reproductions match their published 100-cell numbers to 3-4 digits.

### Q2 Bracket: `e749f71` vs our T2
- **Theirs** (`srm_1d/burn_rate.py:482-500`): keep `r_hi = 20 r0`; if `F(20 r0) < 0` double `r_hi` up to 60 times; if still negative fall back to `r0 + max(ebf*h0, 0)`. Bisection stays at 30 iterations with the `1e-8` break.
  Also `D_hyd <= 0 -> (r0, 0)` and a guard `r_hi <= r_lo -> r_lo + max(ebf*h0, 1e-12)`.
- **Ours** (`burn_rate.py` on `openmotor-frontend`, lines ~660-704): if `F(20 r0) < 0 and r0 > 0` set `r_hi = r0 + ebf*h0` (analytic, guaranteed since `h(r) <= h0`); iteration cap 30 -> 60 with the same `1e-8` break; `r0 == 0` keeps the degenerate bracket and returns `(0, 0)`.
- **States one handles and the other does not:**
  - `D_hyd <= 0`: theirs returns `(r0, 0)`; ours has no guard (unreachable at runtime: burning cells have a positive port; defensive only).
  - `r0 == 0`: ours returns `(0, 0)`; theirs widens to `[0, ebf*h0]` and returns a nonzero erosive rate at zero normal burn (non-physical; unreachable since the 1 kPa floor makes `r0 > 0`).
  - Precision when the root is far above `20 r0`: theirs bisects `[r0, 20*2^k*r0]` with only 30 halvings (width up to 2x the root, so error up to ~root*2^-29; their 150-point envelope closed to 5.73e-8 m/s); ours bisects `[r0, r0+ebf*h0]` with up to 60 halvings (our T2: 8.74e-9 m/s). Both are far below any physical scale.
  - Cost: theirs adds up to 60 extra residual evaluations in clamped states only; ours a single one.
- **Bit-identical on non-clamped states?** Essentially yes for theirs: when `F(20 r0) >= 0` the expansion loop does not execute and the bisection is the original (with the break firing before 30 halvings). One non-bit-identical spot: the `F_lo` sign test now uses the Taylor `h0*(1-beta/2)` branch for `beta < 1e-6` where the old code used `h_lo = h0` (`F_lo` differs by <= 5e-7 relative, sign cannot differ unless `ebf*h0 ~ 0`). Their evidence is the Chunc baseline; I observed identical Chunc B0 digits (headline 1). Ours: 3210/3210 bit-identical (T2 worklog).
- Verdict: equivalent in intent; no port needed beyond the `D_hyd <= 0` guard.

### Q3 Energy: `79867c0` vs our T6
- **Reproduced their numbers** (their tree, MTV stale config, 100 cells, 50 ms; `C:/tmp/sv_cx/q_extract.py`): source ledger `1.12e-8 W` (theirs 1.30e-8), nozzle telescoping `1.86e-8 W` (2.05e-8), `energy_residual` max `2.0e-5 W` (1.79e-5), **storage jump max 2.627 J (theirs 2.63), sum|.| 13.297 kJ (13.30)**.
  So their jump numbers are for a 50 ms startup, not the full burn; our T6 totals (Chunc 57.4 kJ over 1 s; fill 4.46 kJ + spike 1.37 kJ) are on a different window.
  With canonical knobs + library BPNV: max jump **2.726 J** (== T6 "total max single 2.726 J"), sum 12.68 kJ over 50 ms, source ledger 1.12e-8 W.
- **What each measures that the other does not:**
  - Theirs only: nothing numerically new. Their `energy_convention` summary string lists excluded terms (we have `energy_omitted_terms`).
  - Ours only: (i) decomposition of the between-step jump into EOS reconciliation / pressure floor / mixture refresh / geometry refresh / remainder (remainder 1e-5 J), which shows geometry refresh is always positive (+40.8 kJ Chunc) and EOS carry dominates the max and the fill; (ii) the PISO-vs-reconstructed nozzle mismatch (2.2-3.4 kW in fill, a by-design explicit/implicit difference) and mass counterpart; (iii) the omitted kinetic-energy magnitude (up to 80% locally in the fill); (iv) the radiation left-neighbour sink leak (292 W vs 570 W at emissivity 0.05) which their emissivity-0 ledger cannot see; (v) the integrated closure with jumps (1e-15) vs without (0.51%, 30% in the fill).
- **Numerical disagreement:** none. Their "nozzle 2.05e-8 W" is the telescoping flux check (`convective_scalar_flux_power + nozzle_enthalpy_power`, `ignition_diagnostics.py` in `79867c0`), the same one T6 notes as a different comparison from its PISO-vs-STEP-6 nozzle ledger.
- Verdict: COMPARE; T6 is a superset. The ledger channels agree.

### Q4 Limit diagnostics: what `d08bcc5`/`d773b55` track beyond T3
Theirs adds (per cell: activation count, accumulated active time, first/last activation time; plus |correction energy| for T limits): **temperature floor** (`max(1, T_ambient)`), **temperature ceiling**, pressure floor, and per **face** the optional `port_mach_cap`. Later commits add per-cell max pressure deficit and the event trace.
Ours (T3/T6) tracks the pressure floor only, with scalar/per-step activations, mass and energy added by the clamp, min raw P, and first time.
- Not tracked by either: dt limiter (CFL vs source-cap), nozzle state, bisection bracket events, Z-N.
- Theirs lacks: mass/energy the clamp adds (ours has it); EOS/geometry jumps.
- **Port value: yes for the T clips.** The measured clip activity is large at startup: canonical 100 cells/50 ms: temperature ceiling 254,483 cell-steps (8.3% of cell-steps), floor 1,975, 532 J integrated |correction| (their MTV config: 62,215 / 13,944 / 307 J). The per-cell map tells which cells (the fill front) are clipped and when; relevant to thread F (the clip is a second non-conservative term next to the EOS carry) and to the energy audit. Per-cell pressure-floor maps: modest value (the location is already evident: x = 0.817 m, both configs).

### Q5 Pressure-floor probes
- **What they found** (`PRESSURE_FLOOR_EVENT_PROBE.md`, JSON data): at 806 cells (MTV stale config, 50 ms), CFL 0.3: 2 cells, 110 activations, raw min -3882 Pa, first activation 9.44 ms; CFL 0.15: 0 activations, min P 10.097 kPa; CFL 0.075: 0 activations, 12.081 kPa. Peak P 12.4154 / 12.4198 / 12.4219 MPa; t_peak 23.74 / 23.62 / 23.56 ms. Event cell 803 (x = 0.8164 m) just downstream of the bore-to-gap area jump, local source ~0, nozzle in subsonic-outflow state (so not a choked-boundary effect).
  Their fill Mach column (0.4464 at every CFL) is the 3 ms-window metric and carries no information about the later excursion (see Q1).
- **Is it temporal or spatial in OUR configuration? I re-ran** (their tree, canonical knobs, library BPNV, 806 cells, 50 ms; `C:/tmp/sv_cx/q5_run.py`):

| CFL | steps | floor cell-steps / cells | raw deficit | min P | pre-peak max Mach | P_peak MPa |
|---|---|---|---|---|---|---|
| 0.3 | 246,185 | 193 / 2 | 1734.3 Pa (raw P = -734.3 Pa) | 1000 (clamped) | 17.47 | 12.5107 |
| 0.15 | 487,572 | 0 / 0 | - | 8868 Pa | 6.02 | 12.5255 |
| 0.075 | 973,366 | 0 / 0 | - | 10521 Pa | 5.63 | 12.5328 |

  193 cell-steps and -734 Pa reproduce our T3 exactly (event x = 0.8174 m, t = 1.774-1.779 ms). **The floor events are TEMPORAL**: they vanish at CFL 0.15 in both configurations; peak P moves +0.12% / +0.06%, min P keeps rising (8.87 -> 10.52 kPa, not converged), and the pre-peak fill Mach falls 17.5 -> 6.0 -> 5.6 (-66% then -6%).
  At 100 cells the same trend is milder (T7: 3.94 -> 3.38 at CFL 0.15), and with their MTV config at 100 cells the pre-peak Mach is 3.21 / 3.03 / 3.02 at CFL 0.3 / 0.15 / 0.075 (their 3 ms-window metric stays 0.4976 at all three).
  The remaining Mach spread across grids at fixed CFL 0.15 (3.38 at 100 cells (T7) vs 6.02 at 806) shows a spatial component persists in the fill; the floor events themselves are not spatial.
- **Bearing on thread F and T3's -734 Pa:** a pressure undershoot that disappears when dt is halved is the signature of a lagged/explicit coupling error that scales with dt; thread F (isothermal `rho_new_1 = P/(R*T_old)` closure + later T update + EOS re-derivation) predicts exactly that, and T7 already shows fill quantities not converged in dt.
  It does not prove F: the event sits at an area discontinuity next to the nozzle, where an explicit area-jump flux or the nozzle boundary linearization could also scale with dt. The discriminating test is F's test #1 (iterate T within the corrector at fixed CFL 0.3): if the 193 events at 806 cells disappear without lowering CFL, F is the cause. Their event ledger cannot attribute (they say so). Therefore -734 Pa is a transient numerical undershoot with CFL dependence, not a state; min P itself is a poor convergence metric at present.

### Q6 `3a6ff95` registry / locks / runner vs our B0
- B0 (`probes/baseline_gate.py`) is a regression gate: canonical knobs, three motors, full-precision and exact step counts, pass/fail (rtol 1e-12). Their runner records provenance and metrics with no pass/fail and a non-canonical Chunc config; the registry is a data-availability catalogue with checksums.
- Useful: registry quality flags / readiness classification of measurement traces and a loader that refuses to tare/align raw traces silently (relevant to the open validation work); manifest fields (git state, dependency versions, input hashes) are worth printing from B0.
- Harmful/brittle: sha256 pinning of `.ric` (openMotor-owned files change on re-save; 8/35 already mismatch); macOS lock; workspace verifier pinning db0a980 / 92e425f; `verify_workspace.py` assumes a sibling `../openMotor-compatible`.
- **Recommendation:** SKIP lock/verifier; ADAPT (optional, low priority) the manifest provenance into B0 (print git rev, versions, `.ric` sha256 next to results); defer the registry unless we need a validation-data catalogue (then adapt the quality-flag idea, without checksums).

### Q7 Can T8 reuse `run_spatial_extension.py` / `run_fill_configuration_comparison.py`?
Yes, with edits (~half a day, not a rewrite):
- T8-A (grid ladder to 1600): `run_spatial_extension.py` already has the ladder, resume-by-JSON, successive-change screen, Richardson order helper and floor-censoring. Needed edits: read canonical options (not `baseline_configs.json` Chunc), `t_max` 1.0 not 0.05 (T8 spec), replace `_extract_metrics` (its limit/ledger keys come from `numerical_limits`, which we lack until P2; use our `pressure_floor_*` summary keys instead), add plateau and `fill Mach = max(max_mach[t < t_peak])` plus the 3 ms-window value so both definitions are on record.
- T8-B (fill attribution): `run_fill_configuration_comparison.py` `_variants` + `_resolve_options` give the legacy-BPNV construction; we also need our "embedded legacy via `pyrogen=None`" variant. Replace the "historical knobs" set with canonical knobs (including k_solid 0.271); record time of max Mach (my runs show the window effect).
- New information from this review to put in T8 expectations: canonical 806 cells at CFL 0.3 gives Mach 17.5 and floor activations; T8 should run the finest grids at CFL 0.15 as well, since fill Mach and min P at CFL 0.3 are not time-converged.
- Caveat: the 1600-cell headline numbers in their summary (finest-pair 0.151%, extrapolated 12.4478 MPa / Mach 0.4343) have NO committed data in their repo (only the two 806-cell JSON files); they cannot be checked.

## 3. Port plan (ordered)

Constraint: T9/T10 (in flight) edit `simulation.py`; anything touching `_run_time_loop`/`solver.py` tuples should follow their merge. Gate for all solver-touching items: B0 (rtol 1e-12, exact steps) + `trajectory_ab.py` byte A/B (T3 tooling) + full pytest.

| order | task | files | changes default results? | size | gate | expected conflicts |
|---|---|---|---|---|---|---|
| P0 | **Slag unit fix** (`slagCoeff/1e6` in `convert_nozzle`; fix docstrings lines 32/1145, print at ~1666; `nozzle.py:97` comment) + `test_slag_unit_conversion` | `srm_1d/openmotor_adapter.py`, `srm_1d/nozzle.py` (comment), `tests/test_adapter.py` | No for committed motors (slag 0); YES (correct) for any user .ric with slag > 0 - note in DEVNOTES API log | XS (~10 lines) | new test; B0 unchanged; `tests/test_adapter.py`, `test_ric_igniter.py` | none |
| P1 | **Hardening bundle**: validate `burn_update_interval >= 1`; `D_hyd <= 0 -> (r0, 0)` in `burn_rate_cell`; optional `pressure_floor_pa` kwarg with validation (default 1e3 identical) | `simulation.py`, `burn_rate.py`, `solver.py` | No at defaults | S | B0 bit-identical; unit tests (copy theirs: validation, floor 750 deficit 250, `D_hyd=0`); clear Numba cache | `solver.py` PISO signature (T6 tuple); do after T9/T10 |
| P2 | **Limit diagnostics** adapted to our T3/T6 structure: per-cell T-floor/T-ceiling/pressure-floor counts, duration, first/last time, |energy| (T limits), per-cell max deficit; optional port-Mach-cap faces. Add summary view `numerical_limits` + a `tools/ignition_diagnostics` helper | `solver.py`, `simulation.py`, `tools/ignition_diagnostics.py`, `tests/test_pressure_floor.py` or new | No (read-only) | M | B0; byte-identical trajectory A/B at 106/406/806 cells; their single-cell PISO tests ported; cost < 10% steps/s (T6 measured -1.5% for 11 arrays) | T9/T10 `simulation.py` edits; PISO return tuple (currently 21) |
| P3 | **Refresh-delay statistic**: post-hoc `ignition_refresh_delays(result)` (uses `ignition_time_by_cell`, `time`, `burn_refresh_interval_steps`) in `tools/ignition_diagnostics.py`; add the numbers to the T7 probe / worklog correction | `tools/ignition_diagnostics.py`, probe, `T5.md/T7.md` note | No | XS | test on a small run; reproduce 29.8 us-scale values on Chunc | none |
| P4 | **T8 probes adapted** from `run_spatial_extension.py` + `run_fill_configuration_comparison.py` (canonical config, both Mach windows, finest grids also at CFL 0.15) | new files in `docs/solver_verification/probes/` | No | S-M | probe output reproduces B0 digits at 100 cells | none |
| P5 | **Thread-F floor discriminator**: run 806-cell canonical at CFL 0.3 with the F test #1 probe (T iterated in the corrector); count floor events with P2/T3 counters | probe only (+ `solver.py` experiment on a branch) | No (experiment) | M | events at 806 cells: 193 baseline vs F-probe | depends on thread F work |
| P6 (optional) | B0 manifest provenance (git rev, versions, `.ric` sha256) | `probes/baseline_gate.py` | No | XS | B0 output only | none |
| hold | in-kernel floor event trace (`55b19cb`), registry/locks, limiting-case script | - | - | - | - | - |

## 4. Bugs and risks found

In our code:
- **Slag coefficient unit mismatch** (`openmotor_adapter.py:1151` on `openmotor-frontend`): 1e6x for nonzero slag (P0). Neither tests nor committed motors exercise it.
- `burn_update_interval <= 0` unvalidated; `burn_rate_cell` has no `D_hyd <= 0` guard (unreachable at runtime; hardening).
- T5/T7 worklogs say Codex's 29.84 us "not reproduced": it is the per-ignition refresh delay statistic, reproduced here (2.98e-5 s). Amend the notes.
- T7's `fill_mach` ("pre-peak") and Codex's (3 ms window) must not be compared without the window; add both to T8.

In their code:
- Chunc baseline config is stale/non-canonical (`baseline_configs.json:39-46`); all downstream headline numbers (12.4478 MPa, Mach 0.4343, 12.1574 -> 8.1357 MPa, floor probes) belong to it, so none can be compared with our B0 (11.4728 MPa).
- Fill Mach window of 3 ms hides the excursion (3.69 at 3.69 ms with library BPNV; 3.21 at 6.9 ms with MTV); their CFL-probe "fill Mach stable to 0.003%" is blind to it. The legacy-BPNV attribution is thus overstated (1.2-1.3x, not 2.2x, on the full window).
- Their 1600-cell convergence headline has no committed data.
- Energy ledger blind to the left-neighbour radiation sink (emissivity 0 only).
- `normal_rate_preserved` finite-crossflow check is tautological.
- `np.trapz` used throughout the scripts (deprecated in newer NumPy; fine on py3.10.5 here).
- `limit_first_time_s[...]` indexed whenever only `limit_activation_counts` is passed (Optional-None under numba).
- Registry checksums already mismatch our `.ric`; workspace lock pins pre-provider openMotor.

## 5. Open items
- Their pytest suite was not run; test quality judged by reading.
- 50/200-cell fill-Mach points of the factorial were not reproduced (100-cell points match to 3-4 digits).
- The event-trace fields of `55b19cb` were not verified against solver algebra.
- Whether the temporal floor events are due to thread F or to an area-jump/nozzle linearization is undecided (P5).
- Whether a floor-free finer-CFL baseline should replace CFL 0.3 for the T8 fine grids is the user's call.
- Cause of the remaining spatial growth of fill Mach (3.4 at 100 cells vs 6.0 at 806 cells at CFL 0.15) is unexplained.
- Raw run outputs: stdout only (not saved); scripts in `C:/tmp/sv_cx/`.
