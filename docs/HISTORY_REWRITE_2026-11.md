# History rewrite of 2026-11-11

Before the first push since v0.8.1, the 199 local commits made after
`db0a980` (the remote tip, 2026-07-16) were rewritten to keep four reference
PDFs out of the public history (copyright; see `docs/references/README.md`).
Only those PDFs were removed. Every tree is otherwise identical, and
messages, authors and dates are unchanged. Commit `3f11b01` added only the
PDFs, so it is now empty.

Worklogs, reports and commit messages written before the push cite the OLD
hashes. This table maps them to the published ones (topological order; short
hashes; pairs checked by subject and author date):

| old | new | date | subject |
|---|---|---|---|
| e595d01 | e595d01 | 2026-10-05 | feat(burn_rate): erosive-model registry + opt-in Mukunda-Paul ceiling |
| a2e3e8e | a2e3e8e | 2026-10-05 | docs(ignition): cross-motor MP measurement, lit fold-in, registry plan + worklog |
| 3f11b01 | d98fb53 | 2026-10-05 | docs(refs): add Ma 2022, Sivan 2019, Luke 1996, Kuo & Summerfield |
| 73df8c6 | 3aa728d | 2026-10-05 | docs(claude): ignition spike is OPEN and the foremost target; roadmap reorder |
| 0896bb9 | 27e5bb0 | 2026-10-05 | docs(verification): solver-verification package (Codex-review replication) |
| 31c779b | d5de919 | 2026-10-05 | docs(verification): worktree location + mandatory SRM1D_OPENMOTOR_PATH |
| b36063d | 75a186d | 2026-10-05 | fix(adapter): lossless igniter block round-trip + stale-material warning |
| 9b2623b | 1f23f18 | 2026-10-05 | docs(verification): T1 worklog (igniter block round-trip + stale warning) |
| fd36097 | fc71357 | 2026-10-05 | test(verification): T4 tab parity vs openMotor, zero-crossflow, normal-invariance, mode invariants |
| ec34109 | 4deaa61 | 2026-10-05 | docs(verification): T4 worklog |
| f818794 | 6f465d1 | 2026-10-05 | fix(burn_rate): guaranteed upper bracket for the Ma erosive bisection |
| b4b9c4b | c59f678 | 2026-10-05 | docs(verification): T2 worklog (bisection bracket; B0 ULP finding) |
| e0c5b5f | dbba01e | 2026-10-05 | review(verification): r0==0 guard on the T2 bracket widening; B0 tolerance policy |
| 936bbdf | f84a4d8 | 2026-10-05 | feat(solver): account for the 1 kPa pressure-floor clamp (accounting only) |
| 7b7c348 | 85b2eee | 2026-10-05 | docs(verification): T3 worklog + A/B and floor-measurement scripts |
| 9b3ef00 | 16c32ed | 2026-10-05 | docs(verification): fold wave-1 lessons into T6/T7/T8 specs |
| 5b91583 | 6fc7439 | 2026-10-05 | docs(igniter): plugin-provided igniter design + agent tasks |
| 2ddb372 | efe44ca | 2026-10-05 | docs(igniter): I4 - stale GUI solver-config defaults (pre-v0.7.5 knobs) |
| ef0b5eb | 7c7262b | 2026-10-05 | docs(verification): T8 lead - our fill Mach (3.94 library BPNV) vs Codex's ~2.2 |
| 7ada321 | 8013d09 | 2026-10-05 | feat(simulation): split burn/geometry refresh cadence + refresh-gap diagnostics (T5) |
| 56b2412 | c9671b1 | 2026-10-05 | feat(igniter): srm_1d provides the igniter model (hooks, library, superseded-snapshot upgrade) |
| 6145158 | 71c89e7 | 2026-10-05 | test(igniter): offscreen GUI smoke for the provider path; plugin test follows the raw-block contract |
| 5b7db7d | f8d8cc5 | 2026-10-05 | feat(igniter): refresh committed motors/*.ric through openMotor (D7) + refresh tool |
| 54878fc | f60614f | 2026-10-05 | docs(igniter): worklog for I1/I2/I3 |
| 3128d7c | 41bbd04 | 2026-10-05 | feat(energy-audit): ledger completeness, storage-jump decomposition, sensible-enthalpy audit (T6) |
| a3be70d | 4ed126e | 2026-10-05 | docs(verification): T6 worklog |
| ad535c0 | b777c8b | 2026-10-05 | docs: T9/T10 specs (radiation sink, (rho_p - rho_g) mass source) + spike thread F |
| 1da1691 | 342804f | 2026-10-05 | feat(igniter): refresh example_bates/zerox/zerox_LHS (permit config keys at MotorConfig default) |
| 608076c | 7afde26 | 2026-10-05 | feat(igniter): derive GUI solver-config defaults from run_simulation; upgrade stale saved ones (I4) |
| 3212722 | 249e643 | 2026-10-05 | docs/test(igniter): smoke rebuilds the legacy state; worklog I3b/I4 |
| 504e1d5 | e3ba465 | 2026-10-05 | test(verification): T7 cadence/CFL sensitivity probe + worklog |
| e4bdab8 | cf788e9 | 2026-10-05 | docs(verification): R1 - review spec for the friend's Codex commits |
| d5f0daf | aca46a9 | 2026-10-05 | docs(ignition): thread F gets T7's time-convergence evidence |
| 1bc4304 | f575a67 | 2026-10-05 | docs(verification): R1 Codex review report + wave-5 port specs |
| dea31a9 | ede61e0 | 2026-10-05 | fix(simulation): debit the left-neighbour radiation sink from the gas (T9) |
| df35e80 | 813c401 | 2026-10-05 | fix(simulation): debit the gas filling newly opened port volume from the regression source (T10) |
| 8061e8a | 4eda31b | 2026-10-06 | fix(adapter): convert slagCoeff (m*Pa)/s -> (m*MPa)/s (P0) |
| c7e4a52 | e132ee8 | 2026-10-06 | feat(tools): ignition_refresh_delays statistic + T5/T7 correction notes (P3) |
| 7fdad35 | 8606dc5 | 2026-10-06 | chore(b0): print git/version/.ric sha256 provenance in the gate (P6) |
| 5e96c4d | 5e18a45 | 2026-10-06 | docs(verification): PA worklog |
| e74cec2 | 69f75c5 | 2026-10-06 | review(verification): B0 float tolerance 1e-12 -> 1e-9 (post-T10 compile-state noise) |
| 7cbca1d | aa770d2 | 2026-10-06 | fix(simulation): size history from the hottest sound speed; history_capacity override |
| 182e9d0 | 6dfb126 | 2026-10-11 | feat(probes): T8 grid-convergence + fill-Mach attribution probe and raw results |
| be60fea | 5db9064 | 2026-10-11 | docs(verification): T8 worklog + REOPENED section 12 DRAFT (grid convergence, fill-Mach attribution) |
| a611b67 | 602adc8 | 2026-10-11 | docs(verification): P5 spec - energy-coupled PISO iteration experiment (thread F) |
| 52c9a97 | 4de61cc | 2026-10-11 | feat(solver): P1 hardening - pressure_floor_pa kwarg, burn_update_interval validation, D_hyd<=0 guard |
| db286ac | 4295a58 | 2026-10-11 | feat(solver): P2 per-cell numerical-limit diagnostics (T floor/ceiling, pressure floor, Mach cap) |
| 4c65763 | 7011b82 | 2026-10-11 | docs(verification): T11 spec - adopt P5 PISO corrections as defaults |
| 2a54323 | 4666a4c | 2026-10-11 | feat(solver): opt-in energy-coupled PISO iteration (P5) + corr2-storage ablation |
| 45a2915 | ad173fa | 2026-10-11 | docs(verification): P5 worklog, REOPENED section 13 DRAFT, probe results |
| d9796f7 | c983ddd | 2026-10-11 | Merge sv/p5 into sv/t11: energy-coupled PISO iteration + corr2 storage (opt-in), keeping PB limit diagnostics |
| 511a45f | d3557e7 | 2026-10-11 | feat(solver): T11 - make piso_corr2_storage=True and piso_energy_coupling_iters=1 the defaults |
| d54343c | 2dc259f | 2026-10-11 | docs(verification): T11 worklog, DEVNOTES behaviour log + calibration note, REOPENED 13.5 adoption note |
| 9b6d4b3 | 25fe98b | 2026-10-11 | review(verification): re-record B0 for the T11 defaults; t_peak shift ruled expected |
| 8d414fc | 37f58ae | 2026-10-11 | docs(claude): roadmap records the 2026-10-11 PISO default change |
| daf367c | b745061 | 2026-10-11 | review(verification): B0 tolerances re-scaled to the post-T11 noise floor |
| e882134 | 31c51d5 | 2026-10-11 | docs(pyrogen): IG1 diagnosis spec - realistic igniter charges over-spike |
| 9469538 | 24809d2 | 2026-10-11 | docs(pyrogen): IG1 diagnosis report, paywalled-reference list, probes + results |
| cbf52dd | aca8e5d | 2026-10-11 | docs(pyrogen): resolve R6 (Keller 1966) from the repo PDF; geometry caveat |
| e17106b | 55f242c | 2026-10-11 | docs(pyrogen): IG2 spec - option A (ignition heat transfer) + the DeMar question |
| 93058f2 | 97d936b | 2026-10-11 | feat(ignition): opt-in ignition_heat_transfer modes (peretz, keller) + IG2 probes |
| 3234470 | 787e414 | 2026-10-11 | docs(pyrogen): IG2 report, results, paywalled-reference additions; docstring fix |
| 895c723 | c1f7673 | 2026-10-11 | docs(refs): rename Peretz AMS-1100 PDF to peretz1974.pdf; fix citations; record fetched references |
| 2c08c46 | 5740259 | 2026-10-11 | docs(pyrogen): OB1 spec - option B, literature-grounded pyrogen burn model |
| 966261e | b1bea60 | 2026-10-11 | docs(pyrogen): R9 source search - where the Rocketdyne duct data can be found |
| a85493b | 93d033e | 2026-10-11 | docs(pyrogen): IG3 spec - Peretz default, Ap-reference test, drop DeMar, DeMar sizing |
| 9790459 | 0a75606 | 2026-10-11 | IG3 stage 1 WIP: peretz default ignition h, peretz_ap_reference (local|head|aft), plugin choice field, tests,  |
| 1c00408 | f0395f6 | 2026-10-11 | IG3 stage 1: B0 re-record (Peretz default), DEVNOTES, IG3 measurement probes + results |
| bb2f5ab | 1e684d2 | 2026-10-11 | IG3 stage 1: pin two tests to gnielinski ignition h, record suite run, R10 constant check |
| de69787 | da44046 | 2026-10-12 | test(B0): compare t_peak with the step tolerance (argmax over steps) |
| 5a2b28f | e0decb9 | 2026-10-12 | docs(CLAUDE.md): note the Peretz ignition-h default and the pyrogen-realism work in progress |
| ebd593e | c519329 | 2026-10-11 | OB1 WIP: opt-in pyrogen rate plateau + basket area law, probes, first results |
| 4801db8 | 4a86dc5 | 2026-10-11 | OB1: tests, DEVNOTES, bpnv.yaml provenance correction, paywalled refs, probes + results |
| 2a51817 | 143da65 | 2026-10-11 | OB1: report (gates partially open), P_plateau YAML/block round-trip fix + test, B0 evidence |
| 6981fc0 | 20cfc27 | 2026-10-12 | OB1: rebase cleanup (PAYWALLED R9 line) |
| 6da2c92 | cf63ce9 | 2026-10-12 | OB1: exact cylinder form function (plenum 0d + basket default), I-A pellet geometry, SUPERSEDED snapshots (WIP |
| 82b1b64 | 0f223eb | 2026-10-12 | OB1: tests for exact cylinder form function, I-A geometry, superseded snapshots |
| 020212c | b5cb01f | 2026-10-12 | OB1: B0 re-recorded (warm-state values; previous kept in comments), clean+warm evidence |
| 3cedd4a | ebccd0a | 2026-10-12 | OB1: post-rebase measurements (ob1r_*), report rewritten for the new defaults, DEVNOTES |
| 3f34259 | 8d45318 | 2026-10-12 | docs(pyrogen): SM1 spec - solid-heating model comparison; reviewer dead-zone probes |
| 358d83a | 9dd31e3 | 2026-10-12 | docs(pyrogen): Carlson & Seader 1967 obtained - test ranges, Ap definition, head-end correlation form |
| 28d4fa9 | c0b6f69 | 2026-10-12 | IG3 stage 2 WIP: DeMar flux off by default, DeMar sizing, auto-upsize fallback with physical early stop |
| b0eaee3 | 3c69f0a | 2026-10-12 | IG3 stage 2 WIP: grain-gas early-stop criterion, stage-2 tests |
| dc248d1 | 4b8b8bc | 2026-10-12 | IG3 stage 2 WIP: test fixes (impetus in sibling pyrogen, robust wrapper) |
| 39a5eb5 | 2062000 | 2026-10-12 | IG3 stage 2 WIP: pin DeMar-driven tests to heat_delivery_mode demar; history-full is not a no-ignition retry |
| 2e36dc1 | 450db51 | 2026-10-12 | IG3 stage 2 WIP: success = >=99% grain cells lit (ignition_min_lit_fraction), partial-ignition tests |
| 827e5b6 | 2c90e29 | 2026-10-12 | IG3 stage 2 WIP: DEVNOTES entry |
| 5b61dbf | 1b8c2cd | 2026-10-12 | IG3 stage 2: provisional B0 / per-motor results (pre-OB1, not recorded) |
| 101927f | e4c5a45 | 2026-10-12 | IG3 stage 2: rebase onto OB1 (3cedd4a); resize keeps basket_area_law/bed_flame_speed; post-OB1 superseded snap |
| 3028e47 | 910ca19 | 2026-10-12 | IG3 stage 2: B0 re-recorded (DeMar-sized auto charge + fallback, on top of OB1), per-motor results |
| 5e1070b | 0a88eea | 2026-10-12 | IG3 stage 2: DEVNOTES per-motor table |
| a71b887 | b3e0db5 | 2026-10-12 | IG3 stage 2: pin explicit charge in refresh-cadence test; suite run; DEVNOTES |
| 74260a6 | 535f3e3 | 2026-10-12 | IG3 stage 2: success = self-sustaining + warn on never-lit cells (replaces the 99%-lit test) |
| 0f84314 | 9ebb636 | 2026-10-12 | IG3 stage 2: B0 re-recorded for the self-sustaining success test |
| ad5c169 | 805fa79 | 2026-10-12 | IG3 stage 2: DEVNOTES for the self-sustaining success test |
| b597d94 | 23eb56e | 2026-10-12 | IG3 stage 2: final full-suite log (801 passed) |
| 146d3ab | 8a14cfb | 2026-10-12 | SM1 Part 1+2: solid-heating models (goodman|peretz|conduction), exact-solution verification, tests |
| 076d25f | 99e8fb3 | 2026-10-12 | SM1 Part 3+4: sim comparison (dead zones, lighting time, thresholds, cost), SM1_REPORT, HEAD_END_NOTES, DEVNOT |
| b98d90c | 091d712 | 2026-10-13 | SM1: post-rebase gate logs (B0 PASS, 821 passed) |
| 3bd3eba | 3606095 | 2026-10-13 | fix(plugin): GUI snapshot interval floor 5 ms -> 1 ms (resolve the fill/ignition transient) |
| 0eada1e | 9072e0b | 2026-10-13 | docs(pyrogen): SM2 (controlled Peretz-ODE vs conduction comparison) and HE1 (head-end jet h) specs |
| 324956f | d43406b | 2026-10-13 | fix(plugin): 1 ms GUI frames up to Max Simulation Time 20 s |
| c5f4af9 | 1278d70 | 2026-10-13 | docs(solver): EN1 spec - conservative total-energy formulation (opt-in), exact-solution verification |
| 1d17499 | b233539 | 2026-10-13 | fix(adapter): explicit igniter-sizing kwargs now override the .ric igniter block |
| a22b0e1 | 85fd971 | 2026-10-13 | SM2 WIP: ignition-event diagnostics + open-loop forcing recorder, replay tool, case/sweep runners, group A |
| 13e5606 | 43c8158 | 2026-10-13 | SM2 WIP: groups A H D C K G X results |
| 369a6b0 | 2cf675d | 2026-10-13 | SM2: controlled Peretz-ODE vs conduction comparison (open-loop replay, calibration, h/DeMar/charge/grid factor |
| cca685e | 4778128 | 2026-10-21 | feat(tools): compute safeguards - guarded_run (memory cap + timeout + kill-on-close), 2-worker sweep default,  |
| ae86f9b | 1c117c0 | 2026-10-21 | docs(DEVNOTES): API log - sweep default workers 2, guarded_run |
| e768d18 | 380b085 | 2026-10-21 | docs(pyrogen): GS1 spec - grid and contact-front convergence of the startup (picks the LHS grid) |
| 6ca930a | cf48233 | 2026-10-13 | IG4: self-sustaining = grain gas-generation rate exceeds the igniter's peak after the igniter is spent; plugin |
| 8db7a11 | 6cf8014 | 2026-10-13 | IG4: tests (pathfin54 hang-fire, canonical success, plugin alert) |
| 80da565 | 048a600 | 2026-10-21 | IG4: B0 and canonical per-motor result files |
| fb5aa52 | f6630cd | 2026-10-21 | IG4: guarded B0/suite attempts (hit the 4 GB cap; heavy-compute request) |
| 44d34be | 2ce2a8e | 2026-10-21 | MEM1 WIP: no history floor, in-place trim, auto-grow re-run, probe |
| 8145fbb | 3d30c35 | 2026-10-21 | MEM1: tests (no floor, auto-grow re-run, in-place trim), sweep probe |
| aaf9625 | 6485f42 | 2026-10-21 | MEM1: bound the auto history at 3M rows (grow once to 6M); tests |
| ad1e0ab | 421cd69 | 2026-10-21 | MEM1: all-motor memory sweep results |
| 4a625fe | 5883855 | 2026-10-21 | MEM1: DEVNOTES entry; station-viz payload test runs 1 s instead of 6 s (5.3M-step history) |
| d592779 | 0aae4b5 | 2026-10-21 | MEM1: gate logs (B0 PASS, suite 850 passed under the 4 GB guard) |
| 130afc2 | 3463e07 | 2026-10-21 | conduction default: solid_heating_model='conduction', plugin choice field, failed result for PyrogenAutoSizeEr |
| ec05d04 | c7a67d0 | 2026-10-31 | conduction default: B0 baseline re-recorded (previous values kept in comments), per-motor probe results |
| 5cbc300 | f61f08e | 2026-10-31 | conduction default: B0 PASS log |
| 4d68573 | c75629e | 2026-10-31 | conduction default: pin the Goodman-specific diagnostics test, relax the never-lit-warning test, DEVNOTES entr |
| 9c38600 | 56a15b3 | 2026-10-31 | docs(CLAUDE.md): roadmap note - current ignition/pyrogen defaults, open branches, findings |
| cf3b5a0 | 096a400 | 2026-10-31 | GS1: row probe + driver |
| 6792727 | 12d5553 | 2026-10-31 | GS1: probes, raw results, tables |
| 7dc409e | 78a903a | 2026-10-31 | GS1: report (convergence verdicts, contact front, recommended LHS grid + cost, MUSCL recommendation) |
| 912d437 | 798d9ec | 2026-10-13 | EN1 WIP: opt-in conservative total-energy PISO step + conservative momentum predictor, exact-solution harness |
| 0a6c588 | a2c6d85 | 2026-10-13 | EN1 WIP: duct ladder, motor tables, tests, DEVNOTES |
| b24a2c9 | 1cf62e0 | 2026-10-21 | EN1: bounded duct marching (fixes the unbounded-residual hang), duct ladder re-run, Chunc rows on rebased base |
| abf8f25 | cddc846 | 2026-10-21 | EN1: BALLSstick rows re-run on the rebased base (identical); Hasegawa A rows blocked by the 4 GB cap |
| 7aaec6b | 9d77dec | 2026-10-21 | RT1: opt-in ignition driver temperature (adiabatic-wall recovery T_aw), tests, driver flags |
| 11d27ae | 3f4cddd | 2026-10-21 | RT1: ignition driver temperature grid, regime probe, worklog, DEVNOTES |
| bb81394 | 41f750a | 2026-10-31 | EN1/RT1: Hasegawa A rows on the MEM1 base (full-length late-peak grid), docs |
| 394c638 | fa17175 | 2026-10-31 | HE1: opt-in head-end igniter-jet heat transfer (Peretz III-28..30), tests, driver flag |
| 92eb63c | d550346 | 2026-10-31 | HE1: report, measurements (fixed-charge dead zones, auto-size, thresholds), Carlson-Seader figure crops |
| a6a01fe | 9a5d86c | 2026-10-31 | docs(CLAUDE.md): EN1/RT1/HE1 options now on main (opt-in) |
| 2ee2755 | cf2f696 | 2026-10-31 | docs(CLAUDE.md): the startup grid dependence is the 3-cell head gap (MU1), not resolution |
| b77d51b | d6b6697 | 2026-10-31 | docs(solver): G1 spec - physical head/aft ends replacing the 3-cell full-diameter buffers |
| fcf6a90 | eddd93b | 2026-10-31 | D1: EN1/RT1/HE1 options become the run_simulation defaults; plugin choice fields; B0 re-recorded; tests pinned |
| 6ff55f3 | 5e63357 | 2026-10-31 | docs(CLAUDE.md): D1 - total energy / conservative momentum / recovery / head-end jet are the defaults |
| ca27241 | 01eaf9e | 2026-11-10 | docs(pyrogen): Jungowski comment + Carlson reply (1969) - jet regimes, x_im basis; Re_d still undefined |
| d11601d | d9b7ba5 | 2026-11-10 | QF1: total-energy regression-fill debit carries kinetic energy; p dV audit channel |
| fdc3534 | 73c0702 | 2026-11-10 | QF1: re-record the B0 baseline (KE term in the total-energy fill debit) |
| 64d2398 | 04e2415 | 2026-10-31 | G1 WIP: no head/aft gas-buffer cells; head gap + nozzle convergent lumped into cell storage (vol_extra); cone  |
| 93a05f6 | 0d56d96 | 2026-10-31 | G1 WIP: lumped convergent default, BC-reservoir cut for cone cells, stability probe, notes draft |
| 07ea969 | b1c1cd6 | 2026-11-10 | G1 WIP: DEVNOTES entry, stability probe default jet |
| 299e201 | 5197d54 | 2026-11-10 | G1: per-motor head-gap override (run_from_ric head_gap kwarg, igniter-block head_gap_m, plugin chamber field), |
| 921b23b | a21d12a | 2026-11-10 | G1: fix tests for the new ends, B0 re-recorded (D1 values kept in comments) |
| 5738c47 | dd7f3c6 | 2026-11-10 | G1: worklog, results |
| 3770ef5 | 1978482 | 2026-11-10 | G1: final gates (B0 on QF1, suite 938 passed), worklog |
| 2d22beb | baac10f | 2026-11-10 | G1: worklog suite note |
| 331c0b3 | 583f031 | 2026-11-10 | G1: worklog suite history |
| fc7e0de | a11f05d | 2026-11-10 | docs(CLAUDE.md): G1 - physical head/aft ends, grid-independent but more violent startup |
| d66d059 | 0d5defb | 2026-11-10 | docs(fable): prompt documents for the spike audit (task 1) and the blind core audit (task 2) |
| 364d282 | 12a4d0e | 2026-11-10 | AS1 WIP: prompt-start criterion, window early stop, plugin WARNING |
| 2018d36 | 2635806 | 2026-11-10 | AS1: tests, B0 re-record (BALLSstick row), validation probe |
| 9b27ef9 | 6e826eb | 2026-11-11 | AS1: worklog, DEVNOTES entry, validation results |
| e1214e8 | be081f8 | 2026-11-11 | AS1: fix the two suite failures (flame-front window, SimpleNamespace geometry), suite log |
| adc68e1 | d5fd2e4 | 2026-11-11 | AS1: arm the pressure cutoff only after the grain has started (prompt fraction lit or window elapsed) |
| 0efd4a7 | c43045c | 2026-11-11 | AS1: arming follow-up - IG4 hang-fire tests use a charge that genuinely hang-fires (0.25 g), worklog, gate log |
| db86fda | d6177ed | 2026-11-11 | Merge AS1: prompt-start auto-charge criterion + pressure-cutoff arming |
| ea8c8f7 | 0bb47ac | 2026-11-11 | docs(CLAUDE.md): AS1 prompt-start auto charge + cutoff arming |
| 7b56144 | 3d64a65 | 2026-11-11 | docs: NB1 task spec - nozzle boundary uses the true (isentropic) stagnation state |
| 5cf2219 | cabe061 | 2026-11-11 | NB1: nozzle boundary uses the true isentropic stagnation state (default), 'static' kept for A/B |
| 1b06f95 | 30cb0fb | 2026-11-11 | NB1: B0 re-record on the isentropic default (AS1 values kept in comments), ledger and cone-reservoir probes |
| e715296 | f3be8ec | 2026-11-11 | NB1: worklog, plugin schema test for nozzle_bc_stagnation |
| 2232f8d | 4f3620b | 2026-11-11 | NB1 worklog: B0 pass after re-record |
| de419d6 | ec0d169 | 2026-11-10 | G2 step 1 (memo 2.2): cumulative-position rounding of segment lengths + per-segment burning-perimeter scale L_ |
| 32f1cee | 86c9fe2 | 2026-11-11 | G2 step 2 WIP (memo 2.1): inter-grain gaps lumped 50/50 into the neighbouring cells, no dx clamp, x_offset phy |
| ad1a9ea | 4406cf0 | 2026-11-11 | G2: spread every lumped volume over cells (lump_ratio_cap=1): fixes the fine-grid fill collapse of a big conve |
| b2c3bab | 7a05428 | 2026-11-11 | Merge openmotor-frontend (AS1) into g2/snap-lump |
| 41802e2 | 166c359 | 2026-11-11 | G2: B0 re-recorded on top of AS1 (BALLSstick row; AS1 values in comments) |
| 1cb26fc | 38011b6 | 2026-11-11 | G2: worklog, results, memo status update, trace probe --no-lump-gaps |
| 3008961 | c3d9a71 | 2026-11-11 | Merge G2: mass-preserving grain snapping + lumped inter-grain gaps + spread lumps |
| a5fb637 | a09dee4 | 2026-11-11 | docs(CLAUDE.md): G2 grain snapping + lumped gaps |
| 6eff280 | 36d250c | 2026-11-11 | AW1: prompt window default 0.5 -> 0.1 s (user decision); B0 Chunc + Hasegawa A rows re-recorded |
| b4abfdc | fe1ede6 | 2026-11-11 | AW1: gate logs (suite 947 passed, B0 PASS) |
| 53a8d52 | 7e76cab | 2026-11-11 | Merge openmotor-frontend (G2) into aw1/prompt-window |
| 706a24c | 962919b | 2026-11-11 | AW1: B0 Chunc + Hasegawa A rows re-recorded on top of G2 |
| 712080c | 8948ea6 | 2026-11-11 | AW1: merged-tree gate logs (B0 PASS, suite 974 passed) |
| 0b4a16f | 12e16c8 | 2026-11-11 | Merge AW1: prompt window default 0.5 -> 0.1 s (user decision) |
| cd1e946 | d162d7f | 2026-11-11 | docs(CLAUDE.md): AW1 prompt window 0.1 s |
| e341f31 | aabf8fc | 2026-11-11 | Merge openmotor-frontend (G2, AW1) into nb1/stagnation-bc |
| 0f70601 | ed4406c | 2026-11-11 | NB1 step 5: bc_reservoir_area_ratio is a resolution guard (docstring + DEVNOTES), behaviour unchanged |
| 08625f2 | 3fc98e0 | 2026-11-11 | docs: TP1 task spec - thrust/performance from the nozzle-entrance stagnation pressure |
| d7226bf | 85ecea2 | 2026-11-11 | TP1: thrust/CF/Isp/impulse from the nozzle-entrance stagnation pressure P0_nozzle |
| efc8870 | d0f39a1 | 2026-11-11 | TP1: worklog, results, pathfin54 comparison probe/figure |
| 697cbdf | dc2d5ad | 2026-11-11 | NB1: re-record B0 on top of G2 + AW1 (G2+AW1 values kept in comments) |
| 303416f | 7e73181 | 2026-11-11 | NB1 worklog: merged-tree tables, step 5/6 refresh, gates |
| 4e4c796 | c3e8f00 | 2026-11-11 | Merge NB1: nozzle boundary uses the true isentropic stagnation state (default) |
| cc147cb | 2c655b3 | 2026-11-11 | docs: NB1 in CLAUDE.md; Fable spike-audit prompt brought up to G2/AW1/NB1 defaults |
| 32e3742 | 1be4f7b | 2026-11-11 | Merge openmotor-frontend (NB1) into tp1/thrust-p0 |
| 007cfd7 | 602e586 | 2026-11-11 | TP1: worklog results block (suite, post-merge B0) |
| 9ac2c2c | 14653a7 | 2026-11-11 | Merge TP1: thrust/CF/Isp/impulse from the nozzle-entrance stagnation pressure |
| 64f2b82 | 312009f | 2026-11-11 | docs: TP1 in CLAUDE.md and the Fable spike-audit prompt |
| 8296515 | a257a37 | 2026-11-11 | chore(refs): stop distributing reference PDFs (copyright) |
| 47187cd | 3a83e52 | 2026-11-11 | docs: CHANGELOG.md - everything since v0.8.1 (0.9.0-dev) |
