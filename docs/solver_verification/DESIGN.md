# Solver verification — DESIGN

**Goal:** replicate, in our tree, the verification and bug-fix work a
friend's Codex review reported (2026-10-05), prioritised by the user:
**solver timestep/refresh sensitivity, energy audit, and bug-fixing
tests first.** We do not have Codex's diff, only its summary (numbers
below), so each item is re-derived here against our code. If the friend's
fork commits become available, they are ported/cross-checked separately.

This package serves the ignition-spike investigation (the foremost
development target — `docs/v0_7_4/SPIKE_REOPEN_BRIEF.md`): every item either
removes a numerical/bookkeeping confounder from the spike metrics or
corrects a claim in `IGNITION_SPIKE_REOPENED.md`.

Execution: `TASKS.md` (agent-executable, one worktree per task).
Per-task records: `worklog/T<n>.md`. Shared gate: `probes/baseline_gate.py`.

---

## 1. Codex's reported findings (as given; their configs unknown)

| Area | Reported |
|---|---|
| Refresh cadence | default 21-step burn refresh delays up to **29.84 µs**; burn every step → peak P **+0.073%** but transient **minimum** P **+10.18%**; geometry cadence → peak **+0.008%** |
| Erosive bisection | fixed `20*r0` upper bracket fails in small-port/high-Re states; expanding bracket closed a 150-point envelope to **5.73e-8 m/s**; Chunc baseline unchanged |
| Tab lookup | pressure-table boundary behaviour pinned with tests vs the openMotor fork |
| Limiting cases | zero crossflow → **7.5049 mm/s**, zero erosive; crossflow adds **3.4996 mm/s**, normal part unchanged; erosion off → peak **12.1574 → 8.1357 MPa**, impulse **+0.91%**, burn duration **+4.35%**; a 0.6–1.8 s "steady" window **failed** a stationarity test |
| Energy | Cp·T sensible-enthalpy balance, internally tight (source-ledger **1.30e-8 W**, nozzle **2.05e-8 W**, PISO residual **1.79e-5 W**) but NOT total energy (omits KE, formation, solid storage, pressure work, nozzle KE); between-step storage jumps from mixture/geometry refreshes: max **2.63 J**, Σ\|·\| **13.30 kJ** |
| Fill Mach (Chunc, 3 grids) | MTV **0.448/0.498/0.479**; current BPNV **1.987/2.200/2.311**; legacy BPNV **3.271/4.821/11.722** (≈ REOPENED §6's 3.36/4.81/12.04) → the §6 grid-divergent supersonic fill is mainly the pre-`bfc2f3f` BPNV gas |
| Grid convergence | to 1600 cells: finest-pair Δ peak **0.151%**, t_peak **1.373%**, fill Mach **1.389%**; ~first order; extrapolated **12.4478 MPa / Mach 0.4343**; min P unconverged — **1 kPa pressure floor activates at ≥400 cells** |

## 2. Verified against our code (2026-10-05)

- **Bisection clamp — CONFIRMED.** `burn_rate.py` `burn_rate_cell` code-0
  path: `r_hi = r0 * 20.0` with no check of `F(r_hi)`; if the root is above
  20·r0 the bisection silently returns ≈20·r0. Because the transpiration
  factor `β/(e^β−1) ≤ 1` for β ≥ 0, `h(r) ≤ h0`, so
  **`r_hi = r0 + ebf·h0` is a guaranteed bracket** (F ≥ 0 there) — an
  analytic fix, no expansion loop needed.
- **Refresh cadence — CONFIRMED.** One `burn_update_interval =
  max(10, N//5)` (`simulation.py` ~2454) gates BOTH `update_cell_geometry`
  (~1645) and `compute_burn_rates` (~1666). Regression itself advances
  every step.
- **Pressure floor — CONFIRMED unaccounted.** `solver.py`
  `_piso_step_with_energy_diagnostics` STEP 4 (~687): `P_new =
  max(P_new, 1e3)` then `rho_new = P/(R·T)` — mass and energy enter outside
  every ledger, uncounted.
- **EOS reconciliation (new observation).** The energy ledger closes on
  `new_mass = rho_new_1·A·dx` (PISO intermediate density, ~627) but the
  state carried to the next step is `rho_new = P_new/(R·T_new)` (~692). The
  between-step storage jump therefore has ≥4 candidate sources: EOS
  reconciliation, the pressure floor, the Cp refresh
  (`_advect_species` → `_refresh_mixture_arrays`, ~2029-2045), and the
  geometry refresh (`A_port` change). Codex attributed it to the last two
  only — T6 decomposes all four.
- **No ρ_gas debit on the regression mass source** — the newly opened port
  volume is not debited, consistent with a geometry-refresh mass/energy jump.
- **Tab lookup — already matches.** `select_tab_idx` mirrors openMotor's
  `Propellant.getCombustionProperties` (strict containment, then
  closest-boundary). Only regression tests needed.
- **NEW BUG (not in Codex's list) — igniter material.**
  `openmotor_adapter._pyrogen_to_block` / `_block_to_pyrogen` omit
  `gas_mass_fraction` (added to `Pyrogen` by `bfc2f3f`), and EVERY committed
  `.ric` — plus fresh GUI saves (e.g. untracked `L1843mod.ric`, 2026-08-15) —
  embeds pre-`bfc2f3f` BPNV (M 0.030, γ 1.25, T 2800; gas fraction → 1.0).
  `run_from_ric(pyrogen=None)` and the GUI plugin treat that block as
  authoritative → **~2.1× igniter gas on every GUI run.** Validation probes
  pass `pyrogen='bpnv'` and are unaffected.

## 3. Principles (apply to every task)

1. **Defaults stay bit-identical unless a task says otherwise.** Gate B0
   (`probes/baseline_gate.py`) compares Chunc/BALLSstick/Hasegawa A at full
   float precision + exact step counts. A B0 failure is a finding to
   REPORT, never a tolerance to loosen.
2. **Bookkeeping must not perturb the trajectory.** Diagnostics read state;
   they never write it.
3. **No tuning, no smoothing.** Fixes are correctness fixes with a stated
   physical/numerical reason (`feedback_no_unfounded_smoothing`).
4. **`motors/*.ric` are openMotor-owned** — never edit them. Temp copies
   outside `motors/` are fine for measurements.
5. **Numbers, not adjectives.** Every gate and measurement is reported with
   actual values in the task's worklog.

## 4. Task map

| Wave | Task | Scope | Touches |
|---|---|---|---|
| 1 | T1 | igniter block round-trip + stale-material warning | `openmotor_adapter.py`, `tests/test_ric_igniter.py` |
| 1 | T2 | erosive bisection bracket fix | `burn_rate.py`, `tests/test_burn_rate.py` |
| 1 | T3 | pressure-floor accounting | `solver.py`, `simulation.py`, `tools/ignition_diagnostics.py`, new test |
| 1 | T4 | verification tests (tab parity, limiting cases, mode invariants) | new `tests/test_verification_limits.py` only |
| 2 | T5 | split burn/geometry refresh + refresh-gap diagnostics | `simulation.py`, `openmotor_adapter.py` (passthrough), tests |
| 2 | T6 | energy audit (ledgers, jump decomposition, omitted terms) | `solver.py`, `simulation.py`, `tools/ignition_diagnostics.py`, tests |
| 3 | T7 | timestep/cadence sensitivity sweep (measurement) | new probe only |
| 3 | T8 | grid convergence to 1600 cells + fill-Mach attribution + REOPENED §12 draft | new probe + doc |

Wave 1 tasks touch disjoint files and run in parallel. Wave 2 starts after
wave 1 is reviewed and merged (T5/T6 both build on T3). Wave 3 is
measurement-only and needs T3+T5 (counters, split cadence).

## 5. Decisions reserved for the user (tasks measure, never decide)

- How stale `.ric` igniter blocks get refreshed (GUI re-select vs a
  programmatic re-seed) — T1 adds the warning + impact numbers only.
- Whether to change the default refresh cadence — only if T7 shows a
  validation metric moving > 1%.
- Whether to move to a total-energy formulation — only if T6 shows the
  omitted terms matter where the spike lives.
- Whether to change the pressure floor — T3/T8 only measure activations.
