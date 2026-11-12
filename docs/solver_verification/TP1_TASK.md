# TP1: thrust and performance from the nozzle-entrance stagnation pressure

**User decision (2026-11-11):** "Run thrust-pressure fix now", in parallel with
the end of NB1.

## The defect (found by NB1, `worklog/NB1.md` step 7)

`nozzle.compute_thrust_history` and `nozzle.compute_motor_performance` take
`P_head` (cell 0, at the head end) as the chamber pressure of the
thrust-coefficient formula: F = CF · P_c · A_t. Two things differ from the
pressure that actually drives the throat:

- **Static versus stagnation.** At the head end the gas is nearly at rest, so
  static pressure equals stagnation pressure there. That part is fine.
- **Port losses.** The flow loses stagnation pressure along the port
  (friction, and mass addition, which accelerates the gas). The pressure that
  feeds the nozzle is the stagnation pressure at the port exit,
  P0_noz = P_exit · g(M_e), where g(M) = (1 + (γ-1)/2 M²)^(γ/(γ-1)).
  NB1's boundary state already computes M_e.

At the plateau, P_head/P0_noz is 1.030 (Chunc), 1.034 (BALLSstick) and 1.011
(Hasegawa A). Thrust and total impulse are therefore over-predicted by about
that much. During the startup the gap is larger.

`P_head` stays the right quantity to compare with a head-end pressure
transducer. Keep it in the results and plots as the "chamber pressure" trace.
Only the thrust and performance path changes.

## What to do

1. **Record the nozzle-entrance stagnation state per history row**:
   `P0_nozzle` and `T0_nozzle`.
   - Take them from the SAME boundary state the solver used:
     `_nozzle_boundary_flow` with this step's `M_choke`, in either
     `nozzle_bc_stagnation` mode (`'static'` gives P0 = P_exit).
   - Use the same decimation and capacity as the other history arrays (MEM1:
     bounded, no per-step lists).
   - Reverse inflow or no flow: define a sensible value (e.g. P_exit) and
     document it.
2. **Use `P0_nozzle` in the thrust path.**
   - Update `compute_thrust_history`, `compute_motor_performance` and every
     consumer: thrust, CF, Isp, total impulse, average thrust, designation,
     CSV export, `channels` / `build_channels`, and the openMotor plugin's
     thrust and performance channels.
   - The GUI pressure plot should still show the head-end pressure. Find out
     which channel openMotor plots as "pressure" and keep it on `P_head`.
   - If a result lacks `P0_nozzle` (an old pickled result), fall back to
     `P_head` with a warning.
3. **Check for related inconsistencies in the same path** (report, and fix
   only if clearly the same defect):
   - **Throat erosion.** `compute_thrust_history` re-integrates the throat
     diameter from `P_head`. Does the solver erode the throat itself? If it
     does, the post-processor must use the solver's recorded throat history,
     not a second integration with a different pressure.
   - **Gamma.** The CF formula uses a scalar γ. During startup the gas at the
     nozzle is a mixture (air, igniter gas, propellant gas). Report how much
     the mixture γ at the exit cell differs in the first 50 ms. Do not change
     it unless it is clearly wrong.
   - **C\* consistency.** C\*_sim = P0_nozzle · A_t / ṁ_nozzle should match
     the propellant's C\* (or the mixture's) at the plateau to within the
     nozzle-efficiency definitions. With `P_head` it is inflated by
     P_head/P0_noz. Show both.
4. **Validation against measured thrust.** `static_fire_data/pathfin54.csv`
   has a measured thrust curve (`time(s), force(N), pressure(Pa)`). The motor
   is `motors/2025.12.25 54-2800 ST2.0.ric`; use the committed version.
   - Run it at the defaults (auto charge).
   - Compare measured with simulated, before and after the fix: total
     impulse, average thrust, peak thrust, burn time, and thrust-over-pressure
     (the effective CF · A_t) at the plateau. The measured pressure transducer
     location is unknown; say what the comparison can and cannot tell.
   - Also tabulate before and after for Chunc, BALLSstick and Hasegawa A:
     total impulse, average thrust, peak thrust, Isp, and C\*_sim against
     C\*_prop.

## Gates

- Unit tests:
  - static mode reproduces the old thrust exactly when P0 = P_exit;
  - for a no-flow cell (u = 0, no port loss) P0_nozzle = P_head and the thrust
    matches the old formula;
  - `P0_nozzle` ≤ `P_head` on the canonical plateaus;
  - the fallback path works;
  - the history arrays stay bounded.
- B0 (`docs/solver_verification/probes/baseline_gate.py`) is pressure-only and
  must stay UNCHANGED. Run it; if it moves, find out why.
- Run the full pytest suite under the guard with a 60-min timeout (it takes
  about 20 min). Set
  `SRM1D_OPENMOTOR_PATH="C:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/openMotor/openMotor"`.
- Log any signature change of an `@njit` function or of a result or
  performance key in the DEVNOTES "API Breaking Changes Log".

## Branching

- Your branch `tp1/thrust-p0` starts from `nb1/stagnation-bc` at 0f70601,
  which already contains G2, AW1 and NB1.
- NB1 is still finishing on its own branch (B0 re-record, tables, suite) and
  will merge into `openmotor-frontend` soon. When the coordinator tells you it
  has landed, merge `openmotor-frontend` into your branch (merge, not rebase)
  and re-run the gates.

## Rules

- **Compute:** every job through
  `python -m srm_1d.tools.guarded_run --mem-gb 8 --timeout-min 60 -- <cmd>`.
  At most 4 simulation processes. Warm the Numba cache with one small run
  first. Probes are memory-bounded with max-step stops.
- Clear `srm_1d/__pycache__/` after every `@njit` edit.
- No tuned constants. Do not hand-edit `motors/*.ric`. Do not touch the main
  checkout or the NB1 worktree (`C:/tmp/sv/nb1`).
- Commit on `tp1/thrust-p0` in `C:/tmp/sv/tp1`. Do not merge into
  `openmotor-frontend` and do not push.
- Write `docs/solver_verification/worklog/TP1.md`.
- Report:
  - the before/after tables;
  - the pathfin54 measured-versus-simulated comparison;
  - the throat-erosion, γ and C\* findings;
  - the gates;
  - your commits.
