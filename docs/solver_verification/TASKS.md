# Solver verification — TASKS (agent-executable)

Read `DESIGN.md` first (§2 has the verified facts with line numbers, §3 the
principles). Each task is self-contained: it names its files, steps, and
acceptance gates and needs no architectural decisions. **If something is
ambiguous or a step is impossible as written, STOP and report — do not
invent.** Line numbers are approximate (2026-10-05); find the code by the
names given.

## Worktree protocol (every task)

You run in an isolated **git worktree** (at `C:/tmp/sv/t<n>`, branch
`sv/t<n>`) created from local `openmotor-frontend` HEAD. Your worktree path
is given in your launch prompt; `cd` there first and stay there.

- Work ONLY inside your worktree. Never touch the main checkout
  (`C:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d`),
  never switch/modify other branches, never push. Scratch files go in
  `C:/tmp/sv_<task>/` (e.g. `C:/tmp/sv_t1/`).
- **Interpreter:** always `"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe"`
  (system Python lacks the deps). Bash (Git Bash) is available; quote paths
  with spaces.
- **Import-path trap (critical):** `srm_1d` is pip-installed in EDITABLE mode
  pointing at the MAIN checkout. A script run as `python some/dir/x.py`
  imports the MAIN checkout's code, not yours. Therefore:
  - pytest: run from your worktree root (the root `conftest.py` fixes the path).
  - any script: run from your worktree root with `PYTHONPATH="$PWD"` set.
  - **Preflight (record the output in your worklog):**
    `PYTHONPATH="$PWD" "<python>" -c "import srm_1d; print(srm_1d.__file__)"`
    must print a path INSIDE your worktree. If not, STOP.
- **openMotor (mandatory):** a worktree outside the repo tree cannot
  auto-discover the openMotor checkout, so export
  `SRM1D_OPENMOTOR_PATH="C:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/openMotor/openMotor"`
  in EVERY shell command that runs Python (pytest, scripts, gates).
- **Numba cache:** after editing ANY `@njit` function, delete your
  worktree's `srm_1d/__pycache__/` before running anything (stale compiled
  code is the #1 "my fix didn't work" trap). The first run after a clear
  recompiles for 2–4 min; that is normal.
- **Baseline first:** before editing, run the full suite once
  (`"<python>" -m pytest tests/ -q -p no:cacheprovider`) and record the
  pass/fail counts. In a clean worktree the expected baseline is 425 passed,
  0 failed (the 2 failures seen in the main checkout come from an untracked
  `.ric` that is not in your worktree). Record whatever you actually get.
- **Gate B0** (when a task requires it): from the worktree root,
  `"<python>" docs/solver_verification/probes/baseline_gate.py`
  (≈2–4 min warm; it fixes its own import path and prints provenance). It
  must print `B0 RESULT: PASS`: step counts EXACT and floats within the
  default rtol 1e-9. (Differences ≤ ~2e-10 with identical step
  counts are Numba fastmath compile-state noise — see the script docstring —
  not a trajectory change.) A FAIL is a finding: report the full output, do
  not loosen tolerances or "adjust" the baseline.
- **Commits:** commit your work to your worktree branch (one or a few
  focused commits, conventional style e.g. `fix(burn_rate): ...`), each
  message ending with the line
  `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
  Never push.

## Hard rules

- Never edit `motors/*.ric` (openMotor-owned). Temp copies under
  `C:/tmp/sv_<task>/` are fine.
- Do not refactor, reformat, or "improve" code outside your task's spec.
- Diagnostics never write simulation state (bookkeeping must not perturb
  the trajectory — B0 enforces this).
- Do not change defaults or physics unless the task says so.

## Worklog (the review artifact)

Write `docs/solver_verification/worklog/T<n>.md` (your task only, so
parallel tasks never conflict) and commit it. Include: preflight output;
baseline suite counts; every file changed (one line each, what and why);
every gate run with its ACTUAL output numbers (paste tables, not "passed");
every measurement requested; anything that deviated from this spec and
why; open items. Write for a reviewer who has not seen your session.

Canonical knobs for every measurement run (same as the spike probes):
`roughness=32e-6, kappa=0.44, T_ignition=756.0, print_interval=0.0,
verbose=False`, plus `pyrogen='bpnv', pyrogen_mass=None` unless the task
says otherwise. Motors: `motors/machbusterNew.ric` (= "Chunc"),
`motors/BALLSstick.ric`, `motors/hasegawa_a.ric`. Grid resolution is the
`target_propellant_cells` kwarg of `run_from_ric` (default 100). Cadence is
`burn_update_interval`; CFL is `cfl_target` (default 0.3).

---

## WAVE 1 (parallel)

### T1 — igniter material: lossless `.ric` round-trip + stale-material warning

Files: `srm_1d/openmotor_adapter.py`, `tests/test_ric_igniter.py`,
`srm_1d/DEVNOTES.md`. Read DESIGN §2 "NEW BUG".

1. Read the `Pyrogen` dataclass (`srm_1d/propellant.py`), `load_pyrogen`,
   `_pyrogen_to_block`, `_block_to_pyrogen`, `load_igniter`,
   `default_igniter_block` (`openmotor_adapter.py`). Before editing, record
   in the worklog: the list of `Pyrogen` init fields; which are missing from
   `_pyrogen_to_block`; and any field whose absent-key default differs
   between `load_pyrogen` and `_block_to_pyrogen`. Also capture, with a
   scratch script, the `Pyrogen` produced by `load_pyrogen` for every
   `srm_1d/pyrogens/*.yaml` (for the no-change check in step 2).
2. Add one private helper `_pyrogen_from_mapping(d)` that builds a `Pyrogen`
   from a dict using the defaults `load_pyrogen` uses today (the YAML loader
   is the reference for defaults). `load_pyrogen` (keeping its
   required-key check) and `_block_to_pyrogen` both call it. Verify every
   library YAML loads to an identical `Pyrogen` as before the change.
3. `_pyrogen_to_block`: emit EVERY `Pyrogen` init field, derived from
   `dataclasses.fields(Pyrogen)` so a future field cannot be dropped again
   (key names = field names, which matches today's keys). Include only
   `init=True` fields.
4. Staleness warning: add `class StaleIgniterMaterialWarning(UserWarning)`
   in `openmotor_adapter.py`. In `load_igniter`, after building the embedded
   `Pyrogen`, find the library YAML(s) in `srm_1d/pyrogens/` whose `name`
   equals the embedded `name` case-insensitively. If at least one matches
   and the embedded material differs from EVERY matching entry on any field
   other than `name` (floats compared with rel tol 1e-9; None vs a value =
   different; an absent key compares as its default), emit ONE
   `warnings.warn(..., StaleIgniterMaterialWarning)` whose message names the
   motor's material, the library file, and each differing field as
   `field: embedded X vs library Y`, ending with: "The embedded values are
   used. Re-select the igniter material or re-seed the motor's igniter
   block to pick up library corrections." Never modify the values. No
   library match → no warning.
5. Tests (append to `tests/test_ric_igniter.py`):
   a. Round-trip: a `Pyrogen` with EVERY init field set to a non-default
      value survives `_block_to_pyrogen(_pyrogen_to_block(p))` field by
      field, including through `yaml.safe_dump` → `yaml.safe_load`.
   b. The block's key set equals the set of `Pyrogen` init field names.
   c. `load_igniter` on a motor dict whose igniter block is
      `default_igniter_block('bpnv')` emits NO
      `StaleIgniterMaterialWarning`; one whose material has `M=0.03` and no
      `gas_mass_fraction` key DOES, and the message mentions `M` and
      `gas_mass_fraction`.
   Do NOT assert anything about the committed motor files' current state.
6. Measurements (worklog only):
   a. Table: for every committed `motors/*.ric` with an igniter block, the
      fields that differ from the matching library entry.
   b. Chunc impact, `t_max=1.0`: run A = `pyrogen=None` (embedded block —
      the GUI path; same knobs otherwise) vs run B = canonical
      `pyrogen='bpnv', pyrogen_mass=None`. Report peak, t_peak, plateau
      (median P_head over t∈[0.4,0.9] s), ratio for both.
   c. Equivalence: copy the committed Chunc `.ric` to `C:/tmp/sv_t1/`,
      replace its igniter block with `default_igniter_block('bpnv')`
      (yaml load/dump of the TEMP copy), run it with `pyrogen=None`,
      `t_max=0.05`, and compare P_head to run B over the same window
      (report max |ΔP|, and whether step counts match). If not identical,
      explain the cause; do not change sizing logic to force a match.
7. Read-only investigation of the openMotor fork checkout
   (`C:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/openMotor/openMotor`):
   how does its GUI create/edit a motor's igniter block, and where does its
   igniter-material list come from (imports srm_1d? a vendored copy? its
   own file)? Report file paths and whether, after this fix, a newly
   GUI-created motor would embed current library values including
   `gas_mass_fraction`. Do NOT edit the fork.
8. DEVNOTES: one entry in the API Breaking Changes Log describing the
   complete igniter block + the warning.

Gates: B0 with `--only chunc` (the probes pass `pyrogen='bpnv'`, so it must
be identical); full suite = baseline count + your new tests, 0 new failures.

### T2 — erosive bisection: guaranteed bracket

Files: `srm_1d/burn_rate.py`, `tests/test_burn_rate.py`, `srm_1d/DEVNOTES.md`.
Read DESIGN §2 "Bisection clamp".

1. In `burn_rate_cell`'s Ma (code 0) bisection, after the existing
   `F_lo >= 0` early return: evaluate `F_hi` at `r_hi = r0 * 20.0` using
   exactly the same `h(β)` branch logic as the loop (the `< 1e-6` /
   `> 500` / exp branches). If `F_hi < 0.0` (root above 20·r0 — today it is
   silently clamped), set `r_hi = r0 + ebf * h0` (guaranteed bracket:
   `β/(e^β−1) ≤ 1` ⇒ `h ≤ h0` ⇒ `F(r0 + ebf·h0) ≥ 0`). Add a short comment
   stating that reason.
2. Raise the bisection iteration cap from 30 to 60, keeping the existing
   `(r_hi - r_lo) < 1e-8` break unchanged. Comment why defaults stay
   bit-identical: a valid ≤19·r0-wide bracket reaches 1e-8 in < 30 halvings
   for any r0 < 0.5 m/s, so the break fires exactly where it did before.
3. Change nothing else in the code-0 path (the registry rule "code-0 path
   literally identical" is relaxed ONLY for these two edits). Leave the
   `erosive_mp_cap` block, codes 1/2, and all signatures alone.
4. Tests (`tests/test_burn_rate.py`; reuse the propellant constants the
   existing direct `burn_rate_cell` call sites use):
   a. Envelope: D_hyd ∈ {2, 4, 8, 16, 32} mm × Re ∈ 10 log-spaced values
      1e4…1e7 × P ∈ {2, 7, 15} MPa (150 states). For each: result finite,
      `r_total ≥ r0`, and the Ma residual `|F(r_total)| ≤ 1e-7` m/s. To get
      `h0`, `ebf`, `alpha` for the residual, build them in a test helper by
      calling the SAME module-level helper functions `burn_rate_cell` calls,
      in the same order (cite the source lines in a comment). If
      `burn_rate_cell` computes them inline with no reusable helpers, STOP
      and report rather than re-deriving correlations in the test.
   b. Clamp regression: assert at least one envelope state has
      `r_total > 20·r0·(1 + 1e-6)` (the old clamp would have bound). If none
      do, extend to D_hyd = 1 mm and Re up to 3e7; if still none, STOP and
      report (that would contradict Codex's finding).
5. Worklog: number of envelope states with `r_total > 20·r0`, the max
   `r_total/r0`, and the max residual over the envelope.

Gates: B0 FULL (all three motors) must PASS. If it fails, a canonical motor
was being clamped — STOP and report the full B0 output (important finding;
do not proceed to "fix" anything else). Full suite = baseline + new tests.

### T3 — pressure-floor accounting

Files: `srm_1d/solver.py`, `srm_1d/simulation.py`,
`srm_1d/tools/ignition_diagnostics.py`, new `tests/test_pressure_floor.py`,
`srm_1d/DEVNOTES.md`. Read DESIGN §2 "Pressure floor".

1. `solver._piso_step_with_energy_diagnostics`, STEP 4: before the clamp,
   for each cell with `P_new[i] < pressure_floor` accumulate
   `n_floor_cells += 1`,
   `floor_mass_added += (pressure_floor - P_new[i]) / (R_arr[i] * T_new[i]) * A_port[i] * dx`
   (the mass the clamp adds relative to the EOS at the unclamped pressure),
   `floor_energy_added += that_mass * Cp_arr[i] * T_new[i]`, and
   `floor_min_raw_P = min(raw P_new)` over clamped cells (start at +inf).
   Then clamp exactly as today. Do NOT change the floor value or behaviour.
   Append `(n_floor_cells, floor_mass_added, floor_energy_added,
   floor_min_raw_P)` to the END of the return tuple. Update `piso_step`
   (the wrapper below it) and the `_run_time_loop` call site in
   `simulation.py` to unpack the new values.
2. `simulation.py`: per-step history arrays `pressure_floor_cells` (count)
   and `pressure_floor_mass` (kg), following the exact pattern used for
   `min_pressure` (allocation, recording, return, result dict). Summary
   keys: `pressure_floor_activations` (Σ cell-steps),
   `pressure_floor_steps` (steps with ≥1 clamp),
   `pressure_floor_mass_added_kg`, `pressure_floor_energy_added_J`,
   `pressure_floor_first_time_s` (−1.0 if never),
   `pressure_floor_min_raw_P` (NaN if never). Add the two per-step keys to
   `step_diagnostics_timeseries` in `tools/ignition_diagnostics.py`.
3. Tests (`tests/test_pressure_floor.py`): (a) a short Chunc run
   (`t_max=0.05`, default cells) has all summary keys and the per-step
   arrays with history length; if it records zero clamps assert zero,
   otherwise assert internal consistency (activations ≥ steps ≥ 1, mass ≥ 0)
   and record the numbers. (b) If you can construct, in reasonable effort, a
   tiny-N direct call of `_piso_step_with_energy_diagnostics` that forces a
   clamp, assert the counters; if not, skip (b) and say so in the worklog.
4. Measurement (worklog table): Chunc `t_max=0.05` at
   `target_propellant_cells` ∈ {100, 200, 400, 800}, and BALLSstick
   `t_max=0.05` at {100, 400}: actual cell count, activations, steps,
   mass added (kg and as a fraction of `summary['mass_produced']`), energy
   added, first time, min raw P, global min of the `min_pressure` history.
   (Codex: the floor starts activating at ≥400 cells.)
5. DEVNOTES API log: the `_piso_step_with_energy_diagnostics` / `piso_step`
   return tuple grew.

Gates: B0 FULL must PASS (accounting only). Full suite = baseline + new.

### T4 — verification tests (tests only, no source edits)

File: new `tests/test_verification_limits.py`. If a test reveals a real
defect, do NOT fix source — mark the test `xfail(strict=True)` with a reason
string and report it prominently.

a. **Tab-boundary parity vs openMotor.** Import `motorlib` the way
   `tests/test_srm1d_plugin.py` does (`_setup_openmotor_path`; skip the
   module cleanly if unavailable). Build an openMotor propellant with three
   tabs incl. a gap — e.g. [0.1, 2], [2, 6], [8, 15] MPa with distinct
   a/n — and the srm_1d tab arrays FROM THE SAME DATA via the adapter's
   propellant conversion (so units are handled by the adapter, not the
   test). Pressures: every tab min and max exactly, the shared boundary
   2 MPa, in-gap 7 / 6.9 / 7.1 MPa, below all (0.05), above all (20), and
   each tab's midpoint. Assert srm_1d's `select_tab_idx` /
   `saint_robert_from_tabs` selects the same tab (compare a and n) as
   `getCombustionProperties`.
b. **Zero crossflow.** For each registry code (0, 1, 2) and Re ∈ {0, 50, 99}
   (below the `Re < 100` guard): `r_total == saint_robert_from_tabs(P)`
   exactly and `r_erosive == 0`.
c. **Normal component invariance.** Grid P ∈ {1, 5, 12} MPa ×
   Re ∈ {1e4, 1e5, 1e6} × D_hyd ∈ {10, 40} mm: for codes 0 and 1,
   `r_total − r_erosive == r0` (rel 1e-12) and `r_erosive ≥ 0`; for code 2,
   `r_erosive == 0`.
d. **Mode invariants (integration).** Measure warm runtimes first and pick
   the cheapest committed motor that burns out fully (`summary['termination']`
   should indicate burnout / P_cutoff, not `t_max`) in ≤ 30 s/run warm —
   e.g. `example_bates.ric`, `hlbltest.ric`, `L1843.ric`; record the
   candidates' timings. For `erosive_model` ∈ (ma2020, mukunda_paul,
   saint_robert): assert propellant mass delivered (`summary['mass_produced']`)
   agrees within 0.5% across modes and total impulse
   (`compute_motor_performance(...)['total_impulse']`, `srm_1d.nozzle`;
   `run_from_ric` returns the result plus the objects it needs — read its
   return statement) is within 3% of ma2020's. Record actual Δimpulse and
   Δburn-duration per mode (Codex, erosion off: +0.91% / +4.35%).

Gates: no source files changed (`git diff --stat` touches only the new test
file + your worklog); full suite = baseline + new (xfails listed).

---

## WAVE 2 (after wave 1 is merged)

### T5 — split burn / geometry refresh + refresh-gap diagnostics

Files: `srm_1d/simulation.py`, tests (new `tests/test_refresh_cadence.py`),
`srm_1d/DEVNOTES.md`; `openmotor_adapter.run_from_ric` only if its kwargs
do not already pass through (`**sim_overrides` should).

1. `run_simulation` gains `geometry_update_interval=None`; `None` → the
   resolved `burn_update_interval` (so defaults are unchanged). Thread it
   into `_run_time_loop` next to `burn_update_interval`. The
   `update_cell_geometry` block (~1645) uses `geometry_update_interval`; the
   `compute_burn_rates` block (~1666) keeps `burn_update_interval`. List in
   the worklog every OTHER use of `burn_update_interval` in the loop and
   leave them on `burn_update_interval`.
2. Track the simulated time between consecutive refreshes (Σdt since the
   last refresh, per kind). Summary keys:
   `burn_refresh_interval_steps`, `geometry_refresh_interval_steps`,
   `max_burn_refresh_gap_s`, `max_geometry_refresh_gap_s`.
3. Test: Chunc `t_max=0.05`: (a) keys present; (b)
   `burn_update_interval=1, geometry_update_interval=1` vs default — record
   peak and global min P for both; assert peak within 0.5%.
4. Worklog: the two max-gap values for Chunc at default (Codex: 29.84 µs).

Gates: B0 FULL PASS (defaults unchanged). Full suite = baseline + new.

### T6 — energy audit

Files: `srm_1d/solver.py`, `srm_1d/simulation.py`,
`srm_1d/tools/ignition_diagnostics.py`, new `tests/test_energy_audit.py`,
`srm_1d/DEVNOTES.md`. Read DESIGN §2 "EOS reconciliation".

A. **Ledger completeness.** List every term added to `thermal_source` in
   `_run_time_loop`. Give each term that lacks its own per-step power history
   one, at minimum `convective_wall_sink_power` (the v0.7.4 convective
   wall heat-loss sink, ~1245, currently folded into `thermal_source`).
   Record per step `source_ledger_mismatch_W` = Σ(channels) −
   `thermal_source_power`. Compare the PISO's `nozzle_scalar_flux_power`
   with the post-PISO nozzle reconstruction in STEP 6
   (`_nozzle_boundary_flow` on final P/T) → `nozzle_ledger_mismatch_W`
   per step; explain any by-design difference.
B. **Between-step storage jumps (energy J and mass kg per step).** Decompose
   the gap between step n's stored Σm·Cp·T (and Σm) and step n+1's
   pre-PISO value into: (i) EOS reconciliation (`rho_new` from P/(R·T_new)
   vs the `rho_new_1` mass the ledger used), (ii) pressure-floor clamp (T3
   records the per-step MASS `pressure_floor_mass` but only a run-total
   energy — add a per-step `pressure_floor_energy` history [J] from the PISO's
   `floor_energy_added`, same pattern), (iii) mixture refresh (Cp_mix change from
   `_advect_species` + `_refresh_mixture_arrays` at fixed T), (iv) geometry
   refresh (`A_port` change in `update_cell_geometry`), (v) remainder.
   Implement as pure bookkeeping (evaluate Σ before/after each operation;
   never alter state). History keys `storage_jump_{eos,floor,mixture,
   geometry,remainder}_J` and `mass_jump_{...}_kg`.
C. **`sensible_enthalpy_audit(result)`** in `tools/ignition_diagnostics.py`
   returning: the convention string; max/RMS `energy_residual`; max |source
   ledger mismatch|; max |nozzle mismatch|; per jump category max single,
   Σ|·| and Σ signed (J and kg); the integrated run budget (sources in,
   nozzle out, clip, floor, jumps, Δstorage) with closure error as a
   fraction of total source energy; and the same split by phase — fill
   (t < first time the head pressure reaches 50% of its peak), spike
   (± 5 ms around the peak), plateau (the rest) — using `pressure_landmarks`
   where it helps.
D. **Omitted-term magnitude** (post-processing of `result['snapshots']`,
   which carry `Mach`, `u`, `T`, `rho`, `P`, `D_port`): per snapshot, the
   kinetic/sensible ratio `(γ−1)/2·M²` per cell (use the propellant's γ;
   note the approximation) → bore max, a mass-weighted mean
   (weight ρ·D_port²), and the nozzle-inlet cell (N−1) value (= the
   stagnation-vs-static enthalpy outflow fraction). Report each by phase.
   If snapshots lack a needed field, STOP and report which.
E. Summary keys `energy_convention = 'sensible_enthalpy_CpT'` and
   `energy_omitted_terms = ['kinetic_energy', 'pressure_work',
   'formation_energy_implicit_via_T_flame', 'solid_energy_storage',
   'nozzle_kinetic_energy_static_outflow']` (keep the existing
   `energy_residual_convention` key).
F. Cost: report Chunc `t_max=1.0` steps/s before vs after. If the slowdown
   exceeds 10%, put part B behind a new `run_simulation` kwarg
   `energy_audit=False` (default off) and say so; otherwise always on.
G. Tests: audit keys present on a short Chunc run; max |source ledger
   mismatch| ≤ 1e-6 × max thermal_source_power; jumps finite.
H. Measurements (worklog): run the audit on Chunc (`t_max=1.0`) and
   BALLSstick (`t_max=2.0`) — full tables incl. the phase split, next to
   Codex's numbers (source 1.30e-8 W, nozzle 2.05e-8 W, residual 1.79e-5 W,
   max jump 2.63 J, Σ|jump| 13.30 kJ).

Gates: B0 FULL PASS (bookkeeping only). Full suite = baseline + new.

---

## WAVE 3 (measurement only; after T3 + T5 are merged)

### T7 — timestep / refresh sensitivity sweep

File: new `docs/solver_verification/probes/cadence_sensitivity.py` (+ its
results in the worklog). No source edits.

Motors (canonical knobs, ma2020, cap 0): Chunc `t_max=1.0`; BALLSstick
`t_max=2.0`; Hasegawa A `t_max=10.0, P_cutoff=0.05e6`. Configs
(burn_update_interval / geometry_update_interval / cfl_target):
default/default/0.3 (baseline), 1/default/0.3, default/1/0.3, 1/1/0.3,
default/default/0.15, 1/1/0.15 → 18 runs. Support `--motor` / `--config`
flags and append each row to a JSONL so the sweep can be resumed.

Metrics: peak, t_peak, motor plateau/trough/mid and ratios exactly as
defined in `docs/v0_7_4/probes/mp_cap_crossmotor.py`; global min of the
`min_pressure` history + its time AND `summary['pressure_floor_min_raw_P']`
(the `min_pressure` history pins at exactly 1000 Pa once the floor clamps,
hiding the true undershoot — T3); max of `max_mach` for t < t_peak (fill
Mach); T3 floor counters; T5 max refresh gaps; n_steps; wall s; steps/s.
Stationarity of each plateau window: linear-fit slope (%/s of the window
mean) and std/mean.

Worklog: the full table plus each metric's relative change vs the baseline
config; flag every validation metric (Chunc peak/plateau/ratio;
BALLSstick pk/trough, pk/mid; Hasegawa peak) that moves > 1%.

### T8 — grid convergence + fill-Mach attribution + REOPENED §12 draft

**Updated after R1 (`CODEX_REVIEW.md` Q1/Q5/Q7) — read those sections first.**
Runs AFTER T10 merges (T10 changes the mass source). ADAPT the friend's
`scripts/run_spatial_extension.py` (part A) and
`scripts/run_fill_configuration_comparison.py` (part B) from remote
`codexfork` instead of writing probes from scratch (port P4): switch them to
the canonical config (B0 knobs, library BPNV), `t_max=1.0`, and report fill
Mach in BOTH windows — Codex's `max(max_mach[t<=3 ms])` and the full
pre-peak window — plus the time of the max. Run every grid ≥ 400 cells at
CFL 0.3 AND 0.15 (R1: the 806-cell floor events and much of the fill Mach
are temporal). Their 1600-cell headline (12.4478 MPa / Mach 0.4343) came
from a stale config with no committed data — do not compare against it.

Files: new `docs/solver_verification/probes/grid_convergence.py`;
`docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` (append a §12 marked
**DRAFT — for review**). No source edits.

A. Chunc canonical at `target_propellant_cells` ∈ {50, 100, 200, 400, 800,
   1600}, `t_max=1.0` (1600 may take tens of minutes — run it in the
   background, resumable via JSONL like T7). Metrics: actual cell count,
   peak, t_peak, plateau, fill Mach (max `max_mach` for t < t_peak), global
   min P AND `summary['pressure_floor_min_raw_P']` (raw undershoot — the
   history pins at 1000 Pa once clamping; T3 saw −734 Pa at 806 cells),
   floor activations, n_steps, wall s. Observed order from the three
   finest grids `p = ln((f3−f2)/(f2−f1))/ln 2` and Richardson
   `f∞ = f1 + (f1−f2)/(2^p − 1)` (f1 finest), per metric. Compare with
   Codex (12.4478 MPa / Mach 0.4343) and REOPENED §6 (peak climbing
   11.76/12.64/13.20 at 50/100/200, "→ ~14 MPa").
B. Fill-Mach attribution: Chunc `t_max=0.05` at cells ∈ {50, 100, 200} ×
   igniter ∈ {library BPNV (`pyrogen='bpnv'`), embedded legacy BPNV
   (`pyrogen=None` on the committed `.ric` — the GUI path), library MTV
   (`pyrogen='mtv'`)}: fill Mach + peak. Compare with Codex (MTV
   0.448/0.498/0.479; current BPNV 1.987/2.200/2.311; legacy
   3.271/4.821/11.722) and §6 (3.36/4.81/12.04).
   Reviewer data point (2026-10-05, plugin-path triage, 100 cells,
   committed Chunc, canonical knobs): max of the `max_mach` history for
   t ≤ 20 ms = 3.94 @ 2.48 ms with library BPNV and 4.81 @ 1.38 ms with the
   embedded legacy BPNV, both near the aft end (x ≈ 817–833 mm) — well above
   Codex's "current BPNV" 1.987/2.200/2.311. Find and state the reason for
   the difference (metric definition? cell-centre vs face Mach? config?)
   rather than assuming either is right.
   Lead from T1 (worklog/T1.md step 6b): Chunc with the embedded legacy
   igniter (`pyrogen=None`) peaks at 12.654 MPa @ 11.55 ms at 100 cells vs
   11.473 with library BPNV — §6 reports 12.64/12.65 MPa at 100 cells. Check
   the §6 probe scripts in `docs/v0_7_4/probes/` (e.g.
   `chunc_mach_convergence.py`) for whether they passed a `pyrogen` kwarg,
   and state whether §6's runs used the legacy embedded igniter.
C. Draft §12 "Solver verification (Codex-review replication)": the
   grid-convergence result, the fill-Mach attribution, and an explicit list
   of which §6 claims stand and which are corrected. Numbers only — no
   recommendations, no default changes.

---

## WAVE 4 — fixes approved from the T6 energy audit (user, 2026-10-05)

One agent, T9 then T10 (both edit `simulation.py`), in one worktree whose base
INCLUDES T6 (its audit channels are T10's conservation gate). Read
`worklog/T6.md` first.

### T9 — debit the left-neighbour radiation sink

T6 found that with `radiation_emissivity > 0` the left-neighbour radiation
sink is recorded in its power history but never subtracted from the gas
(`thermal_source`) — an energy leak (Chunc, emissivity 0.05: 292 W
source-ledger mismatch vs a 570 W recorded sink). It is pinned by a strict
`xfail` in `tests/test_energy_audit.py`.

1. Find the sink (T6 worklog names it) and debit it from `thermal_source`
   exactly as the existing radiation sink is debited — mirror that code; no
   new physics.
2. Turn the strict `xfail` into a normal passing test.
3. Gates: B0 PASS (canonical emissivity is 0 — must be unchanged); Chunc
   `radiation_emissivity=0.05`, `t_max=0.05`: max |source-ledger mismatch|
   ≤ 1e-6 × max `thermal_source_power` (was 292 W) — report before/after
   P_peak at that emissivity; full suite = baseline + 0 new failures.
4. DEVNOTES API/behaviour log entry. Commit separately from T10.

### T10 — subtract the gas that fills newly opened port volume (physically standard mass source)

T6 found the geometry refresh adds gas from nothing: the solver advances
continuity at fixed `A_port` and then grows `A_port` at fixed ρ, so each
refresh adds ρ_gas·ΔA·dx (Chunc +6.4 g / BALLSstick +9.0 g over the burn,
always positive) plus its sensible enthalpy. The standard fixed-area
formulation is a net source of **(ρ_p − ρ_g)·r·P_b**: the burned propellant
first fills the volume it vacates. User decision: adopt it, accepting that
every result shifts slightly and the B0 baseline is re-recorded.

1. In `_run_time_loop`, wherever the regression mass source is assembled
   (sidewall `r·C_burn`-type term AND the end-face source — both open new
   volume), subtract the local gas that fills the opened volume, using the
   SAME effective regression rate × area as the source (keep the burnout
   ramp `f_active` consistent — DEVNOTES gotcha "Mass conservation"):
   - mass: `−ρ_g · (r · burning area per unit length)`;
   - energy (`thermal_source`): `−ρ_g · Cp_mix · T_gas · (same)` — the
     sensible enthalpy of that gas (local cell mixture Cp and T);
   - species: per species `−ρ_g · Y_k · (same)` (the gas filling the volume is
     the LOCAL mixture, not propellant gas).
   Use the cell state at the start of the step. Momentum: leave unchanged
   and say so in the worklog (not measured by T6).
2. Gates (behaviour CHANGES — report, don't hide):
   a. Conservation: with T6's audit, over Chunc `t_max=1.0` and BALLSstick
      `t_max=2.0`, report Σ geometry storage jump vs Σ new debit (mass and
      energy) — they should now cancel to round-off — and
      `summary['mass_balance_error']` before vs after.
   b. B0 before/after table for all three motors (peak, t_peak, plateau or
      trough/mid, n_steps). Then RE-RECORD the `BASELINE` constants in
      `docs/solver_verification/probes/baseline_gate.py` from this commit,
      with a dated comment naming T10, and confirm B0 PASSes against them.
   c. Validation vs experiment before/after: Chunc peak / plateau / ratio
      (exp plateau 8.66 MPa, ratio 1.027), BALLSstick pk/trough and pk/mid
      (~1.50 / ~1.35), Hasegawa A global peak (~6.5 MPa). Expected order of
      magnitude: ~ρ_g/ρ_p ≈ 0.5% less mass source → roughly −1% on plateau
      pressure. If ANY metric moves more than 3%, STOP and report before
      re-recording anything.
   d. Full suite: list every test whose pinned simulation numbers moved
      (old → new). Update a pinned value only where a test pins a sim
      output, never loosen a tolerance, and flag each in the worklog.
3. DEVNOTES: behaviour-change log entry + the Calibration State note
   (canonical validation values updated). Do NOT edit CLAUDE.md or any
   memory — the reviewer does.

---

## WAVE 5 — ports from the Codex review (R1, `CODEX_REVIEW.md`)

The review's §3 table is the authoritative spec for each port (files,
expected conflicts, gates); read the matching per-commit section in §1 and
the friend's code (remote `codexfork`, read-only) before porting. Ports are
ADAPTATIONS into our tree — never cherry-pick their commits (they predate
the registry, T1–T6, T9/T10 and the igniter provider). Port their tests where
the table says so, adapted to our APIs.

### PA — non-solver ports (P0 slag fix + P3 refresh-delay statistic + P6 B0 provenance)

1. **P0 slag unit fix.** openMotor stores `nozzle.slagCoeff` in (m·Pa)/s and
   uses `slagRate = slagCoeff / P[Pa]` (`motorlib/motor.py` ~452); srm_1d
   passes it through unconverted and divides by P in MPa
   (`simulation.py` `s_rate = slag_coeff / P_throat_MPa`) → 1e6× too fast for
   any nonzero slag. Convert at the adapter boundary (÷1e6 → (m·MPa)/s, the
   srm_1d convention — UNITS stay human-readable per
   `feedback_openmotor_alignment`), fix the docstrings/comments R1 lists, and
   add `test_slag_unit_conversion` asserting srm_1d's throat slag rate equals
   openMotor's formula for the same nozzle + pressure. B0 must PASS (committed
   motors have slag 0). DEVNOTES API-log entry: user `.ric` files with slag > 0
   now behave correctly (previously 1e6× too fast).
2. **P3 refresh-delay statistic.** Post-hoc
   `ignition_refresh_delays(result)` in `tools/ignition_diagnostics.py` (per
   ignited cell: time from ignition to the next burn-rate refresh; summary
   max/median), test on a short run, and reproduce the ~29.8 µs-scale value on
   canonical Chunc. Append a correction note to `worklog/T5.md` and
   `worklog/T7.md`: Codex's "29.84 µs" is THIS statistic, not the max refresh
   gap (R1 reproduced 2.984e-5 s).
3. **P6 B0 provenance.** `probes/baseline_gate.py` prints git HEAD + dirty
   flag, numpy/numba versions, and the sha256 of each `.ric` it runs. Output
   only — gate logic unchanged.

Gates: B0 PASS; full suite = baseline + new tests, 0 new failures.

### PB — solver-touching ports (P1 hardening + P2 limit diagnostics) — AFTER T9/T10 merge

As specified in `CODEX_REVIEW.md` §3 rows P1 and P2. P1: validate
`burn_update_interval >= 1`; `D_hyd <= 0 → (r0, 0)` in `burn_rate_cell`;
optional `pressure_floor_pa` kwarg (default 1e3, bit-identical). P2: per-cell
T-floor / T-ceiling / pressure-floor activation counts, durations, first/last
times, |energy| for the T limits, per-cell max floor deficit, optional
port-Mach-cap face counts; summary view `numerical_limits` + a
`tools/ignition_diagnostics` helper. Gates: B0 (rtol 1e-9,
exact steps); byte-identical trajectory A/B at 106/406/806 cells
(`worklog/T3_scripts/trajectory_ab.py`); steps/s cost < 10%; full suite.

### P5 — thread-F discriminator: energy-coupled PISO iteration (EXPERIMENT)

**Approved by the user after T8 (2026-10-11).** An isolated experiment on its
own branch — it is NOT merged by default; the reviewer and user decide
afterwards from the numbers.

Question (spike thread F, `docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` §8 F and
§12): is the energy–pressure lag — continuity closed on
`rho_new_1 = P/(R·T_old)`, then the state re-derived as `rho = P/(R·T_new)`
(T6: ~30% of fill-phase source energy carried outside the ledgers) — what
drives the fill pathology T8 measured: fill Mach that grows WITHOUT BOUND
with refinement (CFL 0.3: 3.9 → 13.9 → 16.9 → 49.3 at 106→1606 cells; CFL
0.15: 3.4 → 7.25), the temporal pressure-floor events (806 cells @ CFL 0.3,
raw P down to −23 kPa at 1606), and — the real target — does it move the
startup spike (grid-converged ≈ 12.55 MPa vs the 100-cell 11.28)?

1. **Implement an opt-in energy-coupled iteration** in the PISO step, kwarg
   `piso_energy_coupling_iters=0` on `run_simulation` (0 = today's step,
   bit-identical; B0 must PASS with 0). With K > 0: after the energy update
   gives `T_new`, rebuild the density/compressibility the pressure correction
   uses with the updated T (`P/(R·T_new)`), redo the pressure correction(s)
   and the energy update, and repeat up to K times or until
   `max|T_k+1 − T_k|/T_k < 1e-6`. Record per step the iterations used and the
   final relative T change. You decide the exact insertion point after
   reading `_piso_step_with_energy_diagnostics` end to end — document the
   algebra and the choice in the worklog (this is research: implementation
   judgement is expected, but every choice is written down and defaults never
   change).
2. **Measure** (canonical Chunc, library BPNV, post-T10 tree), K ∈ {0, 1, 3}
   (add 10 if K=3 has not converged the EOS jump):
   a. EOS storage jump (T6 `storage_jump_eos_J` / `mass_jump_eos_kg`) and the
      fill-phase closure fraction (`sensible_enthalpy_audit`) — does it go to
      ~0?
   b. Fill Mach (max of `max_mach` pre-peak, and its time), min P, raw floor
      P and floor activations at 106 / 406 / 806 cells × CFL 0.3 / 0.15
      (t_max 0.05 is enough for these).
   c. Startup peak, t_peak, first-ignition time, time until 90% of grain
      cells are ignited (`ignition_time_by_cell`), and plateau (t_max 1.0) at
      106 and 406 cells, CFL 0.3.
   d. BALLSstick at 106 cells, CFL 0.3: peak, pk/trough, pk/mid, fill Mach.
   e. Cost: steps/s and mean iterations per step.
   Run fine-grid rows ONE AT A TIME (history memory ~450 B/step).
3. Deliverables: `worklog/P5.md` with every table, and a §13 in
   `IGNITION_SPIKE_REOPENED.md` marked **DRAFT — for review**, numbers only,
   stating which of {fill Mach bounded?, floor events gone?, spike moved?}
   the iteration changes. No default changes; no recommendations.
