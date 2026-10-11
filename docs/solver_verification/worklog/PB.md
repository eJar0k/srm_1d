# PB - P1 hardening + P2 numerical-limit diagnostics (worklog)

Branch `sv/pb`, worktree `C:/tmp/sv/pb`, based on `a611b67`. Executed 2026-10-11.
Adaptations of Codex `e46562b`/`e749f71` (P1: validation, D_hyd guard),
`86c587f` (pressure_floor_pa), `d08bcc5`/`d773b55` (P2 per-cell limit tallies).
No default, equation or trajectory change.

## Preflight

```
$ PYTHONPATH="$PWD" python -c "import srm_1d; print(srm_1d.__file__)"
C:\tmp\sv\pb\srm_1d\__init__.py
```
`SRM1D_OPENMOTOR_PATH` exported on every Python command. pycache deleted after every @njit edit.

## Suite

Baseline (clean worktree, before any edit): **695 passed, 0 failed** (499 s).
Final: **719 passed, 0 failed** (= 695 + 24 new). The first full run after the P2 code had 2 failures
in my own new P1 test (it passed `pressure_floor_pa` positionally, which P2 shifted); fixed by passing it
as a keyword (in the P2 commit); the two new test files re-run: 24 passed. The 719 figure is 717 passed
in the full run + those 2 re-verified.

## Commits

| commit | content |
|---|---|
| P1 `52c9a97` | `pressure_floor_pa` kwarg (+validation, summary key), `burn_update_interval` validation, `D_hyd <= 0` guard; `tests/test_hardening.py` |
| P2 (next commit) | per-cell limit tallies, `result['numerical_limits']`, `numerical_limits_summary`, `tests/test_numerical_limits.py`, DEVNOTES entry (covers P1+P2), this worklog |

## Files changed

- `srm_1d/solver.py`: `LIMIT_*`/`N_LIMIT_DIAGNOSTICS` constants; `_piso_step_with_energy_diagnostics`
  gains 12 optional trailing args (tally arrays, `step_start_time_s`, `pressure_floor_pa`); return tuple
  (21) UNCHANGED. `pressure_floor = pressure_floor_pa`. Mach-cap face tallies in the cap block; floor
  tallies + max deficit inside the existing T3 accounting loop; T-clip tallies in a SEPARATE loop after the
  energy loop (see "Deviation" below).
- `srm_1d/simulation.py`: `_run_time_loop` gets `pressure_floor_pa` + 10 tally arrays; `run_simulation`
  allocates them, validates `pressure_floor_pa` and `burn_update_interval`, adds `numerical_limits` to the
  result (builder `_build_numerical_limits`) and `summary['pressure_floor_pa']`.
- `srm_1d/burn_rate.py`: `D_hyd <= 0 -> (r0, 0.0)` in `burn_rate_cell` (before the friction factor).
- `srm_1d/tools/ignition_diagnostics.py`: `numerical_limits_summary(result)`.
- `tests/test_hardening.py` (P1, 8 tests), `tests/test_numerical_limits.py` (P2, 16 tests).
- `srm_1d/DEVNOTES.md`: API log entry.

## Design notes

- All tally arrays default to `None` in the PISO kernel (Numba prunes the branches), so existing direct
  callers (`piso_step`, tests) are untouched. Passing `limit_activation_counts` requires the other five
  `limit_*` arrays (Codex's risk note: first/last arrays are indexed whenever counts are given); the Mach
  arrays are optional as a group. `run_simulation` always passes all of them.
- Times: `first` = start time of the first active step; `last` = END of the last active step
  (`step_start + dt`), as in Codex. `-1.0` sentinel in the kernel -> NaN in the result.
- Energy for the T limits = sum |clip correction power| * dt, where the signed power is exactly the term
  added to `clipping_correction_power` (tested: single-step equality with `|out[10]| * dt`).
  The pressure floor has no energy tally here (T3 already carries floor mass/energy).
- A cell hitting the floor and then the ceiling on the same step is labelled ceiling (Codex's ordering,
  which mirrors the clip order).
- `burn_update_interval`: previously `0` meant `step % 0` inside the njit loop. Non-integer floats are
  rejected; integer-valued floats/numpy ints are accepted.

## Deviation (finding): T-clip tallies in a separate loop

First implementation put the tally inside the energy loop (Codex's placement). The byte-level A/B then
showed exactly two arrays differing (`clipping_correction_power`, `energy_residual`; 3446 of 30462 steps
at 106 cells, ~1 ULP, rel ~5e-16); every STATE array (P, T, rho, u, ...) was byte-identical. Fastmath
re-associated the summation of `clipping_correction_power` once the loop body gained stores. Fix: the
energy loop is restored to its original arithmetic with one extra store `T_raw_arr[i] = T_raw`; the tally
re-derives the clip in a separate loop afterwards. Result: all 80 arrays byte-identical (below).
Cost: one `np.empty(N)` per step (inside noise, see steps/s).

## Gates

### B0 (final code; warm cache)
```
[chunc] peak 11.28243546264288 rel 0  t_peak 0.016286760473606336 rel 0  plateau 8.150933877646011 rel 0  n_steps 531883 OK
[balls] peak 8.331331016878805 rel 0  trough 4.343222420159596 rel 0  mid 4.879629094113075 rel 0  n_steps 843274 OK
[haseg] peak 6.108574454139173 rel 0  t_peak 2.391606230911318 rel 0  n_steps 2336140 OK
B0 RESULT: PASS (n_steps exact, floats within rtol 1e-09)
```

### Byte-identical trajectory A/B (`T3_scripts/trajectory_ab.py`, adapted: RIC path; the pre tree is a
`git archive HEAD srm_1d` copy). Chunc `t_max=0.05`, canonical knobs. Like-for-like compile state: 106 + 406
in one fresh-compile process per tree; 806 in a second (cache-loaded) process per tree. 80 common arrays
(the new `numerical_limits` dict is not an ndarray key, so it is outside the set).

| grid (actual) | steps | floor cell-steps (T3) | byte mismatches |
|---|---|---|---|
| 106 | 30462 | 0 | **0 (IDENTICAL)** |
| 406 | 122639 | 26 | **0 (IDENTICAL)** |
| 806 | 246254 | 193 | **0 (IDENTICAL)** |

(The abandoned in-loop variant gave 2 mismatches of 80 at every grid; see Deviation.)

### Cost, Chunc `t_max=1.0` (531883 steps, warm, after an in-process warm-up), steps/s

| tree | reps | mean |
|---|---|---|
| pre (`a611b67`) | 50556 / 49595 / 50820; re-measured alongside post: 49666 / 49904 / 49625 / 49823 | ~49.8k (first batch ~50.3k) |
| post (final) | 49598 / 49106 / 49047 / 48928 | ~49.2k |

Cost about **-1.1%** (< 10% gate). (The abandoned in-loop variant measured 50.5 / 50.3 / 50.0k.)

## Headline: canonical Chunc, 100 target cells (106 actual), library BPNV (`numerical_limits_summary`)

Identical at `t_max = 0.05` and `t_max = 1.0` (all activity ends before 31 ms). Cell-steps = 106 x steps:
3,228,972 (0.05 s) / 56,379,598 (1.0 s).

| limit | activations (cell-steps) | active cells | first t [s] | last t [s] | sum abs energy [J] | frac. of cell-steps (0.05 s / 1.0 s) |
|---|---|---|---|---|---|---|
| T ceiling | 256,570 | 106 | 2.13e-5 | 0.03089 | 525.78 | 7.95% / 0.455% |
| T floor | 1,977 | 106 | 0.0 | 0.002437 | 14.16 | 0.061% / 0.0035% |
| pressure floor | 0 | 0 | n/a | n/a | n/a (max raw deficit 0 Pa) | 0 |
| Mach cap (off) | 0 | 0 | n/a | n/a | n/a | n/a |

Codex measured 254,483 ceiling / 1,975 floor / 532 J on their older tree; ours (with T10) is 256,570 /
1,977 / 540 J (525.8 + 14.2). Where: ceiling activations peak at x = 0.069 m, ceiling energy at
x = 0.020 m (near the head); the T floor fires in every cell during the first 2.4 ms (initial fill), with
the highest count at x = 0.858 m (aft end). The net signed `clipping_correction_power` integral is bounded
by the |energy| tally (tested).

## Tests (24 new)

`tests/test_hardening.py`: default floor bit-identical to explicit 1e3; override 750 clamps to 750;
`pressure_floor_pa` in {0, -1, nan, inf} rejected; summary key; `D_hyd` in {0, -0.01} -> (r0, 0);
`burn_update_interval` in {0, -3, 1.5} rejected.
`tests/test_numerical_limits.py`: ceiling clip single-cell (count/duration/energy/first/last), energy ==
|out[10]|*dt, T-floor single-cell (enthalpy sink), floor override deficit 250 Pa (Codex's case), tallies
leave every PISO output bit-identical (ceiling and floor), Mach-cap face tallies (cap 1e-6), run-level
shapes/NaN/first<=last, enabled Mach cap on a run, floor override reported, summary helper consistent, and
two Chunc 406-cell 50 ms tests tying the per-cell pressure-floor map to the T3 run counters (activations,
first time, min raw P) and bounding the net clip power by the |energy| tally.

## STOP conditions / open items

- None hit. One surprise (fastmath ULP in the in-loop variant), resolved above.
- "Summary view `numerical_limits`" was implemented as the top-level `result['numerical_limits']`
  per-cell dict plus `tools.ignition_diagnostics.numerical_limits_summary(result)` (scalar view); I did
  not add the totals as `summary[...]` scalar keys (only `summary['pressure_floor_pa']`).
- DEVNOTES entry is one combined P1+P2 entry (lives in the P2 commit).
- Scratch (A/B dumps, logs) in `C:/tmp/sv_pb/` (not committed).
