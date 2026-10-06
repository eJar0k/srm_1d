# PA worklog — P0 slag unit fix, P3 refresh-delay statistic, P6 B0 provenance

Branch `sv/pa`, rebased onto `openmotor-frontend` (T9 + T10 merged; tip df35e80).
Their code (`codexfork/main`) was read only; everything is adapted, nothing cherry-picked.

## Preflight
`PYTHONPATH=$PWD python -c "import srm_1d; print(srm_1d.__file__)"` -> `C:\tmp\sv\pa\srm_1d\__init__.py` (inside worktree, before and after rebase).
`srm_1d/__pycache__/` deleted after the rebase (T9/T10 changed @njit code).

## Baselines
- Pre-rebase base (before my edits): **682 passed, 1 xfailed** (478 s).
- Rebased base (T9/T10 agent's measurement): 685 passed, 0 failed.
- Final (rebased + PA): **688 passed, 0 failed** (519 s) = 685 + 3 new tests (slag, 2 refresh-delay).

## P0 — slag unit fix
`slagCoeff` in committed motors: `grep slagCoeff motors/*.ric` -> **21 files, all `slagCoeff: 0.0`** (so B0 cannot move).
Files changed:
- `srm_1d/openmotor_adapter.py`: `convert_nozzle` divides `slagCoeff` by 1e6 ((m*Pa)/s -> (m*MPa)/s); module docstring (line ~32) and `convert_nozzle` docstring corrected; the CLI print (~1666) now labels the raw value (m*Pa)/s and shows the converted value.
- `srm_1d/nozzle.py`: `Nozzle.slag_coeff` doc comment says the .ric value is (m*Pa)/s and the adapter divides by 1e6 (field itself stays (m*MPa)/s).
- `tests/test_adapter.py`: `test_slag_unit_conversion` — 2e-3 -> 2e-9, AND runs `nozzle.compute_thrust_history` (same `slag_coeff / P_MPa` formula as the in-loop step 3b) at constant P = 5 MPa, dt = 1 ms: srm_1d rate = -dD/(2 dt); openMotor formula `slagCoeff / P[Pa]` = 2e-3 / 5e6 = 4.0e-10 m/s; they agree to rel 1e-9 (pytest.approx).
- `srm_1d/DEVNOTES.md`: API-log entry (user .ric with slag > 0 now correct; previously 1e6x too fast). Merge conflict with T10's entry resolved by keeping both.

## P3 — refresh-delay statistic
- `srm_1d/tools/ignition_diagnostics.py`: `ignition_refresh_delays(result, interval=None)` (post-hoc, read-only; interval defaults to `summary['burn_refresh_interval_steps']`; ValueError if < 1). Returns delays_s, n_ignitions, n_with_refresh, max_s, median_s, interval. Same definition as Codex's `_ignition_refresh_delays` (history index used as the step clock, so accurate to one step).
- `tests/test_ignition_diagnostics.py`: a synthetic exact-value test and a short real-run test (`_run_small`, interval 5).
- `docs/solver_verification/worklog/PA_scripts/refresh_delay_chunc.py`: reproduction on canonical Chunc (canonical knobs, BPNV, 100 cells, t_max 0.05, interval 21).
- Result (pre-rebase tree): 100 ignitions, **max 2.87475e-5 s, median 9.947e-6 s**; `max_burn_refresh_gap_s` = 1.438321e-4 s (30457 steps). Matches R1's canonical-BPNV numbers (2.875e-5 / 9.9e-6); Codex's 2.984e-5 was MTV config, same order. Rebased tree: max 2.874755e-5 s, median 9.948085e-6 s (30462 steps; gap unchanged 1.438321e-4 s).
- `worklog/T5.md` and `worklog/T7.md`: correction notes appended ("29.84 us" is this statistic, not the max refresh gap).

## P6 — B0 provenance
`probes/baseline_gate.py`: new `_print_provenance` prints git HEAD + dirty flag (tracked files), numpy/numba/python versions, and sha256 of each .ric the gate runs (resolved via the probe's own committed-copy logic). Output only; pass/fail logic and T10's re-recorded `BASELINE` constants untouched. Example: `git HEAD 71596f3... dirty=False`, numpy 2.2.6, numba 0.63.1; BALLSstick.ric 7fd189d62aabcdc5a5a788611d7b2c3cd94b91fbcc3e111f15e89423074cbd73; hasegawa_a.ric ccb7fccc972fd6074beff0b0099d1fdd2f731549d4628ea0eccbef97ca69185d; machbusterNew.ric a234a1dcf054dfb337e508360170e504228e9be85dfbd984940ea8326819df0e.

## Gate B0 (rebased tree, T10 constants)
Run 1 (cold cache, executed concurrently with the pytest run so both processes were compiling into the same `__pycache__`): steps exact on all three, but plateau rel 1.6e-10 (Chunc) and Hasegawa peak 7.0e-12 / t_peak 4.7e-12 -> reported FAIL (> 1e-12). Run 2 (warm cache): **PASS, every float rel 0.000e+00, steps exact**:
chunc peak 11.28243546264288 / t_peak 0.016286760473606336 / plateau 8.150933877646011 / 531883; balls peak 8.331331016878805 / trough 4.343222420159596 / mid 4.879629094113075 / 843274; haseg peak 6.108574454139173 / t_peak 2.391606230911318 / 2336140.
Interpretation: run 1 was compile-state noise from two processes compiling concurrently (fastmath machine code depends on compile state, per the gate docstring); no PA change touches the solver (adapter unit conversion, tools, tests, gate output). Not a trajectory change. Open item: the noise was larger (1.6e-10) than the 1e-15 quoted in the gate docstring when caches race; flagging for the reviewer.

## Deviations / open items
- Baseline suite was measured while I edited (tests collected before the edits; modules imported before the edits), so it is effectively the pre-edit tree.
- slag test uses the post-processing `compute_thrust_history`, not a full `run_simulation`: both use the identical `slag_coeff / P_MPa` expression.
