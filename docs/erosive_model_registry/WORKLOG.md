# Erosive model registry — WORKLOG

Review artifact for `docs/erosive_model_registry/{DESIGN,TASKS}.md`. One
entry per completed task, written for a reviewer who has not seen the
session.

---

## 2026-08-08 — T1: registry + kernel branch + threading (core)

**Scope:** DESIGN §1-§4, TASKS.md T1. Implemented exactly as specified;
no ambiguities hit that required stopping.

### Files changed

**`srm_1d/burn_rate.py`**
- Added `from dataclasses import dataclass` import.
- New Python-level (NOT `@njit`) registry block inserted right after the
  numba/no-numba shim: `ErosiveModel` frozen dataclass
  (`name`/`code`/`label`/`description`/`params`), `EROSIVE_MODELS` dict
  with `ma2020`→0, `mukunda_paul`→1, `saint_robert`→2 (descriptors per
  DESIGN §1), `resolve_erosive_model(name)` (dict lookup, `ValueError`
  listing sorted registered names on miss), `register_erosive_model(model)`
  (rejects duplicate `name` OR duplicate `code`, in that order, before
  inserting).
- `burn_rate_cell`: appended `erosive_model_code` as the new trailing
  parameter (after `erosive_mp_cap`). Two branch insertions, both
  additive — the existing code-0 (Ma) lines are untouched, not moved,
  not reformatted:
  - Immediately after `r0 = saint_robert_from_tabs(...)`, **before** the
    `Re_local < 100.0` guard: `if erosive_model_code == 2: return r0, 0.0`.
  - Immediately after `ebf = T_diff_gas / (rho_p * Cps * T_diff_prop)`
    (i.e. after the Re<100 / h0<1e-6 / T_diff_prop<1.0 guards have all
    passed), **before** the `alpha = ...` bisection setup line:
    `if erosive_model_code == 1:` computes `G = Re_local * mu_gas / D_hyd`,
    `eta_mp_value = mukunda_paul_eta(G, r0, D_hyd, mu_gas, rho_p)`, returns
    `(r0 * eta_mp_value, r0 * eta_mp_value - r0)` directly (no clamp — always
    ≥0 since `mukunda_paul_eta` floors at 1.0). Falls through to the
    existing bisection for code 0.
  - Docstring: added `erosive_model_code` parameter doc block; no other
    docstring text touched.
- `compute_burn_rates`: appended `erosive_model_code` (after
  `erosive_mp_cap`, before `N`) to both the signature and the internal
  `burn_rate_cell(...)` call; added the parameter's docstring entry.

**`srm_1d/simulation.py`**
- Import line: added `resolve_erosive_model` to the `from .burn_rate
  import (...)` block (converted to multi-line import to fit).
- `_run_time_loop` signature: appended `erosive_model_code,` immediately
  after the existing `erosive_mp_cap,` line in the "spike-reopen probe"
  argument group, with the comment
  `# int — registry code, 0=ma2020, 1=mukunda_paul, 2=saint_robert`
  (exact comment text specified in TASKS.md item 4).
- Inside the time loop, the `compute_burn_rates(...)` call (STEP 2: BURN
  RATES block) now passes `erosive_model_code` through in the matching
  slot.
- `run_simulation`: new kwarg `erosive_model='ma2020'` added right after
  `erosive_mp_cap=0.0,` in the signature, with an inline comment block
  documenting the three registry values and the roughness/kappa
  ignition-gate-vs-burn-rate distinction (DESIGN §2 scope note); matching
  `Parameters` docstring entry added after the existing `erosive_mp_cap`
  entry.
- Validation block: added right after the existing
  `igniter_axial_momentum_fraction` range check (before the
  `pyrogen_heat_flux` / `pyrogen_chamber.pyrogen...` access — i.e. before
  any sim-object attribute access):
  ```python
  _ero_code = resolve_erosive_model(erosive_model)
  if erosive_mp_cap > 0.0 and erosive_model != 'ma2020':
      raise ValueError(...)  # "ma2020-only option"
  ```
- `_run_time_loop(...)` call site: added `int(_ero_code),` immediately
  after the existing `float(erosive_mp_cap),` line (the `# spike-reopen
  probe` comment marks the group; matches the `bool(...)`/`float(...)`
  slot-ordering convention of the surrounding call).
- `summary` dict: added `'erosive_model': erosive_model,` (string
  provenance stamp per DESIGN §4), inserted after `'ambient_temperature'`
  and before `'diagnostic_disable_erosive'`, with a short comment.

**`tests/test_burn_rate.py`**
- Import block: added `mukunda_paul_eta`, `ErosiveModel`, `EROSIVE_MODELS`,
  `resolve_erosive_model`, `register_erosive_model` to the `from
  srm_1d.burn_rate import (...)` line.
- Refactored the Hasegawa-Propellant-1 fixture out of `TestBurnRateCell`
  into a new **non-collected** mixin `_HasegawaPropellant1Fixture` (no
  `Test` prefix, so pytest doesn't gather it) holding the same constants
  and the `_call(...)` helper — now with an `erosive_model_code=0` kwarg
  appended (default 0 = ma2020, preserving every existing call site's
  behavior unchanged). `TestBurnRateCell(_HasegawaPropellant1Fixture)`
  and three new test classes all subclass the mixin so the reference
  propellant is shared without re-running `TestBurnRateCell`'s own tests
  under each subclass (a plain `class Foo(TestBurnRateCell)` would have
  caused pytest to re-collect and re-execute all of `TestBurnRateCell`'s
  inherited test methods under `Foo` too — avoided by extracting the
  fixture instead of subclassing the test class directly). This is a
  **deviation from literal task wording** ("update the two direct call
  sites with the new trailing arg (`0`)") only in that the two direct
  call sites (`_call`'s `burn_rate_cell(...)` call and
  `TestComputeBurnRates`'s `compute_burn_rates(...)` call) both got the
  trailing `0` as specified — the fixture extraction is purely an
  additional refactor to avoid test-duplication, not a change to what
  gets asserted.
- Updated the two direct positional-arg call sites with the trailing `0`
  (default `erosive_model_code`):
  - `_HasegawaPropellant1Fixture._call`'s `burn_rate_cell(...)` call
    (now uses the `erosive_model_code` parameter, default 0).
  - `TestComputeBurnRates.test_non_burning_cells_zero`'s
    `compute_burn_rates(...)` call — literal `0` appended before `N`.
- New test classes added:
  - `TestErosiveModelRegistry` (registry unit tests): resolve known
    names → codes 0/1/2; descriptor sanity; unknown name raises
    `ValueError` matching `"unknown erosive_model"` and listing
    `"ma2020"`; `register_erosive_model` rejects a duplicate name
    (`match="name already registered"`) and a duplicate code
    (`match="code already registered"`, and confirms the new name was
    NOT inserted despite being novel).
  - `TestMukundaPaulModel(_HasegawaPropellant1Fixture)` — **G2**: at
    `P=5e6, Re=5e5, D=0.04` asserts `r_total == pytest.approx(r0 *
    mukunda_paul_eta(...), rel=1e-12)` (exact closed form, no bisection);
    at `Re=50000` (g≈11.4 < 35 threshold) asserts exact `r0`; at `Re=50`
    (< the shared Re<100 guard) asserts exact `r0`.
  - `TestSaintRobertModel(_HasegawaPropellant1Fixture)` — **G2b**: `(r0,
    0.0)` at `Re=50` (low) and `Re=5e6` (very high — plus a same-condition
    ma2020 sanity check confirming ma2020 is NOT bare r0 there, i.e. the
    test is actually exercising the erosive-off branch and not a
    degenerate condition).
  - `TestMPModelMatchesCappedMa(_HasegawaPropellant1Fixture)` — **G3**:
    at `P=5e6, Re=5e5, D=0.04` (the exact condition TASKS.md item 5(d)
    names), confirms the premise `eta_Ma(1.713) > eta_MP(1.621)` holds,
    then asserts `mukunda_paul` mode's `(r_total, r_erosive)` equals
    `ma2020` mode's with `erosive_mp_cap=1.0` **exactly** (`==`, not
    `approx` — both paths converge to `r0 * eta_mp` and the observed
    values were bit-identical in the interactive check, see Gate results).
  - `TestRunSimulationErosiveModelValidation` — TASKS.md item 6
    integration guard: `run_simulation(object(), object(), object(),
    object(), erosive_model='nope')` raises `ValueError` matching
    `"unknown erosive_model"`; same sentinel call with
    `erosive_model='mukunda_paul', erosive_mp_cap=1.0` raises matching
    `"ma2020-only"`; a third confirmatory test with
    `erosive_model='ma2020', erosive_mp_cap=1.0` (the ma2020-only
    validation must NOT fire) instead hits `AttributeError` on the
    sentinel's `.pyrogen` access next — proving the earlier `ValueError`
    genuinely did not trigger, and confirming the whole check fires
    **before** any geo/propellant/nozzle/pyrogen_chamber attribute
    access (no real motor construction, no time loop, per DESIGN §3 /
    TASKS item 6).

**No other files touched.** `run_from_ric` needed no changes (confirmed:
`**sim_overrides` already forwards arbitrary `run_simulation` kwargs
generically — verified by reading `srm_1d/openmotor_adapter.py` around
`run_from_ric`/`run_from_dict`).

### Numba cache

`srm_1d/__pycache__/` deleted before every test/probe run that touched
the edited `@njit` kernels (three times across the session, after each
batch of kernel edits and again before the final full-suite run).

### Gate results (actual numbers)

**G1 — default identity.** Canonical Chunc run against the **committed**
`motors/machbusterNew.ric` (extracted via `git show
HEAD:motors/machbusterNew.ric > C:/tmp/machbusterNew_committed.ric` —
the working-tree copy is modified with different transport per `git
status`, untouched), `roughness=32e-6, kappa=0.44, pyrogen='bpnv',
pyrogen_mass=None, T_ignition=756.0, P_cutoff=0.01e6, t_max=1.0`, all
other kwargs default (`erosive_model` defaults to `'ma2020'`,
`erosive_mp_cap` defaults to `0.0`):
```
erosive_model (summary) = 'ma2020'
P_peak = 11.4728 MPa @ t_peak = 16.45 ms
plateau (median 0.4-0.9s) = 8.2225 MPa
```
Matches the spec target (11.47 MPa @ ~16.5 ms, plateau 8.22 MPa)
exactly. **PASS.**

**G2 — mukunda_paul unit tests.** `TestMukundaPaulModel`, 3/3 passed:
closed-form exact match at high Re; exact `r0` at `g<35`; exact `r0` at
`Re<100`. Confirmed interactively before writing the test:
`P=5e6, Re=5e5, D=0.04` → `r0=4.92979e-3, eta_mp=1.62071,
r_mp=7.98975e-3 == r0*eta_mp` (exact). **PASS.**

**G2b — saint_robert unit tests.** `TestSaintRobertModel`, 2/2 passed:
`(r0, 0.0)` at `Re=50` and `Re=5e6`; interactively confirmed `sr_lo =
sr_hi = r0 = 4.92979e-3` at both ends, vs. ma2020's `0.008446934...` at
the same high-Re condition (confirms the branch is live, not a
no-op). **PASS.**

**G2c — saint_robert == ma2020+diagnostic_disable_erosive integration.**
Short Chunc run (`t_max=0.05`, canonical knobs, committed
`.ric`), once with `erosive_model='saint_robert'`, once with
`erosive_model='ma2020', diagnostic_disable_erosive=True`:
```
n_steps saint_robert = 30890, ma2020+diag = 30890
summary erosive_model (sr)    = 'saint_robert'
summary erosive_model (ma_dd) = 'ma2020'
max |dt| = 4.378699e-10 s
max |dP_head| = 4.230994e+00 Pa (4.230994e-06 MPa)
P_head range (sr): 9.915e+04 .. 4.418e+06 Pa
```
Max |ΔP_head| ≈ 4.23 Pa against a ~4.4 MPa peak (relative ≈1e-6) — float
tolerance, not model divergence. The two code paths are mathematically
equivalent (`saint_robert` short-circuits to `r0` before any flow-state
read; `ma2020+diagnostic_disable_erosive` runs the full bisection and
then the time loop's post-hoc `r_normal = r_total_new[i] -
r_erosive_new[i]` subtraction zeroes the erosive component) but not
bit-identical because they touch different intermediate floats — the
tiny divergence compounds over 30890 adaptive-CFL steps without ever
becoming physically meaningful. **PASS** (not written up as a
`pytest` test per TASKS.md item 7's "script or test, your choice" —
kept as the standalone script `C:/tmp/g2c_check.py`, not committed to
the repo; see Deviations).

**G3 — mukunda_paul == ma2020+cap1.0.**
`TestMPModelMatchesCappedMa`, 1/1 passed. At `P=5e6, Re=5e5, D=0.04`:
`eta_Ma(uncapped) = 1.7134 > eta_MP = 1.6207` (premise confirmed), then
`r_mp == r_ma_capped` and `r_ero_mp == r_ero_ma_capped` bit-exact
(`==`, no tolerance needed). **PASS.**

**`pytest tests/test_burn_rate.py -q`:** `34 passed in 3.13s`.

**`pytest tests/ -q`:** `2 failed, 423 passed in 396.81s`. The 2
failures are the pre-existing, out-of-scope ones named in the task
brief:
- `tests/test_ric_igniter.py::test_all_repo_motors_have_parseable_igniter`
  — `AssertionError: BALLSstick VANILLA.ric has no igniter block`
- `tests/test_ric_transport_migration.py::test_repo_motors_all_migrated_and_loadable`
  — `ValueError: Propellant 'Chase - Energy': frozen gas transport
  (mu/k/Cp) is not provided in the .ric file` (this is the propellant
  embedded inside the untracked `motors/BALLSstick VANILLA.ric`, hit
  first alphabetically by the repo-wide glob scan)

Both are scoped entirely to the new, untracked `motors/BALLSstick
VANILLA.ric` file (visible in `git status` as `??`), not touched by
this task, not caused by this task's edits (verified: isolated re-run
of just these two tests reproduces the identical failures/messages).
**PASS per TASKS.md's explicit "the 2 BALLSstick VANILLA.ric ...
failures are pre-existing" acceptance clause.**

### Deviations from the spec

1. **Test-fixture extraction (`_HasegawaPropellant1Fixture`).** TASKS.md
   item 5 only asked to "update the two direct call sites with the new
   trailing arg (0)" and add new test classes. A naive implementation
   would have had the new test classes subclass `TestBurnRateCell`
   directly to reuse its constants/`_call` helper, but that causes
   pytest to re-collect and re-run every one of `TestBurnRateCell`'s own
   `test_*` methods again under each subclass — redundant test execution
   with no additional coverage. Extracted the shared state into a
   non-`Test`-prefixed mixin instead so it isn't collected on its own.
   Net effect on what's asserted: none — purely a collection-hygiene
   fix. Flagging it since it's a structural change beyond the literal
   task wording.
2. **G2c not committed as a repo test.** TASKS.md item 7 says "script or
   test, your choice — if a test, keep it short." Given the run takes
   ~15-20s (a short Chunc sim, not the full G1 1-second run) and the
   project's existing test suite already runs in ~6.5 minutes, I chose
   the "script" option and did not add it to `tests/`. The script lives
   at `C:/tmp/g2c_check.py` (scratch location per the task brief's
   temp-file guidance) — **it will not survive this session** and is
   not part of the diff. If a permanent regression guard is wanted here,
   that's an open item (see below).
3. **`register_erosive_model` duplicate-check order.** DESIGN/TASKS don't
   specify which duplicate check (name vs. code) runs first when a
   registration violates both. Implemented name-check first, code-check
   second (arbitrary but deterministic); noted here since it's an
   unspecified but observable behavior.

### Open items

- G2c's integration check exists only as an ad hoc script in `C:/tmp/`
  for this session — not a committed regression test. Future work (or a
  T1 follow-up) could promote a short version of it into `tests/` if
  standing coverage of the saint_robert/diagnostic-disable-erosive
  equivalence is wanted.
- T2 (cross-motor validation, model axis in
  `docs/v0_7_4/probes/mp_cap_crossmotor.py`), T3 (GUI leg in
  `srm1d_plugin.py`), and T4 (docs: ARCHITECTURE/DEVNOTES/README/
  IGNITION_SPIKE_REOPENED) are unstarted — out of scope for T1 per
  TASKS.md's stated review order (T1 → T2 → (T3, T4 parallel)).
- Did not touch `motors/*.ric` (hard rule) and did not commit/push (hard
  rule) — working tree still has all T1 changes unstaged at the end of
  this session, plus the pre-existing unrelated modifications noted in
  the initial `git status` (docs/v0_7_4/*, motors/*.ric transport edits,
  untracked reference PDFs and new motor files) which this session did
  not create and did not touch.

---

## 2026-08-08 — T2: cross-motor validation (model axis)

**Scope:** DESIGN §G4, TASKS.md T2. Extended
`docs/v0_7_4/probes/mp_cap_crossmotor.py` to the model axis and ran the
full 15-row sweep (3 motors × 5 `(erosive_model, erosive_mp_cap)` rows).
No ambiguities hit that required stopping.

### File changed: `docs/v0_7_4/probes/mp_cap_crossmotor.py`

The file was untracked (`??` in `git status`, added by the T1 session as
a new probe) so there is no `git diff` against a prior committed version;
described here against the pre-T2 content read at session start (the
version quoted in TASKS.md's own description, 94 lines).

- **Module docstring**: rewritten to a 3-part structure (was 2-part):
  kept (1) the differential test and (2) the cap-validation description
  verbatim in substance, added (3) describing the model-axis sweep
  (`MODEL_AXIS`) and the `saint_robert`-is-recorded-not-gated framing
  from DESIGN §2/G4. The caveat paragraph was extended to describe the
  new `_committed_ric_path()` helper (auto-detects working-tree
  modification via `git status --porcelain` and extracts the committed
  `HEAD` copy — no more manual `git show` step) instead of just warning
  about the working-tree state. Added a runtime note (measured, not the
  stale 2-8 min/row estimate — see Deviations).
- **New import**: `subprocess` (for the git-status/git-show plumbing).
- **New constants**: `TMP_DIR = Path("C:/tmp")`; `CAP_AXIS = (0.0, 1.0,
  1.1)` (kept for callers wanting only the old cap-only sub-axis);
  `MODEL_AXIS` — the 5-tuple of `(erosive_model, erosive_mp_cap)` pairs
  from TASKS.md: `('ma2020', 0.0)`, `('ma2020', 1.0)`, `('ma2020', 1.1)`,
  `('mukunda_paul', 0.0)`, `('saint_robert', 0.0)` (kept the old cap=1.0
  row per TASKS.md's "keeping the old cap=1.0 row is fine — 5 rows").
- **New function `_committed_ric_path(name)`**: runs `git status
  --porcelain -- motors/<name>` from `ROOT`; if the output is non-empty
  (file modified/untracked), runs `git show HEAD:motors/<name>` and
  writes it to `C:/tmp/<stem>_committed.ric`, prints
  `[ric-source] motors/<name>: MODIFIED in working tree -> using
  committed HEAD extract at <path>`, and returns that path; otherwise
  prints `[ric-source] motors/<name>: unmodified -> using working-tree
  copy` and returns `MOTORS / name` unchanged. Applied uniformly to all
  three motors (including `hasegawa_a.ric`, which is unmodified — it
  routes through the same helper for a consistent header note, per
  TASKS.md's example wording, and correctly resolves to the plain
  `motors/` path with the "unmodified" note, confirmed in the run
  output below).
- **`run()`**: signature gained `erosive_model` as a new leading
  parameter (before `cap`); passes both through to `run_from_ric(...,
  erosive_model=erosive_model, erosive_mp_cap=cap, ...)`; now also
  captures and returns `res['summary']` as a third return value (was
  2-tuple `(t, P)`, now 3-tuple `(t, P, summary)`) so callers can read
  the `erosive_model` provenance stamp (DESIGN §4) and `termination` /
  `steps` fields.
- **`chunc_row` / `balls_row` / `haseg_row`**: each gained
  `erosive_model='ma2020'` as a new leading parameter (default preserves
  the old call shape for any external caller still passing only `cap`
  positionally... actually note: call shape changed — see Deviations,
  item 1), call `_committed_ric_path(...)` internally instead of the
  caller passing `MOTORS / "..."` directly (was hardcoded `MOTORS /
  "machbusterNew.ric"` etc. inline in the `run()` call), and each now
  returns 3 additional dict keys: `erosive_model_stamp` (from
  `summary['erosive_model']`, the provenance check), `termination`
  (from `summary['termination']`, e.g. `'t_max reached'` — anomaly
  flag surface), `n_steps` (from `summary['steps']`).
- **`__main__`**: the old `for cap in (0.0, 1.0, 1.1):` loop replaced
  with `for erosive_model, cap in MODEL_AXIS:`; each of the three
  per-motor print blocks gained the `erosive_model` label, the
  `erosive_model_stamp` (printed in the row header, so a stamp mismatch
  is visible without parsing), and the `[termination=..., n_steps=...]`
  suffix. Metrics printed (peak/plateau/ratio for Chunc; peak/pk_tr/
  pk_mid for BALLSstick; global peak/t_peak for Hasegawa A) are
  otherwise unchanged from the pre-T2 script.

**No other files touched** by the probe edit itself. `srm_1d/burn_rate.py`
/ `srm_1d/simulation.py` / `tests/test_burn_rate.py` are T1's changes,
untouched by T2.

### Execution harness (not committed — `C:/tmp` scratch only)

Per the hard-rule "never touch files outside the probe script + WORKLOG +
`C:/tmp` scratch" and the runtime-chunking instruction, the actual 15-row
sweep was **not** run via `python mp_cap_crossmotor.py` directly (which
would be one ~unbounded-length process). Instead:
- `C:/tmp/t2_run_row.py`: a driver that imports `mp_cap_crossmotor` as a
  module, takes `(motor, erosive_model, cap)` on argv, calls the
  matching `*_row()` function, and appends the resulting dict (plus
  `motor`/`erosive_model`/`cap`/`wall_s`) as one JSON line to
  `C:/tmp/t2_results.jsonl`.
- Invoked once per row from Bash, batched 1 (first row, doubling as a
  timing probe) then 4-per-motor (the remaining 12 rows, one Bash call
  per motor covering its other 4 `MODEL_AXIS` rows) — well inside the
  ~9 min/call budget (see Deviations: actual runtimes were far below
  the doc's estimate, so no run came close to the cap).
- `srm_1d/__pycache__/` was cleared once at the start of the T2 session
  (before the first row) even though T2 makes no `@njit` edits, as a
  precaution against any stale cache from the T1 session's last edit;
  not needed again since no kernel code changed during T2.
- Final table below was assembled by reading back `C:/tmp/t2_results.jsonl`
  (15 lines, one per row, verified by count) plus one extra experimental-
  reference computation (Chunc `pl_exp` / peak / ratio from
  `static_fire_data/thomas_chunc_firing.csv`).

### `.ric` provenance (per-run header note, confirmed in every run's stdout)

```
[ric-source] motors/machbusterNew.ric: MODIFIED in working tree -> using committed HEAD extract at C:\tmp\machbusterNew_committed.ric
[ric-source] motors/BALLSstick.ric: MODIFIED in working tree -> using committed HEAD extract at C:\tmp\BALLSstick_committed.ric
[ric-source] motors/hasegawa_a.ric: unmodified -> using working-tree copy
```
Matches `git status --porcelain` at session start (`machbusterNew.ric`
and `BALLSstick.ric` both `M`; `hasegawa_a.ric` absent from the status
list). All 15 runs used these three sources consistently (each row
re-invoked `_committed_ric_path`, re-extracting the same committed copy
every time — confirmed identical extract, cheap no-op cost).

### Full results — all 15 rows (actual numbers)

Chunc metrics in MPa (`plateau` = median P_head, t∈[0.4,0.9]s;
`ratio` = peak/plateau). BALLSstick metrics in MPa (`peak`/`t_peak` are
`argmax` restricted to the t<0.3s window — built for the erosive-spike
shape, see the saint_robert note below for what happens when that
assumption doesn't hold; `trough` = min P_head over t∈(t_peak, 1.2)s,
unrestricted; `mid` = P_head nearest t=1.03s). Hasegawa A: global peak
over the full t_max=10s trace. All rows: `termination='t_max reached'` for all 15 — no
early terminations, no ignition-failure stalls, in any model arm on any
motor (including all three `saint_robert` rows).

| erosive_model | cap | Chunc peak / t_peak(ms) / plateau / ratio | BALLSstick peak / trough / mid / pk_tr / pk_mid | Hasegawa A peak / t_peak(s) | provenance stamp ok |
|---|---|---|---|---|---|
| ma2020 | 0.0 | 11.4728 / 16.45 / 8.2225 / 1.3953 | 8.4235 / 4.3667 / 4.9134 / 1.9290 / 1.7144 | 6.1439 / 2.382 | yes (3/3) |
| ma2020 | 1.0 | 7.8410 / 1000.0* / 7.2999 / 1.0741 | 6.0130 / 3.9613 / 4.3961 / 1.5179 / 1.3678 | 6.0305 / 2.564 | yes (3/3) |
| ma2020 | 1.1 | 9.3566 / 14.50 / 8.2405 / 1.1354 | 7.4347 / 4.3776 / 4.9054 / 1.6984 / 1.5156 | 6.1589 / 2.396 | yes (3/3) |
| mukunda_paul | 0.0 | 7.8410 / 1000.0* / 7.2999 / 1.0741 | 6.0130 / 3.9613 / 4.3961 / 1.5179 / 1.3678 | 6.0307 / 2.564 | yes (3/3) |
| saint_robert | 0.0 | 7.0174 / 1000.0* / 5.9966 / 1.1702 | 3.0248 / 2.7540 / 3.2503 / 1.0983 / 0.9306 | 6.3635 / 4.039 | yes (3/3) |

\* Chunc `t_peak` at (or numerically indistinguishable from) `t_max=1.0s`
for the `ma2020 cap=1.0`, `mukunda_paul`, and `saint_robert` rows means
the trace is still rising (or effectively flat) at cutoff — there is no
early spike in these arms, so `argmax` lands on the last sample rather
than a genuine local maximum. This is expected shape (not a bug): the
whole point of the MP/SR arms on Chunc is that they suppress the
Ma-only early spike (see the ratio column: 1.07-1.17 vs Ma's 1.395),
so the "peak" for these arms is really "value at the end of the
recorded window," not a startup transient. Flagged here so the table
isn't misread as three coincident spikes at t=16ms.

Raw JSONL (15 lines, `C:/tmp/t2_results.jsonl`, not committed):
```
{"peak": 11.472831280605046, "t_peak": 0.016454109339050758, "plateau": 8.222465813350324, "ratio": 1.395303007764082, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 531651, "motor": "chunc", "erosive_model": "ma2020", "cap": 0.0, "wall_s": 10.4}
{"peak": 8.42348927652058, "trough": 4.366702806070186, "mid": 4.913378178011638, "pk_tr": 1.9290273807530538, "pk_mid": 1.7143987235131626, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 842922, "motor": "balls", "erosive_model": "ma2020", "cap": 0.0, "wall_s": 25.6}
{"peak": 6.143912554793205, "t_peak": 2.3816074267208034, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 2335770, "motor": "haseg", "erosive_model": "ma2020", "cap": 0.0, "wall_s": 40.9}
{"peak": 7.8410307141508975, "t_peak": 0.9999987168860677, "plateau": 7.299913238085178, "ratio": 1.074126562661402, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 540585, "motor": "chunc", "erosive_model": "ma2020", "cap": 1.0, "wall_s": 10.6}
{"peak": 9.356565481628524, "t_peak": 0.014503375437932642, "plateau": 8.24050103735905, "ratio": 1.135436478826918, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 534035, "motor": "chunc", "erosive_model": "ma2020", "cap": 1.1, "wall_s": 10.7}
{"peak": 7.8410307141508975, "t_peak": 0.9999987168860677, "plateau": 7.299913238085178, "ratio": 1.074126562661402, "erosive_model_stamp": "mukunda_paul", "termination": "t_max reached", "n_steps": 540585, "motor": "chunc", "erosive_model": "mukunda_paul", "cap": 0.0, "wall_s": 9.2}
{"peak": 7.017426352358137, "t_peak": 0.9999983820583602, "plateau": 5.996614436020544, "ratio": 1.1702313742577422, "erosive_model_stamp": "saint_robert", "termination": "t_max reached", "n_steps": 558124, "motor": "chunc", "erosive_model": "saint_robert", "cap": 0.0, "wall_s": 9.1}
{"peak": 6.013008306594397, "trough": 3.9613394974074554, "mid": 4.396106884883187, "pk_tr": 1.5179229931011164, "pk_mid": 1.3678030275540434, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 855898, "motor": "balls", "erosive_model": "ma2020", "cap": 1.0, "wall_s": 25.8}
{"peak": 7.4346623109211345, "trough": 4.377554697627317, "mid": 4.905366243846018, "pk_tr": 1.6983596606915738, "pk_mid": 1.5156181906393271, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 845064, "motor": "balls", "erosive_model": "ma2020", "cap": 1.1, "wall_s": 25.6}
{"peak": 6.013008306594397, "trough": 3.9613394974074554, "mid": 4.396106884883187, "pk_tr": 1.5179229931011164, "pk_mid": 1.3678030275540434, "erosive_model_stamp": "mukunda_paul", "termination": "t_max reached", "n_steps": 855898, "motor": "balls", "erosive_model": "mukunda_paul", "cap": 0.0, "wall_s": 23.6}
{"peak": 3.0247741780519344, "trough": 2.753989745082887, "mid": 3.2502615099789764, "pk_tr": 1.0983244158597611, "pk_mid": 0.9306248647271153, "erosive_model_stamp": "saint_robert", "termination": "t_max reached", "n_steps": 908437, "motor": "balls", "erosive_model": "saint_robert", "cap": 0.0, "wall_s": 24.4}
{"peak": 6.030535817720915, "t_peak": 2.563807518243538, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 2344248, "motor": "haseg", "erosive_model": "ma2020", "cap": 1.0, "wall_s": 41.6}
{"peak": 6.158893638141223, "t_peak": 2.395829759084376, "erosive_model_stamp": "ma2020", "termination": "t_max reached", "n_steps": 2336505, "motor": "haseg", "erosive_model": "ma2020", "cap": 1.1, "wall_s": 40.8}
{"peak": 6.030667408394668, "t_peak": 2.5637369255059004, "erosive_model_stamp": "mukunda_paul", "termination": "t_max reached", "n_steps": 2344247, "motor": "haseg", "erosive_model": "mukunda_paul", "cap": 0.0, "wall_s": 37.8}
{"peak": 6.363499432521944, "t_peak": 4.039247872000466, "erosive_model_stamp": "saint_robert", "termination": "t_max reached", "n_steps": 2401665, "motor": "haseg", "erosive_model": "saint_robert", "cap": 0.0, "wall_s": 36.7}
```
Provenance check: `erosive_model_stamp == erosive_model` for all 15
rows (scripted check, 0 mismatches) — the `summary['erosive_model']`
stamp from DESIGN §4 is confirmed live end-to-end through
`run_from_ric`.

### G1 re-check (Chunc, `ma2020`/cap=0.0` row)

Target: P_peak ≈ 11.47 MPa @ ~16.5 ms, plateau ≈ 8.22 MPa.
Actual: **P_peak = 11.4728 MPa @ t_peak = 16.45 ms, plateau = 8.2225 MPa.**
Exact match to T1's own G1 numbers (same committed `.ric`, same canonical
knobs). **PASS.**

### G4 gate verdicts (tolerance ±0.03 on ratios, ±0.15 MPa on peaks)

| motor | metric | model/cap | expected | actual | delta | verdict |
|---|---|---|---|---|---|---|
| Chunc | ratio | ma2020, cap=0.0 | ≈1.395 | 1.3953 | +0.0003 | **PASS** |
| Chunc | ratio | ma2020, cap=1.1 | ≈1.135 | 1.1354 | +0.0004 | **PASS** |
| Chunc | ratio | mukunda_paul | ≈1.07 | 1.0741 | +0.0041 | **PASS** |
| BALLSstick | pk/trough | ma2020, cap=0.0 | ≈1.93 | 1.9290 | −0.0010 | **PASS** |
| BALLSstick | pk/trough | ma2020, cap=1.1 | ≈1.70 | 1.6984 | −0.0016 | **PASS** |
| BALLSstick | pk/trough | mukunda_paul | ≈1.52 | 1.5179 | −0.0021 | **PASS** |
| Hasegawa A | global peak (MPa) | ma2020, cap=0.0 | ≈6.14 | 6.1439 | +0.0039 | **PASS** |
| Hasegawa A | global peak (MPa) | ma2020, cap=1.1 | ≈6.16 | 6.1589 | −0.0011 | **PASS** |
| Hasegawa A | global peak (MPa) | mukunda_paul | ≈6.03 | 6.0307 | +0.0007 | **PASS** |

**All 9 gated rows PASS**, every delta at least an order of magnitude
inside its tolerance window.

### cap=1.0 vs mukunda_paul cross-check (DESIGN's "flag if they don't match")

- **Chunc**: bit-identical across all fields (`peak`, `t_peak`, `plateau`,
  `ratio` — every digit matches). 
- **BALLSstick**: bit-identical across all fields.
- **Hasegawa A**: near-identical but **not** bit-identical —
  `peak` 6.030535817720915 (cap=1.0) vs 6.030667408394668
  (`mukunda_paul`), Δ = 1.32e-4 MPa (0.0022%); `t_peak` 2.563807518
  vs 2.563736926 s, Δ = 7.06e-5 s. This is the DESIGN-documented,
  expected case: `erosive_mp_cap=1.0` takes `min(Ma, MP-capped-Ma)` per
  cell/step, so it equals pure `mukunda_paul` mode ONLY where
  `eta_Ma ≥ eta_MP` everywhere active; DESIGN's G4 note explicitly
  anticipates "small deviations where MP > Ma are expected and worth
  recording." The Hasegawa A divergence (absent on Chunc/BALLSstick)
  means somewhere in Hasegawa A's burn `eta_Ma < eta_MP` transiently, so
  the capped-Ma branch and the pure-MP closed form pick different
  values in that cell/step window. Flagged as expected, not a defect —
  magnitude is negligible (< 0.003% of peak P).

### saint_robert decomposition (RECORDED, not gated — DESIGN §2/G4)

**Chunc plateau_SR — the headline number: 5.9966 MPa** (median P_head,
t∈[0.4,0.9]s, pure Saint-Robert / zero erosive contribution).

| source | plateau (MPa) | plateau_SR as % of it |
|---|---|---|
| experimental (`thomas_chunc_firing.csv`) | 8.6611 | 69.2% |
| ma2020 (cap=0.0) | 8.2225 | 72.9% |
| mukunda_paul | 7.2999 | 82.2% |
| saint_robert | 5.9966 | 100% |

Erosive-closure plateau lift over the Saint-Robert floor:
ma2020 adds **+2.226 MPa (+37.1%)** over SR to reach 8.2225 (94.9% of
experimental 8.6611); mukunda_paul adds **+1.303 MPa (+21.7%)** over SR
to reach 7.2999 (84.3% of experimental). Reading for the spike/plateau
decoupling question (REOPENED §10.5/§11.4): a large share of Chunc's
plateau match under `ma2020` is genuinely erosive-closure-supplied, not
a Saint-Robert `a`/`n` calibration artifact — pure Saint-Robert alone
undershoots the experimental plateau by ~31% (5.9966 vs 8.6611). This
doesn't by itself say whether Ma's *peak* overshoot and its *plateau*
adequacy are separable design choices, but it quantifies how much
plateau work would need to move to Saint-Robert `a`/`n` recalibration if
a future closure suppressed the peak further — that gap is ~2.2 MPa
(ma2020 case) or ~1.3 MPa (MP case) at the plateau, not the full
2.66 MPa gap to experimental.

**BALLSstick saint_robert**: peak 3.0248 MPa, trough 2.7540 MPa, mid
(t≈1.03s) 3.2503 MPa → pk/tr = 1.0983 (vs ma2020's 1.9290 — the erosive
closure is responsible for essentially all of BALLSstick's startup
spike shape), pk/mid = **0.9306 (< 1)** — i.e. under saint_robert the
value at t≈1.03s exceeds the local (t<0.3s window) peak. This is not an
error: absent the erosive-driven early transient, BALLSstick's P_head
rises smoothly through the whole window with no distinct early local
max, so the `<0.3s` peak-search window used by `balls_row()` (built for
the erosive-spike shape) finds a smaller local max than the later,
still-climbing trace at 1.03s. Flagged as an expected consequence of the
metric definition, not a simulation anomaly.

**Hasegawa A saint_robert**: global peak **6.3635 MPa @ t=4.039s** —
higher in absolute P than any of the three erosive-active arms
(6.03-6.16 MPa, all peaking at t≈2.4-2.6s) but nearly 1.5-1.7s later.
Consistent with the motor's documented "late progressive" real behavior
(script comment: exp ~6.5 MPa, late progressive): Hasegawa A's global
peak is a geometric/burnback-area feature, not an erosive-startup
feature — removing the erosive closure doesn't remove the peak, it
shifts it later and (slightly) taller, because the early-burn erosive
augmentation in the ma2020/MP arms is apparently pulling mass flow
forward in time rather than adding to the eventual area-driven maximum.

### Anomalies / deviations flagged

1. **Runtimes far below the doc's estimate.** TASKS.md's header estimated
   Chunc ~2-4 min, BALLSstick ~3-6 min, Hasegawa A ~4-8 min per run
   (15 rows ≈ 45-75 min total). Actual: Chunc ≈9-11s, BALLSstick
   ≈24-26s, Hasegawa A ≈37-42s per row; **all 15 rows completed in
   under 6 minutes of total wall time**, comfortably inside a single
   Bash call even before batching. Batched into 3 per-motor calls (after
   1 solo first-row call used to confirm the G1 numbers and get a timing
   read) purely out of caution per the task's explicit chunking
   instruction — not because any call was close to the 9-minute or
   10-minute (tool max) ceiling. Left as an open question for whoever
   maintains the runtime estimate in TASKS.md/DESIGN.md — this machine/
   session measured much faster than documented, cause not
   investigated (could be a faster machine, a warm Numba disk cache
   from T1's own runs, or the original estimate being conservative/
   stale from an earlier, slower architecture).
2. **`(ma2020, cap=1.0)` vs `mukunda_paul` not bit-identical on Hasegawa
   A** — see the dedicated subsection above. Expected per DESIGN, not a
   gate failure (this row is not gated).
3. **BALLSstick `saint_robert` `pk_mid < 1`** — see the saint_robert
   decomposition subsection above. Metric-definition artifact of a
   window built for the erosive-spike shape being applied to a
   non-spiking trace; not a simulation defect.
4. No ignition failures, no early terminations, no anomalous
   `termination_code` in any of the 15 runs (all `'t_max reached'`) —
   explicitly checked because DESIGN/TASKS flagged "a saint_robert run
   failing to ignite or terminating early" as a thing to watch for. All
   three `saint_robert` rows ignited and ran to `t_max` normally.
5. **No gate misses.** All 9 gated G4 rows pass; the two consistency
   checks above (cap=1.0 vs MP; provenance stamps) both come back clean
   except the one documented-as-expected Hasegawa A micro-divergence.

### Deviations from the spec

1. **`chunc_row`/`balls_row`/`haseg_row` call-signature change.** TASKS.md
   only says to extend `__main__`'s sweep; it doesn't specify whether the
   per-motor helper functions gain a new leading parameter or whether
   `erosive_model` is threaded some other way (e.g. a wrapper function
   outside them). Chose to add `erosive_model` as each function's new
   *first* positional parameter (with `cap` shifting to second, both
   defaulted so existing zero-arg / cap-only call patterns aren't
   silently broken at the Python level, but any caller passing `cap`
   positionally as the first arg — e.g. old code calling
   `chunc_row(1.0)` — would now bind `1.0` to `erosive_model` and error
   loudly at `resolve_erosive_model()` rather than silently misrun. No
   such external caller exists in the repo (grepped: only `__main__` in
   this same file called these functions, and it call-site was rewritten
   in the same edit), so this is theoretical, not observed, breakage —
   flagged per the project's "fail loudly on ambiguity" convention.
2. **`_committed_ric_path` re-extracts on every call rather than caching
   within a process.** Each of the 15 row invocations is a separate
   Python process (per the execution-harness design above), so
   in-process caching wouldn't have helped across rows anyway; within a
   single row it's called once. Not a real deviation, noting the design
   choice for a reviewer wondering why there's no memoization.
3. **Execution driver (`C:/tmp/t2_run_row.py`) not part of the probe
   script.** TASKS.md's acceptance is "paste the output table into the
   WORKLOG entry" — it doesn't require the chunking mechanism itself to
   be committed, and the hard rules restrict touched files to "the probe
   script + WORKLOG + `C:/tmp` scratch." Kept the row-runner as a
   throwaway `C:/tmp` script (per the mechanism the T1 worklog also used
   for its own uncommitted G2c script) rather than adding CLI-argument
   plumbing to the committed probe — keeps the reviewed diff focused on
   exactly what TASKS.md asked for (the `__main__` sweep extension).

### Open items

- T2's runtime finding (item 1 above) may be worth a one-line update to
  TASKS.md's stated 45-75 min estimate for future task-planning, but
  that's a docs edit outside this task's touched-files scope (`docs/
  erosive_model_registry/TASKS.md` is not in "the probe script +
  WORKLOG + C:/tmp scratch") — left to a maintainer to decide, not
  edited here.
- T3 (GUI leg, `srm1d_plugin.py`) and T4 (docs: ARCHITECTURE/DEVNOTES/
  README/IGNITION_SPIKE_REOPENED) are unstarted — out of scope for T2
  per TASKS.md's stated review order (T1 → T2 → (T3, T4 parallel)).
- The Chunc/BALLSstick/Hasegawa `_committed_ric_path` extracts live at
  `C:/tmp/machbusterNew_committed.ric` and `C:/tmp/BALLSstick_committed.ric`
  for this session only — scratch, not committed, will not persist.
- Did not touch `motors/*.ric` (hard rule) and did not commit/push (hard
  rule).

---

## 2026-08-08 — T2 review note (main-session reviewer)

**Verdict: APPROVED.** Diff reviewed (probe extension is clean, additive,
self-documenting via `_committed_ric_path`); all 9 G4 gates verified; the
provenance-stamp sweep and the cap1.0≡mukunda_paul cross-check both confirm
the registry end-to-end.

**Runtime open item — RESOLVED.** The 45-75 min estimate assumed cold-start
Numba JIT per run (the earlier investigation sessions cleared
`srm_1d/__pycache__/` before each measurement batch, per the DEVNOTES
stale-cache discipline after kernel edits, forcing ~2-4 min of recompilation
every time). T2 made no kernel edits after T1's final compile, so all 15
rows hit the warm persistent `cache=True` cache and ran at the documented
~45-90k steps/s (10-40 s/row). Both numbers are "correct"; TASKS.md's
estimate is the cold-cache bound. No action needed.

**Interpretation correction (the decomposition — important).** The T2
entry's reading "not a Saint-Robert a/n calibration artifact" is TOO
STRONG. Converting the plateau decomposition to burn-rate space
(P ∝ r^(1/(1-n)), n=0.4 → r-lift = P-lift^0.6):

| | plateau P-lift over SR (8.66/8.22/7.30 ÷ 6.00) | implied avg rate-lift η̄ |
|---|---|---|
| required (experimental) | 1.444 | **≈1.25** |
| ma2020 supplies | 1.371 | ≈1.21 |
| mukunda_paul supplies | 1.217 | ≈1.13 |
| MP's +10% band edge (1.1·η̄_MP) | — | ≈1.24 |

The experimental requirement (≈1.25) sits essentially AT the +10% edge of
MP's own accuracy band (≈1.24), with Ma comfortably inside the band. So the
data admits TWO consistent worlds: (a) current a/n correct → true plateau
erosive lift sits at MP's band edge, where Ma happens to sit → Ma right at
plateau; or (b) Chunc's characterized a/n is ~10% low → MP central + a
corrected a/n reproduces the plateau AND Ma then over-predicts everywhere.
Chunc alone CANNOT break this degeneracy — a ~10% a/n error and a
band-edge erosive lift are observationally equivalent on one motor.
**Discriminator: cross-propellant.** Run the model axis on Zerox
(different propellant, mid-L/D, has static-fire data; igniter calibration
stale but the plateau decomposition is mid-burn and igniter-insensitive).
If MP's plateau deficit recurs across propellants → closure-systematic
(favors Ma-at-plateau). If Chunc-only → a/n debt (favors full MP + a/n
refinement through the legitimate propellant-data channel — the complete
spike/plateau decoupling). Hasegawa A hints the closure supplies genuinely
needed physics (its SR peak timing is grossly late: 4.04 s vs ~2.4 s
experimental/Ma) but its g-regime cannot separate Ma from MP.

**New finding worth promoting (Chunc saint_robert row):** t_peak = t_max —
the pure-Saint-Robert Chunc trace has NO startup spike at all (argmax is
the final sample of a still-rising trace). This is the cleanest, first-class
confirmation yet that the ENTIRE Chunc startup over-pressure enters through
the erosive channel (upgrades the old `diagnostic_disable_erosive`
observation in REOPENED §1 to a registry-mode result). It also kills the
hypothesis that a non-erosive fill overshoot contributes a floor: there is
none.

---

## 2026-08-08 — Zerox discriminator (cross-propellant plateau decomposition)

**Scope:** run the erosive-model registry axis on Zerox (`motors/zerox.ric`,
n=0.456, mid L/D, `static_fire_data/zerox_data.csv`) to try to break the T2
review note's Chunc-only degeneracy — "Ma right at plateau" vs "Chunc a/n
~10% low + MP right" are observationally identical on Chunc alone. Also
measures whether Zerox is usable as a validation motor at all (known-stale
igniter-era calibration).

### File added: `docs/v0_7_4/probes/zerox_model_axis.py`

New probe, mirrors `mp_cap_crossmotor.py`'s conventions: same `CANON` dict
(`roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None,
T_ignition=756.0, print_interval=0.0, verbose=False`), a duplicated
`_committed_ric_path()` helper (self-contained probe, not imported), and the
same git-status-driven provenance note printed per run. Routes
`motors/zerox.ric` through the helper even though it is unmodified in the
working tree (confirmed via `git status --porcelain` — resolves to the plain
`motors/` path), purely for a uniform provenance note across every run.
`motors/machbusterNew.ric` (used for the Chunc cross-check rows) IS modified
in the working tree (still `transportVariant: effective`, the anomaly T2
flagged) so its rows are always served from the `git show HEAD:...` extract
at `C:/tmp/zerox_axis/machbusterNew_committed.ric`.

4-row model axis (`ROWS`): `('ma2020', 0.0)`, `('ma2020', 1.1)`,
`('mukunda_paul', 0.0)`, `('saint_robert', 0.0)` — a trimmed sub-axis of
`mp_cap_crossmotor.py`'s 5-row `MODEL_AXIS` (drops `(ma2020, cap=1.0)`,
already shown bit-identical to `mukunda_paul` on Chunc/BALLSstick in T2).
Provides `zerox_row`/`chunc_row` (each saves its full `(t, P_head_MPa)` trace
as `.npy` under `C:/tmp/zerox_axis/` and returns peak/ref/ratio/duration
metrics + `erosive_model` provenance stamp + termination + n_steps),
`zerox_metrics`/`chunc_metrics`, `compute_verdict`/`print_verdict` (task
brief part C), and four plotting functions (part D) writing to
`artifacts/erosive_model_axis/`. A `--row MOTOR MODEL CAP` CLI flag runs one
row and appends it to a `C:/tmp/zerox_axis/results.jsonl` accumulator (used
for the actual chunked execution below); `--plots-only` regenerates the 4
PNGs from cached `.npy` without re-running sims; no args runs the full
8-row sweep + verdict + plots end to end (the committed reference path).

### Execution harness (not committed — `C:/tmp` scratch only)

Row-by-row via the script's own `--row` flag (no separate driver needed for
the first two rows), then a throwaway 2-row chained driver
(`C:/tmp/zerox_axis/run_remaining.py`, not committed, same pattern as T2's
`t2_run_row.py`) for the last two Zerox rows once the first row's timing
made clear that a single 8-row run would risk the wall-clock budget. Actual
order run: Zerox `(ma2020, 0.0)` → Zerox `(saint_robert, 0.0)` → Zerox
`(mukunda_paul, 0.0)` + Zerox `(ma2020, 1.1)` (chained, one process) → Chunc
all 4 rows (one Bash call, sequential, fast) → 4 plots. Priority order
`ma2020 > saint_robert > mukunda_paul > cap1.1` (set by the coordinator
mid-run) was respected.

### `.ric` provenance (confirmed in every run's stdout)

```
[ric-source] motors/zerox.ric: unmodified -> using working-tree copy
[ric-source] motors/machbusterNew.ric: MODIFIED in working tree -> using committed HEAD extract at C:\tmp\zerox_axis\machbusterNew_committed.ric
```

### Experimental comparators — cross-checked in-script (not just quoted)

Zerox (`zerox_data.csv`, reversed-column recipe `time_col=1, pressure_col=0`,
58 pt, t=[0, 8.113] s):

| metric | task-brief value | in-script recomputation |
|---|---|---|
| peak | 3.994 MPa @ 0.43 s | **3.9938 MPa @ 0.43 s** |
| ref_1_2 (median, t∈[1,2]s) | 3.355 MPa | **3.3551 MPa** |
| startup ratio (peak/ref_1_2) | ≈1.19 | **1.1904** |
| burn-duration proxy (last t, P>0.5MPa) | ≈6.4–6.6 s | **6.56 s** |

All four match the brief to ≤0.03% — confirms the reversed-column recipe is
correct and the brief's numbers are this exact CSV, not a rounded/digitized
copy.

Chunc (`thomas_chunc_firing.csv`, `time_col=0, pressure_col=1`, 346 pt):
peak **8.8939 MPa** @ 0.397 s, plateau (median t∈[0.4,0.9]s) **8.6611 MPa**,
ratio **1.0269** — matches WORKLOG T2's `pl_exp=8.6611` exactly (same CSV,
same window).

### Full results — all 8 rows (actual numbers)

**Zerox** (MPa; `ratio` = startup peak (t<0.8s) / ref_1_2; `burn_duration` =
last t with P>0.5 MPa, seconds). All 4 rows `termination='t_max reached'`,
`erosive_model_stamp` matched the requested model in all 4 cases.

| erosive_model | cap | peak | t_peak (s) | ref_1_2 | ratio | burn_duration (s) | n_steps |
|---|---|---|---|---|---|---|---|
| ma2020 | 0.0 | 7.7835 | 0.0251 | 5.1557 | 1.5097 | 4.8428 | 3,593,660 |
| ma2020 | 1.1 | 7.5091 | 0.0241 | 5.1566 | 1.4562 | 4.8422 | 3,594,942 |
| mukunda_paul | 0.0 | 6.3919 | 0.0228 | 4.9260 | 1.2976 | 5.0457 | 3,610,471 |
| saint_robert | 0.0 | 3.8802 | 0.8000* | 4.1862 | 0.9269* | 5.6744 | 3,705,120 |

\* `saint_robert`'s `t_peak = 0.8000` is the right edge of the `t<0.8s`
startup-search window, not a genuine local max — the trace is still rising
through the whole window (same shape as Chunc's `saint_robert` row in T2:
`t_peak = t_max`). `ratio < 1` follows directly: the "startup peak" the
window finds is below the later `ref_1_2` value. Not an error, an artifact
of applying a peak-search window built for a spiking trace to a
non-spiking one — flagged so the table isn't misread as a sub-unity spike.

**Chunc** (MPa; `plateau` = median P_head, t∈[0.4,0.9]s; `ratio` =
peak/plateau) — identical committed `.ric` + canonical knobs as WORKLOG T2,
included here as the cross-check plot's data source. All 4 rows match T2's
15-row sweep **exactly** (same inputs, same outputs — confirms
determinism/reproducibility across sessions):

| erosive_model | cap | peak | t_peak (ms) | plateau | ratio | n_steps |
|---|---|---|---|---|---|---|
| ma2020 | 0.0 | 11.4728 | 16.45 | 8.2225 | 1.3953 | 531,651 |
| ma2020 | 1.1 | 9.3566 | 14.50 | 8.2405 | 1.1354 | 534,035 |
| mukunda_paul | 0.0 | 7.8410 | 1000.0* | 7.2999 | 1.0741 | 540,585 |
| saint_robert | 0.0 | 7.0174 | 1000.0* | 5.9966 | 1.1702 | 558,124 |

\* Same "trace still rising at `t_max`" shape as T2's Chunc `mukunda_paul`/
`saint_robert` rows — not a coincident early spike, see T2's own footnote.

### FITNESS VERDICT (task brief part C)

Criterion 1 — ma2020(cap=0.0) `ref_1_2` vs experimental 3.355 MPa:
**5.1557 MPa, error = |5.1557−3.355|/3.355 = 53.7%** (threshold ±25%) →
**FAILS**.

Criterion 2 — ma2020(cap=0.0) `burn_duration` vs experimental ~6.5 s:
**4.8428 s, error = |4.8428−6.5|/6.5 = 25.5%** (threshold ±25%) → **FAILS**
(marginally — 0.5 points over the line).

**VERDICT: Zerox is RADICALLY DIVERGENT / UNFIT** by both stated criteria,
one grossly (peak/ref-window pressure is 54% high) and one marginally
(burn-duration 25.5% short, just over the line). Per the task brief, this is
reported and NOT acted on — no fallback-motor selection made here.

### Zerox plateau-decomposition arithmetic (mirrors T2 review note's table)

Propellant `n`: read from `motors/zerox.ric`'s single propellant tab —
**n = 0.456** (`1−n = 0.544`), used below for the P→rate conversion. (Chunc's
`n = 0.4`, `1−n = 0.6`, confirmed from the committed `machbusterNew.ric`
extract's single propellant tab — same value T2's review note used.)

**Absolute values (MPa), both motors, saint_robert / mukunda_paul / ma2020
in the plateau/ref window** (Chunc: median t∈[0.4,0.9]s; Zerox: median
t∈[1,2]s):

| | Chunc plateau | Chunc plateau as % of exp (8.6611) | Zerox ref_1_2 | Zerox ref_1_2 as % of exp (3.355) |
|---|---|---|---|---|
| saint_robert | 5.9966 | 69.2% | 4.1862 | **124.8%** |
| mukunda_paul | 7.2999 | 84.3% | 4.9260 | **146.9%** |
| ma2020 | 8.2225 | 94.9% | 5.1557 | **153.7%** |
| ma2020 (cap=1.1) | 8.2405 | 95.1% | 5.1566 | **153.7%** |

**Zerox's picture is the OPPOSITE of Chunc's.** On Chunc, `saint_robert`
(zero erosive contribution) UNDER-shoots experimental and each closure adds
lift toward it. On Zerox, `saint_robert` ALONE already OVER-shoots
experimental by 24.8% — before any erosive enhancement is added — and every
erosive closure makes the over-shoot worse (mukunda_paul +46.9%, ma2020
+53.7%). `erosive_mp_cap=1.1` changes nothing at the ref window (5.1566 vs
5.1557, +0.02%) — the cap only binds during the early transient (compare the
peak column: 7.51 vs 7.78 MPa), consistent with T2's finding that the cap
converges to uncapped Ma at lower-g mid-burn conditions.

**Rate-space conversion** (P ∝ r^(1/(1−n)) ⇒ rate-lift η̄ = (P-lift over
SR)^(1−n), same identity T2's review note used):

| | Chunc (n=0.4, 1−n=0.6) | | Zerox (n=0.456, 1−n=0.544) | |
|---|---|---|---|---|
| | P-lift over SR | implied η̄ | P-lift over SR | implied η̄ |
| required (experimental) | 1.4443 | 1.2468 | **0.8014** | **0.8866** |
| ma2020 supplies | 1.3712 | 1.2085 | 1.2316 | 1.1200 |
| mukunda_paul supplies | 1.2173 | 1.1252 | 1.1767 | 1.0926 |
| ma2020 (cap=1.1) supplies | — | — | 1.2318 | 1.1201 |

On Chunc, "required" (η̄=1.2468) sits above both closures, closer to ma2020
(1.2085) than mukunda_paul (1.1252) — the T2 finding. On Zerox, "required" is
**below 1** (η̄=0.8866, i.e. Zerox's own experimental ref-window pressure
implies LESS burn rate than pure Saint-Robert alone already gives) while
BOTH closures push further above 1 (1.09–1.12) — every erosive closure
choice moves in the wrong direction relative to what Zerox's data would
need. This is not informative for the Chunc a/n-vs-closure degeneracy in the
way the discriminator was designed to be: it doesn't show "MP's deficit
recurring" (there is no deficit on Zerox — MP already exceeds experimental)
or "MP right, Ma over" in the Chunc sense. It instead shows Zerox's own
Saint-Robert `a`/`n` calibration is very likely too high on its own,
independent of the erosive-closure question — consistent with the
known-stale igniter-era calibration flagged in
`project_zerox_calibration_state` and confirmed quantitatively here for the
first time against the registry.

**One cross-propellant pattern DOES hold, and is worth keeping**: the
*relative* ordering and *magnitude* of the two closures' lift over
Saint-Robert is similar on both propellants — mukunda_paul supplies
+21.7% (Chunc) / +17.7% (Zerox) lift over SR; ma2020 supplies +37.1%
(Chunc) / +23.2% (Zerox); Ma/MP directly is 1.126 (Chunc) vs 1.047 (Zerox).
Both closures add a broadly similar-sized enhancement over SR on a
different propellant/geometry (finocyl+BATES vs Chunc's fore-tapered
grain), even though Zerox's absolute comparison to experiment is unusable.
This is weak, qualitative support for "the closure choice's effect size is
propellant-independent" (consistent with a closure-systematic reading) but
is explicitly NOT the clean discriminator the task was designed to produce,
because Zerox's own a/n calibration confounds the absolute comparison.

### Plots (`artifacts/erosive_model_axis/`, gitignored, dpi=130)

1. `zerox_models_full.png` — experimental + all 4 arms, t∈[0,8.2]s. All four
   sim arms sit well above experimental for the entire mid-burn (5–8 MPa vs
   exp's declining ~4→2 MPa); `ma2020` and `ma2020(cap=1.1)` are visually
   near-indistinguishable (their lines overlap almost everywhere) except a
   small kink ~1.8–2.1s where the cap briefly binds.
2. `zerox_models_startup.png` — same, zoomed t∈[0,1.5]s. Shows the real
   experimental ignition delay (P stays near 0 until ~0.3s, ramps to peak at
   0.43s) that none of the sim arms reproduce (all sim arms rise to their
   peak by t≈0.02–0.03s) — a second, independent sign the Zerox igniter/
   ignition-transient calibration is stale, separate from the plateau
   over-shoot.
3. `chunc_models.png` — same 4 arms + 346-pt experimental, t∈[0,1.0]s.
   Reproduces T2's Chunc picture (ma2020 spikes above experimental then
   undershoots the plateau slightly; mukunda_paul and saint_robert are both
   still rising, well below experimental, at t_max).
4. `plateau_decomposition.png` — grouped bar (Chunc, Zerox) × (saint_robert,
   mukunda_paul, ma2020), each bar annotated with its value, experimental
   reference line per group. Visually confirms the opposite-direction
   finding: Chunc's exp line sits ABOVE all three bars; Zerox's sits BELOW
   all three.

### Anomalies

1. **Row 1 wall-clock throughput was ~44× the normal rate — isolated,
   non-recurring.** Zerox `(ma2020, cap=0.0)`, run as a fresh standalone
   process, took **51m9.4s wall** for 3,593,660 steps (**1,171 steps/s**,
   vs the documented ~45–90k steps/s). Every other row was in the normal
   range: `saint_robert` (fresh process) 67.68s / 3,705,120 steps =
   **54,743 steps/s**; the chained `mukunda_paul` + `ma2020(cap=1.1)` run
   (same warm process, timed via `.npy` write timestamps rather than
   `time`, so approximate) ≈ 45,000 steps/s and **53,393 steps/s**
   (`ma2020(cap=1.1)`'s figure is exact: 3,594,942 steps / 67.34s between
   its and `mukunda_paul`'s `.npy` timestamps). The `ma2020(cap=1.1)` row —
   the SAME erosive-model branch as the anomalous row 1, same motor, same
   canonical knobs, differing only in cap value and in running inside an
   already-warm process — ran at a normal rate. Checked (light,
   non-exhaustive): `fmm_grain.py`'s FMM-table cache (line 737) is an
   **in-process Python dict, not a disk cache** — so it can't explain a
   fresh-process-only slowdown, and it doesn't, since `saint_robert` (also a
   fresh process needing fresh FMM resolution) was fast. Root cause not
   investigated further per the task brief's scope (do not investigate
   why — record and continue); flagged here because it is a genuine outlier
   against every other row measured in this session and in T1/T2, and
   because it burned the majority of this session's wall-clock budget.
2. **`saint_robert`'s and Chunc's `mukunda_paul`/`saint_robert`
   `t_peak`/`ratio` artifacts** — see the footnotes under the results
   tables; both are the known "still-rising trace, peak-search window built
   for a spiking trace" artifact already documented in T2, now also
   observed on Zerox.
3. **Zerox ignition-transient timing mismatch** — see plot 2 above; sim
   ignites near-instantly (t_peak≈0.02–0.03s) vs experimental's ~0.3–0.43s
   ramp. Not scored by the part-C verdict criteria (which target ref_1_2 and
   burn_duration, not ignition delay) but visually obvious in every plot and
   worth flagging as a second, independent symptom of stale calibration.
4. No ignition failures, no early terminations, no `erosive_model_stamp`
   mismatches across all 8 rows.

### Open items

- Zerox is confirmed UNFIT as a validation motor at its current calibration
  (verdict above) — its own Saint-Robert `a`/`n` appears too high
  independent of the erosive-closure question, on top of the already-known
  stale igniter timing. The T2 degeneracy (Chunc a/n-vs-closure) remains
  UNRESOLVED — Zerox could not serve as the clean discriminator it was
  designed to be. A future discriminator needs either a re-calibrated Zerox
  (outside this task's scope) or a different second propellant with usable
  static-fire data.
- The row-1 throughput anomaly (Anomaly 1) is unexplained; if it recurs on
  a future first-row-of-session Zerox invocation, worth a dedicated
  (separate-task) investigation — candidates not checked here: antivirus/OS
  scan of a freshly spawned `python.exe`, disk contention, thermal/power
  state, or something specific to first-touch of the Finocyl/FMM grain path.
- `C:/tmp/zerox_axis/run_remaining.py` (throwaway 2-row chained driver) and
  the full `C:/tmp/zerox_axis/*.npy` / `results.jsonl` accumulator are
  session-scratch only, not committed, will not persist.
- Did not touch `motors/*.ric` (hard rule) and did not commit/push (hard
  rule).
