# Erosive model registry — TASKS (agent-executable)

Read `DESIGN.md` first. Each task is self-contained, names its files and
acceptance gates, and requires no architectural decisions — if something is
ambiguous, STOP and ask rather than inventing. Run everything with the pyenv
interpreter `"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe"`
from the repo root, and delete `srm_1d/__pycache__/` after ANY edit to an
`@njit` function (Numba stale-cache is the #1 "my fix didn't work" trap).

**Work documentation:** on completing a task, append a dated entry to
`WORKLOG.md` in this directory: what was changed (file by file, one line
each), every gate run with its ACTUAL output numbers (not "passed"), anything
that deviated from this spec and why, and anything left open. The worklog is
the review artifact — write it for a reviewer who has not seen your session.

**Hard rules:** never edit `motors/*.ric` (openMotor-owned data); never
commit or push; do not fix the 2 pre-existing `BALLSstick VANILLA.ric` test
failures; do not refactor, reformat, or "improve" code outside the specified
edits — especially not the Ma code-0 path, which must stay literally
identical.

## T1 — registry + kernel branch + threading (core)

Files: `srm_1d/burn_rate.py`, `srm_1d/simulation.py`, `tests/test_burn_rate.py`.

1. `burn_rate.py`: add the `ErosiveModel` dataclass, `EROSIVE_MODELS` dict
   (`ma2020` code 0, `mukunda_paul` code 1, `saint_robert` code 2 —
   descriptors per DESIGN §1), `resolve_erosive_model(name)` (accepts a
   registered name; ValueError listing registered names otherwise) and
   `register_erosive_model(model)` (rejects duplicate name or code). Pure
   Python — NOT @njit.
2. `burn_rate_cell`: append `erosive_model_code` after `erosive_mp_cap`.
   Immediately after `r0` is computed (BEFORE the `Re_local < 100` guard):
   `if erosive_model_code == 2: return r0, 0.0` (pure Saint-Robert — no flow
   dependence). Keep the `Re_local < 100` early-return where it is. After the
   Re/h0/ebf guards pass, branch: code 1 computes
   `G = Re_local * mu_gas / D_hyd` and returns
   `(r0 * eta_mp_value, r0 * eta_mp_value - r0)` with
   `eta_mp_value = mukunda_paul_eta(G, r0, D_hyd, mu_gas, rho_p)`
   per DESIGN §2 — no bisection; code 0 is the existing path UNCHANGED.
   NOTE: the code-0 path must stay literally identical — do not refactor it.
3. `compute_burn_rates`: same appended parameter, passed through.
4. `simulation.py`: `_run_time_loop` gains `erosive_model_code` in the
   signature group right after `erosive_mp_cap` (comment: `# int — registry
   code, 0=ma2020, 1=mukunda_paul, 2=saint_robert`), passed to
   `compute_burn_rates`.
   `run_simulation` gains kwarg `erosive_model='ma2020'` (documented per
   DESIGN); at the top it resolves `_ero_code = resolve_erosive_model(
   erosive_model)` and raises ValueError if `erosive_mp_cap > 0` and
   `erosive_model != 'ma2020'`; the `_run_time_loop` call site passes
   `int(_ero_code)` in the matching slot. Stamp
   `result['summary']['erosive_model'] = <name>` where the summary dict is
   assembled.
5. Tests (`tests/test_burn_rate.py`): update the two direct call sites with the
   new trailing arg (`0`). Add: (a) registry tests — resolve known names,
   ValueError on unknown, duplicate-rejection; (b) G2 — MP mode returns exactly
   `r0 * mukunda_paul_eta(...)` at high Re, `r0` when g<35, `r0` when Re<100;
   (c) G2b — saint_robert returns `(r0, 0.0)` at low AND very high Re;
   (d) G3 — pick P/Re/D where `eta_Ma > eta_MP` (P=5e6, Re=5e5, D=0.04 works)
   and assert MP mode == Ma mode with `erosive_mp_cap=1.0`.
6. Add a `tests/` integration guard for the kwarg validation: `run_simulation`
   raising on `erosive_model='nope'` and on `mukunda_paul`+`erosive_mp_cap=1.0`
   can be tested via `pytest.raises` WITHOUT running a sim (validation happens
   before the time loop) — if it doesn't, move the validation earlier.
7. G2c integration check (script or test, your choice — if a test, keep it
   short): run Chunc (`motors/machbusterNew.ric` — but see the ⚠️ committed-
   vs-working-tree caveat in `docs/v0_7_4/probes/mp_cap_crossmotor.py`; use
   `git show HEAD:motors/machbusterNew.ric` extracted to a temp file) with
   `t_max=0.05`, canonical knobs (roughness=32e-6, kappa=0.44,
   pyrogen='bpnv', T_ignition=756.0), once with
   `erosive_model='saint_robert'` and once with `erosive_model='ma2020',
   diagnostic_disable_erosive=True`. The P_head traces must match to
   float tolerance (report max |ΔP|).

Acceptance: G1 (run `docs/v0_7_4/probes/mp_cap_crossmotor.py` chunc row at
cap 0.0 → peak 11.47 @ ~16.5 ms, plateau 8.22), G2, G2b, G2c, G3, plus
`python -m pytest tests/test_burn_rate.py -q` green and full
`python -m pytest tests/ -q` at the pre-existing pass count (the 2
`BALLSstick VANILLA.ric` repo-scan failures are pre-existing — do not fix).
Write the WORKLOG.md entry (see header) with all gate numbers.

## T2 — cross-motor validation (model axis)

Files: `docs/v0_7_4/probes/mp_cap_crossmotor.py`.

Extend the `__main__` sweep to a model axis: rows for
`(erosive_model='ma2020', cap 0.0)`, `('ma2020', 1.1)`, `('mukunda_paul', 0.0)`,
`('saint_robert', 0.0)` on all three motors (keeping the old cap=1.0 row is
fine — 5 rows). Print the same metrics table. The `saint_robert` rows are
recorded, not gated — they are the plateau-decomposition data for the
spike/plateau decoupling investigation (report Chunc plateau_SR prominently).
Read the ⚠️ caveat in the module docstring: the run must use the COMMITTED
motors/*.ric (if `git status` shows machbusterNew.ric or BALLSstick.ric
modified, extract via `git show HEAD:motors/<name>.ric` to a temp path and
point the runs there — note it in the output header).

Acceptance: G4 windows from DESIGN (±0.03 on ratios, ±0.15 MPa on peaks) for
the gated rows; paste the output table into the WORKLOG entry. Flag, don't
hide, any miss.

## T3 — GUI leg (openMotor fork; separate branch discipline)

Files: `srm_1d/srm1d_plugin.py` (+ verify-only in the openMotor fork checkout).

1. `_transient_config_schema()`: add
   `EnumProperty('Erosive Model', list(EROSIVE_MODELS))` default `'ma2020'`
   and `FloatProperty('MP Erosive Ceiling (0=off)', '', 0.0, 2.0)` seeded 0.0
   mapped to `erosive_mp_cap`. Keys MUST equal the run_simulation kwarg names
   (`erosive_model`, `erosive_mp_cap`) — that is the whole passthrough
   contract. No unit conversion (cf. the roughness um→m special case — not
   needed here).
2. VERIFY (do not assume): the fork's motor-editor config screen renders
   `EnumProperty` generically (check how `uilib/widgets/motorEditor.py`'s
   config controller builds property editors). If it does not, STOP and
   report — do not hand-roll a widget.
3. Round-trip test: build the schema, `setProperties({'erosive_model':
   'mukunda_paul'})`, `getProperties()`, assert survival; and an offscreen
   smoke if the existing GUI smoke harness covers the config screen.

Acceptance: G6 + existing plugin tests green
(`python -m pytest tests/test_srm1d_plugin.py -q`).

## T4 — docs

Files: `srm_1d/ARCHITECTURE.md` (burn_rate.py one-liner gains "erosive-model
registry"), `srm_1d/DEVNOTES.md` (API Breaking Changes Log: burn_rate_cell /
compute_burn_rates / _run_time_loop signature + the new kwargs),
`srm_1d/README.md` (public-API mention: `erosive_model=` with the
default-ma2020 / no-tuning framing and the MP option's tradeoff, 3-4 lines,
pointing at `docs/erosive_model_registry/` + REOPENED §10),
`docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` (one line in §10.5 noting the swap
feature now exists via the registry).

Acceptance: docs mention nothing not actually implemented; tone matches
surrounding files; no claims of "recommended" — model choice is explicitly the
user's tradeoff.
