# Plugin-provided igniter — worklog (I1 / I2 / I3)

Executed 2026-10-05. srm_1d worktree `C:/tmp/sv/ig` (branch `ig/provider`, base
`5b91583`); fork worktree `C:/tmp/sv/ig-om` (branch `igniter/plugin-provider`,
base `92e425f`). Scratch: `C:/tmp/sv_ig/`. Interpreter pyenv 3.10.5.

The first run was interrupted by an API session limit after the I1 edits (not
by a task problem); on resume the I1 diff was reviewed against TASKS I1 and kept
(nothing needed redoing); the fork worktree was still untouched.

## Preflight (env as TASKS.md; re-run after the resume)

```
C:\tmp\sv\ig\srm_1d\__init__.py
C:\tmp\sv\ig-om\motorlib\__init__.py
```

## Baselines (before any edit of the respective repo)

| | result | how |
|---|---|---|
| srm_1d suite (`pytest tests/ -q -p no:cacheprovider`, clean tree `5b91583`) | **604 passed, 0 failed**, 484 s | MEASURED by me before editing (matches the reviewer's verified baseline) |
| B0 (`baseline_gate.py`) | **PASS** (chunc/balls/haseg all rel 0) | MEASURED |
| fork unit suite | **72 tests OK** | MEASURED before any fork edit |

Fork suite command: the CI one is `PYTHONPATH=. python test/unit.py`, which
needs the Cython `mathlib._find_perimeter_cy` (not built in this worktree).
Equivalent run used here: `python C:/tmp/sv_ig/run_fork_unit.py` — puts
`ig-om/test` on `sys.path`, installs srm_1d's existing Cython-replacement shim
(`srm_1d.fmm_grain._install_fake_find_perimeter_cy`), loads `unit` and runs it
with `unittest.TextTestRunner` (same discovery as `test/unit.py`).
`mathlib/_find_perimeter_cy` shim and the pyuic-generated `uilib/views/*_ui.py`
(gitignored build artifacts, needed only by the GUI smoke) are the only
non-source things in the fork worktree.

## I1 — srm_1d provider (commit `0828366`)

Files: `srm_1d/openmotor_adapter.py`, `srm_1d/srm1d_plugin.py`,
`srm_1d/__init__.py` (exports), `srm_1d/DEVNOTES.md`,
`tests/test_igniter_provider.py` (new, 39 tests).

* `SUPERSEDED_PYROGENS['BPNV']` = the fork's former `DEFAULT_IGNITER['pyrogen']`
  (verified by `ast`-extracting the fork literal: identical to
  `DEFAULT_PYROGENS[0]`, 16 keys). `upgrade_igniter_block`, `library_materials`,
  `IgniterMaterialUpgradedWarning`; `load_igniter` upgrades first (one warning),
  then T1's stale check (silent for an upgraded block). `name` is compared
  case-insensitively (the lookup is).
* Plugin: `_PYROGEN_FIELD_META` (17 fields incl. new `gas_mass_fraction`
  0..1; labels/units/ranges ported verbatim), `_PYROGEN_GUI_PASSTHROUGH =
  {Cp_gas}`, `_IGNITER_CHAMBER_META` (9 keys), five hooks on
  `Srm1dTransientSolver`. Schemas are fresh each call and seeded with srm_1d's
  own defaults (chamber: `_IGNITER_SIZING_DEFAULTS`; material: the `Pyrogen`
  dataclass defaults, so a material missing `kappa_jet` shows 8.0 not the
  range minimum).
* **Enum-choice sources (recorded as required):** `injection_topology` =
  `igniter_plenum._VALID_TOPOLOGIES` (srm_1d validates; used directly).
  `burn_law` (`'0d','end_burning'`) and `heat_delivery_mode`
  (`'demar','radiation','none'`) are validated only by if-chains
  (`_burn_law_code`, `simulation._heat_delivery_code`); the lists are kept in the
  fork's order and a test pins every choice as accepted (and a bogus value as
  rejected) by those validators. `form` is not validated anywhere in srm_1d
  (informational) -> the fork's former list `pellets, powder, chunks`.
* Tests (all pass, 39): meta ∪ passthrough == `Pyrogen` init fields (disjoint);
  chamber meta == sizing keys; enum choices vs validators; schema structure;
  range fit for all 4 library materials + default chamber (exact read-back);
  hooks (library, default, fresh objects, upgrade delegation); upgrade (legacy →
  current, chamber untouched, idempotent, `Cp_gas: None` ok, case-insensitive
  name, rel tol 1e-9 boundary, each of the 15 non-name legacy fields altered →
  untouched, missing/extra non-None key → untouched, MTV etc. untouched, input
  not mutated, degenerate blocks, registry well-formed); `load_igniter` (one
  `IgniterMaterialUpgradedWarning`, no stale warning, `gas_mass_fraction==0.781`,
  `M==0.0482`; customised legacy → stale warning only). **User-run addition
  (a):** a legacy dict with the user's GUI float noise (`rho:
  1699.9999999999998`, `particle_diameter_m: 0.0032000000000000006`; values
  copied into the test, user file never read) is recognised and upgraded.
  Mutation checks (in-process): range `gamma.max=1.1` → range-fit test fails;
  dropping `gas_mass_fraction` meta → coverage test fails.
* Gates: B0 PASS (identical). srm_1d suite after I1 (before the 39th test):
  **642 passed, 0 failed** (604 + 38), 427 s. No existing test's numbers moved.
  (After the GUI-noise test: new file 39/39.)

## I2 — fork consumes the provider (commit `9ed4bf8` fork; `e5e1e08` srm_1d)

Fork files: `motorlib/solvers.py` (5 hooks + `get_igniter_provider`),
`motorlib/igniter.py` (hardcoded `Pyrogen`/`Igniter` removed; helpers
`material_collection`, `chamber_collection`, `changed_values`, `merge_edits`),
`motorlib/motor.py` (`igniterData` raw dict; applyDict deep-copies + provider
upgrade; getDict deep-copies, omits the key when `None`), `uilib/fileIO.py`
(`DEFAULT_IGNITER` deleted; 0.6.1→0.7.0 seeds provider default, nothing without a
provider), `uilib/defaults.py` (`DEFAULT_PYROGENS` deleted),
`uilib/pyrogenManager.py` (raw dicts; D6 load), `uilib/widgets/pyrogenMenu.py`,
`solverConfigController.py`, `motorEditor.py`, `mainWindow.py`, `app.py`
(discovery before `PyrogenManager`), `test/unit/igniter.py` (+ registered in
`test/unit/__init__.py`).

Call-site re-grep (`igniterPyrogen|\.igniter\b|motorlib.igniter|DEFAULT_PYROGENS|
DEFAULT_IGNITER|igniterSource`): all hits were in the listed files; none left.

Design details beyond the spec text (all inside D3/D4):
* **Widget baseline.** The editors round-trip values through display units and
  8-decimal spin boxes, so comparing widget output to the pre-widget collection
  snapshot would flag every unit-converted field as "changed" (e.g. `rho` kg/m³ →
  lb/in³). The controller / pyrogen menu therefore also snapshot the WIDGET
  values right after the page is rendered (after the 'auto' sentinel tweaks) and
  only feed widget values that moved off that baseline into the working
  collection; `igniterProps()` is then `changed_values(collection snapshot,
  working collection)`. Verified by the smoke (a no-edit Apply leaves the igniter
  block bit-identical).
* **Library pick** is reported only when it differs from the motor's current
  material name (a no-op re-pick of the same name does not overwrite an embedded
  customised material); a different pick deep-copies the library dict.
* **Motor with no block** (new motor): getDict has no `igniter` key;
  config page shows the provider default; the block is created
  (`merge_edits(None,…)` = provider default + changes) only if the user changes
  something; simulation falls back to library BPNV (`simulate_motor`).
* `SolverConfigController(editor, general, solverConfigs, withIgniter=False,
  igniterBlock=None)` replaces `igniterSource=`.
* `PyrogenManager`: provider-material "absent" test is exact-name; no provider →
  empty list and no file written.
* Known corner: a schema field whose raw value is out of range is held at the
  property minimum in the working collection; editing it to exactly that minimum
  is not detected as a change (any other edit is).
* The existing srm_1d test `test_motor_without_igniter_keeps_defaults` encoded
  the retired "Motor invents a default block" contract; rewritten as
  `test_motor_without_igniter_has_no_block_and_sim_falls_back` (deviation).

Fork unit tests (`test/unit/igniter.py`, 14 tests): provider registered +
base-class defaults; no block → no key; lossless round trip (unknown key, `None`,
out-of-range, nested unknown, deep-copy in/out); legacy BPNV (with the user's GUI
noise) upgraded on applyDict; customised not touched; no provider (patched) →
verbatim; 0.6.1→0.7.0 migration seeds provider default / nothing without one;
`changed_values` / `merge_edits` / `material_collection`; `PyrogenManager` with
temp `getConfigPath`: missing file → seeded (4 materials), **legacy single-BPNV
library → upgraded and every srm_1d material present (BPNV, MTV, MTV_fast,
Cu/Al thermite) — user-run addition (b)**, customised entry kept verbatim,
current file not rewritten, no provider → empty without crashing. The tests
point `platformdirs` at a temp dir while `uilib` imports (its logger opens a log
file at import), so the user's openMotor dirs are untouched.

Gates: fork unit suite **86 tests OK** (72 + 14, 0 skipped, 0 failures);
`tests/test_srm1d_plugin.py` + `test_igniter_provider.py` + `test_ric_igniter.py`
vs the fork worktree: 66 passed.

**Smoke** (`docs/igniter_provider/smoke/gui_smoke.py`, `QT_QPA_PLATFORM=offscreen`):
runs the real `App` (temp `platformdirs`, modal dialogs stubbed), PASS:
solver + provider discovered; GUI library holds BPNV, MTV, MTV_fast, Cu/Al
thermite; temp copy of committed `machbusterNew.ric` embeds legacy BPNV (T_flame
2800) on disk and loads as current BPNV (`gas_mass_fraction` 0.781); no-edit
Apply → igniter block bit-identical (only `solverConfigs[srm_1d]` is added —
openMotor's normal Apply behaviour); edit `basket_fill_fraction` 0.4 and `mass`
0.002 through the widgets → saved → reloaded: material == current BPNV, both
edits persisted (`mass` 0.0019999999999999996 = unit-conversion noise on the
EDITED key), zero other differing key paths; library pick MTV deep-copied, chamber
unchanged; pyrogen-library editor: no-edit Apply bit-identical, one edit changes
exactly one key, `Cp_gas: None` and an unknown key preserved; motors/ untouched.

## I3 — refreshed motor files (commit `75d59b8`)

Tool `srm_1d/tools/refresh_ric_igniters.py` (`python -m
srm_1d.tools.refresh_ric_igniters <paths> [--dry-run]`): fork `loadFile` →
`Motor.applyDict` → `getDict` → diff of the parsed data (tuples normalised) →
`saveFile` only if every changed key path is permitted (the upgraded
`igniter.pyrogen`, plus ADDED `config/taperSlices` and `solverConfigs`); re-reads
the saved file and re-checks. Redirects the openMotor log dir to a temp dir.

Dry run over all 21 committed files, then the 18 clean ones were saved (textual
diff ≈ 6 lines per file: 4 material lines + `taperSlices: 0` + `solverConfigs: {}`).

| file | material | upgraded | changed key paths |
|---|---|---|---|
| 2025.12.25 54-2800 ST2.0, 75_5120ChaseRed_BATES, 76_8000FatBoy…, Aeroblast_flight_2, BALLSstick, BATES Test, ISP_Super_Loki, L1843, L3035, hasegawa_a/b/c, hlbl38, hlbltest, ivanO25k, machbusterNew (16 files) | BPNV | yes | `igniter/pyrogen/{M 0.03→0.0482, T_flame 2800→2719, gamma 1.25→1.17}`, +`igniter/pyrogen/gas_mass_fraction 0.781`, +`config/taperSlices`, +`solverConfigs` |
| M1378 Geo Test, machbusterFAIL (2) | BPNV | yes | same, + `igniter/pyrogen/Cp_gas` (None); already had `solverConfigs` |
| **example_bates** | BPNV | would be | **REFUSED (STOP)**: round trip ADDS `config/{burnoutThrustThres, burnoutWebThres, flowSeparationWarnPercent, mapDim, maxMachNumber, maxMassFlux, maxPressure, minPortThroat, sepPressureRatio, timestep}` (hand-authored file lacking them; openMotor would fill `MotorConfig` minimum-value defaults) — not saved |
| **zerox**, **zerox_LHS** | BPNV | would be | **REFUSED (STOP)**: round trip ADDS `config/maxMachNumber` (default 0.0 would silently change the file's config) — not saved |

Refused files are untouched and still embed the legacy BPNV; they are upgraded
at load by `load_igniter` (one `IgniterMaterialUpgradedWarning`). Decision on
whether to hand-complete their `config` blocks / re-save from the GUI is the
user's (open item).

Gates: **B0 PASS** (run on the committed refreshed `machbusterNew/BALLSstick/
hasegawa_a`: identical, n_steps 531651 / 842922 / 2335770). **Chunc equivalence**
(`run_from_ric` committed refreshed `machbusterNew.ric`, canonical knobs,
`t_max=0.05`; `pyrogen=None` vs `pyrogen='bpnv', pyrogen_mass=None`): 30457 vs
30457 steps, time arrays equal, **max |ΔP_head| = 0.0**, peak 11.4728 MPa.
Warning scan: `load_igniter(load_ric(f))` is **silent for all 18 refreshed files
and for the Chunc run**; only the 3 refused files emit
`IgniterMaterialUpgradedWarning`. Final srm_1d suite (after refresh, against the final fork worktree):
**643 passed, 0 failed** (604 + 39 new), 428 s, 98 warnings (was 106 at I1: the refreshed
files no longer emit stale-material warnings). No existing test's numbers moved.

## Deviations / open items

1. Fork suite run through a shim wrapper (Cython ext not built) — same tests.
2. 3 motor files refused by the refresh safety check (above).
3. `test_motor_without_igniter_keeps_defaults` rewritten (retired contract).
4. Widget-baseline mechanism and "pick only if changed" semantics (above) are
   implementation choices within D4.
5. Out-of-range-field corner case (above). Exact-name rule for the library append.
6. The user's real `%LOCALAPPDATA%\openMotor` files were never read or written
   (all tests/smoke/tool runs use temp dirs); the user's live `pyrogens.yaml` will
   be upgraded + completed with MTV/MTV_fast/thermite on their next GUI start
   with the new fork.
7. Not pushed; no stash used.

## I3b — the three refused files (follow-up)

`refresh_ric_igniters.py` permit rule extended: an ADDED `/config/<key>` is
permitted iff its value equals the fork's `MotorConfig` default for that key
(`motor_config_defaults()` reads `MotorConfig().getProperties()`; nothing
hardcoded). Rationale: srm_1d reads only `ambPressure` from a motor's config
and openMotor fills absent keys with those defaults on load, so this is pure
normalisation. Any other added/removed/changed key is still a STOP. Test:
`tests/test_refresh_ric_igniters.py` (7 tests: at-default permitted;
non-default / changed / removed / unknown key / other section refused; old rules
unchanged).

Tool re-run on the three files (all `[ok]`, saved, RESULT OK):

| file | material | upgraded | changed key paths |
|---|---|---|---|
| example_bates | BPNV | yes | added `config/{burnoutThrustThres, burnoutWebThres, flowSeparationWarnPercent, mapDim, maxMachNumber, maxMassFlux, maxPressure, minPortThroat, sepPressureRatio, taperSlices, timestep}` (all at MotorConfig defaults); `igniter/pyrogen/{M,T_flame,gamma}` changed, `+gas_mass_fraction`; added `solverConfigs` |
| zerox | BPNV | yes | added `config/maxMachNumber`, `config/taperSlices`; pyrogen M/T_flame/gamma/+gas_mass_fraction; added `solverConfigs` |
| zerox_LHS | BPNV | yes | same as zerox |

Sim equivalence (`run_from_ric`, `pyrogen='bpnv', pyrogen_mass=None`, roughness
32e-6, kappa 0.44, T_ignition 756, `t_max=0.05`; before = `git show
HEAD:motors/<f>.ric` in scratch — the files were unmodified at HEAD, so this is
the same content as the requested `HEAD~`):

| file | steps before/after | time equal | max abs dP_head | peak MPa |
|---|---|---|---|---|
| example_bates | 57591 / 57591 | yes | 0.0 | 3.2566 |
| zerox | 27723 / 27723 | yes | 0.0 | 7.7835 |
| zerox_LHS | 27742 / 27742 | yes | 0.0 | 6.7610 |

All 21 committed files now load silently through `load_igniter` (no stale /
upgraded warning).

## I4 — stale GUI solver-config defaults (spec cherry-picked as `02e72f4`)

srm_1d (commit after `a8efc19`): `_transient_config_schema` seeds are derived
from `inspect.signature(run_simulation)` (`_solver_config_defaults`; roughness
m→µm via `round(x*1e6, 9)`): t_max 10, P_cutoff 5e5, cfl 0.3, roughness 32 µm,
kappa 0.44, T_ignition 756. `SUPERSEDED_SOLVER_CONFIG = {roughness:[50.0],
kappa:[0.45], T_ignition:[850.0]}`; `upgrade_solver_config(values) ->
(new_values, upgraded_keys)` (pure, idempotent, rel tol 1e-9) and the plugin
hook. `tests/test_solver_config_upgrade.py` (5 tests): schema defaults ==
`run_simulation` defaults for every key; superseded ≠ current; the user's saved
config (copied into the test) upgrades exactly roughness/kappa/T_ignition and
keeps P_cutoff 344750; customised (roughness 40) untouched; idempotent / noise /
empty. `tests/test_srm1d_plugin.py::test_solver_config_schema` asserted the stale
kappa 0.45 / roughness 50 → updated to 0.44 / 32 (deviation: it encoded the bug).

Fork (`d7ece77`): generic `SolverPlugin.upgrade_solver_config` (default
passthrough) + `solvers.upgrade_solver_configs(configs)`; applied in
`PreferencesManager.loadPreferences` (upgraded keys logged; preferences saved if
changed) and `Motor.applyDict` (per-motor `solverConfigs`; logged via `logging`).
`app.py` now discovers plugins BEFORE the `PreferencesManager` (it was after).
Fork tests +4 (`test/unit/igniter.py`): base passthrough; helper (not mutated,
unknown solver / non-dict pass through); `Motor.applyDict` upgrades per-motor
config, customised kept; preferences load with a temp config dir applies the hook
and saves (reloaded file holds the upgraded values). Fork suite **90 OK, 0 skipped**
(72 + 14 + 4).

Smoke: the committed `machbusterNew.ric` is now refreshed, so the smoke restores
the legacy BPNV snapshot into a temp copy of it first; 21 PASS, still green after
the `app.py` reorder.

Gates: B0 PASS; srm_1d suite **655 passed, 0 failed** (643 + 7 refresh-tool + 5
config-upgrade), 441 s; fork 90 OK. No pre-existing test moved except the one
schema-default assertion above.

Open: a user who deliberately types exactly 50 µm / 0.45 / 850 K in the dialog is
indistinguishable from the stale default and is upgraded on the next load (inherent
to exact-match upgrade, as for materials). The user's saved preferences file is
upgraded on their next GUI start with the new fork.
