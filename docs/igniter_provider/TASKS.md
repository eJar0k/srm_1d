# Plugin-provided igniter — TASKS (agent-executable)

Read `DESIGN.md` first; decisions D1–D7 are settled — do not redesign. If a
step is ambiguous or impossible as written, STOP and report. Execute I1 → I2
→ I3 in order (I2 consumes I1's hooks; I3 needs both).

## Workspaces (two repos, two worktrees)

| repo | worktree | branch | base |
|---|---|---|---|
| srm_1d | `C:/tmp/sv/ig` | `ig/provider` | srm_1d `openmotor-frontend` HEAD |
| openMotor fork | `C:/tmp/sv/ig-om` | `igniter/plugin-provider` | fork `staging` |

Never touch the main checkouts
(`C:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d`
and `.../openMotor/openMotor`), other worktrees, other branches, or the
user's `%LOCALAPPDATA%\openMotor` config. Never push. Never `git stash`.
Scratch: `C:/tmp/sv_ig/`.

**Every shell command that runs Python must export all three:**
```
export SRM1D_OPENMOTOR_PATH="C:/tmp/sv/ig-om"
export OPENMOTOR_SRM1D_PATH="C:/tmp/sv/ig"
export PYTHONPATH="C:/tmp/sv/ig;C:/tmp/sv/ig-om"
```
Interpreter: `"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe"`
(has numba, scikit-fmm, PyQt6 6.7.1). `srm_1d` is pip-installed EDITABLE
pointing at the main checkout, so the env above is what makes your worktrees
win. **Preflight (record in the worklog):**
`python -c "import srm_1d, motorlib; print(srm_1d.__file__); print(motorlib.__file__)"`
must print paths inside `C:/tmp/sv/ig` and `C:/tmp/sv/ig-om`. If not, STOP.
After editing any `@njit` code (you should not need to) delete
`srm_1d/__pycache__/`.

**Baselines first:** before editing, record (a) the srm_1d suite
(`python -m pytest tests/ -q -p no:cacheprovider` from `C:/tmp/sv/ig`;
expect 604 passed, 0 failed), (b) the fork unit suite (find how
`test/unit.py` / `test/unit/` run; record the exact command and counts).

**Gate B0** (srm_1d): from `C:/tmp/sv/ig`,
`python docs/solver_verification/probes/baseline_gate.py` must print
`B0 RESULT: PASS` (n_steps exact, floats within rtol 1e-12).

**Commits:** focused conventional commits on each branch, each message
ending `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
**Worklog:** `docs/igniter_provider/worklog.md` in the srm_1d worktree,
sections I1 / I2 / I3, with actual numbers, every file changed, every
deviation, open items. Commit it.

---

## I1 — srm_1d: igniter provider + superseded-snapshot upgrade

Files: `srm_1d/openmotor_adapter.py`, `srm_1d/srm1d_plugin.py`, new
`tests/test_igniter_provider.py`, `srm_1d/DEVNOTES.md`.

1. `openmotor_adapter.py`:
   - `SUPERSEDED_PYROGENS = {'BPNV': [<legacy dict>]}` where the legacy dict
     is EXACTLY the `pyrogen` dict of the fork's `DEFAULT_IGNITER` in
     `uilib/fileIO.py` (copy it from `C:/tmp/sv/ig-om` BEFORE I2 deletes
     it; comment: "pre-bfc2f3f library BPNV; superseded 2026-06-18").
   - `class IgniterMaterialUpgradedWarning(UserWarning)`.
   - `library_materials()` → `[ _pyrogen_to_block(load_pyrogen(path)) ]`
     for every `srm_1d/pyrogens/*.yaml`, sorted by filename (reuse T1's
     `_library_pyrogen_files` if it fits).
   - `upgrade_igniter_block(block)` → `(new_block, changed)`; pure (deep
     copy; never mutates input). Rule (DESIGN D5): the block's `pyrogen`
     matches a superseded snapshot of the same name (case-insensitive) iff
     every snapshot key is present and equal (floats rel tol 1e-9) and every
     key NOT in the snapshot is absent or `None`. On a match, replace
     `block['pyrogen']` with the current library material of that name;
     leave every other key of the block untouched.
   - `load_igniter`: call `upgrade_igniter_block` first; if changed, emit
     one `IgniterMaterialUpgradedWarning` naming the material and saying the
     embedded value was a superseded library snapshot, the current library
     values are used, and re-saving the motor persists the upgrade. Then run
     T1's stale check on the upgraded material (it must stay silent for an
     upgraded block).
2. `srm1d_plugin.py`:
   - `_PYROGEN_FIELD_META`: per `Pyrogen` init field → label, unit, min, max
     (FloatProperty) or label + choices (EnumProperty). Port label/unit/range
     verbatim from the fork's current `motorlib/igniter.py` (`Pyrogen`) for
     the 16 fields it has; add `gas_mass_fraction` → ('Gas Mass Fraction',
     '', 0, 1). `_PYROGEN_GUI_PASSTHROUGH = {'Cp_gas': '<reason: optional,
     None = derive from gamma/M; FloatProperty cannot hold None>'}`.
   - `_IGNITER_CHAMBER_META`: per `_IGNITER_SIZING_DEFAULTS` key, ported
     verbatim from the fork's `Igniter` class.
   - Enum choices (`form`, `heat_delivery_mode`, `burn_law`,
     `injection_topology`): if srm_1d itself defines/validates the allowed
     values (search `propellant.py`, `igniter_plenum.py`, the adapter), use
     those as the source; otherwise use the fork's lists. Record which.
   - Implement the five DESIGN D2 hooks on `Srm1dTransientSolver`. Schema
     hooks return a FRESH `PropertyCollection` each call (build like
     `_transient_config_schema`); `get_igniter_materials` →
     `library_materials()`; `get_default_igniter` →
     `default_igniter_block('bpnv')`; `upgrade_igniter_block` → the adapter
     function.
3. Tests (`tests/test_igniter_provider.py`; import motorlib the way
   `tests/test_srm1d_plugin.py` does, skip cleanly if unavailable):
   a. Material meta keys ∪ passthrough == `Pyrogen` init fields, disjoint;
      chamber meta keys == `_IGNITER_SIZING_DEFAULTS` keys.
   b. Range fit: for every library material, loading its non-None values
      into the material schema and reading back returns each value exactly
      (catches a range that would make `FloatProperty` silently drop it);
      same for `default_igniter_block('bpnv')`'s chamber values.
   c. Hooks: materials == `library_materials()` with unique names; default
      == `default_igniter_block('bpnv')`; schema hooks return distinct
      objects per call.
   d. Upgrade: legacy BPNV block → material == current BPNV library dict,
      every chamber key unchanged, `changed` True; second call `changed`
      False; legacy with one field altered → untouched, False; legacy plus
      `Cp_gas: None` → upgrades; an MTV block → untouched; input dict not
      mutated.
   e. `load_igniter` on a legacy block → exactly one
      `IgniterMaterialUpgradedWarning`, returned `Pyrogen` has
      `gas_mass_fraction == 0.781` and `M == 0.0482`, and NO
      `StaleIgniterMaterialWarning`.
4. Gates: B0 PASS; srm_1d suite = baseline + new, 0 new failures. Tests that
   simulate a committed `.ric` WITHOUT `pyrogen=` now get the upgraded
   igniter — if any test's numbers move, report test + old/new values (do
   not loosen silently).

## I2 — fork: consume the provider

Files (fork worktree): `motorlib/solvers.py`, `motorlib/igniter.py`,
`motorlib/motor.py`, `uilib/fileIO.py`, `uilib/defaults.py`,
`uilib/pyrogenManager.py`, `uilib/widgets/pyrogenMenu.py`,
`uilib/widgets/mainWindow.py`, `uilib/widgets/motorEditor.py`,
`uilib/widgets/solverConfigController.py`, `app.py`, new fork unit test.

1. `solvers.py`: the five DESIGN D2 hooks on `SolverPlugin` (defaults:
   `None`, `None`, `[]`, `None`, `(block, False)`) + `get_igniter_provider()`.
2. `igniter.py`: remove the hardcoded `Pyrogen` / `Igniter` classes. Provide
   provider-driven helpers: `material_collection(values)` /
   `chamber_collection(values)` (provider schema populated from a dict,
   skipping `None` values; `None` if no provider) and a D4 write-back helper
   (`changed_values(snapshot, collection)` → dict of keys whose value differs
   from the post-load snapshot). Update the module docstring.
3. `motor.py`: `self.igniterData` raw dict (D3); `applyDict` deep-copies the
   block and, if a provider exists, applies `upgrade_igniter_block`;
   `getDict` writes a deep copy when not `None`. Remove `igniter` /
   `igniterPyrogen` and update every caller. Known call sites (2026-10-05;
   re-grep `igniterPyrogen|\.igniter\b|motorlib.igniter|DEFAULT_PYROGENS|DEFAULT_IGNITER`
   to confirm): `mainWindow.py` `applyMotorConfig` (~105–120) and the
   `loadMotorConfig` call (~347); `motorEditor.py` (~516, ~534);
   `solverConfigController.py` `__init__` / `usesIgniter` (~50–90);
   `pyrogenManager.py` (whole); `pyrogenMenu.py` (~13–16, ~110).
4. `fileIO.py`: delete `DEFAULT_IGNITER`; `migrateMotor_0_6_1_to_0_7_0`
   seeds `provider.get_default_igniter()` when a provider exists (adds no
   igniter key otherwise). `defaults.py`: delete `DEFAULT_PYROGENS`.
5. `pyrogenManager.py`: raw-dict materials; D6 load (seed from provider if
   file missing; upgrade each entry; append missing provider materials; save
   if changed); no provider → empty library without crashing.
6. `pyrogenMenu.py`: the "new material" seed is
   `provider.get_default_igniter()['pyrogen']`; the editor edits through the
   provider material schema with D4 write-back.
7. `solverConfigController.py`, `motorEditor.py`, `mainWindow.py`: chamber
   editing via the provider chamber schema with D4 write-back into
   `motor.igniterData`; a library pick deep-copies the library dict into
   `igniterData['pyrogen']`; the current material name is read from
   `igniterData`.
8. `app.py`: call `solvers.discover_external_solvers()` BEFORE constructing
   `PyrogenManager`.
9. Fork unit tests (new module, registered the way the existing
   `test/unit/*` modules are): (a) with the srm_1d provider registered,
   `Motor` applyDict→getDict is lossless for a block containing an unknown
   key, a `None` value and an out-of-range value; a legacy BPNV block is
   upgraded on applyDict; (b) with no provider (patch
   `get_igniter_provider` → `None`) the block passes through verbatim;
   (c) the 0.6.1→0.7.0 migration seeds the provider default; (d)
   `PyrogenManager` with a temp config dir (patch `getConfigPath`): missing
   file → seeded from provider; a legacy `pyrogens.yaml` (old BPNV) →
   upgraded and saved; a customised entry → kept verbatim.
10. Offscreen GUI smoke — `docs/igniter_provider/smoke/gui_smoke.py` in the
    srm_1d worktree, run with `QT_QPA_PLATFORM=offscreen`: start the app the
    way `main.py` does (without entering the event loop); load a temp copy
    of the committed `motors/machbusterNew.ric` (from
    `git -C C:/tmp/sv/ig show HEAD:motors/machbusterNew.ric`); assert the
    material is now the current BPNV; via the motor-config path, change one
    chamber field and apply a library material pick; save to a temp path;
    reload; assert `gas_mass_fraction == 0.781`, the edited chamber field
    persisted, and no other key changed. If the GUI cannot be driven
    headlessly, STOP this step and report exactly why.

Gates: fork unit suite = baseline + new, 0 new failures; srm_1d
`tests/test_srm1d_plugin.py` and `tests/test_igniter_provider.py` green
against the fork worktree; smoke passes.

## I3 — refresh the committed motor files

1. Tool `srm_1d/tools/refresh_ric_igniters.py`, run as
   `python -m srm_1d.tools.refresh_ric_igniters <paths...> [--dry-run]`
   (import the fork lazily with a clear error if unavailable). Per file:
   `uilib.fileIO.loadFile(path, fileTypes.MOTOR)` (applies openMotor
   migrations) → `Motor().applyDict` (applies D5) →
   `uilib.fileIO.saveFile(path, motor.getDict(), fileTypes.MOTOR)`. Print
   per file: material name, upgraded yes/no, any stale-warning fields left.
2. Before saving each file, compare the loaded data dict before and after
   the round-trip (parsed, not text) and print every changed key path.
   Permitted: the `igniter.pyrogen` material values (D5) and keys
   openMotor's save normally adds (e.g. `config.taperSlices`,
   `solverConfigs`) — list them. Any other added/removed/changed key → do
   NOT save that file; report it (STOP for that file).
3. Run it on every committed `motors/*.ric` in `C:/tmp/sv/ig`; commit the
   refreshed files on `ig/provider`.
4. Gates: B0 PASS after the refresh; Chunc equivalence — `run_from_ric`
   with `pyrogen=None` vs `pyrogen='bpnv', pyrogen_mass=None` (canonical
   knobs, `t_max=0.05`) gives identical P_head (max |ΔP| = 0, same steps);
   no `StaleIgniterMaterialWarning` / `IgniterMaterialUpgradedWarning` when
   loading any refreshed committed file; srm_1d suite = I1 final + 0 new
   failures.
