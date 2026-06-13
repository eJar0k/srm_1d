# Upstream openMotor PR(s) for grain tapering — scope / cold-start

> Scope check done 2026-06-12. Records the concrete git facts so a fresh
> session can plan + execute the upstream extraction without re-deriving them.
> Work happens in the **openMotor fork** (`.../openMotor/openMotor`), not srm_1d.
> Line/commit numbers will drift — re-verify with the commands at the bottom.

## Goal

The parametric grain-taper feature is generic (no srm_1d deps), so it fits the
fork strategy (`[[project_openmotor_fork_integration_strategy]]`: upstream PRs =
generic hooks only, NOT the srm_1d vendoring / heavy deps). Upstream it to
`reilleya/openMotor` as 1–2 PRs so the fork can later drop its private copy on
the next sync.

## Base branch — `upstream/staging` (NOT master)

- `upstream/master` is **114 commits behind `upstream/staging`** — stale. The
  earlier plan said "off upstream/master"; that was wrong. Use **`upstream/staging`**
  (it is `upstream/HEAD`, i.e. reilleya's active/default branch).
- The fork's `staging` is **0 behind / 25 ahead** of `upstream/staging` — it
  branched from current staging and added 25 commits (vendoring + file-format +
  plugin + igniter + station-viz + taper). So the taper delta vs `upstream/staging`
  is small and clean (below), NOT the inflated ~1400-line diff vs master.

## The generic taper feature = 3 contiguous fork commits (most recent)

In chronological order (all on fork `staging`, already pushed):

| Commit | Round | What (generic, upstream-able) |
|---|---|---|
| `976e3a9` | 2 | Core: `motorlib/taper.py` + `TaperProperty` (properties.py) + base-`Grain` taper/`isTaperable` (grain.py) + `MotorConfig.taperSlices` + `Motor.runSimulation` expansion (motor.py) + `SimulationResult` snapshot (simResult.py) + `test/unit/taper.py` |
| `359923d` | 3 | GUI bore-taper editor (motorEditor.py) + slice-averaged area graph (grainPreviewWidget.py `averaged_area_curve`) + QS refinements |
| `24d75bb` | 4 | OD/end taper: `taper['od']` schema + `od_diameter_at` + GUI "End taper (OD)" section + `renderGrainLongitudinal` preview + `MainWindow.ui` `verticalLayout_3` warning fix |

**NOT in the PR:** `d75dec0` (this session — "OD-tapered casing in the slice
viewer") is **srm_1d-coupled** (consumes the `srm1d_axial` payload's
`cell_D_outer`) and lives in the fork-only `motorSliceWidget.py`. It stays
fork-only. srm_1d's round-5 transient OD (`cell_D_outer`) is srm_1d-side and is
not upstreamed either; it merely TRACKS the upstream `taper['od']` schema.

## Taper delta vs `upstream/staging` (diffstat)

```
motorlib/taper.py                   | 387 ++  (NEW — clean add)
motorlib/grain.py                   |  44     (isTaperable + base-Grain taper)
motorlib/motor.py                   |  74     (runSimulation expansion + MotorConfig.taperSlices)
motorlib/properties.py              |  19     (TaperProperty)
motorlib/simResult.py               |  13     (snapshot expanded sub-grains for getPortRatio)
uilib/widgets/grainPreviewWidget.py | 112     (averaged_area_curve)
uilib/widgets/motorEditor.py        | 521     (bore + End/OD taper editor + Longitudinal preview)
```
(`MainWindow.ui` `verticalLayout_3` fix is a small additional standalone hunk.)

The small motorlib numbers mean low entanglement with the fork's other 25
commits — the taper hunks should reapply onto vanilla upstream files with few
conflicts.

## PR split

- **PR1 — core motorlib** (clean / low-conflict): `taper.py` + `TaperProperty`
  + base-`Grain` taper/`isTaperable` + `MotorConfig.taperSlices` +
  `Motor.runSimulation` expansion + `SimulationResult` snapshot +
  `test/unit/taper.py` (registered in `test/unit/__init__.py`). Headless,
  testable, no Qt.
- **PR2 — GUI** editor + preview: `motorEditor.py` bore + End-taper editor,
  `grainPreviewWidget.py` `averaged_area_curve`, the Longitudinal preview.
  **PREREQ:** `renderGrainLongitudinal` currently lives in the fork-only
  `uilib/widgets/motorSliceWidget.py` (confirmed absent on `upstream/staging`).
  Relocate it into a non-fork-only file (e.g. `grainPreviewWidget.py` or a small
  preview module) so PR2 doesn't drag in the slice widget.
- **Trivial standalone PR:** `MainWindow.ui` `verticalLayout_3` removal (silences
  the long-standing boot "QLayout … already has a layout" warning). Independent;
  can go first.

## Extraction method

Branch off `upstream/staging`. Either (a) cherry-pick the contiguous range
`976e3a9^..24d75bb` and resolve the (expected few) conflicts from the fork's
earlier config/pyqt6 changes, then SPLIT by file into PR1/PR2; or (b) reapply
the taper hunks additively onto vanilla files guided by those 3 commits.
Verify NO srm_1d/plugin/vendoring leakage rides along (the generic taper has
none, but `motor.py`/`simResult.py` also carry fork-only edits — keep only the
taper hunks).

## Gates

- oM `test/unit` green via the shim runner: `srm_1d.fmm_grain._setup_openmotor_path()`
  (installs the marching-squares find_perimeter shim — the pyenv env lacks the
  Cython build), then `unittest` over `test/unit` (`test/unit/taper.py` = the
  taper suite). Offscreen GUI smoke via `QT_QPA_PLATFORM=offscreen` + `QT_API=pyqt6`.
- **Before opening each PR, run `/code-review ultra` on the PR branch** — it's a
  deep multi-agent review, exactly right for outward-facing upstream diffs.
  (User-triggered/billed; the agent can't launch it.)

## Caveats if accepted

One-time de-dup on the next upstream sync (the fork drops its taper copy for
upstream's). srm_1d then TRACKS upstream's `taper['od']` schema + `motorlib.taper`
API rather than owning it — an API rename in review would ripple to srm_1d's
adapter (`taper_spec_from_props`, `od_ends_from_taper`) + the round-5 transient
OD (`resolve_taper`, `_fill_od_taper`).

## Re-verify commands (fork repo)

```
git fetch upstream
git rev-list --left-right --count upstream/staging...staging      # expect 0  25
git log --oneline upstream/staging..staging                       # the 25 fork commits
git diff --stat upstream/staging staging -- motorlib/taper.py motorlib/grain.py \
    motorlib/motor.py motorlib/properties.py motorlib/simResult.py \
    uilib/widgets/motorEditor.py uilib/widgets/grainPreviewWidget.py
git ls-tree upstream/staging --name-only uilib/widgets/motorSliceWidget.py  # empty = fork-only
```
