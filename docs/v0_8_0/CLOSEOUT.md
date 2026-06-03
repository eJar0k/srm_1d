# v0.8.0 — openMotor Frontend Integration: Close-out

> Companion to [`DESIGN.md`](DESIGN.md) (scoped 2026-06-01) and
> [`TASKS.md`](TASKS.md) (phase plan). Branch `openmotor-frontend`.
> This doc records what v0.8.0 actually delivered and what remains.
> **Not yet tagged** — the tag gate (v0.7.5 base) is still open; see below.

## What shipped

v0.8.0 is the **data-model + channels backbone + plugin contract** that lets
openMotor's GUI run srm_1d's 1-D PISO transient solver alongside its built-in
quasi-steady (QS) solver, on the same `Motor`, rendered in the same widgets.

| Phase | Delivered |
|-------|-----------|
| 1 | **Channel result model** (`srm_1d/channels.py`): `Channel` (scalar/per-grain, openMotor `LogChannel`-aligned), `AxialChannel` (time×N_cells), `SimulationChannels` container. `run_simulation` returns it; it **proxies dict access** so all legacy `result[key]` code is unchanged. |
| 2 | **Unit-aware + generic plotting** (`plot_channels`, flow-snapshot/heatmap over `AxialChannel`); pressure/thrust/summary routed through the generic path. |
| 3 | **Transport-in-propellant** (per-tab `mu` + frozen/effective k/Cp; `transportVariant`, default `frozen`). Transport sidecar **retired** in `run_from_ric`; D7 hard-fault on the sentinel. Migrator ran all repo motors (0 sentinels). |
| 4 | **Igniter as data** (`data.igniter` = embedded pyrogen material + chamber sizing). `run_from_ric` reads the motor's own igniter; explicit `pyrogen=` overrides. |
| 5 | **`motorlib` solver-plugin contract** (`solvers.py` registry + `SolverPlugin`; QS wrapped for coexistence). `srm1d_plugin.py` registers `Srm1dTransientSolver`, maps channels into an openMotor `SimulationResult`. |
| 6 | **GUI integration** (below) + per-grain channels + Motor igniter wiring + docs. |

## Phase 6 — GUI integration (the headline deliverable)

The QS and srm_1d-transient solvers now **coexist in the openMotor GUI**:

- **Solver picker + routing.** `motorlib/solvers.discover_external_solvers()`
  best-effort imports the srm_1d plugin at startup (sibling-checkout walk;
  honors `OPENMOTOR_SRM1D_PATH`; swallows missing-dep failures).
  `SimulationManager` carries an `activeSolverName` and routes `_simThread`
  through the registry (try/except → a solver fault surfaces as an alert
  instead of hanging the progress dialog). `mainWindow` adds a **Simulate →
  Solver** submenu (shown only when >1 solver is registered).
- **Result rendering (coexistence).** The plugin returns a fully-populated
  openMotor `SimulationResult` — scalar **and** per-grain channels
  (`mass/massFlow/massFlux/regression/web/machNumber`, aggregated from the
  per-cell axial snapshots via the `cell_segment_id` map and interpolated onto
  the per-step time grid) plus `volumeLoading` — so the existing graph and
  motor-stats widgets render the transient trace with **no widget changes**.
- **Live progress + Stop button.** The `@njit` time loop is `nogil=True` and
  publishes a progress metric into a shared array each step; a worker-thread
  poller forwards it to openMotor's callback and sets a cancel flag on Stop
  (new `termination_code 5`). Progress metric: web-consumed fraction through
  the burn, then a steady wall-clock fill through the pressure taildown (the
  raw physics metric asymptotes — see Known limitations).
- **Graph-switch performance.** The transient solver emits 10⁵+ per-step
  samples; the GUI `SimulationResult` is **decimated to 5000 points**
  (peak-pressure / peak-thrust / endpoints preserved) so channel switches are
  responsive. srm_1d's own full-resolution channels are untouched.

**User-verified in the running app** (Hasegawa C): solver picker lists both
solvers; transient pressure/thrust + stats render; progress bar advances and
completes; Stop cancels; graph switching is responsive.

## Running openMotor from source (dev / verification)

openMotor is normally a frozen app; GUI changes are tested from the modified
checkout. The pyenv 3.10.5 was provisioned as the single interpreter with
both PyQt6 **and** numba+scikit-fmm:

```
pip install PyQt6==6.7.1 PyQt6-sip platformdirs ezdxf docopt cython wheel
# (numba / scikit-fmm / scipy / matplotlib already present)
python setup.py build_ext --inplace          # compile mathlib._find_perimeter_cy (MSVC)
for f in uilib/views/forms/*.ui:             # generate pyuic6 UI modules
    python -m PyQt6.uic.pyuic -o uilib/views/<stem>_ui.py <f>
python main.py                                # launch
```

Build artifacts (`*_ui.py`, `*.pyd`, the cythonized `.c`) are git-ignored.
**Gotcha:** the plugin only registers if the openMotor runtime Python has
srm_1d's deps; otherwise the Solver submenu silently shows QS only.

## Test status

Full srm_1d suite green — **337 passed**, including 10 plugin tests
(`tests/test_srm1d_plugin.py`): registry selection, per-grain channels,
volumeLoading, igniter round-trip (with/without block), live progress to
completion, cooperative cancel, and decimation peak-preservation. The
progress-metric change touches only `progress_state[0]` (write-only, never
read into physics), so simulation outputs are byte-unchanged.

## Known limitations / deferred

- **Progress-bar taildown is heuristic.** Physics metrics asymptote near the
  end (measured: ~20% of wall time in the last 2% of the bar); the fix fills
  the tail at a steady wall-clock rate (to ~0.96, then snaps to 100% at
  termination). **Deferred TODO:** anchor both the bar and the reported
  burn-time/impulse window on an **industry burn-time definition**
  (tangent-bisector or %-of-max-pressure/thrust per NFPA 1125 / ISO / amateur
  5%–5% · 10%–90%). See TASKS.md Deferred.
- **~2 s post-loop result-build pause at 100%** (per-grain interpolation +
  decimation in the main thread) — separate from the bar, not yet addressed.
- **Capability-gated srm_1d-only panels** (axial viz / igniter editor /
  transport editor surfaced from `capabilities`) deferred to a v0.8.x GUI pass.
- Tapering finocyl grains, RocketCEA transport, rich burnback viz — see
  DESIGN Roadmap / TASKS Deferred.

## Tag gate (still open)

**Do not tag v0.8.0 yet.** Cut it only from a base that already contains
**v0.7.5** (the cross-motor re-LHS), then merge/rebase v0.7.5 into
`openmotor-frontend` so the frontend inherits the recalibrated model. The
re-LHS is STAGED but HELD (user pausing CPU): worktree `../srm_1d-v075-lhs`,
`examples/cross_motor_lhs_v075.py`, N=1000, frozen + F+Z + κ_zn=1. See
TASKS.md "Cross-line sync."
