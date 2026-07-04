# 05 — I/O & openMotor Integration

How motors get **in** (`.ric` → sim), how results get **out** (dict → CSV /
plots / channels), and how srm_1d plugs into the **openMotor** GUI. This is the
"edges" tier of the code map — none of it is in the hot loop, but it's how a
user actually drives the solver.

---

## Part A — Loading a motor: the `.ric` adapter

[`openmotor_adapter.py`](../../srm_1d/openmotor_adapter.py) is the boundary
between openMotor's world and srm_1d's. A **`.ric` file is openMotor's motor
format** (a YAML dict `{nozzle, propellant, grains, config}`). The pipeline:

```
load_ric(path)                       # parse YAML (handles openMotor's Python tags)
convert_propellant(ric_prop, …)      # → srm_1d Propellant  (units converted here)
convert_geometry(ric_grains, …)      # → MotorGeometry     (routes to build_snapped_geometry)
convert_nozzle(ric_nozzle)           # → Nozzle
build_pyrogen_chamber(…)             # → PyrogenChamber (Sutton-default sizing)
  └─ run_simulation(…)               # the transient solve (01–04)
compute_motor_performance(…)         # thrust/Isp (nozzle.py)
```
`run_from_ric(path, …)` wraps all of it — the one-call entry point examples use.

Two conventions you'll see enforced here (they're project rules):

- **Defer to openMotor's data model; units are the exception.** New data
  structures mirror openMotor's field names/semantics, but srm_1d keeps
  human-readable engineering units internally (e.g. `erosion_coeff` in
  μm/(s·MPa)); the adapter converts at this boundary (e.g. openMotor's
  m/(s·Pa) → μm/(s·MPa), MW g/mol → kg/mol).
- **`.ric` files are openMotor's to author** — regenerate them from openMotor,
  never hand-edit. (`feedback_ric_files_openmotor_owned`.)

**Gas transport (a v0.8.0 subtlety):** the combustion-gas transport
(`μ, k, Cp, Pr`, γ) lives **per-tab in the `.ric` propellant** now (selected by
`transportVariant` = `frozen` or `effective`); the old `.transport.yaml`
sidecars are retired. If transport is missing, the adapter **hard-faults rather
than fabricating it** (`feedback_defer_to_thermochem_solvers`) — srm_1d never
invents Cp/k/μ, because the erosive chain (`02`) is sensitive to them. The
fired motors use **effective** RPA transport.

---

## Part B — Results out

`run_simulation` returns a **plain results dict** — this is the source of truth
and the clean inspection surface (`01` §10). Downstream consumers only *read*
it:

- **CSV** — `result_to_csv` / `save_csv` produce an openMotor-compatible CSV
  (`load_openmotor_csv` reads openMotor's back for comparison).
- **Plots** — [`plotting.py`](../../srm_1d/plotting.py) (matplotlib only; it
  takes the dict, no solver deps): `plot_pressure` (with experimental overlay +
  `time_offset` alignment), `plot_thrust`, `plot_flow_snapshot` (2×2 axial
  fields), `plot_summary`. Plus the slice/station axial viewers added in v0.8.0.

---

## Part C — The channels backbone

[`channels.py`](../../srm_1d/channels.py) (v0.8.0) is the migration layer from
"results dict" to openMotor-native **channels**. Two types
(`build_channels` **reshapes the dict without recomputing** — results stay
byte-identical):

- **`Channel`** — a scalar-per-step or per-grain-per-step series, mirroring
  openMotor's `LogChannel` (`name`, `unit`, `getData(unit)`, `getMax/Min/
  Last/Average`). Covers head pressure, thrust, Kn, per-grain mass/regression,
  etc.
- **`AxialChannel`** — the **srm_1d-specific extension**: a per-**cell** axial
  field over time (`time × N_cells`) plus the axial coordinate `x_cells`.
  openMotor has no per-cell concept, so the GUI ignores these until the axial
  visualization work consumes them. This is how srm_1d's unique 1-D spatial
  data is exposed without breaking openMotor's channel API.

Unit conversion defers to openMotor's `motorlib.units` (imported lazily, so the
core stays importable without the openMotor checkout for raw access). The dict
remains the source of truth; consumers migrate onto channels incrementally.

---

## Part D — The openMotor solver plugin

[`srm1d_plugin.py`](../../srm_1d/srm1d_plugin.py) (v0.8.0) makes srm_1d's PISO
transient solver selectable **inside the openMotor GUI**, alongside openMotor's
own quasi-steady solver. Importing the module **registers** the solver against
openMotor's `motorlib.solvers` registry. It:

1. consumes an openMotor `Motor` via `Motor.getDict()` (the same
   `{nozzle, propellant, grains, config}` shape the `.ric` adapter converts,
   now carrying per-tab transport);
2. runs `run_simulation`;
3. maps the results back into an openMotor `SimulationResult` so the GUI renders
   it natively — scalar per-step channels (time/pressure/force/Kn/…) and
   per-grain channels (mass/massFlow/massFlux/regression/web/mach), the latter
   aggregated from srm_1d's per-cell axial snapshots via the cell→segment map.

The igniter is read from the motor's own `data.igniter` block (falling back to
built-in BPNV for older motors). One GUI-specific wrinkle: the transient solver
emits 10⁵⁺ samples per burn, which would hang openMotor's graph widget, so the
plugin **decimates the GUI channels to `GUI_MAX_POINTS` (5000)** while
preserving the peak-P / peak-thrust / endpoint samples — a *display* artifact
only; srm_1d's own full-resolution channels (CSV/analysis) are untouched.

**Architecture context ("one core, two frontends"):** srm_1d is the canonical
core and single source of truth; the openMotor **fork** vendors it (git
subtree) so the same solver runs headless (here) and in the GUI. Generic hooks
(the solver registry, GUI picker) are what go upstream; the vendoring + heavy
deps stay in the fork. See `project_openmotor_fork_integration_strategy` and
`docs/v0_8_0/`.

---

## Part E — Where things live (dev-only siblings)

Not part of the importable package, but part of the repo you'll work in:

| Dir | What |
|---|---|
| `motors/` | Motor data: `<motor>.ric` (transport + igniter embedded) |
| `examples/` | Runnable scripts (`python -m examples.<name>`), incl. the LHS sweeps |
| `tests/` | pytest suite (~406) — run after every `@njit` edit |
| `srm_1d/tools/` | LHS sensitivity sweeps + ignition diagnostics (ships as `srm_1d.tools`) |
| `docs/` | design docs per version + this guide |
| `static_fire_data/` | experimental traces for validation |

---

## You've reached the end of the guide

You now have: the sim core (`01`), the burn-rate + geometry feedback (`02`),
ignition + sources (`03`), the opt-in add-ins (`04`), and I/O + openMotor
(`05`). From here:

- **Function-level reference:** [`../../srm_1d/ARCHITECTURE.md`](../../srm_1d/ARCHITECTURE.md).
- **Gotchas, calibration state, API-break log:** [`../../srm_1d/DEVNOTES.md`](../../srm_1d/DEVNOTES.md).
- **Why each subsystem is the way it is:** the `docs/v0_7_x/`, `docs/v0_8_0/`,
  and `docs/core_loop_opt/` design packages — read the one matching your change.
- **Before a PR:** clear `srm_1d/__pycache__/` after `@njit` edits, run
  `python -m pytest tests/`, and match the conventions in the guide README.
