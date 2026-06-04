# Pyrogen library + igniter editor — design note (v0.8.x)

> Scoped 2026-06-03 from user feedback on the first igniter-editor pass. Mirror
> openMotor's **propellant subsystem** for pyrogens (the `IGNITERS=5` file
> section was reserved for exactly this), add capability-gating, conditional
> fields, and the `burn_law` enum fix. openMotor-fork side.

## Model (mirror propellants exactly)

openMotor **embeds** the propellant in each grain; the propellant *library* is
a **picker** (select → copy into the grain). Mirror that for the igniter:

- **Pyrogen** (material) ↔ Propellant. A reusable **pyrogen library** persisted
  via the `IGNITERS` file section. The motor still **embeds** its
  `igniterPyrogen` (self-contained `.ric`); the library copies a chosen pyrogen
  into it.
- **Igniter** (chamber: mass/throat/volume/topology/...) ↔ per-motor, stays
  embedded like the Nozzle. Edited via the grain-table **Igniter** row.

## Pieces (mirror the propellant files)

| Propellant | Pyrogen (new) |
|---|---|
| `uilib/propellantManager.py` `PropellantManager` | `uilib/pyrogenManager.py` `PyrogenManager` |
| `uilib/widgets/propellantMenu.py` `PropellantMenu` (`.ui`-backed) | `uilib/widgets/pyrogenMenu.py` `PyrogenMenu` (**code-built**, no `.ui`) |
| `DEFAULT_PROPELLANTS` (defaults.py) | `DEFAULT_PYROGENS` (defaults.py) |
| `app.py: self.propellantManager = ...` | `app.py: self.pyrogenManager = ...` |
| Edit → Propellant Editor | Edit → **Pyrogen Library** |
| grain propellant combo (copy from library) | igniter-row **Pyrogen picker** (copy library → `motor.igniterPyrogen`) |

`PyrogenManager`: `pyrogens=[]`, `loadPyrogens`/`savePyrogens`
(`fileTypes.IGNITERS`, `<config>/pyrogens.yaml`, fallback `DEFAULT_PYROGENS`),
`getNames`, `getPyrogenByName`, `showMenu`, `setPreferences`. (Menu created
lazily / `QApplication`-guarded so the data layer is headless-testable —
forced deviation from PropellantManager's eager menu.)

`PyrogenMenu` (code-built QDialog): `QListWidget` of names + a `CollectionEditor`
for the selected `Pyrogen` + New/Delete/Edit buttons; mirrors PropellantMenu's
flow (select→edit, new auto-names, delete, save on change). Gives the rigid
edit-tracking/saving the inline editor lacked.

## Other fixes (from the same feedback)

- **`burn_law` → `EnumProperty(['0d','end_burning'])`** (motorlib/igniter.py).
  DONE 2026-06-03 (valid values per `igniter_plenum._burn_law_code`).
- **Capability-gate the grain-table Igniter row:** add an `igniter` capability
  to the solver plugin (`Srm1dTransientSolver.capabilities`); `updateGrainTable`
  shows the Igniter row only when the active solver declares it (hidden in QS).
- **Conditional fields (the "dropdowns don't propagate to corresponding
  inputs" report):** the enum *values* save correctly (verified); the want is
  **field relevance** — `injection_topology` shows plenum (cartridge/throat/
  volume) vs basket (fill/packing) fields; `form` drives the particle inputs
  (diameter/L-D). Implement as show/hide (or enable/disable) in the igniter
  editor's `propertyUpdate`, keyed on the two dropdowns. (oM's flat
  PropertyCollection editor has no built-in conditional fields, so this is a
  small custom layer in `MotorEditor` for the igniter mode.)

## Build slices

1. **Data layer (this slice):** `DEFAULT_PYROGENS` + `PyrogenManager`
   (load/save/getNames/getByName), persistence round-trip verified headless.
2. **Library UI:** `PyrogenMenu` (code-built) + `app.py` wiring + Edit menu
   action. (`burn_law` enum already done.)
3. **Motor integration:** igniter-row Pyrogen picker (copy library →
   `motor.igniterPyrogen`); capability-gate the row.
4. **Conditional fields** in the igniter editor.

The grain-table igniter editor from the first pass stays functional throughout
(edits the embedded pyrogen inline) until slice 3 swaps the material side to the
library picker.
