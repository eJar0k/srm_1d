# srm_1d v0.7.0 -- Phase Status

The original implementation checklist is preserved in the phase commits:

- `80d00fa` -- v0.7.0 design package and references.
- `c3cad25` -- Phase 1 standalone pyrogen plenum.
- `6e6b367` -- Phase 2 Goodman solid thermal solver.
- `613ae5f` -- Phase 3 main-loop pyrogen + Goodman integration.
- `e58657a` -- Phase 3 status docs and segmented LHS diagnostics.
- `4fc45d1` -- inhibited-interface gap fix plus BALLSstick/L3035 examples.
- `883e1fb` -- snapped-interface cell assignment fix.

This file now tracks the current project state and the remaining work
before a `v0.7.0` tag.

## Phase 1 -- Pyrogen Plenum

Status: complete and committed.

Implemented:

- `srm_1d/propellant.py`: `Pyrogen` dataclass.
- `srm_1d/igniter_plenum.py`: `PyrogenChamber`,
  `initial_plenum_state`, `_step_plenum_ode`, choked/subsonic orifice
  flow, and Sutton sizing helpers through the adapter builder.
- `srm_1d/motors/pyrogens/bpnv.yaml` and `mtv.yaml`.
- `srm_1d/tests/test_igniter_plenum.py`.

## Phase 2 -- Goodman Solid Heating

Status: complete and committed.

Implemented:

- `srm_1d/solid_thermal.py`: Goodman integral surface heating,
  penetration-depth RK4 step, and ignition threshold helper.
- `Propellant.k_solid`, default `0.3 W/(m*K)`.
- `srm_1d/tests/test_solid_thermal.py`.

## Phase 3 -- Main Solver Integration

Status: complete and committed in `613ae5f`.

Implemented decisions:

- The v0.6.0 exponential igniter API is removed. `igniter_mass`,
  `igniter_tau`, `ignition_ramp_tau`, and `P_ignition` are no longer
  accepted.
- `run_simulation` requires `pyrogen_chamber`; `run_from_ric` requires
  `pyrogen=...` or a sibling `<motor>.pyrogen.yaml`.
- The pyrogen plenum is integrated each timestep and injects hot gas
  into cell 0 only.
- Histories include `P_ig`, `T_ig`, `mdot_ig`, and `m_pyrogen`.
- Snapshots include `T_surf` and `is_burning`.
- `piso_step` uses separate `mass_source` and `thermal_source` arrays.
  Propellant/end-face sources use `T_flame`; pyrogen source uses `T_ig`.
- Igniter momentum is deliberately deferred.

Tests:

- Phase 3 integration tests cover pyrogen loading/discovery, missing
  pyrogen errors, removed legacy kwargs, pyrogen-driven ignition,
  thermal-source temperature handling, and snapshot ignition state.

## Phase 4 -- Validation and Diagnostics

Status: in progress.

Completed pre-work:

- Added segmented pressure metrics for Hasegawa A diagnostics:
  spike, post-spike shoulder, plateau, and taildown.
- Added quiet LHS progress modes: `brief`, `verbose`, and `none`.
- Added `BALLSstick` and `L3035` `.ric` example motors with sibling
  transport YAMLs and simple plotting examples.
- Moved current generated outputs under `artifacts/`; `artifacts/` is
  git-ignored. Root-level generated plot/CSV/LHS output files were
  cleaned on 2026-05-10.
- Fixed `.ric` geometry conversion for bonded inhibited interfaces:
  default inter-segment gaps are inserted only when at least one
  interface face is uninhibited.
- Fixed snapped-interface cell assignment: setup now assigns cells to
  the segment with the largest axial overlap, avoiding epsilon-overlap
  misclassification at touching segment boundaries.

Current finding:

- Hasegawa A LHS runs can tune shoulder, plateau, and taildown, but the
  spike segment remains the limiting residual.
- All grain cells become active almost immediately after Goodman surface
  ignition.
- The next model should add post-ignition burn establishment /
  participation before full propellant mass and thermal source are
  applied per cell.

Pending:

- Implement and validate post-ignition burn establishment.
- Re-run segmented Hasegawa A LHS after that model change.
- Re-run Zerox LHS with v0.7.0 pyrogen parameters.
- Update Hasegawa and Zerox calibration tables in `DEVNOTES.md`.
- Revisit L3035/BALLSstick after any geometry or ignition-model change;
  current examples are exploratory and not calibrated predictions.

## Phase 5 -- Release

Status: not started.

Do not tag `v0.7.0` until Phase 4 is complete and pytest is green.

Release checklist:

- Clear `srm_1d/__pycache__/` after any `@njit` edits.
- Run:

```powershell
C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe -m pytest srm_1d/tests/
```

- Confirm docs match the current API and validation state.
- Tag only after approval:

```powershell
git tag -a v0.7.0 -m "Hot-gas plenum igniter model -- replaces v0.6.0 exponential-decay placeholder"
```

Agent memory note: external agent memory is not the source of truth for
this repository. Preserve durable decisions in committed project docs
and git commits. Update external memories separately only when explicitly
requested.

## Deferred Beyond v0.7.0

- Igniter momentum source terms, unless validation shows they are needed.
- Squib stage (electric match to pyrogen).
- Lumped radiation.
- Multi-species/passive-scalar igniter gas transport.
- Head-end primary motor architecture.
- Cavallini-style multi-species Godunov solver.
