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
- Igniter momentum was deliberately deferred in the Phase 3 commit; it
  was reopened and implemented/audited in the Phase 4 work below.

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
- Added ignition-spike diagnostic tooling:
  - opt-in startup controls for ambient initial gas, erosive disable,
    end-face disable, momentum disable, pyrogen surface-heating disable,
    and adjacent-radiation disable.
  - snapshot channels for mass source, thermal source, momentum source,
    pyrogen surface heat flux, and adjacent-radiation heat flux.
  - reducer/classifier plus CLI artifacts under
    `artifacts/ignition_diagnostics/<case>/<variant>/`.

Completed Phase 4 solver work in the current working tree:

- Replaced the temporary nozzle ambient clamp with a signed isentropic
  open-throat boundary. The helper supports subsonic outflow, choked
  outflow, balanced flow, subsonic ambient inflow, and choked ambient
  inflow, and is shared by PISO, energy fluxes, mass-flow history, and
  diagnostics.
- Added `ambient_temperature` to `run_simulation`; `None` defaults to
  `propellant.T_initial` for reverse-flow nozzle inflow.
- Added DeMar pyrogen direct surface heating. Built-in BPNV and MTV
  include heat-flux data; custom pyrogens hard-fault unless surface
  heating is explicitly disabled.
- Added pyrogen axial momentum as an explicit face source plus a
  momentum ledger comparing expected `mdot_ig*v_exit` force to deposited
  force. Hasegawa A baseline and no-momentum runs remain nearly
  identical, so momentum is not the current Hasegawa driver.
- Added gas/solid energy audit histories for pyrogen enthalpy, direct
  surface heat, gas sink, adjacent-radiation heat, nozzle enthalpy,
  thermal-source power, and per-step residual.
- Added adjacent-burning-cell radiation for ignition spread using
  `Propellant.radiation_emissivity` as a material property. The default
  is now 0.0 (opt-in) for all propellants; aluminized `.ric` files no
  longer auto-default to 0.45. Explicit `radiation_emissivity` overrides
  in the .ric (or directly on `Propellant`) are honored.

Current finding:

- Historical hot-fill baseline still ignites essentially instantly and
  retains the spike/erosive-snap behavior.
- The earlier `ambient_initial_gas` + adjacent-radiation combination
  blew up: a `tau_establishment` ramp sweep (2026-05-14) showed
  non-monotonic stability (`tau = 0`, `1 ms` crash to `Mach ~1.9e+05`;
  `0.1`, `0.5`, `5`, `10 ms` "stable" but with interior `Mach 9-50`
  and `P_peak ~14-15 MPa`) -- a discrete PISO/throat resonance, not a
  physical ramp. With radiation now opt-in, the ambient case reduces
  to the previously-stable convective-only spread (~40 ms, ~6.4 MPa)
  observed earlier in `ambient_no_radiation`.
- `ambient_no_surface_heating` remains a no-spread degenerate case,
  confirming the pyrogen-surface-heating path is active.
- Baseline and `no_momentum` Hasegawa A traces are nearly identical in
  the latest smoke run, so momentum is implemented and audited but not a
  dominant effect for this case.
- `tau_establishment` is left in `run_simulation` as an opt-in kwarg
  (default `0.0`, no physical ramp) for diagnostics, but is not used
  for calibration (per `feedback_no_unfounded_smoothing`).

Pending:

- Inspect the new energy/momentum audit CSVs and pressure/x-t plots for
  Hasegawa A before deciding whether a post-ignition burn-establishment
  model is still needed for hot-fill baseline calibration.
- Re-run segmented Hasegawa A LHS after the boundary, direct-heating,
  and adjacent-radiation model changes.
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

- Squib stage (electric match to pyrogen).
- Tuned lumped `C_hc` radiation/heat-transfer multipliers.
- Physical igniter impingement regions unless tied to actual igniter
  basket/jet geometry.
- Multi-species/passive-scalar igniter gas transport.
- Head-end primary motor architecture.
- Cavallini-style multi-species Godunov solver.
