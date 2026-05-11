Read these in order before starting work:

1. `CLAUDE.md` -- session-start orientation and current roadmap.
2. `srm_1d/DEVNOTES.md` -- gotchas, calibration state, and API breaks.
3. `srm_1d/ARCHITECTURE.md` -- function-level map of the current code.
4. `srm_1d/docs/v0_7_0/DESIGN.md` and `TASKS.md` -- implemented
   v0.7.0 design plus remaining validation work.

Current branch is `v0.7.0-phase3`. Phases 1-3 are committed:

- `c3cad25` -- standalone pyrogen plenum.
- `6e6b367` -- Goodman solid thermal solver.
- `613ae5f` -- main-loop pyrogen + Goodman integration.
- `e58657a` -- docs plus segmented LHS diagnostics.
- `4fc45d1` -- L3035/BALLSstick examples and inhibited-interface gap fix.
- `883e1fb` -- snapped-interface cell assignment fix.
- Current working tree adds Phase 4 ignition-transient work: signed
  nozzle open boundary, DeMar pyrogen surface heating, pyrogen momentum
  ledger, gas/solid energy ledger, adjacent-cell radiation, and expanded
  ignition diagnostics.

Workflow: pytest before each commit, never push without approval, never
force-push, never use `--no-verify`. Clear `srm_1d/__pycache__/` after
`@njit` edits to prevent stale Numba cache behavior.

Next modeling work is not "more igniter smoothing." With ambient initial
gas, Hasegawa A now needs both DeMar pyrogen surface heating and
adjacent-cell radiation to spread. Review the new energy/momentum audit
CSVs and x-t plots before deciding whether post-ignition burn
establishment is still needed for the historical hot-fill baseline.

Keep generated plots, CSVs, and LHS artifacts under `artifacts/`.
Root-level generated outputs were cleaned on 2026-05-10. Do not touch
`Zerox Data.xlsx` or local example-script verbosity edits unless asked.
