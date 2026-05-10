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

Workflow: pytest before each commit, never push without approval, never
force-push, never use `--no-verify`. Clear `srm_1d/__pycache__/` after
`@njit` edits to prevent stale Numba cache behavior.

Next modeling work is not "more igniter smoothing." Current Hasegawa
diagnostics point to post-ignition burn establishment: after Goodman
surface ignition, newly ignited cells likely need a finite participation
ramp before they contribute full mass and heat source.
