# R1 — review the friend's Codex commits (agent-executable)

The user's friend ran Codex on srm_1d and pushed 15 commits to
`github.com/israel-shero/srm_1d_Ivans_vibes` (`main`, fetched here as remote
`codexfork`). They are based on our public `db0a980` — i.e. BEFORE our
erosive-model registry, solver-verification T1–T6 and the plugin-provided
igniter, all of which are now on `openmotor-frontend`. We replicated their
work from a summary (`DESIGN.md` §1); now we have the code. **This task is a
REVIEW: it produces a report and a port plan, and changes neither repo.**

## Workspaces (read-only on both repos)

| what | path | note |
|---|---|---|
| their code | `C:/tmp/sv/cx` | detached worktree at `codexfork/main` |
| openMotor fork as THEY expect it | `C:/tmp/sv/cx-om` | detached at fork `92e425f` (pre-igniter-provider) |
| our current code | the main srm_1d checkout | READ ONLY — never edit, never run code there |

Never commit, never push, never `git stash`. Write the report to
`C:/tmp/sv_cx/REVIEW.md` (the reviewer commits it). Scratch: `C:/tmp/sv_cx/`.

To RUN their code: from `C:/tmp/sv/cx`, export
`PYTHONPATH="C:/tmp/sv/cx;C:/tmp/sv/cx-om"`,
`SRM1D_OPENMOTOR_PATH="C:/tmp/sv/cx-om"`, `OPENMOTOR_SRM1D_PATH="C:/tmp/sv/cx"`,
use `"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe"`, and
preflight that `srm_1d` and `motorlib` resolve inside those two worktrees.
Their lockfile targets macOS; note any Windows portability problem you hit
rather than fixing it.

## Read first

- Ours: `docs/solver_verification/DESIGN.md`, `TASKS.md`, and worklogs
  `worklog/T1.md`–`T6.md` (what we built and measured);
  `docs/igniter_provider/DESIGN.md` (stale-igniter finding, which also affects
  their legacy-BPNV numbers); `docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` §6 and
  §8 thread F.
- Theirs: `git -C C:/tmp/sv/cx log --reverse db0a980..HEAD` and each diff.

## Deliverable — `REVIEW.md`

1. **Per commit** (all 15, in order): what it does (≤3 lines); correctness
   review (bugs, numerical risks, any change to DEFAULT simulation behaviour —
   say explicitly whether default results change); test quality; overlap with
   our merged work (name the T-task); verdict **PORT** (adapt into our tree),
   **COMPARE** (we have an equivalent — state concretely how theirs differs and
   whether either is better), or **SKIP** (why).
2. **Answer these questions with evidence** (code citations; short runs only
   where stated):
   - Q1 Fill Mach: we measure max fill Mach 3.94 (library BPNV) / 4.81
     (embedded legacy BPNV) at 100 cells on committed Chunc; their summary
     says "current BPNV 1.987/2.200/2.311", "legacy 3.271/4.821/11.722", "MTV
     0.448/0.498/0.479". From `cab26d0`'s script: what configuration and Mach
     METRIC produce their numbers (cell vs face, which window, which igniter
     source, cells)? You may run their fill-comparison for ONE grid (≤10 min)
     to reproduce one number.
   - Q2 Bracket: `e749f71` (dynamic expansion) vs our T2 (analytic bracket
     `r0 + ebf·h0` when F(20·r0) < 0, plus an `r0 == 0` guard). Any state one
     handles and the other does not? Is theirs bit-identical on non-clamped
     states?
   - Q3 Energy: `79867c0` vs our T6 — what does each measure that the other
     does not? Any numerical disagreement (their source 1.30e-8 W / nozzle
     2.05e-8 W / jump max 2.63 J / Σ 13.30 kJ vs ours in `worklog/T6.md`)?
   - Q4 Limit diagnostics: which numerical limits do `d08bcc5`/`d773b55` track
     beyond our T3 pressure-floor counters (T clip, Mach, dt, …)? Port value?
   - Q5 Pressure-floor probes (`86c587f`, `55b19cb`, `2d664ff`, `32d9230`,
     `9ee1a91`, `docs/v0_8_0/PRESSURE_FLOOR_EVENT_PROBE.md`): what did they
     find — in particular, does the "floor-free finer CFL probe" show the
     806-cell floor events are TEMPORAL (vanish at lower CFL) or spatial? How
     does that bear on our thread F (EOS/energy–pressure lag) and T3's −734 Pa
     raw undershoot?
   - Q6 `3a6ff95` case registry / workspace locks / baseline runner: value vs
     our B0 gate (`docs/solver_verification/probes/baseline_gate.py`)?
     Recommend adopt / adapt / skip.
   - Q7 Can our queued T8 (grid convergence to 1600 cells + fill-Mach
     attribution, `TASKS.md`) reuse `run_spatial_extension.py` /
     `run_fill_configuration_comparison.py` instead of writing new probes?
3. **Port plan:** an ordered list of port tasks for our `openmotor-frontend`
   (which already has the registry, T1–T6, T9/T10 in flight, and the igniter
   provider): for each, files, expected conflicts, whether it changes default
   results, rough size, and a gate. Do not port anything yourself.

Keep every claim checkable (commit + file:line or a command + its output).
If something is ambiguous, say so in the report rather than guessing.
