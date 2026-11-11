# Fable task 2: blind audit of the numerical core

**Run this after task 1, and after AS1 (prompt-start auto charge), G2 (grain
snapping and lumping) and the nozzle-boundary stagnation fix have merged.**
The audit should see the code as it will be, not code that is about to
change.

## Why blind

The auditor gets the code, a neutral statement of the intended physics, and
the primary references. It does not get our docs, worklogs, DEVNOTES,
CLAUDE.md, memory or git history. Those would hand it our assumptions, such as
"the corrector bug is fixed" or "this term is exact", and our sense of what
was already checked. The code's own comments still carry some of that, so its
instructions tell it to treat every comment as an unverified claim.

## What goes in the bundle

`make_core_audit_bundle.py` builds a folder outside the repo, as a fresh git
repo with one commit and no history.

| Path | Contents |
|---|---|
| `srm_1d/` | The package at the chosen git ref, minus every `.md` file, the openMotor plugin and the dev diagnostic tools (`guarded_run` stays) |
| `motors/` | `hasegawa_a.ric`, `machbusterNew.ric`, `example_bates.ric`, all committed versions. They are analytic grains, so no openMotor checkout is needed. BALLSstick is left out because its Finocyl grains need one. |
| `references/` | Ma 2020 (erosive), Gnielinski 2013, Peretz AMS 1100 (1973) and AIAA 72-1119, MIL-P-46994B |
| `SPEC.md` | The intended physics (`docs/fable/core_audit/SPEC.md`) |
| `TASK.md` | The auditor's instructions (`docs/fable/core_audit/TASK.md`) |
| `probes/` | Empty; the auditor's probes go here |

Deliberately left out: static-fire data, the validation papers (Hasegawa) and
every account of the ignition spike. This is verification against physics,
not validation against data. Task 1 covers validation.

## Steps

1. Build the bundle from the repo root:
   ```
   "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" docs/fable/make_core_audit_bundle.py --out C:/tmp/fable_core_audit
   ```
   - `--ref` defaults to `HEAD`; pass the merged commit if you are on another
     branch.
   - The script refuses to write into a folder that is not empty.
2. Open a **new** Claude Code session with `C:/tmp/fable_core_audit` as the
   working folder. Do not open it inside the srm_1d repo: the repo's CLAUDE.md
   and the project memory would load and break the blindness. Choose Fable
   5.2.
3. Say: **"Read TASK.md and SPEC.md, then do the task."**
4. If the session hits a usage limit, start another session in the same
   folder and say: "Read TASK.md, then resume from the checklist in
   AUDIT_REPORT.md."
5. When it is done, copy `AUDIT_REPORT.md` (and any `probes/` worth keeping)
   into `docs/fable/core_audit/` in the repo. Bring it to a normal session
   for triage against the record. Many findings will already be known, but
   the ones that are not are the point.

## Notes

- **Compute budget:** `TASK.md` carries the coexist limits (4 GB, 2
  processes, 30 min, guarded) and tells the auditor to stop and ask before
  heavier work.
- **Fixes:** the auditor may demonstrate fixes only on `fix/<slug>` branches
  inside the bundle repo; it never touches the real repo.
- **Rebuilding:** delete the old bundle folder yourself first, or choose a new
  `--out`.
