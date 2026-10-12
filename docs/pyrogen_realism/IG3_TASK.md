# IG3 — defaults: Peretz ignition h, drop the DeMar flux, DeMar igniter sizing

**User decisions 2026-10-11** (after IG2, `IG2_REPORT.md`):
1. Make `ignition_heat_transfer='peretz'` the default **now**.
2. Port-area reading for Peretz's (Ap/At)^0.4: run the test. Compare the local
   cell area (IG2's choice), the head-end port area and the aft-end port area.
   The conical Chunc is the key test case. The user doubts the aft-end
   reading: an arbitrarily large aft port over a short section of an otherwise
   narrow motor would dominate it. The user picks after seeing the numbers.
3. Drop the DeMar flux, and replace the Sutton auto charge with **DeMar's
   sizing rule at 300 psi** plus an **auto-upsize fallback** (approved).

Two stages. **Stop after Stage 1 and report**; Stage 2 starts on a message.

## Stage 1 — Peretz default + Ap-reference test

1. Default `ignition_heat_transfer='peretz'` (`run_simulation`, and every
   wrapper that carries its own default). Expose it in the openMotor plugin's
   solver config as a choice field: the plugin already has choice fields, e.g.
   `heat_delivery_mode`. Config defaults are generated from
   `inspect.signature(run_simulation)` (I4).
2. New opt-in kwarg `peretz_ap_reference = 'local' (default) | 'head' | 'aft'`.
   It only matters in peretz mode; `'local'` is bit-identical to IG2.
   - `head` = port area of the fore-most grain cell / throat area.
   - `aft` = port area of the aft-most grain cell / throat area.
   - Instantaneous values, like `local`.
   - Report Chunc's head and aft Ap/At.
3. Check the Peretz constant. The thesis (`peretz1974.pdf`) Eq. III-31 gives
   h = 1.56e-3 ... (from a6 = 2.56e-2, p. 40). AIAA Paper 72-1119
   (`peretz1972.pdf`) Eq. 14 prints 1.78e-3 with possibly different
   variables. Read both from page images: render single pages with `fitz` at
   <= 110 dpi, cropped; never Read a whole PDF. Say whether they are the same
   law. **Do not change the constant** without reporting it; the thesis carries
   the derivation.
4. Re-record the B0 `BASELINE` in
   `docs/solver_verification/probes/baseline_gate.py`, keeping the old values in
   comments as before. Update tests pinned to old values, with the reason. Add
   tests for the Ap-reference variants. Add a DEVNOTES API-log entry (default
   change) and a calibration note: the v0.7.5 knobs were fit with the old
   ignition h.
5. Measure local / head / aft:
   - (a) Canonical motors with their usual igniters (Chunc, BALLSstick,
     Hasegawa A; Sutton charge, forward plenum, DeMar step as today): first
     ignition, peak, t_peak, plateau and ratios vs experiment, as in IG2 4.3.
   - (b) Chunc head basket, realistic rows (4.5 / 6 g x `bpnv` / `mtv`,
     `particle_diameter_m=3e-3`, `particle_LD_ratio=1`), with
     `heat_delivery_mode='none'`, as in IG2 4.2.
   - (c) Ignition thresholds, head basket and forward plenum, DeMar off, as in
     IG2 4.1, plus the per-motor forward-plenum thresholds as in the IG2 4.3
     table.
6. Gates: B0 PASS against the new baseline; full suite green. Commit. Report.

## Stage 2 — drop DeMar + DeMar sizing + fallback (on message)

1. DeMar flux off by default. Set `heat_delivery_mode: none` in every library
   YAML in `srm_1d/pyrogens/`; this includes the thermite, where DeMar's own
   conclusion 4 is that thermites do not work as self-contained igniters. Also
   change the `Pyrogen` dataclass and `_pyrogen_from_mapping` defaults. Keep
   `'demar'` selectable.
2. Embedded `.ric` igniter blocks: add the pre-change library entries to
   `SUPERSEDED_PYROGENS`, so old blocks upgrade on load. **Do not modify
   `motors/*.ric`.** The user has uncommitted edits to two of them in the main
   checkout. If a test needs refreshed `.ric` files, report it instead.
3. Sizing: add `demar_pyrogen_mass(V_free_in3, impetus_W, P_target_psi=300.0)`,
   returning kg. DeMar's slides give M[g] = V[in^3] P[psi] / W[psi in^3/g].
   - V = free volume (case minus propellant), the same V the Sutton default
     uses now.
   - W comes from `Pyrogen.impetus_W`. If a library pyrogen lacks it, raise a
     clear error; do not invent a value.
   - Make it the auto default wherever Sutton is used today (`build_pyrogen_chamber`
     and any plugin/GUI auto-charge path). Keep `sutton_pyrogen_mass`.
   - Unit tests: DeMar's three slide examples (2.36 in^3, 300 psi -> 0.14 g;
     17.7 in^3, 200 psi -> 0.71 g; 103.7 in^3, 400 psi -> 8.3 g).
4. Fallback, **auto-sized charges only**:
   - Run. If the motor does not light, stop early and retry at 1.25x the
     charge. Cap the attempts, then raise a clear error.
   - The early-stop test must be physical, not a tuned timeout. For example:
     no grain cell ignited, the pyrogen fully consumed, and every unignited
     surface temperature falling, so nothing can light any more. Propose a
     better first-principles test if you have one.
   - Report the final charge and the attempt count in the result summary, and
     warn naming the charge.
   - An explicitly set charge is never resized. If it fails to light, stop
     early, return the result marked not ignited, and warn.
5. Re-record B0 (the canonical motors now use DeMar-sized charges plus the
   fallback). Report each motor's final charge, attempts and wall-time cost.
6. Gates; commit; report.

## Rules

- No tuned constants. 300 psi and the 1.25 step are the user's approved
  choices; everything else comes from a source or first principles.
- Paywalled-reference rule: list needs in `PAYWALLED_REFERENCES.md`; mark
  dependent conclusions "pending".
- The Ma 2020 erosive path is untouched.
