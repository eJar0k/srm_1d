# AS1 -- prompt-start criterion for the auto igniter charge

Branch `as1/prompt-start` (from `openmotor-frontend` fc7e0de, the G1 geometry). Probe: `docs/pyrogen_realism/probes/as1_valid.py`; raw
results `docs/pyrogen_realism/results/as1/` (`valid.jsonl` = B0-row settings, `valid_cutoff.jsonl` = the GUI / `run_from_ric` default
pressure cutoff 0.5 MPa, `B0_*.txt`, `suite.txt`). Tests: `tests/test_as1_prompt_start.py` (9) + the updated `tests/test_ig4_self_sustain.py`.

## 1. Why

After G1 the DeMar auto charges fell ~25 %, and the IG4 "self-sustaining" test (grain gas generation > igniter peak after burnout) accepts a
hang-fire: BALLSstick's auto 1.268 g lit the grain over ~1 s (t90 971 ms). The user rejected fixed B0 charges ("won't scale against different core
volumes") and asked for a harder auto charge with a prompt-start criterion.

## 2. The criterion and the machinery

- **Success** = at least `ignition_prompt_fraction` = 0.95 of the grain cells ignited within `ignition_prompt_window_s` = 0.5 s of igniter firing (t = 0).
  Both are `run_simulation` kwargs defaulting to the user's numbers (not tuned). A run shorter than the window is judged at its end.
- **Early stop.** For AUTO charges the time loop checks once at the window end (`prompt_early_stop`): fewer than 95 % lit -> termination code 7 ("prompt start not
  met"), the attempt is retried at x1.25. The existing "nothing can light" stop (code 6) is kept. The x1.25 step, the 10-attempt cap and the DeMar 300 psi initial
  guess are unchanged. After the cap: if the grain never sustained combustion `PyrogenAutoSizeError` (as before); if it did light (just not promptly) the last
  attempt is re-run to completion (no early stop) and returned with a warning and `prompt_start = False`.
- **IG4 self-sustaining test: kept as the `summary['ignited']` flag, but it is no longer a retry condition.** Justification: a grain that is >= 95 % lit within 0.5 s
  produces roughly 0.95 of the full-burn gas rate, which exceeds the igniter's peak (a burst of a few ms) on every canonical motor (the two rates stayed ordered
  grain > igniter in all runs below); the IG4 test alone could not catch the slow start, the prompt test subsumes it. It still drives the plugin's no-ignition / hang-fire
  ERROR (a run that never self-sustains), and a warning is issued if the prompt test passes but it fails (never seen).
- **Explicit charges** are never resized and never cut at the window. If one misses the prompt start: `RuntimeWarning` "prompt start not met with the explicitly set ... charge of X g,
  DeMar-equivalent N psi: P % of the grain cells (k of n ever) lit within W s of igniter firing (required 95 %); first light, t90; the grain did light and sustains combustion (a slow
  start) | never reached self-sustaining combustion (a hang-fire). Explicit charges are never resized."; `summary['prompt_start'] = False`, `summary['prompt_lit_fraction']`.
- **Plugin.** `ignited` False (never self-sustaining) stays an ERROR alert; grain lit but not promptly -> a WARNING alert ("slow ignition -- only P % of the grain lit within W s ... a larger
  igniter charge may be needed"), `success` stays True.
- **Summary keys.** `prompt_start`, `prompt_lit_fraction`, `prompt_window_s`, `prompt_window_effective_s`, `prompt_fraction_required`, `ignition_t50_s`, `ignition_t90_s`,
  `ignition_attempt_prompt_fractions`, `ignition_attempt_terminations`, `pyrogen_free_volume_in3`, `pyrogen_demar_equiv_psi` (= M W / V_free of the final charge; the first guess is 300 psi).
  `PyrogenChamber` gained `free_volume_in3`.
- **Pressure cutoff interplay.** An attempt that ends by the pressure cutoff (P < `P_cutoff` after the pyrogen is spent: 0.5 MPa by default in `run_from_ric` and the GUI) before the
  window is judged at its end, so with the default cutoff a slow start fails the attempt as soon as the chamber has blown down. The B0 rows pass `P_cutoff = 0.01 MPa`, which lets a slow start
  run on (Chunc and Hasegawa A below).

## 3. Results (N = 100, current defaults)

B0-row settings (`P_cutoff` 0.01 MPa; Chunc t_max 1 s, BALLSstick 2 s, Hasegawa A 10 s):

| motor | attempts (charges [g]) | lit fraction per attempt | final DeMar-equiv. | lit at 0.5 s | first light | t50 / t90 | startup |
|---|---|---|---|---|---|---|---|
| Chunc | 1 (0.815) | 0.99 | 300 psi | 0.99 | 3.5 ms | 234 / 235 ms | peak 10.51 MPa at 0.259 s, plateau 8.03, pk/plateau 1.31 (exp 1.027): a slow start, the chamber sits near ambient for ~0.23 s |
| BALLSstick | 3 (1.268, 1.585, 1.982) | 0.02, 0.03, 1.00 (the first two cut at 0.5 s) | 469 psi | 1.00 | 1.7 ms | 15.8 / 17.2 ms | peak 8.43 MPa at 29 ms, trough 4.43, mid 4.93: pk/trough 1.90, pk/mid 1.71 (exp ~1.50 / ~1.35): a real startup |
| Hasegawa A | 1 (7.861) | 0.99 | 300 psi | 0.99 | 3.0 ms | 108 / 110 ms | early peak (<0.5 s) 5.86, late peak 6.017 MPa at 2.5 s = 0.93x the 6.44 experiment |
| pathfin54 (GUI path settings) | 4 (0.439, 0.549, 0.686, 0.857) | 0.01, 0.03, 0.05, 1.00 (the first three end by pressure cutoff) | 586 psi | 1.00 | 2.1 ms | 12 / 16 ms | peak 9.24 MPa (early 8.02, late 8.39) |

With the default `run_from_ric` / GUI cutoff (0.5 MPa), the Chunc DeMar charge no longer survives its slow start: Chunc 3 attempts (0.815, 1.019, 1.273 g; fractions 0.03, 0.04, 1.00; 469 psi),
first light 1.7 ms, t50 / t90 7.9 / 8.6 ms, peak 11.82 MPa at 23 ms, plateau 8.22, pk/plateau 1.44 (exp 1.027). BALLSstick, Hasegawa A and pathfin54 are unchanged by the cutoff
(Hasegawa A passes with its first charge, t50 108 ms, and the run ends by pressure cutoff after burnout as before).

Cost: failed attempts are cut at the window or at the blowdown: 5-9 s each for Chunc / BALLSstick / pathfin54 on the 100-cell grid (a full run 12-56 s); the B0 gate takes 166 s (peak 2.5 GB).

## 4. Observations for the user

1. **The criterion passes very slow starts.** Chunc's 0.815 g DeMar charge (B0 settings) has 99 % of the grain lit at 0.235 s, yet the chamber sat at ~0.1 MPa until then: a 0.5 s window admits it;
   Hasegawa A is the same (t50 108 ms). Only the pressure cutoff (default 0.5 MPa) turns this into a retry. A window of ~50-100 ms would be a hard criterion for these two; a 0.5 s window is
   the user's number and was not changed.
2. The DeMar-equivalent pressure of the charges that were actually needed is 469 psi (BALLSstick, Chunc under the default cutoff) and 586 psi (pathfin54) against the 300 psi initial guess:
   raising the initial guess to ~500 psi would save the failed attempts for those three (but would over-charge Chunc / Hasegawa A under the B0 cutoff).
3. BALLSstick's B0 row is a real startup again (G1's hang-fire peak 0.29 / trough 0.076 MPa is gone).

## 5. Tests changed

- `tests/test_ig4_self_sustain.py::TestHangFire::test_explicit_dud_fails_the_test`: the warning text moved from "did not reach self-sustaining combustion" to "prompt start not met ... (a hang-fire)".
- Everything else of IG3 / IG4 (`test_ig3_stage2.py`, the rest of `test_ig4_self_sustain.py`) still holds and passes unchanged: the fallback, never-lit warning, explicit-charge rules and the plugin
  ERROR for a hang-fire are unaffected; the IG4 self-sustaining rate comparison stays as the `ignited` flag.
- New `tests/test_as1_prompt_start.py`: defaults; summary keys and the lit fraction vs the ignition times; DeMar-equivalent pressure; an explicit charge that misses a 10 ms window warns with the lit
  fraction and is neither resized nor cut; an auto charge is cut at the window (code 7), retried at x1.25 until it passes; the cap with a lit grain returns a full run + warning; nothing lit after the cap
  raises; the plugin returns a WARNING (not an ERROR) for a slow ignition.

## 6. Gates

B0 (re-recorded: only the BALLSstick row moved; G1 values kept in comments): PASS, `results/as1/B0_output.txt`. Full suite: see the report.
