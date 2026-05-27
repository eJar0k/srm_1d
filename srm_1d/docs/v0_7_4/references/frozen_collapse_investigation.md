# Phase B.0 cold IC + Sutton-default pyrogen sizing = ignition runaway

**Status (2026-05-27)**: **RESOLVED in v0.7.3.2** (commit pending at
write time). Below is the full discovery → diagnosis → fix narrative.
The investigation was triggered by the v0.7.4 Phase C.2 LHS smoke
test; the root cause turned out to be wider than initially documented
(canonical `hasegawa_motor_a.py` also affected, not just frozen).

## Resolution (v0.7.3.2, 2026-05-27)

Three coupled defaults changed:

1. **Kn-based pyrogen throat sizing** in
   `build_pyrogen_chamber` ([openmotor_adapter.py](../../openmotor_adapter.py)).
   Was `A_throat = 0.01·A_main` (=9 mm² for Hasegawa A, too small).
   Now `A_throat = A_burn / Kn_design` with Kn_design=100 (Sutton 9e
   §14.5 pellet-pyrogen mid-range), bounded to [1, 100] mm². For
   Hasegawa A this gives A_throat ≈ 46 mm² — physically defensible
   and close to the calibrated test config (38.5 mm²).
2. **`cfl_target` default 0.5 → 0.3** in `run_simulation`
   ([simulation.py](../../simulation.py)). Captures the pyrogen-
   plenum-into-cold-cell-0 transient on the µs timescale.
3. **`source_cfl_factor` default 0.10 → 0.05** in `run_simulation`.
   Tightens the source-driven dt cap so the cold-IC transient
   sub-steps properly.

After these defaults change:
- Canonical `hasegawa_motor_a.py` (Sutton-default sizing): clean run,
  P_peak 8.21 MPa, Total Impulse 23009.5 N·s, t_burn=10.0s.
- Frozen transport at same canonical knobs: clean run, P_peak 6.24 MPa.
- All 4 cells of the v0.7.4 Phase C.2 LHS sweep are now active
  (demar/radiation × frozen/effective).
- New pytest gate `test_canonical_hasegawa_motor_a_does_not_collapse`
  prevents silent regression of this fix.
- 277/277 pytest green (272 baseline + 5 new in
  `test_canonical_examples.py`).

The rest of this document is the original investigation log,
preserved for the regression-timeline analysis.

---

## Scope correction (2026-05-27 follow-up)

The original section below described frozen-only collapse. Followup
bisection confirmed:

| Config | Transport | Pyrogen sizing | Result |
|---|---|---|---|
| `hasegawa_motor_a.py` canonical | effective (default) | Sutton-default ~4.2g, 9 mm² throat | **COLLAPSE** (P_peak=1.44 MPa, t_burn=0.001s, term=4) |
| `hasegawa_motor_a.py` + frozen | frozen | Sutton-default | COLLAPSE (same) |
| Test helper `_short_hasegawa_a_run` | effective | EXPLICIT 12.3g, 38.5 mm² throat, 3.2 cm³ vol | RUNS CLEAN (pytest 272/272 green) |

**The discriminator is pyrogen throat area, NOT transport.** Sutton-
default `pyrogen_throat_area = min(max(0.01·A_main, 1e-6), 5e-5)`
gives ~9 mm² for Hasegawa A. The test calibration uses 38.5 mm² —
4× larger. Larger throat → lower plenum P_ig at steady state →
no positive-feedback runaway.

## Regression timeline

- v0.7.0 / v0.7.1 / v0.7.2: hot-bore IC (`T_gas = T_flame_propellant`
  at t=0). Sutton-default sizing worked fine — the bore was already
  at flame T, so pyrogen mass injection didn't drive a giant T
  gradient.
- v0.7.3-phaseB: IC switched to cold-bore (`T_gas = T_ambient`).
  Sutton-default sizing now creates ignition runaway because the
  cold bore + small throat + uncapped Saint-Robert P^n feedback
  blow up plenum P_ig within ~1 ms.
- v0.7.3.1: shipped without re-verifying the canonical example.
  The "23010 N·s impulse" claim in the v0.7.3.1 commit message was
  based on STALE artifact directories from pre-Phase-B runs, not a
  fresh execution.

## Artifact-handling lesson

The `artifact_dir()` helper creates timestamp-stamped subdirs per
run, which is correct. But verification logic that just checks
"PNGs got saved" is unsafe — a collapsed run still saves the
collapsed plots into a fresh timestamp dir. Confirmed by reading
`P_peak=1.44 MPa` directly from the simulation summary line in the
current run, NOT from looking at any saved file.

**Recommended safeguards** (to be added separately):
1. Example scripts should print a CLEAR success/failure banner
   based on `summary['termination_code']` after the run completes.
2. A health-check assertion (e.g., `t_burn > 0.1 s and term_code in
   (1, 2)`) at the top of `if __name__ == '__main__'` would prevent
   silent collapsed-run completion.
3. Verification scripts (LHS, sweep drivers) should NEVER look at
   directory contents as a proxy for "ran clean" — only at fresh
   stdout of the just-completed process.

## Why this is new (original section, retained)

## What was observed

The v0.7.4 Phase C.2 LHS smoke test (5 samples per cell, 4 cells in
the originally-planned 2x2 mode × transport grid) showed that **all
frozen-cell runs tripped numerical collapse (termination_code = 4)
within ~2 ms of t=0**, regardless of LHS knob draw. The collapse
manifests as runaway pyrogen-plenum pressure: `pyrogen_peak_P_MPa`
hits ~190 MPa within a handful of timesteps, far beyond physical
reality (real Hasegawa A peaks at ~6 MPa).

Sample LHS row that collapsed (`demar_frozen`, sample 0):
- roughness = 18.8 µm (near the v0.7.0 calibrated 37.1 µm)
- kappa = 0.451 (essentially v0.7.0 canonical 0.45)
- T_ignition = 875.7 K (near v0.7.0 canonical 850 K)
- k_solid = 0.276 W/(m·K) (just below the literature center 0.30)
- Result: P_peak_sim = NaN, t_burn_sim = 0.0014 s, term_code = 4

This is suspicious because the knobs are *near* the v0.7.0
calibration that produces a clean ignition trace. To rule out the
LHS perturbation as the issue, I ran a control with the canonical
knobs (roughness=37.1 µm, kappa=0.45, T_ign=850, k_solid=default
0.30) but explicitly loaded the FROZEN transport YAML:

```
Canonical knobs + frozen transport: P_peak=1.91 MPa, t_burn=0.002s, code=4 (COLLAPSE)
Canonical knobs + effective (default): runs clean (impulse 23010 N·s)
```

So the collapse is driven by the frozen YAML's lower k_gas under
Phase B.0's cold-bore IC, not by the LHS perturbation.

## Why this is new

The v0.7.1 Phase 5 Task 1 frozen/effective A/B (memory
`project_v0_7_1_phase5_task1_task2_findings`) found that switching
transport changed *spike height* by ~32%, not *stability*. That
analysis was at the **v0.7.0 hot-bore IC** (`T_initial_gas =
T_flame_propellant`). Phase B.0 swapped to a cold-bore IC
(`T_ambient`), which legitimately amplifies the ignition spike. The
combination of:

1. Cold-bore IC → strong T-gradient between bore gas (300 K) and
   pyrogen plenum products (~2800 K)
2. Frozen RPA k_gas (low compared to effective) → bore gas conducts
   less heat away from the impingement zone
3. Phase B.3 pyrogen geometry (Mizushima 3.2 mm pellets, A_burn
   ~46 cm²) → fast pyrogen mass injection at the head end
4. Phase B.4 mode = DeMar (cell-0 lumped flux 69.4 cal/cm²/s) →
   adds to the convective heat at cell 0

…drives the pyrogen plenum into a positive-feedback loop: rising
P_ig → faster Saint-Robert burn rate (a · P^n with n=0.5) → higher
mdot → higher P_ig → higher r_b → …

Effective transport has higher k_gas (literature: k_eff ≈ 0.65 vs
k_frozen ≈ 0.37 for Hasegawa A; memory
`reference_hasegawa_a_effective_rpa_values`). The extra gas
conductivity siphons heat into the cold bore fast enough to break
the feedback before it runs away.

## Plausible fixes (none implemented yet)

1. **Stability knob — source_cfl_factor**: the v0.7.0 design has a
   source-CFL safety factor (default 0.10) that caps dt against
   source magnitudes. Tightening to 0.05 may absorb the high-source
   transient regime. Quick to test.

2. **Adaptive Cp_pyrogen sensible cap**: the cell-0 DeMar sensible
   power `mdot_igniter · Cp_pyrogen · (T_ig - T_surf)` becomes very
   large under cold-bore IC (T_surf at 300 K initially). Cap could
   be made adaptive to current T_surf rather than the IC value.

3. **Soft floor on dt**: numerical collapse currently trips on
   `dt < COLLAPSE_DT_THRESHOLD = 1e-9 s` after
   COLLAPSE_CONSECUTIVE_STEPS = 3 (`simulation.py:1056-1059`). The
   threshold may be too tight for the new dynamics; consider 1e-10.

4. **Pyrogen-plenum P limiter**: cap P_ig at some physical maximum
   (e.g., 30 MPa, the typical Sutton "max sustainable pyrogen
   pressure") to prevent the feedback loop from running away even
   in unstable knob regions.

5. **IC compromise**: preserve cold-bore for uncontained topologies
   (where it's essential — see Phase B.0 motivation) but allow a
   warm-bore option for forward_plenum (e.g.,
   `T_initial_gas = 0.5·T_flame`). Mixed-IC adds complexity but may
   be the cleanest physical compromise.

## Impact on user's note (3)

The user (2026-05-27) asked to A/B-LHS frozen vs effective because
"the earlier v0.6 and v0.7.0ish Hasegawa runs looked great with
frozen values, but they blew up on Chunc, Zerox, etc." The
investigation result extends that finding: at Phase B's IC, **frozen
also blows up Hasegawa A**, the previously-stable canonical case.

This means the user's instinct that frozen vs effective is "genuinely
uncertain" is more correct than the v0.7.1 conclusion suggested. The
Phase 5 task1 finding ("switching k_gas changes spike height +32%
not stability") is **no longer accurate post-Phase-B**.

## Decision driver for v0.7.3.2

When this investigation runs, the architectural choice is:
- **Accept frozen as incompatible** with cold-bore IC; document and
  default to effective everywhere. Risk: motors that need frozen
  thermo (e.g., propellants with strong condensed-phase effects
  that effective captures wrong) lose validation freedom.
- **Add a stability mechanism** so frozen can coexist with cold-bore
  IC. One of the five fixes above. Preferred for architectural
  flexibility, but each fix has knob-calibration debt.

Recommended: try fix #1 (source_cfl_factor tighten) first since it's
a single-line change and the v0.7.0 source-CFL was already designed
exactly for this scenario.

## Related

- Memory `project_v0_7_1_phase5_task1_task2_findings` (old frozen/effective A/B at hot IC)
- Memory `reference_hasegawa_a_effective_rpa_values` (k_eff 0.65 vs k_frozen 0.37)
- Memory `project_v0_7_3_phaseB_state` (cold-IC introduction)
- `srm_1d/docs/v0_7_3/PHASE_B_SCOPE.md` §B.0 (IC change motivation)
- Numerical collapse detector: `simulation.py:1056-1059`
