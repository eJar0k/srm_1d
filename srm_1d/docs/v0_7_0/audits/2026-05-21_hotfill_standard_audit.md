# Hot-fill standard-mode audit (2026-05-21)

## Setup

- Branch: `v0.7.0-phase4` at commit `f82bda3` (final radiation-stability state).
- Motor: Hasegawa A, `t_max = 0.050 s`.
- CLI: `ignition_spike_diagnostic.py --mode standard --t-max 0.05`.
- Artifacts: `artifacts/ignition_diagnostics/hasegawa_a/`.

Purpose: per Step 3 of `continue-with-the-numerical-zippy-dawn.md`,
audit the still-instant ~23 ms hot-fill erosive snap on the current
shipped solver to decide whether a participation/ramp model is
required (Step 5) or whether T_ignition/k_solid tuning (Step 4) alone
can match the experimental Hasegawa A trace.

## Variant summary

| variant                    | P_max  | t_peak | t_firstIg | t_firstEros | t_fullGrainIg | spread 10–90 | max Mach | clip / thermal |
| -------------------------- | ------ | ------ | --------- | ----------- | ------------- | ------------ | -------- | -------------- |
| baseline (hot-fill)        | 12.12  | 24.5   | 2.81      | 3.65        | 4.07          | 0.45         | 1.8      | 0.8 %          |
| ambient_initial_gas        | 12.07  | 24.6   | 2.80      | 3.66        | 5.29          | 1.57         | 48.5     | 23.4 %         |
| ambient_no_surface_heating | 12.07  | 24.6   | 3.38      | 3.73        | 5.33          | 1.53         | 49.5     | 23.5 %         |
| ambient_no_radiation       | 12.07  | 24.6   | 2.80      | 3.66        | 5.29          | 1.57         | 48.5     | 23.4 %         |
| no_erosive                 | 4.87   | 19.8   | 2.81      | --          | 4.12          | 0.49         | 1.4      | 0.6 %          |
| no_endfaces                | 12.12  | 24.5   | 2.81      | 3.65        | 4.07          | 0.45         | 1.8      | 0.8 %          |
| no_momentum                | 12.25  | 24.3   | 2.81      | 3.66        | 4.07          | 0.45         | 1.8      | 0.8 %          |

(Pressures MPa, times ms.)

## Decision-tree branches

**`n_burning` jumps in 1–2 steps after first ignition** → **YES** for
hot-fill (0.45 ms spread, 1.26 ms from first cell to full grain). This
matches the Peretz/Pardue/Cavallini instantaneous-ignition convention
documented in DESIGN.md §"Architectural decisions" point 6: "Sequential
per-cell ignition: flame spread is *emergent* from the per-cell
criterion, not separately modeled."

**Gas-energy clipping > 10 % of thermal source energy** → **NO** for
hot-fill (0.8 %, well under threshold). The chamber absorbs the
propellant gas injection without significant ceiling clipping. Clipping
is significant for ambient_initial_gas variants (~23 %) but those
remain stable due to the 3-cell trailing buffer + source-CFL combo.

**`T_surf − T_ignition` ≈ 0 K at first burn** vs **~50–100 K** → not
directly readable from variant_summary.csv; would need per-cell
inspection from `ignition_times.csv`. Deferred since the first two
branches already resolved the audit.

## Verdict

The hot-fill case is in **branch 1** of the decision tree: ignition
spread itself is fine — the model produces the expected
instantaneous-ignition convention. The spike-shape question is **not**
a participation/ramp problem; it's an **erosive feedback magnitude**
problem.

Direct evidence: turning erosive off (`no_erosive` variant) drops the
peak from 12.12 MPa → 4.87 MPa. The erosive contribution is therefore
~7.2 MPa (60 % of the peak). Hasegawa A's experimental peak is
~6.4 MPa, so we need to either:

1. **Reduce normal burn rate** (smaller `a` or shifted Saint-Robert tab)
2. **Reduce erosive contribution** via `kappa` or `roughness` retuning
3. **Add a Peretz-style burn-establishment ramp** (Step 5)

The plan's pass criterion for Step 4 (T_ignition / k_solid sweep to
MSE < 0.15 MPa²) is therefore the natural next step. Neither
T_ignition nor k_solid directly modifies burn-rate magnitude — they
shift the ignition timing — so Step 4 is mostly about confirming
whether existing roughness / kappa + ignition timing can fit. If it
can't, **Step 5 burn-establishment ramp becomes the v0.7.0-blocker**.

## Other observations

- **First ignition cell = 3** (not 0). The leading 3-cell gas buffer
  pushes the head-end grain start to cell 3, which is correct given
  the geometry preprocessor change committed in `cf488ec`.
- **`baseline`, `no_endfaces`, `no_momentum` are identical to within
  1 %** for hot-fill — end-face injection and pyrogen momentum are not
  driving the spike. This confirms momentum is implemented but not the
  dominant Hasegawa A driver, consistent with the 2026-05-11 finding.
- **`ambient_initial_gas == ambient_no_radiation`** (byte-equivalent
  rows) — confirms that with `radiation_emissivity = 0.0` default,
  radiation contributes nothing, and the `no_radiation` toggle is
  defensible as a diagnostic-only knob.
- **`ambient_no_surface_heating` shifts first ignition from 2.80 ms
  → 3.38 ms** (a 0.58 ms delay) and the first ignition cell moves
  from 3 → 4. The pyrogen direct surface heating is the load-bearing
  mechanism that ignites cell 3 first; without it, ignition propagates
  in from cell 4 onward via convection.
- **Energy residuals close to < 1e-9 relative** across all variants.

## Next step

Proceed to Step 4: implement `--mode ignition-tuning` in
`examples/ignition_spike_diagnostic.py` with the 18-variant Cartesian
sweep over `T_ignition ∈ {650, 750, 850} K`, `k_solid ∈ {0.2, 0.3, 0.5}
W/(m·K)`, `radiation_emissivity ∈ {0.0, 0.45}`. Pass criterion:
≥1 row with MSE vs Hasegawa A < 0.15 MPa². Escalate to Step 5 if
none qualify.
