# Hasegawa A v0.7.0 calibration LHS — 2026-05-21

## Headline

**Step 6 PASS. Step 5 (Peretz participation fraction) is NOT needed.**

The v0.7.0 7-variable LHS rank-1 produces a Hasegawa A pressure trace
with MSE = **0.0968 MPa²** vs experimental, beating both:

- v0.6.0 baseline (MSE = 0.24 MPa², with the now-removed
  `igniter_tau = 127 ms` FSI-cushioning proxy) — **2.5× better**
- The plan's `DESIGN.md` "Validation strategy" target (MSE < 0.15 MPa²)

## Setup

- Branch: `v0.7.0-phase4` at commit `35320c5` (post 7-var LHS
  expansion).
- LHS: `srm_1d/examples/hasegawa_a_lhs.py` with
  `SRM_HASEGAWA_LHS_SAMPLES=500`, `scipy.stats.qmc.LatinHypercube(seed=42)`,
  parallelized across CPU cores.
- Fitness: segmented (spike 0.25 / post 0.35 / plateau 0.20 / tail
  0.20 MSE weights). `mse_all` reported separately is the un-segmented
  full-trace MSE.
- Motor: Hasegawa A, `t_max = 6.0 s`, `P_cutoff = 0.05 MPa`.
- Artifacts: `artifacts/hasegawa_a_lhs/`.

## Rank-1 parameters

| parameter             | v0.6.0 (igniter_tau) | v0.7.0 (pyrogen)        |
| --------------------- | -------------------- | ----------------------- |
| roughness             | 37.1 µm              | **37.5 µm**             |
| kappa                 | 0.45 (locked)        | **0.429**               |
| T_ignition            | (not varied)         | **927 K**               |
| k_solid               | (not varied)         | **0.482 W/(m·K)**       |
| pyrogen_mass          | n/a (exponential)    | **12.3 g**              |
| pyrogen_throat_area   | n/a (exponential)    | **38.5 mm²**            |
| pyrogen_volume        | n/a (exponential)    | **3.2 cm³**             |
| igniter_tau           | 127 ms (FSI proxy)   | **REMOVED**             |
| mse_all (MPa²)        | 0.240                | **0.0968**              |

## Fit quality

- P_peak (sim) = 6.527 MPa vs experimental 6.436 MPa
- **peak_error_pct = +1.41 %**
- trough_error_pct = +3.45 %

## Pass-fraction over the LHS

- N total = 500
- Beat v0.6.0 (mse < 0.24): 4 / 500 (0.8 %)
- Meet plan target (mse < 0.15): 3 / 500 (0.6 %)

The pass band is narrow but does include the v0.6.0-style roughness ≈
37 µm region. Three of the four passing samples have roughness in
24–38 µm, consistent with the physical range for cast AP/HTPB/Al
composites flagged in the v0.6.0 calibration memo. The fourth
(roughness 10.7 µm, mse 0.205) is near-pass and outside the canonical
range — possibly an LHS-sampling artifact at the lower roughness edge.

## What the LHS validated

1. **Pyrogen plenum replaces igniter_tau cleanly.** The v0.6.0
   FSI-cushioning proxy is gone; the v0.7.0 0D pyrogen plenum +
   Goodman ignition produces better fits at a physically grounded
   pyrogen mass (~12 g) within Sutton Eq. 15-4 sizing range.
2. **roughness ≈ 37 µm is reproducible across the API break.** Same
   value, different ignition model, comparable fit quality —
   suggests roughness is a real physical knob for AP/HTPB/Al
   composites, not an FSI-cushion bypass.
3. **kappa = 0.429** (LHS, free) is slightly under v0.6.0's locked
   0.45. Within Ma 2020's reported range.
4. **The 12.0 → 6.5 MPa peak reduction came from kappa + roughness
   joint calibration**, not from any new burn-establishment physics.
   No participation/ramp model was required.

## Why Step 5 isn't needed

The Step 4 audit hypothesized that the residual 12 MPa peak required
a Peretz-aligned participation fraction to bring it down. The LHS
disproves that: with kappa free and roughness in the v0.6.0 range,
the existing Ma erosive-burning model can match the spike directly.
The Step 4 sweep failed only because **kappa was locked at 0.45** for
those variants; opening it up to the LHS solved the problem.

This is a strong outcome for the v0.7.0 calibration philosophy: no
unfounded smoothing / no fitted ramps, and the model still beats v0.6.0.

## Decision

**Step 5 SKIPPED.** Move directly to Step 7 (docs + tag).

The Peretz-aligned participation fraction stays in the v0.7.1 roadmap
as a potential refinement if future motors with different geometries
or aluminum loadings show residuals that the joint roughness/kappa
calibration can't close.
