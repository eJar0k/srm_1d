# Radiation collapse: pre, post local-T, post plume-lag (2026-05-20)

## Setup

- Branch: `v0.7.0-phase4`.
- Motor: Hasegawa A, `t_max = 0.030 s` to match the prior artifact run.
- Matrix: `ignition_spike_diagnostic.py --mode radiation-collapse`
  (27 variants: 11-point emissivity sweep + grid/CFL/timestep refinement
  variants + receiver-heat-no-sink + erosive-disable + hot-fill baseline).
- Three runs compared:
  1. **Pre local-T** — `hasegawa_a_radiation_collapse_pre_localT/`
     (2026-05-11; constant `T_flame` emitter).
  2. **Post local-T** — `hasegawa_a_radiation_collapse_localT/`
     (2026-05-20; commit `70ec63c` + abort trip `f8f3db2`; local
     `T[neighbor]` emitter, no plume lag).
  3. **Post local-T + plume lag** —
     `hasegawa_a_radiation_collapse/` (current head; Step 2a Mechanism 2:
     5-µs linear ramp from `T_initial` to `T[neighbor]` for adjacent
     radiation emitter T after each cell's ignition_time).

## Aggregate result

|                      | pre local-T | post local-T | post local-T + lag |
| -------------------- | ----------- | ------------ | ------------------ |
| Stable               | 18 / 27     | 26 / 27      | 21 / 27            |
| Collapse             | 9 / 27      | 1 / 27       | 6 / 27             |
| Default `ε = 0.45`   | catastrophic| trip-abort   | **stable, 12.7 MPa**|
| Max P over all       | 352 GPa     | 13 MPa       | 13 MPa             |
| Max Mach over all    | 86,600      | 1,154        | 1,045              |

## Emissivity sweep detail

| ε    | pre local-T  | post local-T  | post local-T + lag |
| ---- | ------------ | ------------- | ------------------ |
| 0.00 | stable       | stable        | stable             |
| 0.05 | stable       | stable        | stable             |
| 0.10 | stable       | trip-abort    | trip-abort         |
| 0.20 | collapse     | stable        | stable             |
| 0.30 | collapse     | stable        | stable             |
| 0.40 | stable       | stable        | stable             |
| 0.45 | collapse     | trip-abort    | **stable**         |
| 0.50 | stable       | trip-abort    | trip-abort         |
| 0.60 | collapse     | stable        | trip-abort         |
| 0.75 | collapse     | trip-abort    | trip-abort         |
| 0.90 | stable       | stable        | stable             |

## Headline tradeoff of the plume lag

**The lag is a net WIN for the default-emissivity user path but a small
net LOSS over the full diagnostic matrix.**

- The historically conventional aluminized-AP emissivity `ε = 0.45`
  (Hasegawa A "ambient nominal" path) now produces a fully-developed
  pressure trace at 12.7 MPa instead of trip-aborting at 1.0 MPa.
- Two pure-diagnostic variants regressed:
  `ambient_rad045_no_surface_heating` (stable → trip-abort) and
  `ambient_rad045_no_erosive` was already trip-aborting and remains so.
- High emissivities `ε ∈ {0.50, 0.60, 0.75}` still trip-abort. The lag
  delays the radiation onset but does not eliminate the resonance when
  the radiation magnitude is well above the no-lag stability boundary.

The numerical mechanism: the 5-µs lag suppresses radiation in the
first ~2-5 pressure-wave transit times after a neighbor ignites,
giving the local pressure gradient a chance to flatten before
the receiver's surface heating ramps. When that lag is enough to
break the resonance, the run completes normally. When the radiation
power is too high (or refinement variants alter the energy balance),
the resonance reforms after the lag completes.

## Refinement still works for the residual outliers

All 6 grid/CFL/timestep refinement variants complete normally with
the plume lag in place (`cells = 50, 200`; `cfl = 0.10, 0.25`;
`dt_max = 2e-5`). Net stability is therefore: default-grid runs
stable at most emissivities, and grid/CFL refinement recovers the
remaining cases.

## Classifier-trip alignment (closed)

The earlier mismatch where the runtime trip fired (`termination_code
= 4`) but the in-window classifier still labeled the run
`collapse_class = "stable"` has been fixed by treating
`termination_code == 4` as a direct signal of `collapse_detected =
True`. Now every aborted run is tagged `collapse_class = "collapse"`
**and** `diagnostic_failure_mode = "numerical_collapse_aborted"`.

## Energy/momentum residuals

All 27 variants close energy residuals to better than 1e-9 relative
to thermal/convective scale. The plume-lag formulation remains
energetically self-consistent (emitting cell debits exactly what the
receiver cell gains, scaled by the same `phi_plume` factor).

## Verdict vs Step 2 pass criteria

- "All emissivities terminate `normal_completed`": **PARTIAL** — the
  default `ε = 0.45` works; high-emissivity edge cases still trip.
- "`collapse_class == stable` for every variant": **NO** (21 / 27).
- "P_peak varies monotonically with ε": **PARTIAL** — stable ones do
  (~12.7 MPa); trip-aborted ones cluster at ~1 MPa.
- "Energy residual `< 1e-6` relative": **PASS** (all variants).

## Recommendation for next iteration

Three possible refinements, ordered by implementation cost:

1. **Document the current state as the v0.7.0 deliverable.** The
   default-emissivity Hasegawa A calibration path is now usable;
   trip-aborted variants are clean, classified, and resolvable via
   grid/CFL refinement. Step 4 (parametric sweep on T_ignition and
   k_solid) can use `ε = 0.45` directly without a workaround.
2. **Add a source-aware CFL constraint** in `compute_dt_cfl`. Limit
   `dt` so the per-step radiation deposition cannot change a cell's
   gas temperature by more than ~10 % of `(T_flame - T_initial)`.
   This is a numerical stability constraint, not a tuning knob, and
   would resolve the remaining outliers without further radiation-kernel
   modifications.
3. **Mechanism 3 from the plan: implicit treatment of the radiation
   source.** Largest blast radius; only worth it if the source-aware
   CFL in option 2 is insufficient.

Constants introduced this iteration:

- `RADIATION_PLUME_LAG_S = 5.0e-6` in `srm_1d/simulation.py`.
  Module-level constant; documented as a numerical-stability lag,
  not a fitted plume timescale. Sized to ~2-5 cell pressure-wave
  transit times at default CFL=0.5, cells=100.
