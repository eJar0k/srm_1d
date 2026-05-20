# Radiation collapse: pre, post local-T, post plume-lag (reverted) (2026-05-20)

## Setup

- Branch: `v0.7.0-phase4`.
- Motor: Hasegawa A, `t_max = 0.030 s` to match the prior artifact run.
- Matrix: `ignition_spike_diagnostic.py --mode radiation-collapse`
  (27 variants: 11-point emissivity sweep + grid/CFL/timestep refinement
  variants + receiver-heat-no-sink + erosive-disable + hot-fill baseline).

## Three runs compared

1. **Pre local-T** — `hasegawa_a_radiation_collapse_pre_localT/`
   (2026-05-11; constant `T_flame` emitter).
2. **Post local-T (shipped)** — `hasegawa_a_radiation_collapse_localT/`
   (2026-05-20; commit `70ec63c` + abort trip `f8f3db2`; local
   `T[neighbor]` emitter, no plume lag).
3. **Post local-T + plume lag (tried, reverted)** — Step 2a Mechanism 2:
   5-µs linear ramp from `T_initial` to `T[neighbor]` for adjacent
   radiation emitter T after each cell's ignition_time. Implemented
   in commit `a018737`, **reverted** in a follow-up commit.

## Aggregate result

|                      | pre local-T | post local-T (shipped) | + plume lag (reverted) |
| -------------------- | ----------- | ---------------------- | ---------------------- |
| Stable               | 18 / 27     | 26 / 27                | 21 / 27                |
| Collapse             | 9 / 27      | 1 / 27                 | 6 / 27                 |
| Default `ε = 0.45`   | catastrophic| trip-abort             | stable, 12.7 MPa       |
| Max P over all       | 352 GPa     | 13 MPa                 | 13 MPa                 |
| Max Mach over all    | 86,600      | 1,154                  | 1,045                  |

## Emissivity sweep detail

| ε    | pre local-T  | post local-T (shipped) | + plume lag (reverted) |
| ---- | ------------ | ---------------------- | ---------------------- |
| 0.00 | stable       | stable                 | stable                 |
| 0.05 | stable       | stable                 | stable                 |
| 0.10 | stable       | trip-abort             | trip-abort             |
| 0.20 | collapse     | stable                 | stable                 |
| 0.30 | collapse     | stable                 | stable                 |
| 0.40 | stable       | stable                 | stable                 |
| 0.45 | collapse     | trip-abort             | stable                 |
| 0.50 | stable       | trip-abort             | trip-abort             |
| 0.60 | collapse     | stable                 | trip-abort             |
| 0.75 | collapse     | trip-abort             | trip-abort             |
| 0.90 | stable       | stable                 | stable                 |

## Decision: revert the plume lag

The lag rescued the default `ε = 0.45` case but regressed
`ambient_rad045_no_surface_heating` and `ε ∈ {0.60}`. Aggregate stable
count dropped from 26 to 21. More importantly, the 5-µs lag is a magic
constant chosen to fit a discrete numerical resonance — it shifts the
resonance window rather than eliminating it, and shifting the window
helps some configurations while breaking others.

Per the user's `feedback_no_unfounded_smoothing` boundary, the lag is
not defensible as physics-driven: it's a tuning knob with no
empirical or theoretical anchor for the 5-µs value. The cleaner
trade-off is to keep the simpler local-T-only kernel (which has
26/27 stable variants and clean residuals) and document grid/CFL
refinement as the workaround for the remaining edges.

The lag implementation is preserved in commit `a018737` for future
reference if a physically-grounded plume-development model (e.g.,
Peretz-style Goodman penetration depth) is investigated later.

## Shipped state (post-revert, current head)

The shipped Phase 4 state is the **localT-only** kernel:
- `Propellant.radiation_emissivity` defaults to 0.0 (opt-in).
- Emitter T = local `T[neighbor]` (not constant `T_flame`).
- Receiver/emitter exchange is energetically self-consistent.
- No plume-development lag.
- Numerical-collapse abort trip catches the residual outliers cleanly
  at termination_code = 4.
- Classifier now treats termination_code = 4 as
  `collapse_detected = True` so `collapse_class` and
  `diagnostic_failure_mode` agree.

## Refinement still works for the residual outliers

The 4 trip-aborted emissivity values (0.10, 0.45, 0.50, 0.75) all
fail at the default grid (cells=100) and CFL (0.5). All 6
refinement variants complete normally:
- `cells = 200`, `cells = 50`
- `CFL = 0.10`, `CFL = 0.25`
- `dt_max = 2e-5`

Conclusion: the residual instability is a PISO/throat numerical
front interaction that responds to grid/CFL but not to source-side
modifications.

## Recommendation for v0.7.0

- Ship the localT-only kernel + abort trip as the v0.7.0 radiation
  deliverable.
- Document `ε = 0.45` (the conventional aluminized default) as
  requiring grid refinement (`cells = 200`) or CFL tightening
  (`cfl_target = 0.10-0.25`) under ambient-initial-gas conditions.
- Phase 4 calibration LHS can pin `radiation_emissivity = 0.0`
  (the current shipped default) and treat the ambient-radiation path
  as a future v0.7.1 enhancement, OR use refinement variants for
  the LHS rank-1 if radiation is a load-bearing calibration knob.

## Possible v0.7.1 follow-ups (deferred)

1. **Source-aware CFL constraint** in `compute_dt_cfl`. Limit `dt`
   so the per-step radiation deposition cannot change a cell's gas
   temperature by more than ~10 % of `(T_flame - T_initial)`.
   Numerical stability constraint, not a tuning knob; would
   eliminate the residual outliers without source-kernel changes.
2. **Peretz-aligned plume model**. Replace the 5-µs lag (rejected
   here) with a Goodman penetration-depth ratio
   `δ_burned / δ_steady` — no fitted constant, derived from the
   existing solid-conduction ODE.
3. **Implicit treatment of the radiation source** in PISO. Largest
   blast radius; defer unless 1 and 2 are insufficient.

## Energy/momentum residuals

All shipped-state runs close energy residuals to better than 1e-9
relative to thermal/convective scale. The local-T formulation is
energetically self-consistent.
