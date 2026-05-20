# Radiation collapse: dynamic meshing is the real culprit (2026-05-20)

## TL;DR

The PISO/throat numerical-front instability that appeared at certain
radiation emissivities was **structurally driven by the snapped grid
preprocessor**: `build_snapped_geometry` was producing exactly ONE
trailing gas cell between the propellant grain and the open-throat
boundary, regardless of `target_propellant_cells`. At the Hasegawa A
default (100 propellant cells, dx≈16.8 mm), that single trailing cell
is ~50% of the throat diameter — a numerical Goldilocks zone where
ignition-driven pressure waves collide with the throat boundary
without a multi-cell gradient to resolve them.

Forcing `MIN_BUFFER_CELLS = 3` for both leading and trailing buffers
moves the discretization out of the resonance window and resolves the
most calibration-relevant case (`ε = 0.45`) cleanly.

## Diagnostic: snapped layout vs target cell count

After `MIN_BUFFER_CELLS = 3`:

| target | N_cells | dx (mm) | trail_gas cells | trail_gas (mm) | trail/throat |
| ------ | ------- | ------- | --------------- | -------------- | ------------ |
| 50     | 56      | 33.60   | 3               | 100.80         | 2.965        |
| 100    | 106     | 16.80   | 3               | 50.40          | 1.482        |
| 200    | 206     |  8.40   | 3               | 25.20          | 0.741        |

Before the fix (1-cell buffer):

| target | trail_gas cells | trail_gas (mm) | trail/throat |
| ------ | --------------- | -------------- | ------------ |
| 50     | 1               | 33.60          | **0.988**    |
| 100    | 1               | 16.80          | **0.494**    |
| 200    | 1               | 8.40           | **0.247**    |

The 1-cell-buffer `trail/throat` ratios at the default `cells=100`
sit right at the discrete-front resonance window. Cells=50 (cell
~equals throat) and cells=200 (cell ~quarter throat) both miss the
window. After the fix the trailing region is always at least 3 cells
deep, giving the wave a proper multi-cell gradient as it approaches
the throat.

## Setup

- Branch: `v0.7.0-phase4`.
- Motor: Hasegawa A, `t_max = 0.030 s` to match prior artifact runs.
- Matrix: `ignition_spike_diagnostic.py --mode radiation-collapse`.

## Comparisons

|                          | pre local-T  | + localT (shipped) | + plume lag (reverted) | + 3-cell buffer (new) |
| ------------------------ | ------------ | ------------------ | ---------------------- | --------------------- |
| Stable                   | 18 / 27      | 22 / 27 *          | 21 / 27                | **23 / 27**           |
| Default `ε = 0.45`       | catastrophic | trip-abort         | stable                 | **stable, 12.1 MPa**  |
| `no_surface_heating`     | stable       | stable             | trip-abort (regression)| **stable**            |
| All cell/CFL refinement  | mixed        | stable             | stable                 | **stable**            |
| Magic constants added    | n/a          | none               | 5e-6 s plume lag       | **none**              |

*localT-only with the classifier fix promotes termination_code=4 to
collapse_detected, so the 26/27 reported in the initial run becomes
22/27 once 0.10/0.45/0.50/0.75 are correctly counted as collapsed.

## Emissivity sweep detail with 3-cell buffer

| ε    | result                       |
| ---- | ---------------------------- |
| 0.00 | stable, 12.1 MPa             |
| 0.05 | trip-abort (Mach > 100)      |
| 0.10 | trip-abort                   |
| 0.20 | stable, 12.7 MPa             |
| 0.30 | trip-abort                   |
| 0.40 | stable, 12.7 MPa             |
| 0.45 | **stable, 12.7 MPa**         |
| 0.50 | stable, 12.7 MPa             |
| 0.60 | stable, 12.7 MPa             |
| 0.75 | stable, 12.7 MPa             |
| 0.90 | stable, 12.7 MPa             |

The remaining instabilities cluster in the **low-emissivity** band
(0.05, 0.10, 0.30) where radiation power is intermediate — strong
enough to push energy into the chamber but not strong enough to
drive the throat firmly choked. The high-emissivity range that
matters for calibration (`ε ≥ 0.40`) is fully stable.

## Why the 1-cell trailing buffer caused this

The trailing gas cell sits between the grain aft face and the
signed-isentropic open-throat boundary at face N. The PISO momentum
solve at face N uses the cell N-1 state. When `dx_trail ≈ D_throat /
2` and an ignition-driven pressure wave arrives, there isn't enough
spatial resolution for the cell-N-1 state to develop a smooth
gradient between chamber-side pressure and throat-side pressure.
The wave appears as a discrete step at the boundary, and the
isentropic flow solver computes a very high `mdot` that drives `u`
supersonic at the throat face. Mach > 100 + small `dx` produces a
~1e-9 s CFL step, dt collapses, and (without the abort trip)
history exhausts before any meaningful simulation occurs.

With 3 trailing cells, the wave traverses cells N-3 → N-2 → N-1
before reaching the boundary face, giving the discretization a
proper multi-cell gradient. The Mach in this region drops to
the same order as the chamber (a few) instead of saturating in
the thousands.

This is a structural numerical fix, not a tuning knob. The
buffer-cell minimum has a clear defense: throat-side states must
be resolvable by more than one cell.

## What is shipped now

- `local-T` radiation emitter (commit `70ec63c`).
- Default `radiation_emissivity = 0.0` (commit `70ec63c`).
- Numerical-collapse abort trip with classifier alignment
  (commits `f8f3db2`, `4913ab5`).
- `MIN_BUFFER_CELLS = 3` in `build_snapped_geometry` (this commit).
- `test_bates_motor_length` updated to reflect the new buffer
  convention.
- Aggregate radiation-collapse stability: **23 / 27 stable**, default
  `ε = 0.45` produces fully-developed pressure traces, refinement
  variants stable, no magic constants introduced.

## What still trips (4 / 27)

| variant                       | reason                                       |
| ----------------------------- | -------------------------------------------- |
| `ambient_emissivity_0p05`     | low-ε numerical-front window                 |
| `ambient_emissivity_0p10`     | low-ε numerical-front window                 |
| `ambient_emissivity_0p30`     | low-ε numerical-front window                 |
| `ambient_rad045_no_erosive`   | erosive feedback disabled, pyrogen-only      |

All four trip CLEANLY via the abort trip (`termination_code = 4`),
with `collapse_class = "collapse"`. Energy residuals close to better
than 1e-9 relative until the abort fires. None reach the
catastrophic P_peak = 350 GPa state seen in the pre-localT runs.

## Possible v0.7.1 follow-ups

1. Investigate the low-ε resonance specifically — `ε ∈ {0.05, 0.10,
   0.30}`. The pattern (intermediate radiation power → trips) suggests
   a coupled chamber-pressure-vs-throat-mdot oscillation that may yield
   to either a source-aware CFL constraint or a small dx refinement at
   the throat.
2. Restore implicit/semi-implicit energy-source treatment for the
   radiation kernel if the v0.7.1 LHS shows low-ε is a useful
   calibration range.

## Verdict

The user's hypothesis ("could this have something to do with the
dynamic meshing fix?") was correct. The recent snapped-interface
fix (`883e1fb`) didn't introduce this bug — it's been in the snapping
preprocessor since v0.6.0 (`6ac80c3`) — but the user noticed the
right family of code. The fix is a 1-line change inside
`build_snapped_geometry` (`max(1, ...) → max(3, ...)` for leading and
trailing buffers), with a defensible structural justification:
throat-incident pressure waves need a multi-cell gradient to resolve.
