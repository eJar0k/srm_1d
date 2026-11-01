# G1 — physical head and aft ends: replace the 3-cell full-diameter gas buffers

**User decisions 2026-10-31.**

## Background

`grain_geometry.build_snapped_geometry` adds `MIN_BUFFER_CELLS = 3` gas cells at
the head and at the aft end of every motor. Gas-only cells take the **full case
diameter**. So the head end is a case-diameter pocket 3·dx long, and its size
changes with the grid: for Chunc it is 24 mm at N=100 and 3 mm at N=800,
about 43 cm³ (≈15 % of the free volume) at N=100. The igniter injects into
cell 0, inside that pocket.

MU1 (`docs/solver_verification/worklog/MU1.md`, branch `mu1/muscl`) showed this
pocket is the whole cause of the startup's grid dependence. With a fixed
physical buffer the startup converges at N=100. The startup peak depends
strongly on the gap: Chunc forward plenum 10 g gives 23.6 MPa at 24 mm and
33 MPa at 3 mm. The aft buffer moves results about 1 %.

Real motors:
- The canonical motors (Chunc, BALLSstick, Hasegawa A) have a **zero head gap**:
  the forward closure and insulator are bonded to the grain face. All their
  grains are `inhibitedEnds: Both`, and the adapter already treats a
  both-inhibited interface as bonded (gap 0).
- Many other motors have head gaps of a few mm.

## Part 1 — head end: option A, lumped head volume

1. **No head buffer cells.** The first grain cell becomes cell 0, against the
   closure.
2. **The head gap is toggled by the head face's inhibition:**
   - head face inhibited (`Top` or `Both`): bonded, zero gap, nothing added;
   - head face uninhibited: the face must have room to burn, so the gap
     defaults to the inter-grain default `max(3 mm, 5 %·D_outer)`.
   - The gap's volume `A_case · L_gap` is **lumped into cell 0's storage**. No
     extra cells, so the result is grid-independent by construction.
3. **Condition (user):** do the lumping only if it does NOT require significant
   refactoring of the solver storage terms.
   - First, map where `A_port·dx`, or a cell mass or volume, enters the PISO
     continuity, pressure matrix, energy (sensible and total), species,
     momentum and the ledgers.
   - If a per-cell volume for cell 0 can be added cleanly in a few places, do it.
   - **If not, use the user's FALLBACK:** remove the head cells altogether and
     keep no lumped volume. The uninhibited head-face end burning (the
     existing end-face mass-injection hat function) and the igniter injection
     then start at the first cell. Say which path you took and why.
4. **Re-check the head-end paths with cell 0 as a grain cell:**
   - forward-plenum igniter injection: mass, enthalpy, momentum into cell 0 /
     face 1, and the `L_jet` distribution;
   - `head_grain_cell`;
   - the HE1 jet h (`x_from_head` is measured from the closure, which is now the
     grain face);
   - the head-basket and aft-basket cartridge cell ranges;
   - DeMar and radiation paths (opt-in);
   - the history and snapshot cell indexing, and the GUI station / slice payload
     (cell ranges, `is_grain`).

## Part 2 — aft end: remove the full-diameter aft buffer

The user's direction: remove it, or make it a function of the nozzle geometry.
Remove the default full-diameter cells and add only cells whose width follows the
nozzle's converging section, set by the convergence angle, `convAngle` in the
`.ric` nozzle. Assume the convergent starts at the case ID and ends at the throat.

**Stability is the main risk.** In the 2026-05-20 radiation-collapse audit, a
1-cell trailing buffer tripped a numerical front instability (Mach > 100) when
the trailing cell width was comparable to the throat diameter (see the
`MIN_BUFFER_CELLS` comment). Also, the nozzle BC treats cell N−1 as a
low-Mach reservoir: it uses the cell T as the stagnation T (EN1 open item).

Candidates, to test in this order:
1. **Cone cells:** cells following the convergent cone, truncated where the
   isentropic Mach at that area ratio would invalidate the reservoir treatment.
   Justify the cut from the BC's own assumption (for example the error in
   T/T0), not by tuning. Lump or neglect the small remaining volume, and say
   which.
2. **No aft buffer:** the last grain cell sits next to the nozzle BC.

Measure stability on every committed motor at its default settings, plus the
B0 rows.

## Part 3 — design study only, no implementation: inter-grain gap lumping and length snapping

Today, any inter-grain gap caps dx at the gap size, so BATES motors run on
about 3 mm cells (3-5x the cost). Segment lengths are rounded to whole cells
(up to dx/2 each), which moves propellant mass.

Write a short design memo, `docs/solver_verification/G1_NOTES.md`, covering:
- more elegant options, for example cut-cell / partial-cell finite volume,
  volume-porosity factors with uniform dx, and mass-preserving area scaling for
  the snapping;
- their stability implications for this PISO, given the history above;
- the literature you rely on, with DOIs marked verified or unverified;
- a recommendation.

## Gates and deliverables

- **B0:** this changes every startup, so B0 is re-recorded with the previous
  values in comments.
- **Full suite:** run it in chunks if needed to stay within the 4 GB guard.
- **Stability sweep** of all committed motors, each job guarded.
- **Report:** per canonical motor, before and after (head gap, aft treatment):
  first light, dark cells, startup peak and ratios vs experiment, Hasegawa A
  late peak, auto charges.
- **Specific check:** with the head pocket gone, do the head-end dead zones
  disappear even without HE1? Run with `head_end_jet_h` off and on.
- **Worklog:** `docs/solver_verification/worklog/G1.md`.
- Follow the CLAUDE.md compute budget (coexist mode).
