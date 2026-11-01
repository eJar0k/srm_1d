# G1 Part 3 -- inter-grain gap lumping and length snapping (design study, no implementation)

Status: memo only. Nothing here is implemented. What G1 DID implement (head gap and nozzle convergent lumped into a cell's
gas storage, `A_vol = A_port + vol_extra/dx`) is the mechanism this memo builds on; its measured behaviour is in
`worklog/G1.md`.

## 1. The two problems today

**(a) A gap caps dx.** `build_snapped_geometry` sets `dx = min(L_prop / N_target, min_gap)` when any inter-grain gap exists
(the "Nyquist clamp", `max_allowed_dx = min_gap / 1.0`). The default inter-grain gap is `max(3 mm, 5 % D_outer)`, so a BATES motor
with a 3 mm gap runs on 3 mm cells: a 0.8 m propellant column becomes about 270 cells instead of 100, and the step count
scales with the cell count (CFL) times the cell count again, so the cost is about 7x (the acoustic time step shrinks 2.7x and there are
2.7x the cells). The gap itself holds a volume `A_case * L_gap` of a few cm^3.

**(b) Lengths are rounded to whole cells.** `n_seg = max(1, round(L/dx))`, `snapped_length = n_seg * dx`. Each segment end moves by
up to dx/2 (per segment: up to dx), and the errors accumulate along the motor. Propellant mass, burning area and the gas volume
all follow the rounded length. At dx = 8 mm a 4-segment motor can carry a 1-3 % propellant-mass error that the cell count
(not physics) decides.

## 2. Options

### 2.1 Lump the gap volume into the neighbouring cells (volume porosity with uniform dx)
Remove the gap from the axis: the grain faces on either side of it become neighbouring cell edges, and the gap's gas volume
`V_gap = A_case L_gap` is added to the storage of the adjacent cells (split 50/50 so reverse flow during ignition is
symmetric). In the PISO this is exactly the `vol_extra` mechanism G1 already installed: storage terms use `A_vol dx`, the flow
areas stay `A_port`. Properties:

* conservative by construction (it is a storage term, no flux is changed);
* no dx clamp, so a BATES motor runs at the target cell count;
* what it drops: the gap's axial length (x of downstream grains shifts by `L_gap` per gap, harmless to the solver, but the
  Peretz x-dependence, ignition-front distances and the viz axis need an `x_offset` array if they should stay physical), the gap's momentum
  (a 3 mm region: negligible), and the sudden-expansion / contraction losses at the gap, which the model does not carry today either;
* the end-face (uninhibited) burning faces then sit at cell edges against a lumped volume; the existing hat function plus the wall clauses
  (verified in `tests/test_g1_ends.py::TestEndFaceAtWalls` for the wall case) already conserves face mass, and a face that regresses into the
  lump simply moves the lump boundary. The growth of the gap volume with the receding face (`A_case dL`) must be added to `vol_extra` per step
  (a small, local extension; G1 holds `vol_extra` constant).

### 2.2 Mass-preserving area scaling for the length rounding
Keep the cell grid; keep `n_seg dx` as the modelled length `L'` but scale the burning perimeter of that segment's cells by `f = L / L'`
(and the end-face area is untouched). Then the burning surface `C L'` equals the true `C L`; the swept volume `C L dr` over the web,
hence the propellant mass burned `rho * integral`, equals the true value; the port (flow) area is unchanged. The residual error is only the
gas volume `(L' - L) A_port` (<= dx/2 per segment end) and the axial position (<= dx/2). Rounding the cumulative positions (not the
lengths) keeps the position error below dx/2 for the whole motor rather than accumulating. No new state, no stability implication (a
multiplicative factor on the existing `C_burn`, which already carries `grain_frac * f_active`).
The only subtlety is that the mass flux in the erosive closure uses the generated flow, which is now `f`-consistent, so nothing else changes.

### 2.3 Cut-cell / partial-cell finite volume
Let the cell edges stay uniform and let a cell that contains a segment end or a gap be a partial cell with a volume fraction
`alpha = V_gas/(A dx)` and a face area fraction (the classical Cartesian embedded-boundary method; Colella et al. 2006, Quirk 1994).
In 1-D this collapses to the porosity/lumping form above, with `alpha` varying per cell and per face, plus the known small-cell problem:
an explicit scheme's stable step scales with `alpha dx`. Remedies in the literature are flux redistribution, cell merging and the h-box method
(Berger, Helzel and LeVeque 2003; Berger and Helzel 2012). In one dimension cell merging IS the lumping of 2.1 (merge the tiny cell into its
neighbour and keep the combined storage), so the heavy machinery is not needed here. Full cut-cell with flux redistribution is not recommended.

### 2.4 Non-uniform grid (dx per segment)
Choose `dx_k = L_k / n_k` per segment so every segment is exact and the gap is one cell of its own width. It removes both problems
but touches every `dx` in the PISO, the advection, the ledgers, the CFL and the geometry kernels (`dx` is a scalar in about 300 places) and
introduces cell-size jumps at the interfaces (first-order accuracy loss unless smoothed). Rejected: large refactor, and the "numerical
Goldilocks" history (MU1; the 2026-05-20 radiation-collapse audit) shows the solver is sensitive to the size ratio of neighbouring cells.

## 3. Stability implications for this PISO

* The pressure equation's transient coefficient is `a_t = A_vol dx / (R T dt)`. Adding storage to a cell increases `a_t` only: it makes the
  pressure matrix MORE diagonally dominant, the acoustic response of that cell slower. It cannot create a new instability by itself.
  G1 measured this at the head (cell 0) and the aft (last cell) of the three canonical motors and a 20-motor stability sweep (`worklog/G1.md`).
* The explicit parts (upwind energy / species advection) see a cell of volume `V = A_vol dx`; a lump enlarges V, which only relaxes
  their step restriction (a lump in a cell never shrinks it, unlike a cut cell with alpha << 1).
* What a mid-motor lump changes that the ends do not: a flow passes THROUGH the cell, so the velocity ratio between the neighbouring faces
  changes by the lumped fraction; the face areas stay `A_port`, so continuity is exact and the momentum predictor sees the unchanged face
  areas. A gap that is 5 % of the case diameter is about 1/3 to 1 cell worth of volume at the target grid: far from the "small cell"
  regime of cut cells.
* The rounding fix (2.2) touches no solver term.
* History to respect: a single trailing cell comparable to the throat diameter tripped a numerical front (2026-05-20 audit; the
  `MIN_BUFFER_CELLS` comment). The remedy here is never to make a cell narrow relative to its neighbours; both 2.1 and 2.2 keep dx uniform.

## 4. Recommendation

1. **Do 2.2 first** (mass-preserving perimeter scaling + cumulative-position rounding): purely geometric, no stability risk, removes the
   dx-dependent propellant-mass / burn-area error. Gate: propellant mass and burn area equal the .ric values to round-off at every N; B0 unchanged
   for motors whose lengths are already multiples of dx.
2. **Then 2.1** for inter-grain gaps (drop the dx clamp, lump `V_gap` 50/50 into the two adjacent cells, add the gap-growth term for uninhibited
   faces, keep an `x_offset` array for the physical axis). Gate: BATES test motors (`BATES Test.ric`, `example_bates.ric`, `75_5120ChaseRed_BATES.ric`,
   `Aeroblast_flight_2.ric`) at N = 100 against the current clamped N ~ 270 results (pressure trace within the grid-convergence band of GS1/MU1),
   the stability sweep of G1 repeated, and the end-face conservation test extended to a face regressing through a lump. Expected cost drop: 3-5x on gapped motors.
3. Do not do 2.3 (full cut cell) or 2.4 (non-uniform grid).

## 5. Literature

DOIs were checked against Crossref (title, journal, volume, pages) where marked verified; "unverified" = recalled, not checked. None of the numbers
of any paper are used here, only the method names.

| reference | DOI | status |
|---|---|---|
| Colella, Graves, Keen, Modiano, "A Cartesian grid embedded boundary method for hyperbolic conservation laws", J. Comput. Phys. 211(1), 347-366 (2006) | 10.1016/j.jcp.2005.05.026 | verified |
| Quirk, "An alternative to unstructured grids for computing gas dynamic flows around arbitrarily complex two-dimensional bodies", Comput. Fluids 23(1), 125-142 (1994) | 10.1016/0045-7930(94)90031-0 | verified |
| Berger, Helzel, LeVeque, "h-box methods for the approximation of hyperbolic conservation laws on irregular grids", SIAM J. Numer. Anal. 41(3), 893-918 (2003) | 10.1137/s0036142902405394 | verified |
| Berger, Helzel, "A simplified h-box method for embedded boundary grids", SIAM J. Sci. Comput. 34(2), A861-A888 (2012) | 10.1137/110829398 | verified |
| Issa, "Solution of the implicitly discretised fluid flow equations by operator-splitting", J. Comput. Phys. 62(1), 40-65 (1986) (PISO) | 10.1016/0021-9991(86)90099-9 | unverified (Crossref rate-limited) |
| Chern, Colella, "A conservative front tracking method for hyperbolic conservation laws", LLNL report UCRL-97200 (1987) (flux redistribution / cell merging origin) | none (report) | not checked; not needed |
| LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge (2002) (flux-form conservation, small-cell CFL) | 10.1017/CBO9780511791253 | unverified |

Needs for the user to fetch (paywalled or not checked) if the 2.1 implementation goes ahead: none required for the method itself; the
two items worth reading before implementing are the h-box paper (to confirm that the 1-D merging argument above matches the small-cell CFL result) and Issa
for the PISO storage-term form.
