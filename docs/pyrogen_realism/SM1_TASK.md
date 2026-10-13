# SM1 — solid-heating model comparison (Goodman vs Peretz ODE vs exact conduction)

**User go-ahead 2026-10-12.**

## Why

A probe of BALLSstick (4.8 g BPNV, forward plenum, DeMar off, main `3cedd4a`)
found grain cells 3-11 never light during the 2 s burn, although the gas over
them is 2200-3000 K. Part of the cause is real flow: the closed head end
stagnates once the downstream grain lights, and h drops 5-15x.

Part of it is the solid model. `solid_thermal._compute_T_surf` is Goodman's
constant-h heat-balance-integral solution,
`T_s = (3k T_i + h delta T_g) / (3k + h delta)`. It is algebraic in the
CURRENT h at the current penetration depth delta. When h fell 14x, cell 9 went
746 -> 349 K in 2 ms, discarding ~85 % of the heat stored in the solid. Exact
conduction would cool far more slowly (reviewer's estimate ~570 K after 2 ms).
The same closure makes T_s jump UP instantly when h rises. That may feed the
whole-grain-lights-in-1-2-ms cascade IG1 found behind the startup spike.

Probe scripts: the reviewer's scratch probes are reproduced in Part 3.

## Part 1 — the models (opt-in; default stays Goodman, bit-identical)

New kwarg `solid_heating_model = 'goodman' (default) | 'peretz' | 'conduction'`.

- **goodman**: today's kernel.
- **peretz**: Peretz's integral-method surface-temperature ODE, carrying T_s as a
  state. It is in `peretz1974.pdf` (thesis, Ch. III) and `peretz1972.pdf`
  (AIAA 72-1119, Eq. 15-17). The text layer is garbled, so read the equation
  from page images (single pages, <= 110 dpi, cropped). Record exactly what it
  assumes about time-varying h.
- **conduction**: an exact 1-D transient conduction solve in the solid, per
  unignited grain cell.
  - Implicit (backward Euler or Crank-Nicolson), tridiagonal.
  - A graded grid near the surface.
  - Domain depth several thermal penetration depths for the ignition window.
  - Boundary condition: the same total surface flux the kernel applies today
    (convective h(T_g - T_s) plus any pyrogen / radiation surface flux). Check
    how those extra fluxes enter the Goodman path now and treat them the same.
  - Properties: the same k_solid and alpha the sim uses. No new parameters.
  - Choose the node count by a convergence study, and report it.

## Part 2 — verification against exact solutions (all three models)

1. Constant h, constant T_g: the exact semi-infinite Robin solution (erfc-exp,
   Carslaw & Jaeger).
2. Constant flux: the exact sqrt(t) solution.
3. Step-DOWN in h (x1/14 at t1, as at cell 9) and step-UP in h (x10): reference
   = the converged conduction solve, cross-checked by Duhamel superposition.

Report T_s(t) error and time-to-756 K error per model. The conduction model
must match the exact solutions to the stated tolerance before it is used in
Part 3.

## Part 3 — effect in the sim

Use main's current defaults: Peretz ignition h (local), OB1 pellets and basket.
Run each solid model on:

1. **BALLSstick dead zone**: explicit 4.8 g BPNV, forward plenum,
   `heat_delivery_mode='none'`. Report which grain cells never light, and the
   ignition time per cell. Reviewer's numbers to reproduce with goodman: cells
   3-11 never lit; cells 12-53 light at 5.1-5.45 ms.
2. **Whole-grain lighting time** (first light -> 90 % lit) and peak, on:
   - Chunc realistic head basket (4.5 / 6 g x bpnv / mtv);
   - the three canonical motors with their usual igniters, as the B0 rows.
3. **Ignition thresholds** (head basket and forward plenum, DeMar off), as in
   IG2 section 4.1.
4. Wall-time cost of `conduction` vs `goodman`.

Calibration caveat to state: T_ignition 756 K and k_solid 0.271 were fitted
with Goodman + Gnielinski ignition h (v0.7.5 LHS).

## Part 4 — head-end note (no code)

The user points out that a 1-D model cannot represent the turbulent mixing and
recirculation at a real head end, which likely deliver hot gas to the surface
far faster. While reading Peretz:

- Extract **Eqs. III-28, III-29, III-30** verbatim from the page images. These are
  the Cohen, Derr & Price correlations of the Rocketdyne igniter-duct data:
  impingement region, upstream and downstream of the igniter-jet impingement
  point (x_im about one port diameter).
- Note how Peretz used them (Tig, the effective igniter-gas temperature) and what
  inputs they need (igniter nozzle diameter? jet Re?).
- Write them into `docs/pyrogen_realism/HEAD_END_NOTES.md` for a follow-up task.
  The fuller source is Carlson & Seader, AIAA J. 5(7) 1967, DOI 10.2514/3.4183
  (paywalled; R9).

## Gates and deliverables

- B0 PASS with defaults. Full suite green.
- `docs/pyrogen_realism/SM1_REPORT.md`: the models, the verification tables, the
  sim tables, the cost, and a recommendation (the user decides the default).
- `HEAD_END_NOTES.md`.
- Opt-in implementation, tests (exact-solution checks, default identity) and
  probes on your branch.
- Paywalled-reference rule applies.
