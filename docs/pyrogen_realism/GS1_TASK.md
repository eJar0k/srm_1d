# GS1 — grid and contact-front convergence of the startup (picks the LHS grid)

**User decision 2026-10-21.** The order is 1. RT1 (recovery-temperature ignition
driver), 2. this study, 3. HE1 (head-end jet heat transfer), 4. the LHS re-fit.
The user triggers the LHS later, overnight on an idle machine. GS1 picks the
grid it runs on.

## Why

- **Knife-edge results.** The startup results sit on a knife edge. Head-end
  dead cells end 1-12 K below the 756 K ignition temperature, and the dark-cell
  count swings by about 6 cells between formulations. Auto charges sit at the
  ignition threshold.
- **The 100-cell grid is too coarse.** T8 found it under-resolves the startup
  spike by about 10 %. In SM2 the Chunc head-basket peak went from 15.5 MPa at
  100 cells to 18.5 MPa at 200.
- **A suspected missing limiter for ignition travel: the contact surface.**
  This is the front between hot igniter or combustion gas and the cold gas
  ahead. First-order upwind advection smears a contact front without bound; a
  shock, by contrast, self-steepens. Hot gas may therefore leak ahead of the
  real front and preheat cells. That would fit the downstream grain lighting
  almost at once (first light to 90 % lit in about 2-3 ms) while the stagnant
  head end stays dark.
- **Grid is not an LHS axis.** If it were, T_ignition and k_solid would absorb
  discretization error. Converge the startup first, then fit.

## Configurations

- Use main's defaults at the time: the conduction solid model (it becomes the
  default before GS1), Peretz ignition h, DeMar off, MIL-P-46994B I-A pellets.
- If RT1 has merged, add `energy_formulation='total'`,
  `momentum_formulation='conservative'` and the recovery driver as a second
  configuration.
- **Clean rows**, charges well above threshold, from SM2's high-charge set:
  Chunc forward plenum 10 g, BALLSstick 10 g, Hasegawa A 30 g, plus Chunc head
  basket 4.5 g (the real igniter).
- **Threshold rows**: the fixed IG3 charges (Chunc 2.62 g, BALLSstick 4.83 g,
  Hasegawa A 14.98 g), only to show how the knife edge behaves under refinement.
- **Grids**: 100, 200 and 400 propellant cells. An 800-cell grid needs a
  HEAVY-COMPUTE REQUEST.
- **Startup metrics** use a short `t_max` (0.1-0.3 s). Full burns are needed
  only for plateau and late-peak values: at 100 and 200 cells, and at 400 cells
  only if it fits the budget.

## Metrics, per row and grid

- First light, t50 and t90.
- Ignition-front speed: the speed of the lit boundary, in m/s versus position.
- Contact-front width vs time: the distance between 10 % and 90 % of the gas-T
  rise, and of the propellant or pyrogen species mass fraction, in mm and in
  cells. Compare its growth with the first-order expectation, width ~ sqrt(dx·u·t).
- Startup peak, t_peak, plateau and peak/plateau. For BALLSstick, the
  trough/mid ratios. For Hasegawa A, the late peak vs 6.44 MPa.
- The dark-zone end as a POSITION (mm), not a cell count.
- Observed convergence order and Richardson-extrapolated values, wherever three
  grids allow it.

## Deliverables

- `docs/pyrogen_realism/GS1_REPORT.md`, with:
  - the tables;
  - the convergence verdict per metric;
  - a recommended LHS grid. Give the smallest grid within the tolerances you
    state for the startup metrics; the user picks the tolerance. Include its
    cost per run (wall time and peak memory from the guard), so the overnight
    LHS can be sized.
  - whether the contact smearing is real and how much of the ignition-front
    speed it explains;
  - a recommendation on a higher-order (MUSCL/TVD) advection for T and species
    as the follow-up fix. GS1 measures only; it does not implement it.
- Probes and raw results on your branch.
- No default changes.

## Budget

Follow the CLAUDE.md compute budget: coexist mode, 4 GB, 2 processes, 30 min
per job, every job through `srm_1d.tools.guarded_run`. Measure one 400-cell run
first. Anything over the budget is a HEAVY-COMPUTE REQUEST.
