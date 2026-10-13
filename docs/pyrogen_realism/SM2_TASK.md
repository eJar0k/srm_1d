# SM2 — controlled comparison: Peretz surface-temperature ODE vs exact conduction

**User decision 2026-10-12.** Goodman is rejected as the eventual default. Before
choosing between `solid_heating_model='peretz'` (Peretz's integral-method ODE,
T_s as state) and `'conduction'` (exact 1-D conduction), the user wants them
compared on **lighting time** and **ignition-spike magnitude**, with the
confounding factors isolated. In the user's words: "There are confounding
variables here that don't isolate purely to one or the other ignition/heat
transfer criterion."

Naming: "Peretz ODE" means `solid_heating_model='peretz'`. "Peretz bore h"
means `ignition_heat_transfer='peretz'` (the default ignition h). Keep the two
apart in every table.

Read first: `SM1_REPORT.md` (models, verification, first results),
`IG3_TASK.md`, `OB1_REPORT.md` (current pyrogen defaults).

## Confounders to isolate

1. **Flow feedback (the main one).** In a full run, the solid model changes
   when cells light. That changes the flow, which changes h at every other
   cell. Separate the solid model from that feedback with an **open-loop
   replay**:
   - record per cell the surface forcing (`h_total(t)`, `T_driver(t)`, and the
     surface flux) from a reference run;
   - replay that identical forcing through each solid model offline;
   - compare per-cell ignition times.

   Record the forcing from a Goodman run, a Peretz-ODE run and a conduction run,
   so the forcing itself does not favour one model. Then compare the open-loop
   differences with the closed-loop ones; the gap between them is the feedback
   share.
2. **Calibration.** T_ignition 756 K and k_solid 0.271 were fitted under
   Goodman (v0.7.5 LHS). Run a small sensitivity on both models:
   - k_solid in {0.26, 0.271, 0.32} W/(m K). Published AP/HTPB(+Al) cold
     values span 0.20-0.40, centred on 0.25-0.30, and the project's tightened
     LHS range is 0.26-0.32;
   - T_ignition at three levels spanning the v0.7.5 LHS range (state the range
     you use and where it comes from).

   Ask: is the Peretz-ODE vs conduction difference larger or smaller than the
   spread from re-calibrating either model?
3. **Ignition-h closure:** Peretz bore h (default) vs `gnielinski`.
4. **Igniter charge.** Auto-sizing picks a DIFFERENT charge per model, so use
   **fixed explicit charges** in every comparison row: Chunc 2.62 g, BALLSstick
   4.83 g, Hasegawa A 14.98 g (IG3's auto charges under Goodman, forward
   plenum, library BPNV), plus Chunc 4.5 g in the head basket. Report the auto
   charges separately, as their own row.
5. **DeMar on / off** (default off since IG3).
6. **The ignition criterion itself.** Report, per cell, the surface heat flux
   and T_s at the moment of ignition. That shows whether a model lights a cell
   by a sustained flux or by holding a temperature (the Peretz-ODE failure mode
   SM1 found after an h drop).
7. **Grid.** For one case per motor, repeat at 200 cells. T8 found the
   100-cell grid under-resolves the startup peak by ~10%; state whether the
   model difference survives refinement.

## Metrics (every row)

First light, t50 and t90 (time to 50 % / 90 % of grain cells lit), whole-grain
lighting duration, cells never lit, startup peak, t_peak, plateau,
peak/plateau (BALLSstick: peak/trough, peak/mid), and Hasegawa A's late peak
vs 6.44 MPa.

## Deliverable

`docs/pyrogen_realism/SM2_REPORT.md`:
- the open-loop and closed-loop tables;
- the factor-by-factor attribution: how much of each metric's difference comes
  from the solid model, from feedback, calibration, h closure, charge, DeMar,
  grid;
- a recommendation between Peretz ODE and conduction. The user decides.

Probes and raw results go on your branch. No default changes.
