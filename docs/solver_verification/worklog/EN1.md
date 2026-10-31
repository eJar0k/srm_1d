# EN1 -- conservative total-energy formulation (opt-in) -- worklog

Branch `en1/total-energy`, worktree `C:/tmp/sv/en1`, rebased onto `openmotor-frontend` @ `ae86f9b`.
Executed 2026-10-13 and 2026-10-21 (resumed after a session limit).
Probes: `docs/solver_verification/probes/en1_verification.py` (harness + exact solutions), `en1_sod_table.py`,
`en1_duct_ladder.py` + `en1_duct_table.py`, `en1_motor_tables.py` + `en1_motor_summary.py`. Raw rows: `worklog/EN1_results/`.
Tests: `tests/test_total_energy.py` (19).

## STATUS / what is NOT done

* **B0 PASSED on the rebased base d592779 (MEM1) and the full suite is pending a single gate after the conduction-default flip** (see the report). Earlier evidence
  that the defaults are untouched: B0 PASSED on `c5f4af9` + EN1 (fresh compile, rel diff <= 1.6e-11, `n_steps` exact;
  Chunc/BALLS/Hasegawa A: 3.9e-15 / 3.2e-15 / 1.6e-11), and the Chunc fixed-charge sensible row on the rebased base
  reproduces the B0 Chunc row to the digit (8.698445703737, 533513 steps). Full suite before/after on `c5f4af9`: baseline
  runs were lost to Windows temp-file/cache races (no usable count); `tests/test_total_energy.py` 19 passed,
  `test_solver` + `test_energy_coupling` + `test_pressure_floor` 34 passed on the rebased base. **HEAVY-COMPUTE REQUEST** is
  in the final report.
* **Rebased-base motor rows (2026-10-21 / 10-31, guarded):** all BALLSstick rows and ALL Hasegawa A rows (auto and fixed charge,
  4 forms) were re-run on the rebased base (BALLSstick 2.3 GB peak; Hasegawa A 2.4 GB with the MEM1 history and 10 ms snapshots)
  and are IDENTICAL to the pre-rebase rows (same peaks, step counts, first light, never-lit counts to the digit). SM2/MEM1 changes
  are diagnostics/memory only. Every fixed-charge row has `pyrogen_mass_auto = False`.
* Supersonic-exit isentropic nozzle: not verifiable with the existing nozzle boundary (see section 2.2).

## 1. Formulation (opt-in kwargs, defaults bit-identical)

`run_simulation(..., energy_formulation='sensible'|'total', momentum_formulation='legacy'|'conservative')`;
kernel args `total_energy=False`, `conservative_momentum=False`. Defaults are the shipped step (B0 PASS).

### 1.1 Total energy (`_total_energy_update`, replaces STEP 3b when `total_energy`)

Per cell (between faces i and i+1):

    E = m (Cv T + KE),   Cv = Cp - R  (per-cell mixture),   m_new = rho_new_1 A dx   (the continuity mass, as before)
    E_new = E_old + dt (F_i - F_{i+1} + S_i dx)
    F_j  = mdot_j (Cp T + KE)_{upwind cell, start-of-step state},    mdot_j as in the sensible step (new u, rho_new_1 face avg)
    S_i dx = thermal_source[i] dx    (enthalpy of injected gas; injected at rest -> zero axial KE; wall sinks unchanged)
    Cv T_new = E_new/m_new - KE(u_new),  then the unchanged T floor / ceiling clips

Design choices and why:

* **Cell KE from the face velocities: `KE_i = (u_i^2 + u_{i+1}^2)/4`.** It is the mean of the face squares, i.e. exactly
  the share of the staggered (face-centred) kinetic energy `1/2 rho_f u^2` that belongs to the cell; positive-definite,
  exact for uniform flow, no new interpolation of the velocity. The square of the mean, `((u_i+u_{i+1})/2)^2/2`, would
  under-count the energy of a velocity gradient.
* **Flux = total enthalpy of the UPWIND CELL (`Cp T + KE_cell`)**, not the face velocity's `u_j^2/2`. This is the same
  rule the sensible step already follows (advect the cell's conserved specific quantity, as for species and Cp).
  Consequence: in a steady flow `Cp T_i + KE_i` is constant along the duct, so `T_i` is centred at the cell
  (`T_i = (H - KE_i)/Cp`) instead of being shifted half a cell by a face-velocity flux.
* **T recovery: `Cv T = E/m - KE(u_new)`** with the corrected velocities of the same PISO pass (after the Mach limiter and
  the nozzle-face velocity). The pressure work is then carried by the conserved variable, not dropped.
* **Coupling to the pressure correction and P5.** Unchanged: the pass loop is the same, the step still closes mass on
  `rho_new_1 = P/(R T_comp)` and `T_new` (now from E) feeds the next pass as `T_comp`. In a closed vessel `T_new`
  does not depend on P, so the single K=1 pass makes the EOS carry exact (EOS jump 0, table 2.4). `energy_coupling_iters`
  and `corr2_storage` keep their defaults.
* **Nozzle boundary.** Outflow carries the cell's total enthalpy `Cp T_up + KE_cell`; ambient inflow carries `Cp T_ambient`
  (reservoir at rest). The boundary mass-flow function (`_nozzle_boundary_flow`) is unchanged and still treats the cell T as
  the reservoir stagnation temperature (error `(g-1)/2 M^2` ~ 0.1 % at the Chunc nozzle-inlet cell, T6).
* **Ledgers.** Clip ledger `(T_clipped - T_raw) m Cv/dt`; energy residual closes to round-off (1e-6 W). `eos_energy_jump`
  and `floor_energy_added` are weighted by `Cv T_new + KE_new` in this mode (so the ledger identity holds); the simulation
  ledger helpers `_gas_mass_energy_total` / `_ledger_mass_energy` use the same energy. History channel names keep "sensible".

### 1.2 Conservative momentum (`conservative_momentum`), added beyond the task text

**Finding that forced it.** With total energy alone Sod is NOT exact (table 2.1). The shipped momentum predictor solves
`rho_f (u* - u^n)/dt = RHS`, i.e. `rho u_t = -(rho u u)_x - P_x`. That is the conservative law minus `u rho_t`: exact for
steady flow (and it contains the right dilution `-u S` of injected mass), but in a transient it is neither the conservative nor the
primitive form, so a captured shock moves at the wrong speed with wrong post-shock states on every grid (Sod S = 727 m/s vs 554,
independent of N). The opt-in predictor is `rho_f^{n+1,p} u* - rho_f^n u^n = dt RHS`, with the face density predicted from the
old-state discrete continuity (same face densities and nozzle boundary flow the correctors use). Steady flow is unchanged
(tables 2.2/2.3: identical to 3 digits).

## 2. Verification against exact solutions

BC mapping: head = wall (`u[0]=0`) + mass/enthalpy source in the first cell for the stagnation reservoir (exactly how
propellant burning enters); nozzle end = the motor's own `_nozzle_boundary_flow`: closed (`A_throat = 0`) for Sod and vessel;
for the ducts the subsonic-orifice branch with `A_throat = 3 A_exit` acts as a pressure outlet, with `P_ambient` chosen so
that the boundary passes exactly the exact mass flow from the exact last-cell state. Runs use the shipped defaults K=1,
corr2_storage, CFL 0.3. "sens" = sensible, "tot" = total, "+c" = conservative momentum.

### 2.1 Sod (gamma 1.4, P/rho 1e5/1 | 1e4/0.125, t = 0.2 L/sqrt(P_L/rho_L)); exact shock x = 0.8504, S = 554.1 m/s

L1 errors (scaled) rho / u / P / T and star-state values (numerical/exact; star-L: rho 0.4263, u 293.3, P 30313, T 247.7;
star-R: P 30313, T 397.7):

| form | N | shock x | S [m/s] | L1 rho | L1 u | L1 P | L1 T | star-L u / P / T | star-R T |
|---|---|---|---|---|---|---|---|---|---|
| sens | 100 | 0.957 | 742 | 0.058 | 0.307 | 0.062 | 0.118 | 398 / 26668 / 348.1 | 305 |
| sens | 400 | 0.966 | 742 | 0.062 | 0.324 | 0.064 | 0.117 | 408 / 25908 / 348.4 | 303 |
| sens | 1600 | 0.968 | 742 | 0.064 | 0.329 | 0.064 | 0.116 | 409 / 25811 / 348.4 | 301 |
| sens+c | 400 | 0.810 | 491 | 0.046 | 0.107 | 0.045 | 0.123 | 344 / 33620 / 348.4 | 280 |
| tot | 100 | 0.949 | 728 | 0.046 | 0.252 | 0.059 | 0.081 | 381 / 19062 / 222.4 | 247 |
| tot | 400 | 0.957 | 728 | 0.048 | 0.266 | 0.062 | 0.089 | 390 / 18195 / 218.0 | 226 |
| tot | 1600 | 0.959 | 727 | 0.049 | 0.269 | 0.063 | 0.090 | 391 / 18080 / 217.2 | 218 |
| tot+c | 100 | 0.842 | 542 | 0.0137 | 0.026 | 0.0120 | 0.0170 | 289.3 / 30937 / 252.5 | 400.0 |
| tot+c | 200 | 0.842 | 541 | 0.0092 | 0.019 | 0.0083 | 0.0125 | 289.6 / 30836 / 250.4 | 401.1 |
| tot+c | 400 | 0.842 | 541 | 0.0067 | 0.016 | 0.0062 | 0.0099 | 289.7 / 30813 / 249.7 | 401.0 |
| tot+c | 800 | 0.842 | 541 | 0.0051 | 0.0145 | 0.0050 | 0.0082 | 289.7 / 30805 / 249.3 | 400.9 |
| tot+c | 1600 | 0.842 | 541 | 0.0043 | 0.0143 | 0.0044 | 0.0071 | 289.7 / 30800 / 249.1 | 400.8 |

* Sensible form (either momentum): T stays at 348.4 K in the expansion fan (exact 247.7): the missing pressure work.
  Wrong jump conditions, grid-independent (L1 does not decrease with N; observed order ~0).
* Total energy with legacy momentum: still wrong (S 728 vs 554), grid-independent. The error is the momentum form.
* Total + conservative momentum: L1 falls with N (observed orders 0.2-0.6, the usual 0.5-1 for a shock/contact in L1).
  The residual ~1.6 % plateau error is a TIME-step error: at N = 400, CFL 0.3 / 0.15 / 0.075 / 0.0375, star-L P = 30813 /
  30639 / 30503 / 30415 (exact 30313) and S = 541.0 / 545.5 / 548.8 / 551.6 (exact 554.1): first order in dt
  (the P5 K=1 splitting), no spatial inconsistency.
* Conservation: mass drift (K=1) -5.8e-4 at CFL 0.3, -1.0e-4 at CFL 0.0375; with K=3 it is -8e-7 (the pressure-based
  splitting's EOS carry, pre-existing, not EN1). Total energy drift -1.4e-3 (K=1) / -6e-4 (K=3), vs +8.3 % to +8.5 % for the sensible form.

### 2.2 Steady isentropic quasi-1-D flow

Duct: smooth converging-diverging, A/A* = 8 inlet, 2 exit, throat at x = 0.4, stagnation reservoir via source in cell 0
(P0 = 5e5, T0 = 1000 K). Run to steady state (t = 0.06 s, residual <= 2e-4 per 5 ms), started from the exact solution.

*Subsonic venturi, throat Mach 0.8* (T/T0 = 1/(1 + 0.2 M^2)); L1 normalised by T0/P0/max M:

| form | N | L1 M | L1 T | L1 P | T/T0 throat num / exact | T/T0 exit num / exact |
|---|---|---|---|---|---|---|
| sens | 100 / 200 / 400 | 0.118 / 0.130 / 0.135 | 0.0452 / 0.0452 / 0.0452 | 0.092 / 0.096 / 0.098 | 1.0000 / 0.8872 | 1.0000 / 0.9828 |
| tot | 100 | 0.0072 | 0.0011 | 0.0110 | 0.8977 / 0.8872 | 0.9826 / 0.9828 |
| tot | 200 | 0.0035 | 0.0005 | 0.0055 | 0.8912 / 0.8867 | 0.9829 / 0.9829 |
| tot | 400 | 0.0017 | 0.0002 | 0.0027 | 0.8879 / 0.8865 | 0.9830 / 0.9830 |
| tot+c | 100 / 200 / 400 | 0.0072 / 0.0035 / 0.0017 | 0.0011 / 0.0005 / 0.0002 | 0.0110 / 0.0055 / 0.0027 | same | same |

Total energy converges at first order and the momentum form is irrelevant in steady flow. The sensible form never cools
(T = T0 everywhere) and its Mach / pressure errors do NOT decrease with refinement.

*Supersonic part* is verified in 2.3 (the isentropic region upstream of the shock, M up to 1.8, `L1up` columns).
**A supersonic-exit nozzle (pressure outlet at M = 2.19) could not be verified**: the pinned-pressure boundary is ill-posed for
supersonic outflow (no area makes the orifice branch pass the flow; with a large area the total-energy exit cell, whose KE
uses the boundary-set face velocity, goes unstable within 7 steps from the EXACT initial state). The motor never has a
supersonic last cell, so this is a limitation of the test boundary, not of the interior scheme; it is not a verdict on either form.
Also: with an area of 30 A_exit the same boundary destabilises the total-energy exit cell during a back-pressure ramp
(`T` of the last cell collapses); 3 A_exit is stable. Worth knowing if the nozzle boundary is ever made stiffer.

### 2.3 Normal shock in the divergent section (exact station A/A* = 1.5: x_s = 0.8243, M1 = 1.854, M2 = 0.605, p02/p01 = 0.7884)

Initial state: exact solution with the shock at A/A* = 1.25 (weaker), back pressure then ramped smoothly to the target so
the shock must travel; t = 0.06 s. Post-shock quantities are read 8 cells downstream of the shock vs the exact value there.

| form | N | shock x num / exact | post-shock P/P0 | post-shock T/T0 | post-shock M | upstream (isentropic) L1 M / T / P | residual (5 ms) |
|---|---|---|---|---|---|---|---|
| sens | 100 | 0.943 / 0.824 | 0.228 / 0.666 | 1.000 / 0.953 | 1.53 / 0.50 | 0.080 / 0.137 / 0.148 | 2e-5 |
| sens | 400 | 0.954 / 0.824 | 0.256 / 0.632 | 1.000 / 0.939 | 1.46 / 0.57 | 0.086 / 0.156 / 0.141 | 7e-5 |
| tot | 100 | 0.8102 / 0.8243 | 0.6663 / 0.6662 | 0.9530 / 0.9530 | 0.496 / 0.496 | 0.0101 / 0.0049 / 0.0112 | 6e-7 |
| tot | 200 | 0.8174 | 0.6453 / 0.6455 | 0.9444 / 0.9445 | 0.542 / 0.542 | 0.0055 / 0.0027 / 0.0062 | 2e-6 |
| tot | 400 | 0.8209 | 0.6316 / 0.6318 | 0.9387 / 0.9387 | 0.572 / 0.571 | 0.0030 / 0.0015 / 0.0034 | 2e-6 |
| tot+c | 100 / 200 / 400 | 0.8102 / 0.8174 / 0.8209 | identical to tot to 4 digits | | | identical | 3e-9 |

Total energy: shock position error 0.0141 / 0.0069 / 0.0034 (first order), post-shock P, T and Mach equal the exact
Rankine-Hugoniot values to 3-4 digits, the supersonic isentropic region converges at first order. The sensible form puts the
shock at the wrong station, keeps T = T0 (no cooling, so post-shock T is the stagnation T), and its post-shock pressure is
off by 3x; none of it improves with N.

**Resolved: the "non-converging total + legacy shock case" was a harness bug, not physics.** `duct_case` looped
`while tube.t < t_end` around an inner `advance_to` that stopped at `t_end - 1e-12`; when a chunk ended a sliver short the outer
loop spun forever appending a residual each pass (the 2.4 h / 30 GB job of 2026-10-13). Fixed: chunk count fixed up front,
preallocated residual array, exact landing on `t_end`, plus a step cap, wall-clock cap and dt floor in `Tube.advance_to`
(a stalled/diverged case is returned with `stall_reason`). With the fix the same case runs in 6.6 s / 1.1 GB at N = 100,
t = 0.06, and converges (residual 2e-6 per 5 ms). My earlier "denormals slow the kernel" explanation was wrong.

### 2.4 Closed vessel (A = 1e-3, L = 1, 20 cells; P0 = 1e5, T0 = 300 K; inject T_in = 1200 K, m_in = 2 m0 over 10 ms)

Exact (adiabatic, uniform): T = (m0 Cv T0 + m_in Cp T_in)/(m Cv) = 1220 K, P = 1.22 MPa. The sensible form's limit
(`T = mass-weighted mean of T0, T_in`) is 900 K / 0.90 MPa.

| form | dt | T num / exact | P num / exact | mass num / exact |
|---|---|---|---|---|
| sens | 4e-5, 1e-5, 2.5e-6 | 900.000 / 1220.000 (all dt) | 0.900000 / 1.220000 MPa | exact |
| tot | 4e-5, 1e-5, 2.5e-6 | 1220.000 / 1220.000 (1e-13 relative, all dt) | 1.220000 / 1.220000 MPa | exact |

The total form is exact to round-off at every dt (T_new is independent of P inside the step, so the K=1 pass closes
the EOS carry); the sensible form misses the pressure-work heating, a factor 1220/900 = gamma-ish in the rise.

## 3. Effect on the motors (canonical knobs, library BPNV, Peretz ignition h default, DeMar off, IG3 auto charge)

Charges: Chunc 2.6226 g, BALLSstick 4.8292 g, Hasegawa A 14.9765 g (the auto-size result of the sensible form), passed
explicitly with `pyrogen='bpnv'`; every fixed row has `pyrogen_mass_auto = False` and `pyrogen_mass_initial` = the charge
(driver asserts it). "total" = total energy + legacy momentum; "tot+c" adds conservative momentum; "sens+c" = conservative
momentum alone.

### 3.1 Startup peak / plateau / ignition (B0 definitions)

| motor | form | peak [MPa] @ ms | plateau / trough / mid | ratio | first light / t90 [ms] | cells never lit |
|---|---|---|---|---|---|---|
| Chunc | sens (B0) | 8.698 @ 24.5 | plateau 8.066 | 1.078 | 4.44 / 16.3 | 0 |
| Chunc | total | 9.034 @ 24.8 | 7.160 | 1.262 | 5.21 / 5.92 | 7 (cells 3-9) |
| Chunc | tot+c | 9.080 @ 24.3 | 7.160 | 1.268 | 5.23 / 5.83 | 7 |
| Chunc | sens+c | 7.880 @ 22.6 | 6.385 | 1.234 | 4.40 / 5.19 | 13 |
| BALLS | sens (B0) | 7.193 @ 22.4 | trough 3.522, mid 3.888 | pk/tr 2.042, pk/mid 1.850 | 5.02 / 5.68 | 9 |
| BALLS | total | 7.703 @ 22.0 | 3.689 / 4.089 | 2.088 / 1.884 | 5.09 / 5.78 | 7 |
| BALLS | tot+c | 7.803 @ 21.6 | 3.689 / 4.090 | 2.115 / 1.908 | 5.03 / 5.84 | 7 |
| BALLS | sens+c | 7.282 @ 22.1 | 3.522 / 3.889 | 2.067 / 1.873 | 4.97 / 5.75 | 9 |
| Hasegawa A | sens (B0) | 5.362 @ 2.68 s | | | 4.20 / 5.75 | 8 |
| Hasegawa A | total, 14.98 g | **does not ignite** (0/106 lit, peak 1.2 MPa) | | | n/a | all |
| Hasegawa A | total, auto 18.72 g (4 attempts) | 6.140 @ 22.1 ms | | | 2.29 / 3.51 | 1 |
| Hasegawa A | tot+c, 14.98 g | 5.010 @ 2.84 s | | | 8.82 / 10.04 | 12 |
| Hasegawa A | sens+c, 14.98 g | 6.106 @ 2.40 s | | | 3.97 / 5.40 | 0 |

Experiment (ad hoc tree, not committed): driving the wall heat transfer with the stagnation temperature `T + u^2/2Cp`
instead of the static T in the total form changed Chunc first light 5.21 -> 4.17 ms but left 10 cells unlit, so the dead zones
are not simply a lower driving temperature. In the near-miss cells T_surf peaks at 726-757 K against the 756 K threshold
at 4-5 ms, i.e. whether a cell lights is a marginal event, and which cells light varies by +-6 cells between ANY two
formulations (sens vs sens+c differ as much as sens vs total). Treat dead-zone counts as unconverged, not as a discriminator.

Plateau: the total form lowers the Chunc plateau 8.07 -> 7.16 MPa (-11 %). That is the 7 unlit head-end cells (the grain
that never burns), not a nozzle-BC effect (`T6`: nozzle-inlet KE ~ 0.1 %). Experimental Chunc: peak 8.66, plateau ~8.43, ratio
1.027; so the plateau moves away from experiment for the same reason.

### 3.2 Fill (Chunc, 1 ms snapshots, 100 cells; max |Mach| / x of max [mm] / supersonic cells / downstream end of the supersonic region [mm])

| t [ms] | sensible | total | tot+c |
|---|---|---|---|
| 1 | 1.04 / 386 / 5 / 406 | 1.18 / 410 / 14 / 528 | 1.22 / 410 / 10 / 463 |
| 2 | 1.70 / 825 / 29 / 829 | 1.42 / 825 / 12 / 829 | 1.94 / 825 / 36 / 837 |
| 3 | 1.28 / 264 / 28 / 366 | 1.31 / 158 / 24 / 333 | 1.39 / 386 / 51 / 553 |
| 4 | 1.34 / 386 / 42 / 480 | 1.31 / 158 / 34 / 415 | 1.44 / 337 / 27 / 350 |
| 5 | 0.78 / - / 0 | 1.33 / 289 / 21 / 309 | 1.33 / 224 / 15 / 252 |
| 6 | 0.36 / - / 0 | 0.94 / 817 / 0 | 0.77 / 809 / 0 |
| >= 7 | subsonic (0.35-0.39) | subsonic | subsonic |

Whole-fill max Mach (before the peak): sensible 1.95, total 2.12, tot+c 2.67 (BALLS 1.84 / 1.93 / 2.48; Hasegawa A 3.05 /
3.15 / 4.16). So total energy does NOT remove the supersonic fill; the supersonic region lives about the same time
(through ~5 ms), is a bit longer-lived in the total form, and ends in the same place (the contact moves to the nozzle end).
Bore-gas T at the first-lit cell (cell 11-13): sensible 2261 K (4.0 ms) -> 2458 K (5.0 ms); total 2340 K (5.0 ms) -> 2757 K (6.0 ms);
tot+c 2336 -> 2817 K. The total form's gas is NOT colder at the lit cells; the first light is ~0.8 ms later (5.2 vs 4.4 ms).
The 1 ms snapshot spacing is coarse for a 3 ms fill; the fill-Mach trend above is a trend, not a shock-position measurement.

### 3.3 T6 energy ledger, Chunc fill phase (P < 50 % of peak)

| | sensible | total | tot+c |
|---|---|---|---|
| EOS jump, sum abs [J] | 12.4 | 10.8 | 10.0 |
| EOS jump / fill source energy (17.9 / 18.7 / 18.0 kJ) | 0.07 % | 0.06 % | 0.06 % |
| EOS jump, whole run, sum abs [J] | 151 | **21** | 19 |
| all storage jumps / fill source | 6.8 % | 5.9 % | 6.3 % |
| of which mixture refresh [J] | 1178 | 1077 | 1095 |
| T-clip (T ceiling) [J] | -1693 | -2524 | -2518 |
| energy residual max [W] | 5e-6 | 5e-6 | 5e-6 |

The T6 "30 % of the fill energy is the EOS jump" figure was the pre-T11 state; after the T11 corrector fix the EOS jump is already
0.07 % of the fill source, and the total form changes that little in the fill (the T lag is now removed by K=1 already). Over
the whole run the EOS jump shrinks 7x (151 -> 21 J). The remaining ledger gap in the fill is the species/mixture refresh (Cp
changes at fixed rho, T: 6 %) and the T-ceiling clip (-9.4 % of the fill source for sensible, -13.5 % for total): neither
is touched by EN1.

### 3.4 Cost (Chunc, 1 s, fixed charge, cache-loaded, sequential, 3 repetitions, machine shared with ~6 other busy processes)

| form | steps | steps/s |
|---|---|---|
| sensible | 533513 | 36.5k, 38.0k, 37.4k |
| total | 535037 | 30.6k, 33.0k, 31.6k (-15 %) |
| tot+c | 535066 | 30.6k, 32.5k, 32.5k (-15 %) |

## 4. Recommendation

1. **Keep `energy_formulation='sensible'` as the default for now.** The total form is the correct conservation law and the
   verification is unambiguous (vessel exact, Sod/shock/venturi converge, sensible does not), but on the motors it
   re-opens the head-end ignition (dead zones 0 -> 7 cells on Chunc, Hasegawa A fails to light at the sensible charge) and the
   Chunc plateau falls 11 %: the validated knobs (roughness, kappa, T_ignition, k_solid) were fit with the sensible form, so a flip needs the LHS re-fit
   and a decision on the ignition-side gas temperature (static vs recovery), which is HE1's territory. Do not flip before that.
2. **If the user wants shocks/fill modelled correctly, adopt both opt-ins together** (`'total'` + `'conservative'`): the
   verification matrix shows the energy form alone leaves the transient shock speed wrong by 31 %, and the momentum form alone
   leaves T wrong by 40 %. `'conservative'` momentum alone is also a cheap diagnostic of how much of the motor's transient is
   the legacy predictor.
3. The cost is -15 % steps/s. The EOS-jump motivation of thread F is already spent by T11; the case for the total form is the
   physics of the fast fill (static T, shocks), not the ledger.
4. Not changed and worth a decision: the wall-heat-transfer driving temperature in the total form (static T; a recovery
   temperature would be the physically standard choice) and the nozzle boundary's use of the cell T as stagnation T.

## 5. Files

`srm_1d/solver.py` (`_total_energy_update`, args `total_energy`, `conservative_momentum`), `srm_1d/simulation.py` (kwargs,
`_gas_mass_energy_total`, `_ledger_mass_energy`, summary stamps), `srm_1d/DEVNOTES.md` (API log), `tests/test_total_energy.py`,
the probes above, `worklog/EN1_results/` (sod.json, duct_*.json, motors/*.json).
Process note: all runs after 2026-10-21 went through `guarded_run`; the motor driver sets `history_capacity` per motor
(1.7 GB peak for Chunc instead of 3.9 GB).
