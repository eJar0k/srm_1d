# NB1 -- the nozzle boundary uses the true (isentropic) stagnation state

Branch `nb1/stagnation-bc` (from `openmotor-frontend` 7b56144). Spec: `docs/solver_verification/NB1_TASK.md`.
Probes: `probes/nb1_motor_ab.py`, `nb1_ledger.py`, `nb1_cone_reservoir.py`, `nb1_dump_loss.py`; tests `tests/test_nb1_nozzle_boundary.py`.

## What changed

* `solver._nozzle_boundary_flow(P, T, A_throat, gamma, R, P_amb, T_amb, A_exit, M_choke)`.
  The last cell's (P, T) are the STATIC state at the port exit (flow area `A_face[N] = A_port[N-1]`).
  * Choked: `M_e` = subsonic area-Mach root of `A_exit/A_throat`; `P0 = P g(M_e)`, `T0 = T(1 + (g-1)/2 M_e^2)`,
    `mdot = Gamma A_t P0 / sqrt(R T0)`; `dmdp = mdot/P` exactly (P0 is linear in P at fixed `M_e`).
  * Subsonic outflow: root of `G(M) = A_e P M sqrt(gamma/RT) - mdot_isentropic(P0(M), T0(M) -> P_amb) = 0` on the bracket
    `[0, M_choke]` (`G(0) < 0 < G(M_choke)`), safeguarded Newton with bisection fallback, <= 100 iterations. `dmdp` is the
    implicit derivative `mdot/P - c_e G_P/G_M` (closed form, no finite differences). The choked/unchoked switch compares
    `P_amb` with `crit * P0` (P0 of the moving gas), so it is continuous (mdot and dmdp both).
  * `A_exit <= A_throat`: no subsonic root; the port exit is the sonic section (`M_e = 1`, controlling area `A_exit`,
    `mdot = rho a A_e`). Because `crit * g(1) = 1`, every outflow (`P > P_amb`) is then choked. No NaN, no crash.
    This case is stiff by nature (flow switches from 0 to sonic at `P = P_amb`); it only arises for a port exit smaller than the throat.
  * Reverse inflow (reservoir at rest): unchanged, bit for bit.
* `M_choke` is computed ONCE per PISO step (`_nozzle_choke_mach`; 0.0 = the old static boundary) and passed to every call
  site: momentum-predictor density prediction, both pressure correctors (flow and `dmdp`), the exit-velocity update, the
  energy update (sensible and `_total_energy_update`), the accounting call, and the post-PISO bookkeeping call in
  `simulation.py`. One boundary state everywhere.
* Energy flux trap: the function still returns the cell's STATIC T as `upstream_temperature`. The callers' outflow flux
  `mdot (Cp T_static + KE)` is already the stagnation enthalpy; returning T0 would have double-counted the kinetic energy.
* `run_simulation(nozzle_bc_stagnation='isentropic')` (default) | `'static'` (the pre-NB1 boundary, bit-identical to AS1:
  B0 rows reproduce the AS1 values to all printed digits). Threaded through `run_from_ric(**sim_overrides)`, the plugin
  schema (choice field, default from the signature), the result summary. `piso_step` / `_piso_step_with_energy_diagnostics`
  gained `nozzle_isentropic=True`. DEVNOTES API log updated.
* The EN1 duct harness (`probes/en1_verification.py`) uses the static orifice branch as a pressure outlet
  (`A_throat = 3 A`), so it now passes `nozzle_isentropic=False` and `M_choke=0.0` explicitly.

## Verification against exact answers

1. Unit tests (`tests/test_nb1_nozzle_boundary.py`, 73 tests, gamma 1.1/1.2/1.4, area ratios 1.05..1000):
   * area-Mach root vs independent bisection (1e-12 on `A/A*`);
   * choked mdot vs `Gamma A_t P0/sqrt(R T0)` to 1e-10, equal to `rho u A_e` at the port exit, implied `P0/P = g(M_e)` to 1e-10;
   * subsonic point self-consistent (exit continuity = throat continuity, 1e-9);
   * mdot continuous across the switch to 1e-7, dmdp to 1e-4 (two-sided 1e-9 offset);
   * dmdp vs central finite difference on both branches (1e-5);
   * the switch uses `P_amb/P0` not `P_amb/P` (a pressure between the two thresholds is choked new / subsonic old);
   * `A_exit <= A_throat` (1.0, 0.9, 0.3 x A_t): finite, sonic, `mdot = rho a A_e`; continuous through `A_e = A_t`;
   * reverse flow identical in both modes; `'static'` mode equals the old formulas to 1e-14; `A_exit -> inf` recovers static.
2. Steady duct with uniform mass addition (exact: `P + rho u^2 = P_head`, `T0 = const`, choked throat; A/A_t = 1.6,
   M_e = 0.4046, gamma 1.2, mdot 0.8 kg/s, L 0.5 m; `_run_duct`, total energy + conservative momentum, driven to steady state):

   | N | P_head error, isentropic | P_head error, static | exit-face Mach new / old (exact 0.4046) |
   |---|---|---|---|
   | 25  | -1.106 % | +5.56 % | 0.4046 / 0.3700 |
   | 50  | -0.563 % | +6.04 % | " |
   | 100 | -0.285 % | +6.29 % | " |
   | 200 | -0.143 % | +6.42 % | " |

   New: first-order (halves per refinement; the exit cell centre is dx/2 short of the exit). Old: a constant ~6.4 % miss, as
   predicted (`f^(k-1/2)` = 1.095 on the exit static pressure). The new exit-face Mach equals the exact one to 4 digits.
3. Conservation (`probes/nb1_ledger.py`, Chunc 0.1 s, D1 defaults): the closure measures are unchanged by the new boundary:

   | measure | static | isentropic |
   |---|---|---|
   | PISO energy residual / peak source | 1.2e-12 | 1.2e-12 |
   | source ledger mismatch / peak source | 1.2e-15 | 1.3e-15 |
   | nozzle enthalpy ledger mismatch / peak nozzle power | 6.2e-5 | 6.9e-5 |
   | nozzle mass mismatch / peak mdot | 2.5e-4 | 2.7e-4 |
   | run mass-balance error | 1.06e-3 | 0.99e-3 |

   (The nozzle mismatches are the pre-existing PISO-vs-post-PISO sampling difference, same size in both modes.)

## Canonical motors, N = 100, committed .ric, canonical knobs (static vs isentropic; probe `nb1_motor_ab.py`)

| motor / metric | static (= AS1 B0) | isentropic | change |
|---|---|---|---|
| Chunc peak [MPa] | 10.5104 | 9.8263 | -6.5 % |
| Chunc plateau (0.4-0.9 s) | 8.0348 | 7.9173 | -1.5 % |
| Chunc peak/plateau (exp 1.027) | 1.308 | 1.241 | |
| BALLSstick peak | 8.4329 | 8.0817 | -4.2 % |
| BALLSstick trough | 4.4314 | 4.3310 | -2.3 % |
| BALLSstick mid (1.03 s) | 4.9252 | 4.8380 | -1.8 % |
| BALLSstick peak/trough, peak/mid | 1.903, 1.712 | 1.866, 1.670 | |
| Hasegawa A late peak (exp 6.44) | 6.0172 @ 2.50 s (0.934x) | 5.9609 @ 2.53 s (0.926x) | -0.9 % |
| Hasegawa A early peak (t < 0.5 s) | 5.8626 @ 0.169 s | 5.3051 @ 0.172 s | -9.5 % |

Grain-exit Mach and stagnation ratio at the plateau (M_e from the run's own mass flow and exit-cell state; window Chunc 0.4-0.9 s,
BALLSstick 1.0-1.1 s, Hasegawa A 2-3 s):

| motor | A_e/A_t | M_e static | P0/P static | M_e isentropic | P0/P isentropic | P_head / P0_exit (isentropic) |
|---|---|---|---|---|---|---|
| Chunc | 3.83 | 0.154 | 1.0143 | 0.157 | 1.0148 | 1.0295 |
| BALLSstick | 4.66 | 0.126 | 1.0096 | 0.128 | 1.0099 | 1.0343 |
| Hasegawa A | 5.27 | 0.112 | 1.0074 | 0.113 | 1.0076 | 1.0114 |

Explanation of the moves. The task text quoted M_e ~ 0.35 (-> 7 %); the real grain-exit Mach at the plateau is 0.11-0.16
(port/throat 3.8-5.3), so the steady-state correction is 0.8-1.5 %, and the plateau moves by exactly that (Chunc -1.5 %).
The old boundary said `P0 = P`, so for the gas generated the chamber had to sit g(M_e) higher. The moves are bigger where the
port Mach is higher (the startup fill): Chunc peak -6.5 %, Hasegawa A early peak -9.5 %, BALLSstick peak -4.2 %. The effect is
in the direction the task predicted (lower static pressure for the same generation, largest in the transient) and is a
correction of the nozzle model, not a re-tune: no constant was touched. Auto igniter charges are unchanged (0.815 / 1.982 / 7.861 g).
The Hasegawa A late peak moves from 0.934x to 0.926x of experiment (the plateau-level correction is -0.9 %).

## Step 5: G1's reservoir assumption (`bc_reservoir_area_ratio`, `BC_RESERVOIR_P0_TOL`)

Probe `nb1_cone_reservoir.py` (Chunc, `aft_end='cone'`, plateau 0.4-0.9 s, MPa; `aft_end='lump'` isentropic reference 7.917):

| boundary | cut A_BC/A_t | cone cells | plateau | peak |
|---|---|---|---|---|
| static | 4.64 (G1 default, 1 % P0) | 1 | 8.0239 (lump 8.0348, -0.14 %) | 10.066 |
| isentropic | 4.64 | 1 | 7.9378 (lump 7.9173, +0.26 %) | 9.995 |
| isentropic | 1.5 | 2 | 7.5114 (-5.1 % vs lump) | 9.515 |

Finding. The cut was written as "keep the BC's static-for-stagnation error below 1 %". With the isentropic boundary that error is
gone, but the cut does NOT become unnecessary: loosening it to 1.5 (second cone cell, area ratio 2.25, D 40.8 -> 27.2 mm in ONE
cell at the 45-degree cone) moves the plateau by -5 %, the mirror image of the +5 % G1 measured with the static boundary at N = 200.
The remaining error is the resolution of a steep convergent by one or two cells (area changes by 2x across a cell; the BC
places the cell-centre state at the cell-centre area while mass leaves through the face), a first-order O(dx dA/dx) discretisation
error, not the stagnation assumption. At the default cut both boundaries agree with the lump reference to 0.3 %.
Recommendation: KEEP the function and the default 1 % cut, but REPURPOSE the rationale (a resolution guard on the area drop per
cone cell, not a reservoir-assumption bound) and reword its docstring/DEVNOTES; do not loosen it. A principled replacement
would bound the area ratio between neighbouring cone cells (not a Mach tolerance) -- a design choice for the user. Nothing deleted.

## Step 6 (report only): aft dump loss (Borda-Carnot)

Port exit (A_e, M_e) dumping into the case-ID plenum (A_c = pi/4 D_case^2, the grain OD at the aft cell). Control-volume momentum
balance with the step-face pressure P_e (no friction, adiabatic, no constants):
`P_c(1 + g M_c^2) = P_e(1 + g (A_e/A_c) M_e^2)`, mass, and conserved T0; solved for M_c; loss = 1 - P0_c/P0_e
(`probes/nb1_dump_loss.py`; the incompressible Borda-Carnot `(g/2) M_e^2 (1 - A_e/A_c)^2` agrees to 2 %).

| motor | M_e | A_e/A_c | isentropic correction P0/P - 1 | dump loss dP0/P0 |
|---|---|---|---|---|
| Chunc | 0.157 | 0.557 | 1.48 % | 0.29 % |
| BALLSstick | 0.128 | 0.668 | 0.99 % | 0.11 % |
| Hasegawa A | 0.113 | 0.953 | 0.76 % | 0.002 % |

The dump loss is 5-20 % of the isentropic correction on Chunc/BALLSstick and nil on Hasegawa A (port nearly case ID): it raises the
plateau by the same fraction (0.1-0.3 %), below the validation precision, so it does not matter at the plateau for these motors.
It scales with M_e^2: it would matter for a high-Mach port into a much larger plenum (high-L/D, small port/throat).

## Step 7 (report only): thrust and performance pressure

`nozzle.compute_motor_performance` and `compute_thrust_history` take `P_head` (cell 0) as the chamber pressure:
`thrust = C_F(P) * P * A_t`. The head-end static equals the head-end stagnation pressure (u = 0), but the stagnation pressure at
the nozzle entrance is lower, by the momentum-flux loss along the port (`P + rho u^2 = const`) -- it is the one that feeds the
throat. From the table above `P_head / P0_exit` = 1.030 (Chunc), 1.034 (BALLSstick), 1.011 (Hasegawa A): the thrust path
over-predicts thrust and impulse by that much at the plateau (the same static-as-stagnation family of error, from the head end
instead of the exit cell). Consistent fix (not implemented): record a history channel `P0_nozzle = P_exit g(M_e)` from the
boundary solve (cheap: `M_e` is known per step) and let `compute_motor_performance` use it when present (fall back to `P_head`).
P_head stays the right quantity to compare with a head-end transducer. C* and Isp (which use C_F and c*) are unaffected.

## Gates

* New + touched tests: `tests/test_nb1_nozzle_boundary.py` (73), `tests/test_solver.py`, `tests/test_solver_config_upgrade.py` pass.
* Full suite: see the results block at the end of this file.
* B0 re-recorded in `probes/baseline_gate.py` (AS1 values kept in comments); the same code reproduced the recorded values in
  the A/B probe to all digits.

## Surprising / worth the user's attention

* Plateau shift is 0.8-1.5 %, not the 7 % the task text implied: the grain-exit Mach is 0.11-0.16, not 0.35.
* The startup moves a lot more than the plateau (-6.5 % Chunc peak, -9.5 % Hasegawa A early peak).
* The EN1 harness had to be pinned to the static orifice (it is a pressure-outlet device, `A_e < A_t`).
* The nozzle-convergent resolution sensitivity (step 5) is first order and survives the BC fix.

## Results block (on base 7b56144 + NB1 commits)

* Full pytest suite under the guard: 687 passed in the first chunk (stopped by `-x` at `test_solver_config_upgrade.py`, which
  asserts the exact set of plugin config keys; updated for the new `nozzle_bc_stagnation` choice field), then the remaining
  files from `test_solver_config_upgrade.py` on: 333 passed. Total 1020 passed, 0 failed (22 + 7.4 min).
* B0 after the re-record: see the final line below (re-run on the committed tree).
