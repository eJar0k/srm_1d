# TP1 -- thrust and performance from the nozzle-entrance stagnation pressure

Branch `tp1/thrust-p0` (from `nb1/stagnation-bc` 0f70601 + the task commit). Spec: `docs/solver_verification/TP1_TASK.md`.
Probes: `probes/tp1_motor.py`, `probes/tp1_pathfin_compare.py`; tests `tests/test_tp1_thrust_pressure.py`;
results: `worklog/TP1_results/*.json`, `worklog/tp1_pathfin54.png`.

## What changed

* **Recorded state.** `solver._nozzle_entrance_stagnation(mdot, P_cell, T_cell, gamma, R, A_exit, M_choke) -> (P0, T0, M_e)`.
  It is called right after the post-PISO `_nozzle_boundary_flow` of each step, with that call's own `mdot` and `M_choke`
  (hoisted into one local so both calls share it); there is no second root solve. For outflow the boundary's mass flow is
  `mdot = rho_e u_e A_e = c_e M_e`, `c_e = A_e P sqrt(gamma/(R T))` (exact in the subsonic branch; the choked branch has
  `M_e = M_choke` to round-off, clamped), so `M_e = mdot / c_e`, `P0 = P g(M_e)`, `T0 = T (1 + (gamma-1)/2 M_e^2)`.
  Defined values elsewhere: `'static'` boundary (`M_choke = 0`), no flow, and reverse inflow (reservoir at rest, nothing feeds the
  thrust path through the port) all return the static cell state `(P_exit, T_exit)`, i.e. `'static'` gives `P0 = P_exit`
  exactly as the spec asks. `A_exit <= A_throat` gives `M_e = 1`.
* **History.** `P0_nozzle`, `T0_nozzle`, plus the exit-cell mixture `gamma_nozzle`, `R_nozzle` that the boundary used
  (needed for the gamma and C* findings): four more `np.empty(max_hist)` arrays, same capacity, same trim, same
  `history_capacity` re-run logic as every other per-step array (MEM1: nothing per-step in a Python list; +32 B per row).
  `_run_time_loop` gets four new positional arguments (DEVNOTES log).
* **Thrust path.** `compute_motor_performance` takes `P0_nozzle` (fallback to `P_head` with a `RuntimeWarning` if absent),
  and always the solver's own `D_throat` history when present (it did before, but only after a wasted full pass with the
  old formula). New `compute_thrust_series(P0, D_throat, ...)` (thrust from a GIVEN throat history);
  `compute_thrust_history` keeps its signature and is now only the fallback integrator. New perf keys `P0_nozzle`,
  `c_star_sim`, `c_star_sim_avg`.
* **Consumers.** Everything downstream of `perf` follows automatically: plotting (`plot_thrust`, summary, comparison),
  the openMotor plugin's `force` channel (`perf['thrust']`), designation, average thrust, burn time, the CSV `Force`/`Isp`
  columns. New channels `P0_nozzle`, `T0_nozzle`, `gamma_nozzle`, `R_nozzle`; CSV gains a last performance column
  `Nozzle Stagnation Pressure (Pa)` (the earlier column positions are unchanged; `load_openmotor_csv` reads by name).
* **What stays on P_head.** openMotor's `pressure` channel (`sr.channels['pressure']` in `_result_to_om_simresult`, the one
  the GUI plots), the CSV `Pressure (Pa)` column, `plot_pressure`/summary pressure traces, the Kn-equilibrium column
  (`compute_grain_metrics`), all of `tools/` fitness (the comparisons with head-end transducers). openMotor derives its own
  impulse / max pressure statistics from the `force` and `pressure` channels, so they follow the thrust fix with no change.

## Results (Chunc, BALLSstick, Hasegawa A: committed .ric, canonical knobs, library BPNV auto charge, N = 100)

One run per motor: TP1 changes no dynamics (B0 below), so "before" is the old formula (`compute_thrust_series` fed with `P_head`
and the same throat history) and "after" is `compute_motor_performance`. Plateau windows: Chunc 0.4-0.9 s, BALLSstick 1.0-1.1 s,
Hasegawa A 2-3 s, pathfin54 1-2 s.

| motor | metric | before (P_head) | after (P0_nozzle) | change |
|---|---|---|---|---|
| Chunc | total impulse [N s] | 3337.5 | 3226.3 | -3.33 % |
| Chunc | average thrust [N] | 3363.8 | 3254.4 | -3.25 % |
| Chunc | peak thrust [N] | 3918.8 | 3795.3 | -3.15 % |
| Chunc | average Isp [s] | 243.08 | 242.74 | -0.14 % |
| BALLSstick | total impulse [N s] | 7169.1 | 6752.4 | -5.81 % |
| BALLSstick | average thrust [N] | 3605.4 | 3380.0 | -6.25 % |
| BALLSstick | peak thrust [N] | 5898.5 | 4593.2 | -22.1 % |
| BALLSstick | average Isp [s] | 223.65 | 222.51 | -0.51 % |
| Hasegawa A | total impulse [N s] | 22192.8 | 21569.9 | -2.81 % |
| Hasegawa A | average thrust [N] | 5174.4 | 5028.7 | -2.81 % |
| Hasegawa A | peak thrust [N] | 7219.3 | 7136.8 | -1.14 % |
| Hasegawa A | average Isp [s] | 205.08 | 204.89 | -0.09 % |
| pathfin54 | total impulse [N s] | 2542.3 | 2528.3 | -0.55 % |
| pathfin54 | average thrust [N] | 865.6 | 860.9 | -0.54 % |
| pathfin54 | peak thrust [N] | 1110.1 | 1105.2 | -0.44 % |
| pathfin54 | average Isp [s] | 174.39 | 174.33 | -0.03 % |

`P_head` / `P0_nozzle` (median): plateau Chunc 1.0183, BALLSstick 1.0345, Hasegawa A 1.0110, pathfin54 1.0020. The task text quoted
1.030 / 1.034 / 1.011 for the first three (NB1's measurement); BALLSstick and Hasegawa A reproduce it, Chunc is lower here
(1.018): this branch carries G2 + AW1 (NB1's table was on the earlier base), and the ratio also falls over the burn
(Chunc 1.024 at 0.4-0.6 s, 1.012 at 0.9-1.0 s). An independent reconstruction (NB1's `exit_state` method: `M_e` from the run's
own mass flow, the snapshot's exit-cell T, port area, gamma, R) gives `P0/P_exit = 1.00951` against the recorded `1.00957`
(6e-5), `M_e = 0.126`, `A_e/A_t = 4.75` (the last cell's port area, the one the boundary uses).
The startup gap is much larger: `P_head/P0` is up to 1.2 for Chunc after 10 ms (3.4 inside the first 10 ms while the head
pressure waves are still crossing), up to 4.4 for BALLSstick (with P_head > 1 MPa), i.e. most of the startup pressure hump is a head-end
phenomenon that never reaches the nozzle. That is why BALLSstick's peak thrust falls 22 %.
`P0 <= P_head` holds on every plateau. It does NOT hold at every instant: in the first ~4 ms (P_head < 0.4 MPa, head cell
still expanding, waves not yet crossed) the exit-cell stagnation pressure exceeds the head pressure by up to 1.9x. That is
unsteady-flow physics, not an accounting error (1218 of 531 719 Chunc rows; the unit test checks the plateau only).

### C* consistency (plateau median of `P A_t / mdot`)

| motor | C*_prop (rep. tab) | C*_sim with P0 | with P_head | T0 vs T_flame |
|---|---|---|---|---|
| Chunc | 1572.3 | 1570.5 (-0.11 %) | 1599.3 (+1.72 %) | 3105.3 / 3105 |
| BALLSstick | 1535.8 | 1536.2 (+0.02 %) | 1589.1 (+3.47 %) | 2977.7 / 2977 |
| Hasegawa A | 1543.1 | 1543.2 (+0.01 %) | 1560.1 (+1.10 %) | 3039.4 / 3041 |
| pathfin54 | 1317.4 | 1308.2 (-0.70 %) | 1310.9 (-0.49 %) | 2202.1 / 2202 |

With the stagnation pressure the simulated C* equals the propellant's to 0.1 % on three motors; with the head-end pressure it
is inflated by exactly `P_head/P0`. That is the cleanest evidence that `P0_nozzle` (not `P_head`) is the pressure that
drives the throat: the mixture C* `sqrt(R T0)/Gamma(gamma_mix)` at the exit cell equals the P0-based C*_sim to 1e-7 (it is
the boundary identity), and T0 equals the flame temperature to 0.1 %. pathfin54's -0.7 % is its scalar-vs-mixture gamma
(below). The Isp definition (`CF c*_prop / g0`) is unchanged and only moves through CF(P) (-0.03 to -0.5 %): CF depends on
the exit-pressure ratio, not on the absolute pressure, so the large effect is in thrust, not Isp.

## Findings (report only, nothing changed)

* **Throat erosion.** The solver DOES erode the throat itself: `_run_time_loop` step 3b integrates
  `dD/dt = 2 (erosion_coeff * P[N-1] - slag_coeff / P[N-1])` with the port-exit STATIC pressure each step, feeds the new
  `A_throat` to the boundary and records `D_throat_hist`. The post-processor was already using that recorded history
  whenever present (always, for `run_simulation` results); the old code first ran `compute_thrust_history` on `P_head`
  (a second, different integration) and discarded it. That dead second pass is gone: `compute_motor_performance` now calls
  `compute_thrust_series` with the solver's history, and `compute_thrust_history` is only the fallback for results with no
  `D_throat`. So there was no actual thrust inconsistency, but the second integration used a different pressure and was
  misleading. Remaining family-resemblance: the in-loop erosion uses the exit-cell static pressure, the stagnation pressure
  is 0.1-0.9 % higher on these motors; with openMotor's convention (rate proportional to the chamber pressure) the erosion
  would be ~1 % faster. Chunc is the only one of the four with erosion (184 um/(s MPa): 17.02 -> 19.90 mm in 1 s); the
  effect would be ~0.03 mm. Not changed: it would alter the flow solution (B0 must stay unchanged), and it is a modelling
  choice about which pressure drives wall erosion, not an obvious defect.
* **Gamma.** The CF formula keeps the propellant's scalar gamma. The exit-cell mixture gamma (recorded in `gamma_nozzle`)
  starts at the air value (1.40) and is above the propellant value by more than 3 % for only 3.3-6.0 ms (Chunc 3.5 ms,
  BALLSstick 3.4, Hasegawa A 6.0, pathfin54 3.5); first-50-ms median is within 0.4 % of the plateau value on the three
  large motors. Replacing the scalar by the per-step mixture gamma changes the impulse by -0.11 % (Chunc), +0.011 %
  (BALLSstick), +0.006 % (Hasegawa A) and -0.58 % (pathfin54), and the first-50-ms impulse by < 0.1 %. The only visible
  difference is a plateau one: pathfin54's mixture gamma at the exit is 1.2648 against the tab's 1.24 (the tab gamma
  and the Cp/M of its species disagree), worth -0.58 % in thrust; the other motors agree to 0.3 %. Not clearly wrong,
  not changed. Instantaneous CF differences up to 22-28 % exist in the first 2-3 ms, at thrusts of ~50 N (1 % of peak).
* **C*.** See the table: C*_sim(P0) = C*_prop; C*_sim(P_head) is inflated by `P_head/P0`.

## pathfin54 (`motors/2025.12.25 54-2800 ST2.0.ric`, committed; `static_fire_data/pathfin54.csv`)

Run: the motor's own igniter block (`BPNV`, auto charge 0.857 g), canonical defaults, `t_max` 6 s, `P_cutoff` 0.1 MPa
(ran to t_max = 6 s, 3.9 M steps; thrust is below 5 % of peak by ~2.95 s). Measured pressure at t = 0 is 90 kPa (absolute), F = 80 N: the trace starts mid-ignition.
The measured trace ends at 2.765 s with F = 82 N and 1.08 MPa (truncated tail: thrust already down to 7 % of peak); the sim thrust
reaches zero at ~3.15 s. The sim is aligned to the measured rising edge (20 % of own peak) and the integrals are also compared over the
same measured span (`probes/tp1_pathfin_compare.py`; figure `tp1_pathfin54.png`):

| quantity | measured | sim before | sim after | after vs measured |
|---|---|---|---|---|
| total impulse, same span [N s] | 2668.4 | 2510.7 | 2496.9 | -6.4 % |
| total impulse, whole sim [N s] | (2668.4 truncated) | 2542.3 | 2528.3 | -5.2 % |
| average thrust, same span [N] | 965.1 | 910.4 | 905.3 | -6.2 % |
| peak thrust [N] | 1175.7 | 1110.1 | 1105.2 | -6.0 % |
| burn time, thrust > 5 % of peak [s] | 2.765 (truncated) | 2.937 | 2.937 | |
| peak pressure [MPa] | 10.12 | 9.24 (P_head) | 9.24 (P_head); P0 9.20 | -8.6 % |
| thrust / pressure at the plateau [N/MPa] | 117.4 | 119.2 (/P_head) | 119.0 (/P_head), 119.2 (/P0) | +1.4 % (/P_head) |

What it can and cannot tell:

* The time-resolved comparison does not discriminate the fix. The sim under-predicts the PRESSURE (peak -8.6 %, and a
  different shape: it is high during the first 0.85 s and low during 0.9-2.2 s), so thrust is low by the same ~6 %. That
  gap is in the burn-rate / grain-burnback side (a finocyl) and is untouched by TP1; the fix moves the sim 0.5 % further
  from the measurement (-5.9 % to -6.4 % over the same span), which is below that gap.
* Thrust over pressure is the quantity that tests the nozzle path independently of the burn rate: measured 117.4 N/MPa
  (median over F > 0.8 F_max), sim 119.0 against the head-end pressure and 119.2 against `P0_nozzle` (CF A_t = 1.318-1.320 against
  1.301 measured). The sim is 1.4-1.6 % high, i.e. the whole CF chain (ideal CF with the propellant gamma, the 0.9 user
  efficiency, divergence/throat/skin factors) is consistent with the measurement to ~1.5 %, with either pressure.
* The transducer location is unknown, so the measured F/P cannot say which simulated pressure it matches. This motor has a
  wide port (`P_head/P0 = 1.002`), so the P_head-versus-P0 choice moves F/P by only 0.2 % here, far inside the
  1.5 % agreement: pathfin54 neither supports nor contradicts the fix. The discriminating evidence is the C* identity
  above (and the physics), which is what the fix rests on. A high port-Mach motor with measured thrust AND a head-end
  pressure (Chunc-like, `P_head/P0 = 1.02-1.03`) would be the test; the existing Chunc trace has no thrust.
* The measured trace is gauge-vs-absolute consistent with the sim ambient (90 kPa at t = 0), so there is no offset to
  correct.

## Gates

* **Unit tests** (`tests/test_tp1_thrust_pressure.py`): static mode gives `P0 = P_exit`; choked `P0/P = g(M_choke)` to
  1e-10 and equals the boundary's own `Gamma A_t P0 / sqrt(R T0)` flow (gamma 1.1/1.2/1.4, area ratios 1.05..20); subsonic
  outflow recovers `M_e` and the throat flow from `(P0, T0)` to 1e-8; `A_exit <= A_throat` gives `M_e = 1`; no-flow / reverse
  flow return the static state; P0 = P_head reproduces the old thrust and CF point by point; thrust follows `P0_nozzle`, not
  `P_head`; the fallback warns and uses `P_head`; the fallback throat integrator runs on the series it is given;
  `compute_thrust_series` equals `compute_thrust_history` for the integrated throat; a real short run records all four
  arrays at every step (`len == n_steps <= history_capacity`, finite), `P0 >= P_exit`, `P0 == P_exit` where there is no
  outflow, perf thrust equals the P0-based recomputation, and the CSV keeps `Pressure (Pa)` as column 3 with the new
  stagnation column.
* **B0** (`probes/baseline_gate.py`): unchanged. The same B0 on the base (0f70601, copied out with `git archive`) and on this
  branch gives identical values to every printed digit (Chunc peak 11.09137373410309, t_peak 0.02284186258807952, plateau
  8.156904979898151, 531719 steps; BALLSstick 8.410209348265113 / 4.333772377796931 / 4.844986906954673, 846027;
  Hasegawa A 5.957335057263664 at 2.426343433577292 s, 2337265 steps). The `BASELINE` recorded on this branch is the
  pre-NB1 (G2 + AW1, static boundary) one, so the gate as shipped here reports MISMATCH; NB1's own re-record
  (`697cbdf`: 11.091373734103097 / 8.156904979898155 / 531719; 8.410209348265097 / 4.3337723777969295 / 4.8449869069546745 /
  846027; 5.957335057263658 / 2.4263434335772907 / 2336970) equals these values to within 1e-14 (Hasegawa A steps
  2337265 vs 2336970: 1.3e-4, the documented compile-state noise, inside the 5e-4 gate).
* **Full suite:** see the results block below.

## Results block

* Full suite (pre-merge, 0f70601 + TP1, clean run under the guard): **1081 passed**, 0 failed, 19 min, peak 4.1 GB.
* Merged `openmotor-frontend` (NB1 landed, 4e4c796; the merge brings only docs and NB1's B0 re-record, no code): clean
  `__pycache__`, **B0 PASS** against NB1's new BASELINE (floats within 2e-15; Hasegawa A n_steps 2337265 vs 2336970, 1.3e-4,
  inside the 5e-4 gate). TP1 + NB1 + channels + adapter + plugin + config-upgrade tests after the merge: 183 passed.
* The motor tables above are from the pre-merge tree; the merge changes no code and B0 reproduces the same values, so they stand.
