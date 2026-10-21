# RT1 -- opt-in ignition driver temperature (adiabatic-wall recovery temperature) -- worklog

Branch `en1/total-energy` (on `openmotor-frontend` @ ae86f9b), 2026-10-21. Probes: `probes/en1_motor_tables.py --driver --solid`,
`probes/rt1_summary.py` (grid table), `probes/rt1_regime.py` (regime at ignition). Raw rows + `summary.txt`: `worklog/RT1_results/`.
Tests: `tests/test_rt1_recovery.py`.

## Implementation

`run_simulation(..., ignition_driver_temperature='static'|'recovery')`; the default `'static'` is bit-for-bit the shipped ignition heating.
In the ignition kernel the convective driver was the cell's static T in three places: the driver numerator `h_c T`, the `h_driver` handed to
the solid-heating model, and the wall-loss sink `h_c (T - T_s) C_burn` (plus the SM2 recorder `rec_Td`). All now use

    T_driver = T                                  (static)
    T_driver = T_aw = T + r u^2 / (2 Cp)          (recovery),    u = cell velocity, Cp = cell mixture Cp

with the textbook recovery factor (Schlichting, Boundary-Layer Theory; White, Viscous Fluid Flow; Kays and Crawford, Convective Heat
and Mass Transfer): `r = Pr^(1/3)` in a turbulent boundary layer, `r = Pr^(1/2)` in a laminar one.

* Regime per cell: turbulent when the bulk Reynolds number `Re = rho |u| D_h / mu` of the cell exceeds the pipe-flow transition value 2300,
  laminar otherwise (the `Re` array the Gnielinski path already uses; no new knob).
* Pr: the one the active ignition closure uses. Peretz (default): `Pr = gamma/(1.77 gamma - 0.45)` with the cell's mixture gamma (Eq. III-25).
  Gnielinski and Keller: the run's propellant-gas Pr.
* Only the driver temperature changes. The heat-transfer coefficient is still evaluated at the static / film state of each closure. The
  DeMar/pyrogen flux, adjacent-cell radiation and the energy-equation sinks are untouched. The Ma 2020 erosive path is untouched.
* `ignition_diag` (SM2) is now (N, 14): columns 10-13 = driver T, r, bulk Re, turbulent flag at the cell's ignition step.

Flow regime in practice (`rt1_regime.py`, 30 ms runs, fixed charges): EVERY cell that lights is turbulent (Re 1.9e4 to 1.0e6 at ignition) and
r = 0.894-0.897 (Pr ~ 0.72), so the laminar branch is never active at ignition in Chunc, BALLSstick or Hasegawa A. The driver rise `T_aw - T` at
ignition is 0-570 K (median: Chunc sensible 216 K, Chunc total 89 K, BALLS 304-313 K, Hasegawa A about 100 K); ignition Mach number up to
1.35-1.69 on Chunc/BALLS and about 1.0 on Hasegawa A.

## Grid (fixed charges Chunc 2.6226 g, BALLS 4.8292 g, Hasegawa A 14.9765 g; and the auto charge)

Forms: `sensible` = sensible energy + legacy momentum (shipped); `total_cons` = total energy + conservative momentum. Solid model: `goodman` (the
current default) and `conduction`. `auto` = the IG3 auto-size + upsize fallback (final charge shown when it differs). Hasegawa A rows are the first
0.3 s only (startup metrics). Full table: `worklog/RT1_results/summary.txt`.

### Chunc (exp: peak 8.66 MPa, plateau about 8.43, ratio 1.027)

| solid | charge | form | driver | first / t90 [ms] | never lit | peak [MPa] @ ms | plateau | ratio |
|---|---|---|---|---|---|---|---|---|
| goodman | 2.62 fixed | sensible | static | 4.44 / 16.31 | 0 | 8.698 @ 24.5 | 8.066 | 1.078 |
| goodman | 2.62 fixed | sensible | recovery | 3.47 / 12.14 | 1 | 10.138 @ 26.4 | 7.990 | 1.269 |
| goodman | 2.62 fixed | total_cons | static | 5.23 / 5.83 | 7 | 9.080 @ 24.3 | 7.160 | 1.268 |
| goodman | 2.62 fixed | total_cons | recovery | 4.24 / 4.94 | 10 | 8.572 @ 22.3 | 6.766 | 1.267 |
| goodman | auto | sensible | recovery | 5.22 / 7.43 (2.10 g) | 0 | 10.351 @ 27.7 | 8.138 | 1.272 |
| conduction | 2.62 fixed | sensible | static | 5.31 / 7.01 | 10 | 8.211 @ 26.3 | 6.766 | 1.214 |
| conduction | 2.62 fixed | sensible | recovery | 4.34 / 6.04 | 12 | 7.951 @ 24.4 | 6.511 | 1.221 |
| conduction | 2.62 fixed | total_cons | static | 5.56 / 7.84 | 7 | 8.807 @ 25.3 | 7.158 | 1.230 |
| conduction | 2.62 fixed | total_cons | recovery | 4.91 / 6.81 | 9 | 8.529 @ 25.4 | 6.894 | 1.237 |
| conduction | auto | sensible | static | 8.34 / 11.48 (2.10 g) | 1 | 9.745 @ 35.4 | 7.991 | 1.219 |
| conduction | auto | sensible | recovery | 6.16 / 9.10 (2.10 g) | 11 | 7.396 @ 31.2 | 6.632 | 1.115 |
| conduction | auto | total_cons | static | 5.56 / 7.84 (2.62 g) | 7 | 8.807 @ 25.3 | 7.158 | 1.230 |
| conduction | auto | total_cons | recovery | 8.48 / 11.75 (2.10 g) | 1 | 9.805 @ 33.3 | 7.988 | 1.227 |

(Other auto rows: goodman auto sensible/static is the B0 row at 2.62 g; goodman auto total_cons rows equal the fixed-charge rows.)

### BALLSstick (exp pk/trough about 1.50, pk/mid about 1.35)

| solid | charge | form | driver | first / t90 [ms] | never lit | peak [MPa] | trough / mid | pk/trough |
|---|---|---|---|---|---|---|---|---|
| goodman | 4.83 fixed | sensible | static | 5.02 / 5.68 | 9 | 7.193 | 3.522 / 3.888 | 2.042 |
| goodman | 4.83 fixed | sensible | recovery | 4.18 / 4.79 | 11 | 7.069 | 3.360 / 3.695 | 2.104 |
| goodman | 4.83 fixed | total_cons | static | 5.03 / 5.84 | 7 | 7.803 | 3.689 / 4.090 | 2.115 |
| goodman | 4.83 fixed | total_cons | recovery | 4.16 / 4.93 | 8 | 7.922 | 3.604 / 3.990 | 2.198 |
| goodman | auto | sensible | recovery | 8.25 / 8.96 (3.86 g) | 0 | 6.230 | 4.137 / 4.696 | 1.506 |
| conduction | 4.83 fixed | sensible | static | 6.02 / 8.08 | 9 | 6.812 | 3.523 / 3.884 | 1.933 |
| conduction | 4.83 fixed | sensible | recovery | 5.34 / 6.78 | 10 | 6.855 | 3.441 / 3.788 | 1.992 |
| conduction | 4.83 fixed | total_cons | static | 5.81 / 8.61 | 7 | 7.340 | 3.690 / 4.086 | 1.989 |
| conduction | 4.83 fixed | total_cons | recovery | 5.15 / 7.22 | 8 | 7.442 | 3.604 / 3.986 | 2.065 |
| conduction | auto | sensible | static | 9.59 / 12.93 (3.86 g) | 8 | 5.815 | 3.607 / 3.971 | 1.612 |
| conduction | auto | sensible | recovery | 8.00 / 10.02 (3.86 g) | 9 | 5.931 | 3.523 / 3.875 | 1.684 |
| conduction | auto | total_cons | static | 12.33 / 18.76 (3.86 g) | 1 | 6.861 | 4.242 / 4.723 | 1.617 |
| conduction | auto | total_cons | recovery | 8.44 / 11.69 (3.86 g) | 7 | 6.213 | 3.690 / 4.075 | 1.684 |

### Hasegawa A, first 0.3 s (exp late peak 6.44 MPa; see "Not done")

| solid | form | driver | first / t90 [ms] | never lit | early peak [MPa] @ ms |
|---|---|---|---|---|---|
| goodman | sensible | static | 4.20 / 5.75 | 8 | 4.654 @ 31.1 |
| goodman | sensible | recovery | 3.87 / 5.49 | 0 | 5.441 @ 31.9 |
| goodman | total_cons | static | 8.82 / 10.04 | 12 | 4.003 @ 60.7 |
| goodman | total_cons | recovery | 3.91 / 5.15 | 8 | 4.787 @ 28.5 |
| conduction | sensible | static | 5.88 / 10.05 | 8 | 4.465 @ 61.0 |
| conduction | sensible | recovery | 5.42 / 9.32 | 8 | 4.478 @ 55.9 |
| conduction | total_cons | static | 5.93 / 10.36 | 7 | 4.574 @ 61.3 |
| conduction | total_cons | recovery | 5.31 / 9.55 | 8 | 4.472 @ 54.4 |

(The Hasegawa A auto charge is 14.98 g in all four form/solid combinations, identical to the fixed rows.)

## Reading

1. **The recovery driver advances first light by 0.3-1.0 ms in every motor/form/solid combination** (shipped form, goodman: Chunc 4.44 -> 3.47,
   BALLS 5.02 -> 4.18, Hasegawa A 4.20 -> 3.87) and shortens t90. For the total-energy form it removes the delay the lower static T had opened:
   first light for total_cons goes 5.23 / 5.03 / 8.82 ms (static) -> 4.24 / 4.16 / 3.91 ms (recovery), against the shipped sensible/static
   4.44 / 5.02 / 4.20 ms. On Hasegawa A the total_cons form at the sensible charge lights at 3.9 ms instead of 8.8 ms (early peak 4.0 -> 4.8 MPa).
2. **It does NOT fix the head-end dead zones.** Never-lit counts move by 0-3 cells either way (Chunc sensible/goodman 0 -> 1, total_cons 7 -> 10; BALLS
   9 -> 11 and 7 -> 8). The one large drop is Hasegawa A sensible/goodman 8 -> 0 (total_cons 12 -> 8). With the conduction solid the dead zones are 7-12
   cells in every column. Consistent with the EN1 finding that those cells sit just under the ignition threshold: a hotter driver is not what is missing.
3. **The startup spike is not tamed.** Chunc sensible/goodman recovery gives a 10.14 MPa peak (ratio 1.27, exp 1.027): earlier first light starts more of
   the grain together. With the conduction solid the recovery rows are lower (7.95 sensible, 8.53 total_cons). BALLS pk/trough stays 1.9-2.2 at fixed charges
   (exp about 1.5); only the auto rows with the SMALLER 3.86 g charge reach 1.5-1.7.
4. **With the recovery driver the auto-size fallback finds a smaller charge** (shipped form, goodman: Chunc 2.62 -> 2.10 g, BALLS 4.83 -> 3.86 g, about -20 %).
   With the conduction solid the auto charges are already 2.10 g / 3.86 g (SM2).
5. The conduction solid delays first light by 1-4 ms and lowers the plateau (more unlit cells), as in SM2. The driver temperature stacks on top of it
   (about -0.5 to -1 ms) without changing that picture.

## Not done / decisions

* **Hasegawa A late peak vs 6.44 MPa** needs the full 10 s run (2.3M steps). That run exceeds the 4 GB cap (MemoryError at 4.0 GB at the trimmed-history stage, even with
  an explicit `history_capacity`), so it is blocked on MEM1: HEAVY-COMPUTE REQUEST. The 0.3 s rows above are the startup part only.
* **Nozzle boundary treating the cell T as stagnation T (optional): NOT done.** A clean fix needs the stagnation (P0, T0) of the cell built from (P, T, u) in every
  boundary evaluation (6 call sites in the PISO step plus the implicit `dmdp`, which then depends on u), for about 0.1 % effect (T6: nozzle-inlet KE is 0.1 % of
  the enthalpy at Chunc). Not clean for its size; left as an open item.
* B0 and the full suite are not run (MEM1). The default path is unchanged by construction (`T_driver = T` when `static`), and the Chunc goodman/static/sensible/fixed
  row above reproduces the B0 Chunc row to all digits (8.698445703737 MPa, 533513 steps).
