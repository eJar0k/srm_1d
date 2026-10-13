# SM1 report -- solid-heating model comparison (Goodman vs Peretz ODE vs exact conduction)

Branch `sm1/solid-model` (based on `3f34259`). All sim runs: Peretz ignition h (local), OB1 pellets/basket, canonical knobs
(roughness 32 um, kappa 0.44, T_ignition 756 K, k_solid 0.271), `numba` warm cache. Data: `docs/pyrogen_realism/results/sm1_*`,
probes `docs/pyrogen_realism/probes/sm1_*.py`, kernels `srm_1d/solid_thermal.py`, verification helpers `srm_1d/tools/solid_verify.py`.
**The default stays `solid_heating_model='goodman'` (bit-identical path; B0 PASS). The user decides the default.**

## 0. Headline

1. **The head-end dead zones are NOT a Goodman-closure artifact.** With exact conduction the same cells stay dark: BALLSstick
   cells 3-11 (4.8 g and 4.83 g) and Hasegawa A cells 3-10 (14.98 g), identical to Goodman. Their peak surface temperatures
   differ by only 0-9 K between Goodman and exact conduction (cells sit 1-30 K short of 756 K, section 3.1). Hasegawa A late peak
   5.36 MPa = 0.83x the 6.44 MPa experiment under both. The dead zone is the **gas side**: the 1-D head end stagnates once the
   downstream grain lights, h falls 5-15x, and the surface never reaches T_ignition. (Peretz's ODE lights them, section 3.1, but only
   because it never gives heat back: a non-conservative artifact, verified in section 2.)
2. **Goodman's closure is wrong for varying h, and that matters for the ignition CASCADE, not the dead zone.** Against exact
   solutions Goodman errs by -282 K (step down) and +269 K (step up) right after an h step, and time-to-756 K by 44 % after a x10 h
   rise (section 2.3). In the sim, whole-grain lighting time (first light -> 90 % lit) on the canonical Chunc row goes 0.5 ms
   (Goodman) -> 15.8 ms (exact conduction) -> 35 ms (Peretz); on BALLSstick 0.5 -> 30.5 ms (Peretz 70); on Hasegawa A 1.4 -> 44.3 ms (Peretz 84).
   The Goodman "whole grain lights in 1-2 ms" cascade IG1 found is therefore substantially a closure artifact. The startup peak barely moves
   (Chunc 9.98 -> 9.90 MPa, -0.8 %, t_peak 35 -> 57 ms), plateau and Hasegawa A peak do not move.
3. Exact conduction costs 7-20 % wall time; Peretz is free.

## 1. The three models

All three see the **same combined surface flux** the Goodman kernel applies today: the pyrogen flux and adjacent-cell radiation are
folded into one effective `(h_total, T_driver)` before the solid step (`_goodman_ignition_sources_and_mass`: `h_total = h_conv + h_ig + h_rad`,
`T_driver = sum(h_i T_i)/h_total`), i.e. flux `q = h_total (T_driver - T_s)`. The conduction model uses exactly this as a Robin boundary
condition. The convective wall-loss sink to the gas is unchanged (`h_conv (T_gas - T_s)` at the pre-step T_s).

**goodman (default).** Goodman's cubic-profile heat-balance integral with constant-h assumption. State `delta`; `T_s = (3 k T_i + h delta T_g)/(3 k + h delta)`
is **algebraic in the current h and delta**: a change in h changes T_s instantly.

**peretz.** Peretz's integral-method ODE with T_s as the state. As read from the page image (thesis `peretz1974.pdf` p. 46, Eq. III-50; the same
equation is Eq. 19 of AIAA 72-1119 `peretz1972.pdf` p. 7 -- the task spec said 15-17, the paper's numbering is 17 heat equation, 18 surface BC, 19 the ODE):

```
dT_ps/dt = 4 alpha_pr h_c^2 (T - T_ps)^3 / [ 3 lambda_pr^2 (T_ps - T_pi) (2T - T_ps - T_pi) ]        (III-50)
T_ps(0) = T_pi + eps, eps = 0.1 K                                                                     (III-51)
```
derived from the same cubic profile, the surface BC `h_c (T - T_ps) = -lambda dT/dy|0` and the heat-balance integral
(III-47..III-49), with the explicit assumption (his words): "the rates of change of the heat-transfer coefficient, h_c, and gas
temperature, T, with respect to time, are neglected as compared with the rate of change of T_ps, which is completely justified for
the entire induction interval". He checked it only against a **linearly varying** h_c (max 5 % difference from an explicit finite-difference solution, heat-up 298 -> 700 K), not
steps. For constant h and T it is algebraically the same solution as Goodman (checked: same delta = 3 lambda (T_ps - T_pi)/(h_c (T - T_ps))).
Implementation: RK4 with sub-steps limiting the relative change of (T_s - T_i) and (T_g - T_s) to 10 %; the ODE is derived for heating, so when
T_gas <= T_s (or h = 0) T_s is held. **Consequence of the assumption:** T_s is carried but the implied penetration depth `delta = 3 k (T_s-T_i)/(h (T-T_s))` is recomputed from the
new h, so heat stored in the solid is neither conserved nor allowed to diffuse inward: when h falls, T_s is **held** (it cannot cool by conduction into the solid), and
when h rises it responds sluggishly (section 2.3).

**conduction.** Exact 1-D transient conduction in the solid per unignited grain cell: finite-volume on a geometrically graded grid (y_0 = 0 = surface, half control
volume there), backward Euler with sub-steps <= 0.0125 x heating age, Thomas algorithm, active depth 7 sqrt(alpha (age + dt)) (erfc 3.5 = 7e-7) with the far
boundary at 8 sqrt(alpha t_max). Surface flux implicit: `(V0/dt + alpha/d0 + alpha h/k) th0 - alpha/d0 th1 = V0/dt th0^n + (alpha h/k)(T_drv - T_i)`.
Properties: only `k_solid` and `alpha = k_solid/(rho_p Cps)` (the sim's own). **No new parameters.** Resolution (dy0 = 5e-8 m, growth ratio 1.1, sub-step fraction 0.0125 -> 95 nodes for
t_max = 2 s) was fixed by the convergence study of section 2.1, not tuned to sim results.

Calibration caveat: T_ignition 756 K and k_solid 0.271 were fitted (v0.7.5 LHS) with Goodman + the Gnielinski ignition h (and before the T11 PISO fix). Any other
solid model needs them re-fitted; the whole-grain lighting-time change (section 3.2) is large.


## 2. Verification against exact solutions (all three models)

Setup: k = 0.271 W/(m K), rho = 1800, Cps = 1500 (the sim's fixed Cps) -> alpha = 1.004e-7 m2/s, T_i = 293 K, T_g = 2800 K. Kernels driven with a prescribed h(t) at fixed dt
(`srm_1d/tools/solid_verify.py`); full output `results/sm1_verify.md` (reproduce: `python docs/pyrogen_realism/probes/sm1_verify.py --convergence`).
Reference solutions: exact semi-infinite Robin (`1 - erfcx(beta)`, beta = h sqrt(alpha t)/k), exact constant-flux `2 q sqrt(alpha t/pi)/k`, and for varying h an
**independent Volterra / Duhamel flux-superposition solve** (does not use the conduction grid; Richardson-extrapolated; it reproduces the exact Robin solution to 0.0002 / 0.005 / 0.07 K at
h = 1000 / 5000 / 20000, and agrees with the converged conduction solve to ~0.3 K at the sampled times, 2.4 K max at the step instant). Both reference solves are for the semi-infinite solid; the conduction domain
(3.9 mm, insulated far end, active depth cut at 7 sqrt(alpha t)) is far deeper than any ignition window.

**Stated tolerance for the conduction model before it was used in Part 3:** T_s error <= 1 K over the ignition window (T_s <= 900 K) and time-to-756 K error <= 0.5 % against the exact Robin solution for
h <= 20000 W/(m2 K), at dt = 2 us (the sim's mean step is 2.3 us) and at 20 us. Met (table 2.1: 0.28-0.79 K, 0.12-0.34 %). The step cases are a harder test (the step falls inside one coarse time step):
max error 0.3 K at dt = 2 us, 2.3-2.9 K at dt = 20 us.

### 2.1 Constant h, T_g = 2800 K (exact Robin)

| h [W/m2K] | t756 exact [ms] | dt [us] | goodman: max err K / t756 % | peretz | conduction |
|---|---|---|---|---|---|
| 1000 | 27.036 | 2 | 5.13 / -2.59 % | 5.12 / -2.58 % | 0.28 / +0.12 % |
| 1000 | 27.036 | 20 | 5.13 / -2.59 % | 5.12 / -2.58 % | 0.31 / +0.14 % |
| 5000 | 1.081 | 2 | 5.28 / -2.66 % | 5.12 / -2.58 % | 0.37 / +0.18 % |
| 5000 | 1.081 | 20 | 5.28 / -2.66 % | 5.12 / -2.58 % | 0.70 / +0.30 % |
| 20000 | 0.068 | 2 | 18.71 / -3.80 % | 5.12 / -2.59 % | 0.75 / +0.33 % |
| 20000 | 0.068 | 20 | 9.40 / -3.90 % | 5.12 / -2.63 % | 0.79 / +0.34 % |

(Error = max over T_exact <= 900 K. Goodman and Peretz are the same integral-method solution at constant h: ~2.6 % early in time-to-756 K, consistent with the ~2 % in T_s - T_i that Goodman/Peretz quote.
Goodman's extra error at h = 20000, dt = 2 us is its `delta_0 = 1e-6 m` start-up seed, which Peretz's eps = 0.1 K does not have.)

### 2.2 Constant flux q = 1e7 W/m2 (h = 1, T_drv = T_i + q/h; exact `2 q sqrt(alpha t/pi)/k`)

| t [ms] | exact [K] | goodman err | peretz err | conduction err |
|---|---|---|---|---|
| 0.1 | 424.9 | +3.63 | +3.08 | -0.19 |
| 0.3 | 521.5 | +5.65 | +5.33 | -0.28 |
| 0.6 | 616.1 | +7.76 | +7.53 | -0.32 |
| 1.0 | 710.1 | +9.89 | +9.72 | -0.37 |
| 1.5 | 803.9 | +12.04 | +11.89 | -0.41 |

Time to 756 K: exact 1.2319 ms; goodman -4.57 %, peretz -4.50 %, conduction +0.17 %.

### 2.3 Step changes in h (T_g const., dt = 2 us; reference = Volterra; the converged conduction solve agrees within 0.1-0.3 K at the sampled times)

**Step DOWN, h 5400 -> 5400/14 at 0.8 ms** (cell-9-like):

| t [ms] | reference [K] | goodman err | peretz err | conduction err |
|---|---|---|---|---|
| 0.80 (at the step) | 727.8 | +5.3 | +5.1 | -0.4 |
| 0.90 | 614.7 | **-282.0** | +118.3 | +0.1 |
| 1.30 | 532.5 | -192.3 | +201.0 | -0.1 |
| 1.80 | 498.7 | -150.6 | +235.3 | -0.1 |
| 4.80 | 462.5 | -81.2 | +275.1 | -0.1 |
| 9.00 | 468.0 | -55.8 | +274.4 | -0.1 |

Max |error| after the step: goodman 378 K, peretz 277 K, conduction 1.9 K. **Goodman drops the surface temperature at once to a value far below the true one (it discards the stored heat); Peretz holds T_s at ~740 K for ever (no inward conduction): about +275 K too hot.** Only conduction reproduces the slow relaxation.

**Step UP, h 500 -> 5000 (x10) at 1.0 ms:**

| t [ms] | reference [K] | goodman err | peretz err | conduction err |
|---|---|---|---|---|
| 1.00 (at the step) | 344.4 | +1.1 | +1.1 | 0.0 |
| 1.10 | 485.3 | **+269.2** | -24.8 | -0.5 |
| 1.50 | 646.4 | +174.5 | -14.9 | -0.4 |
| 2.00 | 757.7 | +130.3 | -9.8 | -0.4 |
| 5.00 | 1079.7 | +56.6 | -4.0 | -0.4 |

Time to 756 K: reference 1.991 ms; **goodman 1.108 ms (-44 %)**, peretz 2.044 ms (+2.7 %), conduction 1.993 ms. Goodman's T_s jumps up instantly when h rises (the ignition-cascade mechanism).

### 2.4 Convergence study of the conduction grid (what fixed the node count)

Varying dy0 / growth ratio / sub-step fraction against exact Robin (h = 1000, 5000, 20000) and against the converged solve (dy0 = 1.25e-8, ratio 1.03, 0.003 x age; 308 nodes) on the step cases
(T_s <= 900 K window). Full tables: `results/sm1_verify.md` section 4. Summary, worst of the three h / two step cases, dt = 2 us / 20 us:

| dy0 [m] | ratio | sub-step | nodes | worst Robin err [K] | worst t756 err % | worst step err [K] |
|---|---|---|---|---|---|---|
| 1e-7 | 1.15 | 0.0125 | 63 | 1.19 / 1.22 | 0.53 / 0.54 | 0.61 / 3.0 |
| **5e-8 (chosen)** | **1.1** | **0.0125** | **95** | **0.75 / 0.79** | **0.33 / 0.34** | **0.28 / 2.9** |
| 2.5e-8 | 1.05 | 0.0125 | 183 | 0.51 / 0.55 | 0.22 / 0.22 | 0.14 / 2.8 |
| 2.5e-8 | 1.05 | 0.00625 | 183 | 0.32 / 0.32 | 0.13 / 0.13 | 0.09 / 0.87 |
| 1.25e-8 | 1.03 | 0.003125 | 308 | 0.15 / 0.15 | 0.06 / 0.06 | 0 (reference) |

The chosen setting is the cheapest meeting the stated tolerance; halving the spacing again improves it by <= 0.25 K at 2x the nodes. The residual error is first-order backward Euler in time in the first steps and the surface half-cell.

## 3. Effect in the sim

Conditions: main's current defaults (Peretz ignition h local, OB1 I-A pellets and exact-cylinder basket), motors from the committed `.ric` files in this worktree, canonical knobs. "never" = not lit by the end of the run.
Raw data: `results/sm1_dz/*.json`, `sm1_realistic.jsonl`, `sm1_canon.jsonl`, `sm1_threshold.jsonl`, `sm1_tsmax_*.txt`.

### 3.1 Head-end dead zones (explicit BPNV charge, forward plenum, DeMar off), per model

(`probes/sm1_deadzone_model.py`; BALLSstick 2 s, Hasegawa A 10 s. The Goodman row reproduces the reviewer's numbers: cells 3-11 never lit, 12-53 at 5.1-5.45 ms. Hasegawa A is the 14.98 g IG3 auto-size.)

#### BALLSstick 4.8 g: explicit BPNV, forward plenum, heat_delivery_mode=none

| model | cells never lit (within the run) | x range of dead zone [mm] | P_max [MPa] @ t [s] | steps |
|---|---|---|---|---|
| goodman | 3-11 (9 cells) | 36-117 | 7.15 @ 0.023 | 850277 |
| peretz | 3-3 (1 cells) | 36-36 | 7.54 @ 0.026 | 847250 |
| conduction | 3-11 (9 cells) | 36-117 | 6.78 @ 0.025 | 851577 |

Ignition time per grain cell [ms] ("never" = not lit by t_end):

| cell | x [mm] | goodman | peretz | conduction |
|---|---|---|---|---|
| 3 | 36 | never | never | never |
| 4 | 46 | never | 1608.13 | never |
| 5 | 56 | never | 231.43 | never |
| 6 | 66 | never | 33.23 | never |
| 7 | 76 | never | 22.06 | never |
| 8 | 87 | never | 8.12 | never |
| 9 | 97 | never | 6.65 | never |
| 10 | 107 | never | 6.20 | never |
| 11 | 117 | never | 6.04 | never |
| 12 | 127 | 5.11 | 5.96 | 6.10 |
| 13 | 138 | 5.11 | 5.93 | 6.08 |
| 14 | 148 | 5.13 | 5.94 | 6.10 |
| 15 | 158 | 5.15 | 5.98 | 6.12 |
| 16 | 168 | 5.15 | 6.02 | 6.14 |
| 17 | 178 | 5.16 | 6.07 | 6.16 |
| 18 | 188 | 5.17 | 6.11 | 6.18 |
| 19 | 199 | 5.18 | 6.16 | 6.20 |
| 20 | 209 | 5.19 | 6.21 | 6.22 |
| 21 | 219 | 5.19 | 6.26 | 6.24 |
| 22 | 229 | 5.20 | 6.31 | 6.26 |
| 23 | 239 | 5.21 | 6.36 | 6.29 |
| 33 | 341 | 5.28 | 6.89 | 6.53 |
| 43 | 443 | 5.37 | 7.34 | 6.79 |
| 53 | 545 | 5.45 | 7.70 | 7.03 |
| 63 | 647 | 5.52 | 8.00 | 7.26 |
| 73 | 749 | 5.59 | 8.27 | 7.47 |
| 83 | 851 | 5.66 | 8.67 | 7.78 |
| 93 | 952 | 5.76 | 9.13 | 8.15 |

#### BALLSstick 4.83 g (IG3 auto size): explicit BPNV, forward plenum, heat_delivery_mode=none

| model | cells never lit (within the run) | x range of dead zone [mm] | P_max [MPa] @ t [s] | steps |
|---|---|---|---|---|
| goodman | 3-11 (9 cells) | 36-117 | 7.19 @ 0.022 | 850245 |
| peretz | 3-3 (1 cells) | 36-36 | 7.57 @ 0.026 | 847161 |
| conduction | 3-11 (9 cells) | 36-117 | 6.81 @ 0.024 | 851541 |

Ignition time per grain cell [ms] ("never" = not lit by t_end):

| cell | x [mm] | goodman | peretz | conduction |
|---|---|---|---|---|
| 3 | 36 | never | never | never |
| 4 | 46 | never | 1580.07 | never |
| 5 | 56 | never | 199.39 | never |
| 6 | 66 | never | 32.69 | never |
| 7 | 76 | never | 21.87 | never |
| 8 | 87 | never | 8.06 | never |
| 9 | 97 | never | 6.59 | never |
| 10 | 107 | never | 6.15 | never |
| 11 | 117 | never | 5.98 | never |
| 12 | 127 | 5.02 | 5.91 | 6.03 |
| 13 | 138 | 5.02 | 5.88 | 6.01 |
| 14 | 148 | 5.04 | 5.89 | 6.04 |
| 15 | 158 | 5.06 | 5.93 | 6.06 |
| 16 | 168 | 5.07 | 5.97 | 6.08 |
| 17 | 178 | 5.08 | 6.01 | 6.10 |
| 18 | 188 | 5.09 | 6.06 | 6.12 |
| 19 | 199 | 5.09 | 6.10 | 6.14 |
| 20 | 209 | 5.10 | 6.15 | 6.16 |
| 21 | 219 | 5.11 | 6.21 | 6.18 |
| 22 | 229 | 5.12 | 6.26 | 6.21 |
| 23 | 239 | 5.12 | 6.31 | 6.23 |
| 33 | 341 | 5.20 | 6.83 | 6.47 |
| 43 | 443 | 5.29 | 7.28 | 6.72 |
| 53 | 545 | 5.37 | 7.63 | 6.97 |
| 63 | 647 | 5.44 | 7.93 | 7.20 |
| 73 | 749 | 5.51 | 8.20 | 7.40 |
| 83 | 851 | 5.58 | 8.59 | 7.71 |
| 93 | 952 | 5.68 | 9.05 | 8.08 |

#### Hasegawa A 14.98 g (IG3 auto size): explicit BPNV, forward plenum, heat_delivery_mode=none

| model | cells never lit (within the run) | x range of dead zone [mm] | P_max [MPa] @ t [s] | late peak / 6.44 | steps |
|---|---|---|---|---|---|
| goodman | 3-10 (8 cells) | 59-176 | 5.36 @ 2.682 | 5.36 MPa = 0.83x | 2347426 |
| peretz | 3-4 (2 cells) | 59-76 | 5.92 @ 2.475 | 5.92 MPa = 0.92x | 2338290 |
| conduction | 3-10 (8 cells) | 59-176 | 5.36 @ 2.688 | 5.36 MPa = 0.83x | 2347930 |

Ignition time per grain cell [ms] ("never" = not lit by t_end):

| cell | x [mm] | goodman | peretz | conduction |
|---|---|---|---|---|
| 3 | 59 | never | never | never |
| 4 | 76 | never | never | never |
| 5 | 92 | never | 14.53 | never |
| 6 | 109 | never | 7.64 | never |
| 7 | 126 | never | 6.62 | never |
| 8 | 143 | never | 6.01 | never |
| 9 | 160 | never | 5.77 | never |
| 10 | 176 | never | 5.66 | never |
| 11 | 193 | 4.24 | 5.61 | 5.89 |
| 12 | 210 | 4.20 | 5.60 | 5.87 |
| 13 | 227 | 4.20 | 5.63 | 5.91 |
| 14 | 244 | 4.22 | 5.68 | 5.95 |
| 15 | 260 | 4.25 | 5.74 | 6.01 |
| 16 | 277 | 4.28 | 5.82 | 6.06 |
| 17 | 294 | 4.29 | 5.89 | 6.11 |
| 18 | 311 | 4.31 | 5.97 | 6.16 |
| 19 | 328 | 4.33 | 6.05 | 6.22 |
| 20 | 344 | 4.35 | 6.13 | 6.28 |
| 21 | 361 | 4.36 | 6.21 | 6.33 |
| 22 | 378 | 4.38 | 6.29 | 6.39 |
| 23 | 395 | 4.39 | 6.38 | 6.44 |
| 33 | 563 | 4.54 | 7.16 | 7.06 |
| 43 | 731 | 4.68 | 7.76 | 7.56 |
| 53 | 899 | 4.83 | 8.47 | 8.27 |
| 63 | 1067 | 4.99 | 9.12 | 8.87 |
| 73 | 1235 | 5.21 | 9.58 | 9.31 |
| 83 | 1403 | 5.45 | 9.97 | 9.69 |
| 93 | 1571 | 5.76 | 10.37 | 10.05 |


**Why the dead zone survives exact conduction: peak surface temperature reached by each dead cell (first 50 ms, snapshots every 0.5 ms, T_ignition = 756 K).**

BALLSstick 4.83 g (`results/sm1_tsmax_balls.txt`):

| cell | goodman peak T_s [K] | peretz | conduction |
|---|---|---|---|
| 3 | 543 | 584 | 530 |
| 5 | 698 | 746 | 687 |
| 7 | 731 | 756 (lit) | 722 |
| 9 | 747 | 756 (lit) | 743 |
| 10 | 752 | 756 (lit) | 749 |
| 11 | 755 | 756 (lit) | 753 |

Hasegawa A 14.98 g (`results/sm1_tsmax_haseg.txt`), cells 3-10, goodman / conduction peaks: 595/587, 688/683, 708/709, 719/726, 730/739, 738/747, 744/752, 748/755 K.
At t = 50 ms the dead cells' surfaces are at 308-320 K (goodman, which has forgotten the heat) vs 395-464 K (conduction) and are flat or cooling.

**The call (dead zone: closure artifact or real stagnation?): real stagnation on the gas side of the 1-D model, not the Goodman closure.** Exact conduction and Goodman leave the same cells dark and reach the same
peaks to within 0-9 K (conduction slightly hotter in the deeper dead cells, slightly cooler in the first ones). The surface never reaches 756 K because once the downstream grain lights, u and h at those
cells collapse (reviewer's probe: u 0.3-6 m/s, h down 5-15x). Two caveats: (i) the dead-zone edge is **marginal**: cells 9-11 end 1-12 K short of T_ignition, so a few-percent change in h or in T_ignition
moves the boundary by several cells; (ii) the 1-D flow cannot represent real head-end mixing and recirculation (`HEAD_END_NOTES.md`), so "real stagnation" means real *within the 1-D model*.
Peretz's ODE lights the dead cells only because it holds T_s (section 2.3: it overstates T_s by ~275 K after an h drop); its lighting times (6 ms ... 1.6 s) should not be trusted.
Hasegawa A late peak: 5.36 MPa (0.83x the 6.44 MPa experiment) with Goodman and with conduction, 5.92 MPa (0.92x) with Peretz (more grain lit).


### 3.2 Whole-grain lighting time and peak

Tables A and B (`results/sm1_tables.md`). t90 - first = first light -> 90 % of the grain lit.

Canonical motors, library BPNV igniter (the B0 rows; Sutton mass, forward plenum, DeMar as shipped):

| motor | model | first light [ms] | t90 [ms] | t90 - first [ms] | last cell [ms] | peak [MPa] | t_peak [ms] | plateau / trough / mid | ratios |
|---|---|---|---|---|---|---|---|---|---|
| Chunc | goodman | 5.86 | 6.37 | **0.5** | 6.46 | 9.977 | 35 | 8.142 | pk/pl 1.225 |
| Chunc | peretz | 7.43 | 42.31 | 34.9 | 43.08 | 9.908 | 77 | 8.116 | 1.221 |
| Chunc | conduction | 6.99 | 22.82 | **15.8** | 23.33 | 9.899 | 57 | 8.130 | 1.218 |
| BALLSstick | goodman | 7.27 | 7.80 | **0.5** | 7.91 | 7.061 | 37 | 4.343 / 4.856 | pk/tr 1.63, pk/mid 1.45 |
| BALLSstick | peretz | 9.25 | 78.97 | 69.7 | 80.63 | 6.969 | 114 | 4.357 / 4.760 | 1.60 / 1.46 |
| BALLSstick | conduction | 8.64 | 39.17 | **30.5** | 40.09 | 6.945 | 76 | 4.350 / 4.812 | 1.60 / 1.44 |
| Hasegawa A | goodman | 7.65 | 9.07 | **1.4** | 9.33 | 6.107 | 2411 ms | peak/6.44 = 0.95 | |
| Hasegawa A | peretz | 8.89 | 92.86 | 84.0 | 94.91 | 6.112 | 2493 ms | 0.95 | |
| Hasegawa A | conduction | 8.61 | 52.88 | **44.3** | 53.68 | 6.110 | 2455 ms | 0.95 | |

First light moves +1.0-1.4 ms with conduction (the exact solid is slower than Goodman's cubic early on), and the lighting cascade stretches 30-60x (Peretz 70-170x). The startup peak is unchanged to within 0.8-1.6 %
(the 10 MPa spike is set by the ignited grain, not by how fast it lit: IG1/IG2), the plateau to 0.2 %, Hasegawa A's late peak to 0.05 %. **The Goodman 0.5-1.4 ms whole-grain cascade
is a closure artifact; exact conduction gives 16-44 ms.**

Chunc head basket, OB1 default (3.175 mm I-A pellet, exact-cylinder law), t_max 0.2 s, plateau = median 0.1-0.2 s (24 rows, DeMar on and off in `results/sm1_tables.md` table A; DeMar on shown):

| pyrogen | g | model | first [ms] | t90 [ms] | t90 - first | peak [MPa] | t_peak [ms] | plateau | peak/plateau |
|---|---|---|---|---|---|---|---|---|---|
| bpnv | 4.5 | goodman | 2.51 | 3.84 | 1.3 | 16.72 | 18.3 | 8.44 | 1.98 |
| bpnv | 4.5 | peretz | 3.47 | 6.88 | 3.4 | 15.66 | 20.9 | 8.47 | 1.85 |
| bpnv | 4.5 | conduction | 3.14 | 5.75 | 2.6 | 16.05 | 19.9 | 8.46 | 1.90 |
| bpnv | 6.0 | goodman | 1.63 | 2.81 | 1.2 | 20.89 | 15.4 | 8.39 | 2.49 |
| bpnv | 6.0 | peretz | 2.26 | 4.69 | 2.4 | 19.90 | 16.8 | 8.42 | 2.36 |
| bpnv | 6.0 | conduction | 2.04 | 4.04 | 2.0 | 20.27 | 16.3 | 8.41 | 2.41 |
| mtv | 4.5 | goodman | 2.49 | 12.06 | 9.6 | 14.33 | 40.2 | 8.46 | 1.69 |
| mtv | 4.5 | peretz | 3.26 | 18.79 | 15.5 | 14.11 | 45.4 | 8.51 | 1.66 |
| mtv | 4.5 | conduction | 3.01 | 16.38 | 13.4 | 14.19 | 43.3 | 8.49 | 1.67 |
| mtv | 6.0 | goodman | 1.93 | 5.15 | 3.2 | 16.45 | 30.8 | 8.36 | 1.97 |
| mtv | 6.0 | peretz | 2.54 | 10.85 | 8.3 | 16.16 | 35.8 | 8.40 | 1.92 |
| mtv | 6.0 | conduction | 2.33 | 8.83 | 6.5 | 16.26 | 34.1 | 8.39 | 1.94 |

With a realistic (4.5-6 g) head-basket charge the effect is smaller (the pyrogen itself lights the grain over a few ms): lighting spread 1.3 -> 2.6 ms, peak -4 % (BPNV 4.5 g, conduction), -3 % at 6 g,
-1 % MTV; plateau +0.02-0.04 MPa. The over-spike stays (peak/plateau 1.7-2.4): it is not a solid-model problem.

### 3.3 Ignition thresholds (Chunc, DeMar off, BPNV, 0.2 s; bisected to 3 %; `results/sm1_tables.md` table C)

| topology | criterion | goodman [g] | peretz [g] | conduction [g] |
|---|---|---|---|---|
| head basket | >= 99 % of grain lit | 3.00 | 1.86 | 2.58 |
| head basket | all cells lit | 3.20 | 2.81 | 3.34 |
| forward plenum | >= 99 % of grain lit | 2.27 | 1.38 | 5.40 |
| forward plenum | all cells lit | 3.58 | 7.16 | 8.70 |

Head basket: the exact solid needs ~4 % more charge than Goodman to light every cell (3.34 vs 3.20 g; 2.58 vs 3.00 g on the 99 % criterion, where it is 1 cell). Forward plenum, DeMar off:
2.4x more (8.7 vs 3.6 g all cells lit). The 99 % criterion is brittle there: the head-end cells 3-10 are marginal and `frac lit` is non-monotonic in charge (conduction 4.5 g: 0.97, 5.4 g: 0.99, 9.1 g: 1.00;
Peretz 0.99 up to 4.5 g, all lit only from 7.2 g), so the "all cells lit" row is the reliable one. Both columns are lower bounds on a realistic charge: the real 4.5-6 g head-basket charge is above
all of them (1.3-1.8x the conduction head-basket threshold), while the forward plenum, DeMar off, needs 8.7 g with conduction against the library Sutton size 0.9 g. Per-charge fractions: end of table C.

### 3.4 Wall-time cost (serial, warm cache, best of 2; `results/sm1_cost.txt`)

| case | goodman | peretz | conduction | conduction / goodman |
|---|---|---|---|---|
| BALLSstick dead zone 4.83 g, 2 s (850 k steps) | 31.4 s (27.1 k steps/s) | 31.2 s | 37.7 s (22.6 k steps/s) | **1.20x** |
| Chunc head basket 4.5 g, 0.2 s (115 k steps) | 2.74 s | 2.81 s | 2.94 s | 1.07x |
| Chunc canonical, 1 s (532 k steps) | 12.16 s | 12.68 s | 13.22 s | 1.09x |

(The 1.20x case keeps ~10 cells unlit for the whole 2 s, each carrying a 95-node profile; the others stop paying once the cell ignites. Memory: N x 95 doubles. Peretz is free.)

## 4. Gates

- **B0** (`baseline_gate.py`, defaults, from the worktree root): **PASS** (`results/sm1_B0_output.txt`). Chunc peak 9.977393 (ref 9.977396, 3.0e-7), t_peak identical to 7.6e-9, plateau 7e-10, n_steps 532348 exact;
  BALLSstick peak 2.9e-7, n_steps 843680 exact; Hasegawa A peak 5e-11, n_steps 2336319 vs 2336318 (4e-7). Differences are at the documented compile-state noise floor (1e-7..4e-7); the default path
  executes the pre-SM1 Goodman code unchanged (the new branches are not taken).
- **Full suite** `python -m pytest tests/ -q -p no:cacheprovider`: base `3f34259`: 775 passed; this branch: 795 passed, 0 failed (775 + the 20 new tests). 20 new tests in `tests/test_solid_heating_models.py` (exact Robin / flux / step against Volterra, default
  identity = explicit goodman bitwise, models run). The direct callers of `_goodman_ignition_sources_and_mass` in `tests/test_ignition_heat_transfer.py` and `tests/test_simulation_phase3.py` were extended with the 7 new trailing arguments (internal hard API break,
  logged in `DEVNOTES.md`).

## 5. Part 4

`HEAD_END_NOTES.md`: Peretz Eqs. III-28..III-30 (Cohen, Derr & Price's igniter-duct correlations, Nu_d = 0.0615 Re^0.8 (Ap/At)^0.4 at impingement; 0.7 Re^0.8 (Ap/At)^0.4 (0.022 + 0.066 x/x_im) upstream;
0.7 Re^0.8 (Ap/At)^0.4 [0.031 + 0.057 exp(-0.5 x_r/d)] downstream; x_im ~ one port diameter), how Peretz used them (calibrating a6 and Ap/At^0.4 of III-31 up to x/d = 10, and the effective igniter-gas T_ig, 25-35 % below adiabatic), the comparison with
Carlson & Seader's figure-based correlation, and the follow-up proposal: the dead-zone cells (x/d ~ 3-9) lie in the data range, and an igniter-flow-based Re_d there is the only data-based route to a gas-side fix.

## 6. Recommendation (the user decides the default)

1. **Make `conduction` the default for runs where ignition-front timing matters; keep Goodman only for bit-compatibility with the v0.7.5 calibration.** It is the only model that meets the exact solutions (0.3-0.8 K; Goodman errs by hundreds of K when h changes, which it does constantly in the starting
   transient), it has no new parameters and costs 7-20 % wall time. It removes the unphysical 0.5-1.4 ms whole-grain cascade (-> 16-44 ms) while leaving the startup peak, plateau and Hasegawa A peak within ~1 %.
2. **Do not adopt `peretz` as a default.** It is the same integral method as Goodman at constant h (identical errors), but its carried T_s cannot relax: +275 K error after an h drop, and it lights the head-end dead zones only through that. It is the one model that
   changes the dead-zone answer, and it does so wrongly.
3. **Switching requires re-fitting T_ignition 756 K and k_solid 0.271** (fitted with Goodman, Gnielinski h and the pre-T11 PISO): conduction first light is 1.0-1.4 ms later and the front 30-60x slower than Goodman's. That re-fit is the user's decision (the plateau and Hasegawa A peak are insensitive, so the
   fit targets would be the startup shape).
4. **The dead zones need a gas-side/closure fix, not a solid-model fix** (section 3.1). Candidates: the head-end h of Peretz III-28..III-30 with igniter-flow Re_d (HEAD_END_NOTES.md; needs Cohen/Derr/Price for the Re_d definition), or a realistic igniter that keeps the head-end flow alive. Note the margin: cells 9-11
   are 1-12 K short, so any such term would also shift where the boundary lies.

## 7. Not done / limits

- The Volterra reference and the exact Robin solution are for a semi-infinite solid; no check against a finite web or a curved (cylindrical) port surface (negligible at ignition depths: 0.1-1 mm vs 6-10 mm bore radius).
- The dead-zone runs use one igniter size per motor (4.8 / 4.83 g and 14.98 g); no charge sweep and no h-scaling sensitivity of the 1-12 K marginality was run.
- Threshold bisection was to 3 % and the 99 % criterion is non-monotonic (section 3.3); only the "all cells lit" row is monotone enough to quote.
- The T_surf of the new models after `heat_cell == False` (flame-front gate) is frozen like Goodman's (no diffusion while gated off); the conduction state for cells that are not heated is not advanced.
- Calibration (T_ignition, k_solid) was not re-fitted for any model.
- Part 4 is a note only; I did not re-read the Carlson & Seader PDF (coordinator's instruction) and the question whether III-28..30 are the analytical fit of its figures is left open with the reason (no Dp/Di parameter in III-28..30).
