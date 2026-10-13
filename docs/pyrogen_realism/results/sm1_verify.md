k = 0.271 W/(m K), alpha = 1.0037e-07 m^2/s, T_i = 293.0 K, conduction grid: dy0 = 5.0e-08 m, ratio 1.1, M = 95 nodes to 3.89 mm, BE sub-step <= 0.0125 x heating age, active depth 7 sqrt(alpha t)

### 1. Constant h, T_g = 2800 K: exact Robin (erfcx) solution

Error = max |T_s - exact| over the window T_exact <= 900 K [K]; dt90% = time-to-756 K error [% of exact t756].

| h [W/m2K] | t756 exact [ms] | dt [us] | goodman err K / t756 % | peretz err K / t756 % | conduction err K / t756 % |
|---|---|---|---|---|---|
| 1000 | 27.036 | 2 | 5.13 / -2.59% | 5.12 / -2.58% | 0.28 / +0.12% |
| 1000 | 27.036 | 20 | 5.13 / -2.59% | 5.12 / -2.58% | 0.31 / +0.14% |
| 5000 | 1.081 | 2 | 5.28 / -2.66% | 5.12 / -2.58% | 0.37 / +0.18% |
| 5000 | 1.081 | 20 | 5.28 / -2.66% | 5.12 / -2.58% | 0.70 / +0.30% |
| 20000 | 0.068 | 2 | 18.71 / -3.80% | 5.12 / -2.59% | 0.75 / +0.33% |
| 20000 | 0.068 | 20 | 9.40 / -3.90% | 5.12 / -2.63% | 0.79 / +0.34% |

### 2. Constant surface flux: exact T_s = T_i + 2 q sqrt(alpha t / pi) / k

Realised as h = 1 W/m2K, T_drv = T_i + q/h (flux droop = (T_s - T_i) / (q/h) < 5e-5). q = 1e7 W/m2.

| t [ms] | exact T_s [K] | goodman err [K] | peretz err [K] | conduction err [K] |
|---|---|---|---|---|
| 0.1 | 424.9 | +3.63 | +3.08 | -0.19 |
| 0.3 | 521.5 | +5.65 | +5.33 | -0.28 |
| 0.6 | 616.1 | +7.76 | +7.53 | -0.32 |
| 1.0 | 710.1 | +9.89 | +9.72 | -0.37 |
| 1.5 | 803.9 | +12.04 | +11.89 | -0.41 |

Time to 756 K: exact 1.2319 ms; error goodman -4.57%, peretz -4.50%, conduction +0.17%

### 3. Step changes in h (T_g = 2800 K const.), dt = 2 us

Reference = Duhamel/Volterra flux-superposition solve (independent of the conduction grid, Richardson-extrapolated from dt = 2 us and 1 us); the converged conduction solve (dy0 = 1.25e-8, ratio 1.03, sub-step 0.003 x age) is quoted as a cross-check.

**step DOWN  h 5400 -> 5400/14 at 0.8 ms (cell-9 like)**

| t [ms] | Volterra ref [K] | converged cond. [K] | goodman err [K] | peretz err [K] | conduction err [K] |
|---|---|---|---|---|---|
| 0.50 | 647.2 | 647.1 | +5.2 | +4.9 | -0.4 |
| 0.80 | 727.8 | 727.7 | +5.3 | +5.1 | -0.4 |
| 0.90 | 614.7 | 614.9 | -282.0 | +118.3 | +0.1 |
| 1.30 | 532.5 | 532.5 | -192.3 | +201.0 | -0.1 |
| 1.80 | 498.7 | 498.8 | -150.6 | +235.3 | -0.1 |
| 2.80 | 473.9 | 473.9 | -112.8 | +261.4 | -0.1 |
| 4.80 | 462.5 | 462.5 | -81.2 | +275.1 | -0.1 |
| 9.00 | 468.0 | 468.0 | -55.8 | +274.4 | -0.1 |

max |T_s error| after the step [K]: goodman 378.3, peretz 276.7, conduction 1.9; converged-conduction vs Volterra 2.11
T_s never reaches 756 K (Volterra max 728 K); goodman max 733 K, peretz max 744 K, conduction max 727 K

**step UP    h 500 -> 5000 (x10) at 1.0 ms**

| t [ms] | Volterra ref [K] | converged cond. [K] | goodman err [K] | peretz err [K] | conduction err [K] |
|---|---|---|---|---|---|
| 0.50 | 329.6 | 329.5 | +0.8 | +0.8 | -0.0 |
| 1.00 | 344.4 | 344.4 | +1.1 | +1.1 | -0.0 |
| 1.10 | 485.3 | 485.0 | +269.2 | -24.8 | -0.5 |
| 1.50 | 646.4 | 646.2 | +174.5 | -14.9 | -0.4 |
| 2.00 | 757.7 | 757.5 | +130.3 | -9.8 | -0.4 |
| 3.00 | 901.4 | 901.3 | +90.0 | -5.8 | -0.4 |
| 5.00 | 1079.7 | 1079.6 | +56.6 | -4.0 | -0.4 |
| 9.00 | 1289.7 | 1289.6 | +30.2 | -5.2 | -0.4 |

max |T_s error| after the step [K]: goodman 370.9, peretz 26.6, conduction 2.4; converged-conduction vs Volterra 2.35
time to 756 K: Volterra 1.991 ms; goodman 1.108 ms, peretz 2.044 ms, conduction 1.993 ms

### 4. Conduction-grid convergence study (T_s <= 900 K; dt = 2 us / 20 us)

max |T_s - exact Robin| [K] / time-to-756 error [%] for h = 1000, 5000, 20000; 'step' = max |conduction - converged conduction| [K] on the step-down and step-up cases (T_s <= 900 K window).

dt = 2 us

| dy0 [m] | ratio | sub-step frac | nodes | h=1000 | h=5000 | h=20000 | step-down | step-up |
|---|---|---|---|---|---|---|---|---|
| 1.00e-07 | 1.15 | 0.05 | 63 | 0.60 K / +0.25% | 0.71 K / +0.33% | 1.88 K / +0.95% | 0.56 K | 0.61 K |
| 1.00e-07 | 1.15 | 0.0125 | 63 | 0.60 K / +0.25% | 0.71 K / +0.32% | 1.19 K / +0.53% | 0.55 K | 0.61 K |
| 5.00e-08 | 1.1 | 0.05 | 95 | 0.28 K / +0.12% | 0.41 K / +0.18% | 1.47 K / +0.74% | 0.40 K | 0.26 K |
| 5.00e-08 | 1.1 | 0.025 | 95 | 0.28 K / +0.12% | 0.37 K / +0.18% | 1.12 K / +0.46% | 0.32 K | 0.26 K |
| 5.00e-08 | 1.1 | 0.0125 | 95 | 0.28 K / +0.12% | 0.37 K / +0.18% | 0.75 K / +0.33% | 0.28 K | 0.26 K | (chosen)
| 5.00e-08 | 1.1 | 0.00625 | 95 | 0.28 K / +0.12% | 0.36 K / +0.17% | 0.57 K / +0.24% | 0.25 K | 0.26 K |
| 2.50e-08 | 1.05 | 0.0125 | 183 | 0.08 K / +0.03% | 0.20 K / +0.09% | 0.51 K / +0.22% | 0.14 K | 0.05 K |
| 2.50e-08 | 1.05 | 0.00625 | 183 | 0.08 K / +0.03% | 0.16 K / +0.08% | 0.32 K / +0.13% | 0.09 K | 0.05 K |
| 1.25e-08 | 1.03 | 0.003125 | 308 | 0.03 K / +0.01% | 0.10 K / +0.05% | 0.15 K / +0.06% | 0.00 K | 0.00 K |

dt = 20 us

| dy0 [m] | ratio | sub-step frac | nodes | h=1000 | h=5000 | h=20000 | step-down | step-up |
|---|---|---|---|---|---|---|---|---|
| 1.00e-07 | 1.15 | 0.05 | 63 | 0.63 K / +0.27% | 1.43 K / +0.73% | 2.57 K / +1.09% | 4.75 K | 6.56 K |
| 1.00e-07 | 1.15 | 0.0125 | 63 | 0.62 K / +0.27% | 1.05 K / +0.44% | 1.22 K / +0.54% | 2.15 K | 3.00 K |
| 5.00e-08 | 1.1 | 0.05 | 95 | 0.31 K / +0.14% | 1.16 K / +0.58% | 2.14 K / +0.88% | 4.94 K | 6.43 K |
| 5.00e-08 | 1.1 | 0.025 | 95 | 0.31 K / +0.14% | 0.92 K / +0.44% | 1.24 K / +0.52% | 5.42 K | 6.39 K |
| 5.00e-08 | 1.1 | 0.0125 | 95 | 0.31 K / +0.14% | 0.70 K / +0.30% | 0.79 K / +0.34% | 2.34 K | 2.86 K | (chosen)
| 5.00e-08 | 1.1 | 0.00625 | 95 | 0.31 K / +0.14% | 0.52 K / +0.21% | 0.56 K / +0.25% | 0.69 K | 0.94 K |
| 2.50e-08 | 1.05 | 0.0125 | 183 | 0.12 K / +0.05% | 0.49 K / +0.20% | 0.55 K / +0.22% | 2.48 K | 2.80 K |
| 2.50e-08 | 1.05 | 0.00625 | 183 | 0.11 K / +0.05% | 0.30 K / +0.12% | 0.32 K / +0.13% | 0.82 K | 0.87 K |
| 1.25e-08 | 1.03 | 0.003125 | 308 | 0.06 K / +0.03% | 0.15 K / +0.06% | 0.15 K / +0.06% | 0.00 K | 0.00 K |

