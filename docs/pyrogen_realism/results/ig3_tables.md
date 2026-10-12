## (c) Ignition threshold, Chunc, DeMar off, Peretz h (min charge with >=99 % of grain cells lit within 0.2 s; bisected to 2 %)

| topology | Ap reading | min charge [g] |
|---|---|---|
| forward_plenum | aft | 1.195 |
| forward_plenum | head | 2.205 |
| forward_plenum | local | 2.205 |
| head_basket | aft | 0.889 |
| head_basket | head | 1.433 |
| head_basket | local | 1.433 |

## (c) Minimum forward-plenum BPNV charge that lights each canonical motor (>=99 % lit within 0.25 s, bisected to 3 %)

| motor | DeMar | Ap reading | min charge [g] |
|---|---|---|---|
| balls | off | aft | 1.447 |
| balls | off | head | 7.508 |
| balls | off | local | 7.190 |
| balls | step | aft | 0.055 |
| balls | step | head | 0.055 |
| balls | step | local | 0.055 |
| chunc | off | aft | 1.198 |
| chunc | off | head | 3.245 |
| chunc | off | local | 2.245 |
| chunc | step | aft | 0.027 |
| chunc | step | head | 0.027 |
| chunc | step | local | 0.027 |
| haseg | off | aft | 14.381 |
| haseg | off | head | 14.381 |
| haseg | off | local | 14.381 |
| haseg | step | aft | 0.140 |
| haseg | step | head | 0.140 |
| haseg | step | local | 0.140 |

## (b) Realistic charges (Chunc head basket, 3 mm cubes L/D 1, DeMar off, Peretz h; t_max 0.2 s; plateau = median 0.1-0.2 s)

| pyrogen | g | Ap reading | first ign [ms] | t90 [ms] | peak [MPa] | t_peak [ms] | plateau [MPa] | peak/plateau |
|---|---|---|---|---|---|---|---|---|
| bpnv | 4.5 | aft | 1.72 | 2.54 | 20.89 | 9.63 | 8.43 | 2.48 |
| bpnv | 4.5 | head | 2.51 | 3.23 | 20.16 | 10.12 | 8.44 | 2.39 |
| bpnv | 4.5 | local | 2.51 | 3.19 | 20.17 | 10.12 | 8.44 | 2.39 |
| bpnv | 6.0 | aft | 1.23 | 2.09 | 26.87 | 8.43 | 8.39 | 3.20 |
| bpnv | 6.0 | head | 1.82 | 2.52 | 26.25 | 8.73 | 8.39 | 3.13 |
| bpnv | 6.0 | local | 1.82 | 2.48 | 26.26 | 8.73 | 8.39 | 3.13 |
| mtv | 4.5 | aft | 3.13 | 3.97 | 18.56 | 37.38 | 8.39 | 2.21 |
| mtv | 4.5 | head | 8.34 | 9.06 | 18.54 | 42.93 | 8.42 | 2.20 |
| mtv | 4.5 | local | 8.34 | 8.94 | 18.54 | 42.85 | 8.42 | 2.20 |
| mtv | 6.0 | aft | 2.37 | 3.11 | 21.94 | 34.33 | 8.34 | 2.63 |
| mtv | 6.0 | head | 4.14 | 4.96 | 21.92 | 36.57 | 8.35 | 2.62 |
| mtv | 6.0 | local | 4.14 | 4.83 | 21.92 | 36.50 | 8.35 | 2.62 |

## (a) Canonical motors, library BPNV Sutton-mass forward-plenum igniter, DeMar step, Peretz h (B0 rows)

| motor | Ap reading | igniter [g] | first ign [ms] | peak [MPa] | t_peak [ms] | plateau / trough / mid | ratios | steps |
|---|---|---|---|---|---|---|---|---|
| chunc | aft | 0.90 | 4.11 | 9.990 | 32.6 | plateau 8.144 | pk/pl 1.227 (exp 1.027) | 532059 |
| chunc | head | 0.90 | 5.81 | 9.950 | 36.8 | plateau 8.142 | pk/pl 1.222 (exp 1.027) | 532388 |
| chunc | local | 0.90 | 5.81 | 9.951 | 36.7 | plateau 8.142 | pk/pl 1.222 (exp 1.027) | 532381 |
| balls | aft | 1.39 | 4.68 | 7.105 | 33.1 | trough 4.343 / mid 4.861 | pk/tr 1.64, pk/mid 1.46 (exp ~1.50 / ~1.35) | 843437 |
| balls | head | 1.39 | 7.23 | 7.018 | 39.4 | trough 4.343 / mid 4.855 | pk/tr 1.62, pk/mid 1.45 (exp ~1.50 / ~1.35) | 843792 |
| balls | local | 1.39 | 7.22 | 7.020 | 39.2 | trough 4.343 / mid 4.855 | pk/tr 1.62, pk/mid 1.45 (exp ~1.50 / ~1.35) | 843726 |
| haseg | aft | 4.18 | 7.69 | 6.107 | 2411.0 | - | peak 0.95x exp (6.44) | 2336269 |
| haseg | head | 4.18 | 7.69 | 6.107 | 2411.0 | - | peak 0.95x exp (6.44) | 2336269 |
| haseg | local | 4.18 | 7.69 | 6.107 | 2411.0 | - | peak 0.95x exp (6.44) | 2336271 |
