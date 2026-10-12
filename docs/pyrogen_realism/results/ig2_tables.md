## Ignition threshold (min charge with >=99 % of grain cells lit within 0.2 s; bisected to 2 %)

| topology | ignition h | DeMar | min charge [g] |
|---|---|---|---|
| forward_plenum | gnielinski | off | 0.245 |
| forward_plenum | gnielinski | step | 0.028 |
| forward_plenum | keller | off | 1.738 |
| forward_plenum | keller | step | 0.027 |
| forward_plenum | peretz | off | 2.205 |
| forward_plenum | peretz | step | 0.027 |
| head_basket | gnielinski | off | 0.352 |
| head_basket | gnielinski | step | 0.352 |
| head_basket | keller | off | 1.511 |
| head_basket | keller | step | 1.511 |
| head_basket | peretz | off | 1.433 |
| head_basket | peretz | step | 1.433 |

## Realistic charges (head basket, 3 mm cubes L/D 1, Chunc; t_max 0.2 s; plateau = median 0.1-0.2 s)

| pyrogen | g | ignition h | DeMar | first ign [ms] | t90 [ms] | peak [MPa] | t_peak [ms] | plateau [MPa] | peak/plateau |
|---|---|---|---|---|---|---|---|---|---|
| bpnv | 4.5 | gnielinski | off | 0.68 | 1.79 | 21.60 | 9.07 | 8.42 | 2.57 |
| bpnv | 4.5 | gnielinski | step | 0.66 | 1.76 | 21.62 | 9.05 | 8.42 | 2.57 |
| bpnv | 4.5 | keller | off | 2.18 | 2.89 | 20.47 | 9.92 | 8.43 | 2.43 |
| bpnv | 4.5 | keller | step | 1.96 | 2.77 | 20.66 | 9.80 | 8.43 | 2.45 |
| bpnv | 4.5 | peretz | off | 2.51 | 3.19 | 20.17 | 10.12 | 8.44 | 2.39 |
| bpnv | 4.5 | peretz | step | 2.07 | 3.04 | 20.38 | 9.99 | 8.44 | 2.42 |
| bpnv | 6.0 | gnielinski | off | 0.42 | 1.54 | 27.40 | 8.10 | 8.38 | 3.27 |
| bpnv | 6.0 | gnielinski | step | 0.40 | 1.54 | 27.41 | 8.09 | 8.38 | 3.27 |
| bpnv | 6.0 | keller | off | 1.55 | 2.26 | 26.57 | 8.58 | 8.39 | 3.17 |
| bpnv | 6.0 | keller | step | 1.22 | 2.21 | 26.71 | 8.52 | 8.39 | 3.18 |
| bpnv | 6.0 | peretz | off | 1.82 | 2.48 | 26.26 | 8.73 | 8.39 | 3.13 |
| bpnv | 6.0 | peretz | step | 1.35 | 2.36 | 26.49 | 8.62 | 8.39 | 3.16 |
| mtv | 4.5 | gnielinski | off | 1.15 | 2.23 | 18.58 | 34.97 | 8.38 | 2.22 |
| mtv | 4.5 | gnielinski | step | 0.96 | 2.16 | 18.59 | 34.78 | 8.37 | 2.22 |
| mtv | 4.5 | keller | off | 11.26 | 12.01 | 18.54 | 45.70 | 8.44 | 2.20 |
| mtv | 4.5 | keller | step | 2.35 | 10.02 | 18.55 | 43.79 | 8.43 | 2.20 |
| mtv | 4.5 | peretz | off | 8.34 | 8.94 | 18.54 | 42.85 | 8.42 | 2.20 |
| mtv | 4.5 | peretz | step | 2.22 | 7.16 | 18.56 | 40.76 | 8.41 | 2.21 |
| mtv | 6.0 | gnielinski | off | 0.86 | 1.94 | 21.97 | 32.64 | 8.33 | 2.64 |
| mtv | 6.0 | gnielinski | step | 0.62 | 1.91 | 21.98 | 32.52 | 8.33 | 2.64 |
| mtv | 6.0 | keller | off | 4.71 | 5.52 | 21.91 | 37.38 | 8.36 | 2.62 |
| mtv | 6.0 | keller | step | 1.71 | 4.44 | 21.93 | 35.73 | 8.35 | 2.63 |
| mtv | 6.0 | peretz | off | 4.14 | 4.83 | 21.92 | 36.50 | 8.35 | 2.62 |
| mtv | 6.0 | peretz | step | 1.69 | 3.97 | 21.94 | 35.18 | 8.35 | 2.63 |

## Realistic charges at t_max 1.0 s (plateau = median 0.4-0.9 s; Chunc exp: peak 8.89, plateau 8.66, ratio 1.027)

| tag | peak [MPa] | t_peak [ms] | plateau | peak/plateau | first ign [ms] |
|---|---|---|---|---|---|
| cr_bpnv4.5_gn_step | 21.62 | 9.05 | 8.17 | 2.65 | 0.66 |
| cr_bpnv4.5_pe_off | 20.17 | 10.12 | 8.16 | 2.47 | 2.51 |
| cr_bpnv6_gn_step | 27.41 | 8.09 | 8.18 | 3.35 | 0.40 |
| cr_bpnv6_pe_off | 26.26 | 8.73 | 8.17 | 3.21 | 1.82 |
| cr_mtv4.5_gn_step | 18.59 | 34.78 | 8.17 | 2.27 | 0.96 |
| cr_mtv4.5_pe_off | 18.54 | 42.85 | 8.16 | 2.27 | 8.34 |
| cr_mtv6_gn_step | 21.98 | 32.52 | 8.18 | 2.69 | 0.62 |
| cr_mtv6_pe_off | 21.92 | 36.50 | 8.18 | 2.68 | 4.14 |

## Validation, canonical motors, library BPNV Sutton-mass forward-plenum igniter (B0 rows)

| motor | ignition h | DeMar | igniter [g] | first ign [ms] | peak [MPa] | t_peak [s] | plateau / trough / mid | ratios | steps | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|
| balls | gnielinski | off | 1.39 | 1.77 | 7.256 | 0.027 | trough 4.34 / mid 4.87 | pk/tr 1.67, pk/mid 1.49 | 843086 | 54 |
| balls | gnielinski | step | 1.39 | 1.75 | 7.258 | 0.027 | trough 4.34 / mid 4.87 | pk/tr 1.67, pk/mid 1.49 | 843077 | 54 |
| balls | peretz | off | 1.39 | - | 0.859 | 0.004 | trough 0.10 / mid 0.10 | pk/tr 8.53, pk/mid 8.48 | 491605 | 51 |
| balls | peretz | step | 1.39 | 7.22 | 7.020 | 0.039 | trough 4.34 / mid 4.86 | pk/tr 1.62, pk/mid 1.45 | 843726 | 52 |
| chunc | gnielinski | off | 0.90 | 1.72 | 10.083 | 0.028 | plateau 8.15 | pk/pl 1.238 | 531704 | 24 |
| chunc | gnielinski | step | 0.90 | 1.72 | 10.083 | 0.028 | plateau 8.15 | pk/pl 1.238 | 531701 | 24 |
| chunc | peretz | off | 0.90 | - | 0.734 | 0.002 | plateau 0.10 | pk/pl 7.340 | 312181 | 22 |
| chunc | peretz | step | 0.90 | 5.81 | 9.951 | 0.037 | plateau 8.14 | pk/pl 1.222 | 532381 | 24 |
| haseg | gnielinski | off | 4.18 | 4.63 | 6.108 | 2.409 | - | - | 2336353 | 70 |
| haseg | gnielinski | step | 4.18 | 2.77 | 6.107 | 2.406 | - | - | 2336674 | 71 |
| haseg | peretz | off | 4.18 | - | 0.444 | 0.002 | - | - | 1405595 | 70 |
| haseg | peretz | step | 4.18 | 7.69 | 6.107 | 2.411 | - | - | 2336268 | 71 |
