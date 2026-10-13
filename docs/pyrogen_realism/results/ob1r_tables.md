## T1 rate laws (mm/s)

| P [MPa] | sim law 2e-5 P^0.5 | Nakano 71.1 P^0.589 (valid 0.005-0.2 MPa) | Zhu mean r0 (valid 2.8-50 MPa) | MIL-P-46994B (closed-bomb average) |
|---|---|---|---|---|
| 0.1 | 6.3 | 18.3 | (below range) |  |
| 0.2 | 8.9 | 27.6 | (below range) |  |
| 1 | 20.0 | (71, out of range) | (below range) |  |
| 2.8 | 33.5 | (130, out of range) | 17.0 | 31.8-44.5 |
| 7 | 52.9 | (224, out of range) | 42.5 | 31.8-44.5 |
| 20 | 89.4 | (415, out of range) | 42.5 | 31.8-44.5 |

## T3 realistic rows (Chunc head basket, t_max 0.2 s, plateau = median 0.1-0.2 s; Peretz ignition h = main default)

| pyrogen | g | DeMar | model | first ign ms | t90 ms | peak MPa | t_peak ms | plateau | peak/plateau | pyro t50 | t90 | duration ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bpnv | 4.5 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 2.19 | 3.23 | 20.26 | 10.7 | 8.44 | 2.40 | 7.3 | 10.1 | 10.7 |
| bpnv | 4.5 | on | I-A pellet, constant-area basket law (opt-in) | 2.44 | 3.65 | 20.03 | 12.2 | 8.44 | 2.37 | 8.4 | 11.5 | 12.2 |
| bpnv | 4.5 | on | **new default**: I-A pellet, exact cylinder area law | 2.51 | 3.84 | 16.72 | 18.3 | 8.44 | 1.98 | 9.8 | 17.3 | 25.4 |
| bpnv | 4.5 | on | new default + rate capped at the Zhu plateau (opt-in) | 2.51 | 3.84 | 15.05 | 22.2 | 8.44 | 1.78 | 11.0 | 24.5 | 39.7 |
| bpnv | 4.5 | on | 3 mm cube, constant-area law | 2.07 | 3.04 | 20.38 | 10.0 | 8.44 | 2.42 | 6.8 | 9.4 | 10.0 |
| bpnv | 4.5 | on | 3 mm cube, exact area law (L/D 1) | 2.13 | 3.15 | 17.17 | 15.6 | 8.44 | 2.03 | 8.0 | 14.5 | 23.0 |
| bpnv | 4.5 | on | 3 mm cube, exact area law + Zhu cap | 2.13 | 3.15 | 15.49 | 19.4 | 8.44 | 1.84 | 9.1 | 20.7 | 37.0 |
| bpnv | 4.5 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 2.73 | 3.43 | 20.03 | 10.9 | 8.44 | 2.37 | 7.5 | 10.3 | 10.9 |
| bpnv | 4.5 | off | I-A pellet, constant-area basket law (opt-in) | 3.19 | 3.91 | 19.80 | 12.4 | 8.44 | 2.34 | 8.5 | 11.7 | 12.4 |
| bpnv | 4.5 | off | **new default**: I-A pellet, exact cylinder area law | 3.41 | 4.16 | 16.52 | 18.8 | 8.44 | 1.96 | 10.1 | 17.7 | 25.7 |
| bpnv | 4.5 | off | new default + rate capped at the Zhu plateau (opt-in) | 3.41 | 4.16 | 14.85 | 22.7 | 8.45 | 1.76 | 11.2 | 24.7 | 39.9 |
| bpnv | 4.5 | off | 3 mm cube, constant-area law | 2.51 | 3.19 | 20.17 | 10.1 | 8.44 | 2.39 | 6.9 | 9.6 | 10.1 |
| bpnv | 4.5 | off | 3 mm cube, exact area law (L/D 1) | 2.65 | 3.34 | 16.98 | 15.9 | 8.44 | 2.01 | 8.2 | 14.7 | 23.3 |
| bpnv | 4.5 | off | 3 mm cube, exact area law + Zhu cap | 2.65 | 3.34 | 15.30 | 19.7 | 8.44 | 1.81 | 9.1 | 20.8 | 37.1 |
| bpnv | 6.0 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 1.43 | 2.48 | 26.29 | 9.3 | 8.39 | 3.13 | 6.3 | 8.7 | 9.3 |
| bpnv | 6.0 | on | I-A pellet, constant-area basket law (opt-in) | 1.60 | 2.75 | 25.90 | 10.5 | 8.40 | 3.08 | 7.1 | 9.9 | 10.5 |
| bpnv | 6.0 | on | **new default**: I-A pellet, exact cylinder area law | 1.63 | 2.81 | 20.89 | 15.4 | 8.39 | 2.49 | 8.3 | 15.0 | 22.3 |
| bpnv | 6.0 | on | new default + rate capped at the Zhu plateau (opt-in) | 1.63 | 2.81 | 17.84 | 19.6 | 8.39 | 2.13 | 10.2 | 23.7 | 38.9 |
| bpnv | 6.0 | on | 3 mm cube, constant-area law | 1.35 | 2.36 | 26.49 | 8.6 | 8.39 | 3.16 | 5.8 | 8.1 | 8.6 |
| bpnv | 6.0 | on | 3 mm cube, exact area law (L/D 1) | 1.37 | 2.40 | 21.61 | 13.3 | 8.39 | 2.58 | 6.9 | 12.6 | 20.3 |
| bpnv | 6.0 | on | 3 mm cube, exact area law + Zhu cap | 1.37 | 2.40 | 18.55 | 17.2 | 8.39 | 2.21 | 8.4 | 20.1 | 36.5 |
| bpnv | 6.0 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 1.97 | 2.64 | 26.04 | 9.4 | 8.39 | 3.10 | 6.4 | 8.8 | 9.4 |
| bpnv | 6.0 | off | I-A pellet, constant-area basket law (opt-in) | 2.28 | 2.95 | 25.62 | 10.6 | 8.40 | 3.05 | 7.3 | 10.0 | 10.6 |
| bpnv | 6.0 | off | **new default**: I-A pellet, exact cylinder area law | 2.37 | 3.05 | 20.61 | 15.8 | 8.40 | 2.45 | 8.5 | 15.3 | 22.6 |
| bpnv | 6.0 | off | new default + rate capped at the Zhu plateau (opt-in) | 2.37 | 3.05 | 17.59 | 20.1 | 8.40 | 2.10 | 10.3 | 23.8 | 39.1 |
| bpnv | 6.0 | off | 3 mm cube, constant-area law | 1.82 | 2.48 | 26.26 | 8.7 | 8.39 | 3.13 | 5.9 | 8.2 | 8.7 |
| bpnv | 6.0 | off | 3 mm cube, exact area law (L/D 1) | 1.88 | 2.56 | 21.38 | 13.5 | 8.39 | 2.55 | 7.0 | 12.7 | 20.5 |
| bpnv | 6.0 | off | 3 mm cube, exact area law + Zhu cap | 1.88 | 2.56 | 18.32 | 17.6 | 8.39 | 2.18 | 8.5 | 20.1 | 36.6 |
| mtv | 4.5 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 2.31 | 8.36 | 17.92 | 42.4 | 8.42 | 2.13 | 27.8 | 42.0 | 45.5 |
| mtv | 4.5 | on | I-A pellet, constant-area basket law (opt-in) | 2.48 | 10.80 | 16.95 | 45.4 | 8.44 | 2.01 | 32.2 | 48.3 | 52.3 |
| mtv | 4.5 | on | **new default**: I-A pellet, exact cylinder area law | 2.49 | 12.06 | 14.33 | 40.2 | 8.46 | 1.69 | 38.5 | 81.3 | 134.4 |
| mtv | 4.5 | on | 3 mm cube, constant-area law | 2.22 | 7.16 | 18.56 | 40.8 | 8.41 | 2.21 | 25.5 | 38.8 | 42.0 |
| mtv | 4.5 | on | 3 mm cube, exact area law (L/D 1) | 2.23 | 7.83 | 15.05 | 33.7 | 8.42 | 1.79 | 30.6 | 66.5 | 123.6 |
| mtv | 4.5 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 9.88 | 10.43 | 17.91 | 44.8 | 8.43 | 2.12 | 29.5 | 43.7 | 47.2 |
| mtv | 4.5 | off | I-A pellet, constant-area basket law (opt-in) | 13.07 | 13.58 | 16.93 | 48.6 | 8.45 | 2.00 | 34.4 | 50.6 | 54.6 |
| mtv | 4.5 | off | **new default**: I-A pellet, exact cylinder area law | 15.35 | 15.85 | 14.20 | 44.4 | 8.49 | 1.67 | 41.5 | 84.3 | 137.4 |
| mtv | 4.5 | off | 3 mm cube, constant-area law | 8.34 | 8.94 | 18.54 | 42.9 | 8.42 | 2.20 | 27.0 | 40.2 | 43.5 |
| mtv | 4.5 | off | 3 mm cube, exact area law (L/D 1) | 9.53 | 10.09 | 14.91 | 36.4 | 8.44 | 1.77 | 32.5 | 68.4 | 125.4 |
| mtv | 6.0 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 1.77 | 4.31 | 21.07 | 36.0 | 8.35 | 2.52 | 23.7 | 37.1 | 40.4 |
| mtv | 6.0 | on | I-A pellet, constant-area basket law (opt-in) | 1.92 | 5.03 | 19.71 | 37.6 | 8.35 | 2.36 | 26.7 | 41.9 | 45.7 |
| mtv | 6.0 | on | **new default**: I-A pellet, exact cylinder area law | 1.93 | 5.15 | 16.45 | 30.8 | 8.36 | 1.97 | 31.9 | 73.4 | 126.5 |
| mtv | 6.0 | on | 3 mm cube, constant-area law | 1.69 | 3.97 | 21.94 | 35.2 | 8.35 | 2.63 | 22.1 | 34.5 | 37.6 |
| mtv | 6.0 | on | 3 mm cube, exact area law (L/D 1) | 1.69 | 4.05 | 17.50 | 27.7 | 8.35 | 2.10 | 26.5 | 61.1 | 118.0 |
| mtv | 6.0 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 4.74 | 5.43 | 21.04 | 37.7 | 8.36 | 2.52 | 24.7 | 38.1 | 41.4 |
| mtv | 6.0 | off | I-A pellet, constant-area basket law (opt-in) | 6.32 | 6.98 | 19.69 | 40.2 | 8.37 | 2.35 | 28.4 | 43.6 | 47.4 |
| mtv | 6.0 | off | **new default**: I-A pellet, exact cylinder area law | 6.88 | 7.53 | 16.24 | 34.1 | 8.38 | 1.94 | 34.0 | 75.5 | 128.6 |
| mtv | 6.0 | off | 3 mm cube, constant-area law | 4.14 | 4.83 | 21.92 | 36.5 | 8.35 | 2.62 | 22.9 | 35.3 | 38.4 |
| mtv | 6.0 | off | 3 mm cube, exact area law (L/D 1) | 4.32 | 5.00 | 17.34 | 29.2 | 8.36 | 2.08 | 27.4 | 62.0 | 118.9 |

## T4 t_max 1.0 s, plateau = median 0.4-0.9 s (Chunc experiment: peak 8.89, plateau 8.66, ratio 1.027)

| pyrogen | g | DeMar | model | peak MPa | t_peak ms | plateau | peak/plateau | pyro duration ms |
|---|---|---|---|---|---|---|---|---|
| bpnv | 4.5 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 20.26 | 10.7 | 8.16 | 2.48 | 10.7 |
| bpnv | 4.5 | on | **new default**: I-A pellet, exact cylinder area law | 16.72 | 18.3 | 8.16 | 2.05 | 25.4 |
| bpnv | 4.5 | on | 3 mm cube, exact area law (L/D 1) | 17.17 | 15.6 | 8.16 | 2.10 | 23.0 |
| bpnv | 4.5 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 20.03 | 10.9 | 8.16 | 2.45 | 10.9 |
| bpnv | 4.5 | off | **new default**: I-A pellet, exact cylinder area law | 16.52 | 18.8 | 8.16 | 2.02 | 25.7 |
| bpnv | 4.5 | off | 3 mm cube, exact area law (L/D 1) | 16.98 | 15.9 | 8.16 | 2.08 | 23.3 |
| bpnv | 6.0 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 26.29 | 9.3 | 8.17 | 3.22 | 9.3 |
| bpnv | 6.0 | on | **new default**: I-A pellet, exact cylinder area law | 20.89 | 15.4 | 8.17 | 2.56 | 22.3 |
| bpnv | 6.0 | on | 3 mm cube, exact area law (L/D 1) | 21.61 | 13.3 | 8.17 | 2.64 | 20.3 |
| bpnv | 6.0 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 26.04 | 9.4 | 8.17 | 3.19 | 9.4 |
| bpnv | 6.0 | off | **new default**: I-A pellet, exact cylinder area law | 20.61 | 15.8 | 8.17 | 2.52 | 22.6 |
| bpnv | 6.0 | off | 3 mm cube, exact area law (L/D 1) | 21.38 | 13.5 | 8.17 | 2.62 | 20.5 |
| mtv | 4.5 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 17.92 | 42.4 | 8.16 | 2.20 | 45.5 |
| mtv | 4.5 | on | **new default**: I-A pellet, exact cylinder area law | 14.33 | 40.2 | 8.15 | 1.76 | 134.4 |
| mtv | 4.5 | on | 3 mm cube, exact area law (L/D 1) | 15.05 | 33.7 | 8.16 | 1.85 | 123.6 |
| mtv | 4.5 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 17.91 | 44.8 | 8.16 | 2.20 | 47.2 |
| mtv | 4.5 | off | **new default**: I-A pellet, exact cylinder area law | 14.20 | 44.4 | 8.14 | 1.74 | 137.4 |
| mtv | 4.5 | off | 3 mm cube, exact area law (L/D 1) | 14.91 | 36.4 | 8.15 | 1.83 | 125.4 |
| mtv | 6.0 | on | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 21.07 | 36.0 | 8.18 | 2.58 | 40.4 |
| mtv | 6.0 | on | **new default**: I-A pellet, exact cylinder area law | 16.45 | 30.8 | 8.17 | 2.01 | 126.5 |
| mtv | 6.0 | on | 3 mm cube, exact area law (L/D 1) | 17.50 | 27.7 | 8.17 | 2.14 | 118.0 |
| mtv | 6.0 | off | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 21.04 | 37.7 | 8.17 | 2.57 | 41.4 |
| mtv | 6.0 | off | **new default**: I-A pellet, exact cylinder area law | 16.24 | 34.1 | 8.16 | 1.99 | 128.6 |
| mtv | 6.0 | off | 3 mm cube, exact area law (L/D 1) | 17.34 | 29.2 | 8.17 | 2.12 | 118.9 |

## T5 ignition thresholds (min BPNV charge with >= 99 % of the 100 grain cells lit within 0.2 s; grid scan then bisection to 2 %)

| topology | model | DeMar on [g] | DeMar off [g] |
|---|---|---|---|
| head_basket | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 1.43 | 1.43 |
| head_basket | I-A pellet, constant-area basket law (opt-in) | 1.51 | 1.51 |
| head_basket | **new default**: I-A pellet, exact cylinder area law | 3.01 | 3.01 |
| head_basket | 3 mm cube, exact area law (L/D 1) | 2.71 | 2.71 |
| forward_plenum | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 0.0268 | 2.21 |
| forward_plenum | **new default**: I-A pellet, exact cylinder area law | 0.0289 | 2.27 |
| forward_plenum | 3 mm cube, exact area law (L/D 1) | 0.0258 | 2.17 |

## T6 peak vs charge (BPNV head basket; peak MPa (peak/plateau); NL = not lit within 0.2 s)

| DeMar | model | 0.6 g | 0.9 g | 1.36 g | 1.8 g | 2.27 g | 2.93 g | 3.6 g | 4.54 g | 6.3 g | 9.0 g |
|---|---|---|---|---|---|---|---|---|---|---|---|
| on | pre_ob1 | NL | NL | NL | 9.9 (1.13) | 10.2 (1.19) | 13.4 (1.57) | 16.3 (1.93) | 20.4 (2.42) | 27.5 (3.28) | 37.5 (4.51) |
| on | lib | NL | NL | NL | NL | NL | NL | 14.0 (1.65) | 16.8 (2.00) | 21.8 (2.60) | 29.0 (3.48) |
| on | cube3 | NL | NL | NL | NL | NL | 12.1 (1.46) | 14.4 (1.70) | 17.3 (2.05) | 22.5 (2.69) | 30.1 (3.61) |
| off | pre_ob1 | NL | NL | NL | 9.9 (1.13) | 10.2 (1.19) | 13.4 (1.57) | 16.3 (1.93) | 20.2 (2.39) | 27.3 (3.25) | 37.3 (4.48) |
| off | lib | NL | NL | NL | NL | NL | NL | 14.0 (1.65) | 16.6 (1.97) | 21.5 (2.56) | 28.7 (3.45) |
| off | cube3 | NL | NL | NL | NL | NL | 12.1 (1.46) | 14.4 (1.70) | 17.1 (2.03) | 22.3 (2.66) | 29.9 (3.59) |

## T7 canonical motors, library BPNV, Sutton-mass forward plenum, Peretz ignition h (main default); old = pre-OB1 library pellet (3.2 mm, L/D 1), new = I-A pellet (3.175 mm, L/D 1.504)

| motor | pellet | DeMar | rate law | lit | first ign ms | peak MPa | t_peak | plateau / trough / mid | ratios |
|---|---|---|---|---|---|---|---|---|---|
| chunc | old 3.2 mm L/D 1 | on | library | yes | 5.81 | 9.9509 | 36.7 ms | 8.142 | pk/pl 1.222 |
| chunc | new I-A | on | library | yes | 5.86 | 9.9774 | 35.2 ms | 8.142 | pk/pl 1.225 |
| chunc | new I-A | on | capped (opt-in) | yes | 5.99 | 10.0459 | 35.4 ms | 8.142 | pk/pl 1.234 |
| chunc | old 3.2 mm L/D 1 | off | library | NO | - | 0.73 | | | |
| chunc | new I-A | off | library | NO | - | 0.69 | | | |
| balls | old 3.2 mm L/D 1 | on | library | yes | 7.22 | 7.0203 | 39.2 ms | 4.343 / 4.855 | 1.62 / 1.45 |
| balls | new I-A | on | library | yes | 7.27 | 7.0609 | 37.4 ms | 4.343 / 4.856 | 1.63 / 1.45 |
| balls | new I-A | on | capped (opt-in) | yes | 7.37 | 7.0980 | 37.2 ms | 4.343 / 4.856 | 1.63 / 1.46 |
| balls | old 3.2 mm L/D 1 | off | library | NO | - | 0.86 | | | |
| balls | new I-A | off | library | NO | - | 0.81 | | | |
| haseg | old 3.2 mm L/D 1 | on | library | yes | 7.69 | 6.1073 | 2.411 s |  |  |
| haseg | new I-A | on | library | yes | 7.65 | 6.1073 | 2.411 s |  |  |
| haseg | new I-A | on | capped (opt-in) | yes | 7.47 | 6.1073 | 2.411 s |  |  |
| haseg | old 3.2 mm L/D 1 | off | library | NO | - | 0.44 | | | |
| haseg | new I-A | off | library | NO | - | 0.42 | | | |

## T8 attribution and ablation (BPNV 4.5 g head basket, DeMar on, t_max 0.2 s)

| model | variant | first ign ms | peak MPa | plateau (0.1-0.2 s) | peak/plateau | pyro duration ms |
|---|---|---|---|---|---|---|
| pre_ob1 | Ma 2020 erosive (default) | 2.19 | 20.26 | 8.44 | 2.40 | 10.7 |
| pre_ob1 | erosive = saint_robert | 2.19 | 7.06 | 4.50 | 1.57 | 13.8 |
| pre_ob1 | erosive = mukunda_paul | 2.19 | 14.08 | 6.70 | 2.10 | 11.7 |
| pre_ob1 | grain cannot ignite (pyrogen gas only) | - | 2.04 | - | - | - |
| lib | Ma 2020 erosive (default) | 2.51 | 16.72 | 8.44 | 1.98 | 25.4 |
| lib | erosive = saint_robert | 2.51 | 5.62 | 4.50 | 1.25 | 39.1 |
| lib | erosive = mukunda_paul | 2.51 | 10.86 | 6.69 | 1.62 | 29.8 |
| lib | grain cannot ignite (pyrogen gas only) | - | 1.42 | - | - | - |
| lib_capZ | Ma 2020 erosive (default) | 2.51 | 15.05 | 8.44 | 1.78 | 39.7 |
| lib_capZ | erosive = saint_robert | 2.51 | 5.53 | 4.50 | 1.23 | 40.5 |
| lib_capZ | erosive = mukunda_paul | 2.51 | 10.02 | 6.69 | 1.50 | 39.9 |
| lib_capZ | grain cannot ignite (pyrogen gas only) | - | 1.42 | - | - | - |
| lib | rate flat above 7 MPa (library law) | 2.51 | 15.64 | 8.44 | 1.85 | 33.2 |
| lib | Zhu law, r0 = 5.271e-3 | NOT LIT | 0.17 | | | |
| lib | Zhu law, r0 = 6.074e-3 (mean) | NOT LIT | 0.19 | | | |
| lib | Zhu law, r0 = 6.848e-3 | NOT LIT | 0.21 | | | |
| lib | diagnostic bed flame speed 3 m/s | 8.13 | 16.62 | 8.47 | 1.96 | 31.7 |
| lib | diagnostic bed flame speed 10 m/s | 3.08 | 16.95 | 8.44 | 2.01 | 26.5 |
| lib | diagnostic bed flame speed 30 m/s | 2.63 | 16.84 | 8.44 | 1.99 | 25.6 |
