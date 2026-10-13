## T1 rate laws (mm/s)

| P [MPa] | sim law 2e-5 P^0.5 | Nakano 71.1 P^0.589 (valid 0.005-0.2 MPa) | YAML "candidate" 3.5e-5 P^0.589 | Zhu mean r0 (valid 2.8-50 MPa) | MIL-P-46994B (closed-bomb average) |
|---|---|---|---|---|---|
| 0.1 | 6.3 | 18.3 | 30.8 | (below range) |  |
| 0.2 | 8.9 | 27.6 | 46.4 | (below range) |  |
| 1 | 20.0 | (71, out of range) | 119.7 | (below range) |  |
| 2.8 | 33.5 | (130, out of range) | 219.5 | 17.0 | 31.8-44.5 |
| 7 | 52.9 | (224, out of range) | 376.6 | 42.5 | 31.8-44.5 |
| 20 | 89.4 | (415, out of range) | 698.8 | 42.5 | 31.8-44.5 |

## T2 ranking probes (BPNV 4.5 g, 3 mm cubes, head basket, t_max 0.2 s, plateau = median 0.1-0.2 s)

### combo = current

| option | first ign ms | t90 ms | peak MPa | d peak vs base | t_peak ms | peak/plateau | pyro t50 ms | pyro t90 ms | pyro duration ms |
|---|---|---|---|---|---|---|---|---|---|
| base: current pyrogen model (library law, constant area) | 0.66 | 1.76 | 21.62 | +0.0 % | 9.0 | 2.57 | 6.0 | 8.5 | 9.1 |
| area23: (b) shrinking area A0 (m/m0)^(2/3), library law | 0.66 | 1.77 | 18.10 | -16.3 % | 14.0 | 2.15 | 7.0 | 13.2 | 21.6 |
| cap7: library law, rate flat above 7 MPa (IG1 probe) | 0.66 | 1.76 | 20.61 | -4.7 % | 11.1 | 2.45 | 6.3 | 10.1 | 11.1 |
| capZ: (a) library law capped at Zhu plateau rate 42.5 mm/s (P_cap 4.5 MPa) | 0.66 | 1.76 | 20.05 | -7.3 % | 12.9 | 2.38 | 7.0 | 11.7 | 12.9 |
| zhu: (a) Zhu piecewise law r = r0 P, flat above 7 MPa, constant area | 5.79 | 14.88 | 20.20 | -6.6 % | 28.9 | 2.38 | 23.0 | 27.7 | 28.9 |
| zhu+area23: (a)+(b) Zhu law + shrinking area | 5.94 | 14.97 | 16.27 | -24.7 % | 33.9 | 1.91 | 24.6 | 36.2 | 52.6 |
| capZ+area23: (a)+(b) capped library law + shrinking area | 0.66 | 1.77 | 16.13 | -25.4 % | 17.8 | 1.92 | 8.5 | 20.1 | 36.4 |
| area23+bed3: (b) + DIAGNOSTIC bed flame speed 3 m/s (no sourced value) | 2.82 | 4.23 | 18.17 | -16.0 % | 18.2 | 2.15 | 11.3 | 17.5 | 27.3 |
| area23+bed10: (b) + DIAGNOSTIC bed flame speed 10 m/s (no sourced value) | 1.87 | 3.04 | 17.84 | -17.5 % | 15.6 | 2.12 | 8.5 | 14.8 | 23.4 |
| area23+bed30: (b) + DIAGNOSTIC bed flame speed 30 m/s (no sourced value) | 1.03 | 2.19 | 18.03 | -16.6 % | 14.5 | 2.14 | 7.5 | 13.7 | 22.2 |

### combo = forward

| option | first ign ms | t90 ms | peak MPa | d peak vs base | t_peak ms | peak/plateau | pyro t50 ms | pyro t90 ms | pyro duration ms |
|---|---|---|---|---|---|---|---|---|---|
| base: current pyrogen model (library law, constant area) | 2.51 | 3.19 | 20.17 | +0.0 % | 10.1 | 2.39 | 6.9 | 9.6 | 10.1 |
| area23: (b) shrinking area A0 (m/m0)^(2/3), library law | 2.65 | 3.34 | 16.98 | -15.8 % | 15.9 | 2.01 | 8.2 | 14.7 | 23.3 |
| cap7: library law, rate flat above 7 MPa (IG1 probe) | 2.51 | 3.19 | 19.40 | -3.8 % | 11.8 | 2.30 | 7.1 | 10.9 | 11.8 |
| capZ: (a) library law capped at Zhu plateau rate 42.5 mm/s (P_cap 4.5 MPa) | 2.51 | 3.19 | 18.94 | -6.1 % | 13.5 | 2.25 | 7.6 | 12.3 | 13.5 |
| zhu: (a) Zhu piecewise law r = r0 P, flat above 7 MPa, constant area | NOT LIT | - | 0.23 | - | - | - | - | - | - |
| zhu+area23: (a)+(b) Zhu law + shrinking area | NOT LIT | - | 0.23 | - | - | - | - | - | - |
| capZ+area23: (a)+(b) capped library law + shrinking area | 2.65 | 3.34 | 15.30 | -24.1 % | 19.7 | 1.81 | 9.1 | 20.8 | 37.1 |
| area23+bed3: (b) + DIAGNOSTIC bed flame speed 3 m/s (no sourced value) | 6.73 | 7.51 | 17.03 | -15.6 % | 21.0 | 2.01 | 13.3 | 19.8 | 29.2 |
| area23+bed10: (b) + DIAGNOSTIC bed flame speed 10 m/s (no sourced value) | 3.37 | 4.17 | 17.21 | -14.7 % | 16.9 | 2.04 | 9.4 | 15.8 | 24.6 |
| area23+bed30: (b) + DIAGNOSTIC bed flame speed 30 m/s (no sourced value) | 2.78 | 3.52 | 17.10 | -15.2 % | 16.1 | 2.03 | 8.5 | 14.9 | 23.6 |

## T3 realistic rows (Chunc head basket, 3 mm cubes L/D 1, t_max 0.2 s, plateau = median 0.1-0.2 s)

| pyrogen | g | combo | option | first ign ms | t90 ms | peak MPa | t_peak ms | plateau | peak/plateau | pyro t50 | t90 | duration ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bpnv | 4.5 | current | base | 0.66 | 1.76 | 21.62 | 9.0 | 8.42 | 2.57 | 6.0 | 8.5 | 9.1 |
| bpnv | 4.5 | current | area23 | 0.66 | 1.77 | 18.10 | 14.0 | 8.42 | 2.15 | 7.0 | 13.2 | 21.6 |
| bpnv | 4.5 | current | cap7 | 0.66 | 1.76 | 20.61 | 11.1 | 8.42 | 2.45 | 6.3 | 10.1 | 11.1 |
| bpnv | 4.5 | current | capZ | 0.66 | 1.76 | 20.05 | 12.9 | 8.42 | 2.38 | 7.0 | 11.7 | 12.9 |
| bpnv | 4.5 | current | zhu | 5.79 | 14.88 | 20.20 | 28.9 | 8.50 | 2.38 | 23.0 | 27.7 | 28.9 |
| bpnv | 4.5 | current | zhu+area23 | 5.94 | 14.97 | 16.27 | 33.9 | 8.50 | 1.91 | 24.6 | 36.2 | 52.6 |
| bpnv | 4.5 | current | capZ+area23 | 0.66 | 1.77 | 16.13 | 17.8 | 8.42 | 1.92 | 8.5 | 20.1 | 36.4 |
| bpnv | 4.5 | forward | base | 2.51 | 3.19 | 20.17 | 10.1 | 8.44 | 2.39 | 6.9 | 9.6 | 10.1 |
| bpnv | 4.5 | forward | area23 | 2.65 | 3.34 | 16.98 | 15.9 | 8.44 | 2.01 | 8.2 | 14.7 | 23.3 |
| bpnv | 4.5 | forward | cap7 | 2.51 | 3.19 | 19.40 | 11.8 | 8.44 | 2.30 | 7.1 | 10.9 | 11.8 |
| bpnv | 4.5 | forward | capZ | 2.51 | 3.19 | 18.94 | 13.5 | 8.44 | 2.25 | 7.6 | 12.3 | 13.5 |
| bpnv | 4.5 | forward | zhu | NOT LIT (peak 0.23) | | | | | | | | |
| bpnv | 4.5 | forward | zhu+area23 | NOT LIT (peak 0.23) | | | | | | | | |
| bpnv | 4.5 | forward | capZ+area23 | 2.65 | 3.34 | 15.30 | 19.7 | 8.44 | 1.81 | 9.1 | 20.8 | 37.1 |
| bpnv | 6.0 | current | base | 0.40 | 1.54 | 27.41 | 8.1 | 8.38 | 3.27 | 5.4 | 7.6 | 8.1 |
| bpnv | 6.0 | current | area23 | 0.40 | 1.54 | 22.33 | 12.4 | 8.38 | 2.66 | 6.3 | 11.9 | 19.5 |
| bpnv | 6.0 | current | cap7 | 0.40 | 1.54 | 25.45 | 10.7 | 8.38 | 3.04 | 6.0 | 9.7 | 10.7 |
| bpnv | 6.0 | current | capZ | 0.40 | 1.54 | 24.51 | 12.6 | 8.38 | 2.92 | 6.8 | 11.5 | 12.6 |
| bpnv | 6.0 | current | zhu | 3.81 | 6.03 | 24.64 | 19.7 | 8.41 | 2.93 | 13.8 | 18.5 | 19.7 |
| bpnv | 6.0 | current | zhu+area23 | 3.83 | 6.04 | 19.12 | 23.3 | 8.41 | 2.27 | 15.2 | 26.9 | 43.2 |
| bpnv | 6.0 | current | capZ+area23 | 0.40 | 1.54 | 19.01 | 16.4 | 8.38 | 2.27 | 8.2 | 19.8 | 36.2 |
| bpnv | 6.0 | forward | base | 1.82 | 2.48 | 26.26 | 8.7 | 8.39 | 3.13 | 5.9 | 8.2 | 8.7 |
| bpnv | 6.0 | forward | area23 | 1.88 | 2.56 | 21.38 | 13.5 | 8.39 | 2.55 | 7.0 | 12.7 | 20.5 |
| bpnv | 6.0 | forward | cap7 | 1.82 | 2.48 | 24.54 | 11.1 | 8.39 | 2.92 | 6.4 | 10.1 | 11.1 |
| bpnv | 6.0 | forward | capZ | 1.82 | 2.48 | 23.69 | 12.9 | 8.39 | 2.82 | 7.0 | 11.7 | 12.9 |
| bpnv | 6.0 | forward | zhu | NOT LIT (peak 0.28) | | | | | | | | |
| bpnv | 6.0 | forward | zhu+area23 | NOT LIT (peak 0.28) | | | | | | | | |
| bpnv | 6.0 | forward | capZ+area23 | 1.88 | 2.56 | 18.32 | 17.6 | 8.39 | 2.18 | 8.5 | 20.1 | 36.6 |
| mtv | 4.5 | current | base | 0.96 | 2.16 | 18.59 | 34.8 | 8.37 | 2.22 | 21.7 | 34.8 | 38.0 |
| mtv | 4.5 | current | area23 | 0.96 | 2.17 | 15.43 | 26.9 | 8.38 | 1.84 | 26.1 | 62.0 | 119.2 |
| mtv | 4.5 | forward | base | 8.34 | 8.94 | 18.54 | 42.9 | 8.42 | 2.20 | 27.0 | 40.2 | 43.5 |
| mtv | 4.5 | forward | area23 | 9.53 | 10.09 | 14.91 | 36.4 | 8.44 | 1.77 | 32.5 | 68.4 | 125.4 |
| mtv | 6.0 | current | base | 0.62 | 1.91 | 21.98 | 32.5 | 8.33 | 2.64 | 20.5 | 32.9 | 35.9 |
| mtv | 6.0 | current | area23 | 0.62 | 1.91 | 17.79 | 24.9 | 8.33 | 2.14 | 24.8 | 59.4 | 116.3 |
| mtv | 6.0 | forward | base | 4.14 | 4.83 | 21.92 | 36.5 | 8.35 | 2.62 | 22.9 | 35.3 | 38.4 |
| mtv | 6.0 | forward | area23 | 4.32 | 5.00 | 17.34 | 29.2 | 8.36 | 2.08 | 27.4 | 62.0 | 118.9 |

## T4 t_max 1.0 s, plateau = median 0.4-0.9 s (Chunc experiment: peak 8.89, plateau 8.66, ratio 1.027)

| pyrogen | g | combo | option | peak MPa | t_peak ms | plateau | peak/plateau | pyro duration ms |
|---|---|---|---|---|---|---|---|---|
| bpnv | 4.5 | current | base | 21.62 | 9.0 | 8.17 | 2.65 | 9.1 |
| bpnv | 4.5 | current | area23 | 18.10 | 14.0 | 8.17 | 2.22 | 21.6 |
| bpnv | 4.5 | current | zhu+area23 | 16.27 | 33.9 | 8.15 | 2.00 | 52.6 |
| bpnv | 4.5 | current | capZ+area23 | 16.13 | 17.8 | 8.17 | 1.98 | 36.4 |
| bpnv | 4.5 | forward | base | 20.17 | 10.1 | 8.16 | 2.47 | 10.1 |
| bpnv | 4.5 | forward | area23 | 16.98 | 15.9 | 8.16 | 2.08 | 23.3 |
| bpnv | 4.5 | forward | zhu+area23 | NOT LIT | | | | |
| bpnv | 4.5 | forward | capZ+area23 | 15.30 | 19.7 | 8.16 | 1.88 | 37.1 |
| bpnv | 6.0 | current | base | 27.41 | 8.1 | 8.18 | 3.35 | 8.1 |
| bpnv | 6.0 | current | area23 | 22.33 | 12.4 | 8.18 | 2.73 | 19.5 |
| bpnv | 6.0 | current | zhu+area23 | 19.12 | 23.3 | 8.17 | 2.34 | 43.2 |
| bpnv | 6.0 | current | capZ+area23 | 19.01 | 16.4 | 8.17 | 2.33 | 36.2 |
| bpnv | 6.0 | forward | base | 26.26 | 8.7 | 8.17 | 3.21 | 8.7 |
| bpnv | 6.0 | forward | area23 | 21.38 | 13.5 | 8.17 | 2.62 | 20.5 |
| bpnv | 6.0 | forward | zhu+area23 | NOT LIT | | | | |
| bpnv | 6.0 | forward | capZ+area23 | 18.32 | 17.6 | 8.17 | 2.24 | 36.6 |
| mtv | 4.5 | current | base | 18.59 | 34.8 | 8.17 | 2.27 | 38.0 |
| mtv | 4.5 | current | area23 | 15.43 | 26.9 | 8.17 | 1.89 | 119.2 |
| mtv | 4.5 | forward | base | 18.54 | 42.9 | 8.16 | 2.27 | 43.5 |
| mtv | 4.5 | forward | area23 | 14.91 | 36.4 | 8.15 | 1.83 | 125.4 |
| mtv | 6.0 | current | base | 21.98 | 32.5 | 8.18 | 2.69 | 35.9 |
| mtv | 6.0 | current | area23 | 17.79 | 24.9 | 8.18 | 2.18 | 116.3 |
| mtv | 6.0 | forward | base | 21.92 | 36.5 | 8.18 | 2.68 | 38.4 |
| mtv | 6.0 | forward | area23 | 17.34 | 29.2 | 8.17 | 2.12 | 118.9 |

## T5 head-basket ignition threshold (min BPNV charge with >= 99 % of 100 grain cells lit within 0.2 s; grid then bisection to 2 %)

| option | current combo [g] | forward combo [g] |
|---|---|---|
| base | 0.332 | 1.38 |
| area23 | 0.52 | 2.71 |
| capZ+area23 | 0.52 | 2.71 |
| zhu+area23 | 1.62 | 10.6 |

## T6 peak vs charge (BPNV head basket; peak MPa (peak/plateau); NL = not lit within 0.2 s)

| combo | option | 0.6 g | 0.9 g | 1.36 g | 1.8 g | 2.27 g | 2.93 g | 3.6 g | 4.54 g | 6.3 g | 9.0 g |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current | base | 9.9 (1.15) | 10.0 (1.17) | 10.3 (1.21) | 11.3 (1.33) | 12.8 (1.51) | 15.3 (1.81) | 17.9 (2.12) | 21.8 (2.59) | 28.6 (3.42) | 38.5 (4.63) |
| current | area23 | 9.9 (1.15) | 10.0 (1.17) | 10.5 (1.23) | 11.3 (1.33) | 12.2 (1.44) | 13.7 (1.62) | 15.5 (1.83) | 18.2 (2.16) | 23.2 (2.77) | 30.6 (3.68) |
| current | capZ+area23 | 10.0 (1.15) | 10.0 (1.17) | 10.5 (1.23) | 11.2 (1.31) | 11.9 (1.40) | 13.0 (1.54) | 14.3 (1.69) | 16.2 (1.93) | 19.6 (2.34) | 24.3 (2.92) |
| forward | base | NL | NL | NL | 9.9 (1.14) | 10.2 (1.20) | 13.5 (1.59) | 16.4 (1.94) | 20.3 (2.41) | 27.5 (3.28) | 37.7 (4.52) |
| forward | area23 | NL | NL | NL | NL | NL | 12.1 (1.46) | 14.4 (1.70) | 17.1 (2.03) | 22.3 (2.66) | 29.9 (3.59) |
| forward | capZ+area23 | NL | NL | NL | NL | NL | 11.7 (1.41) | 13.2 (1.59) | 15.4 (1.82) | 19.0 (2.26) | 23.8 (2.86) |

## T7 canonical motors, library BPNV, Sutton-mass forward plenum, current ignition combo (law = rate law; area law not applicable to the plenum)

| motor | law | first ign ms | peak MPa | d peak vs lib | t_peak | plateau / trough / mid | ratios |
|---|---|---|---|---|---|---|---|
| chunc | lib | 1.72 | 10.0827 | +0.00 % | 27.9 ms | 8.147 | pk/pl 1.238 |
| chunc | cap7_lib | 1.72 | 10.1038 | +0.21 % | 27.9 ms | 8.147 | pk/pl 1.240 |
| chunc | capZ_lib | 2.03 | 10.1604 | +0.77 % | 29.2 ms | 8.146 | pk/pl 1.247 |
| chunc | zhu_mean | 11.96 | 10.3661 | +2.81 % | 39.8 ms | 8.135 | pk/pl 1.274 |
| balls | lib | 1.75 | 7.2577 | +0.00 % | 27.2 ms | 4.343 / 4.868 | 1.67 / 1.49 |
| balls | cap7_lib | 1.75 | 7.2572 | -0.01 % | 27.2 ms | 4.343 / 4.868 | 1.67 / 1.49 |
| balls | capZ_lib | 2.26 | 7.2742 | +0.23 % | 29.1 ms | 4.343 / 4.867 | 1.67 / 1.49 |
| balls | zhu_mean | 13.71 | 7.5668 | +4.26 % | 42.8 ms | 4.343 / 4.848 | 1.74 / 1.56 |
| haseg | lib | 2.77 | 6.1073 | +0.00 % | 2405.6 ms |  |  |
| haseg | cap7_lib | 2.77 | 6.1073 | +0.00 % | 2405.6 ms |  |  |
| haseg | capZ_lib | 3.10 | 6.1073 | +0.00 % | 2406.0 ms |  |  |
| haseg | zhu_mean | 12.54 | 6.1071 | -0.00 % | 2416.7 ms |  |  |

## T8 attribution (BPNV 4.5 g head basket, current combo, t_max 0.2 s)

| model | variant | peak MPa | plateau (0.1-0.2 s) | peak/plateau |
|---|---|---|---|---|
| base | Ma 2020 erosive (default) | 21.62 | 8.42 | 2.57 |
| base | erosive = saint_robert | 7.39 | 4.50 | 1.64 |
| base | erosive = mukunda_paul | 14.95 | 6.70 | 2.23 |
| base | grain cannot ignite (pyrogen gas only) | 2.28 | - | - |
| capZ+area23 | Ma 2020 erosive (default) | 16.13 | 8.42 | 1.92 |
| capZ+area23 | erosive = saint_robert | 5.77 | 4.50 | 1.28 |
| capZ+area23 | erosive = mukunda_paul | 10.71 | 6.70 | 1.60 |
| capZ+area23 | grain cannot ignite (pyrogen gas only) | 1.76 | - | - |
| zhu_lo+area23 | Ma 2020 (r0 = 5.271e-3 m/s/MPa) | 15.87 | 8.54 | 1.86 |
| zhu_hi+area23 | Ma 2020 (r0 = 6.848e-3 m/s/MPa) | 16.56 | 8.50 | 1.95 |
