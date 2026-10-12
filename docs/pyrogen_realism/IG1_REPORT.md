# IG1 report -- why realistic igniter charges over-spike (diagnosis only)

Branch `ig1/diagnosis`, worktree `C:/tmp/sv/ig1`, base `e882134` (main `daf367c` defaults, T11 solver).
Committed `motors/machbusterNew.ric` (Chunc), `pyrogen='bpnv'` (library), roughness 32 um, kappa 0.44,
T_ign 756 K, 106 cells, `ma2020`. **No defaults changed.** Probe-only instrumentation was added to
`srm_1d/simulation.py` (all default-off; see section 9). 1 lb = 453.59 g; 0.01 lb = 4.54 g.

Preflight: `import srm_1d` -> `C:\tmp\sv\ig1\srm_1d\__init__.py` (inside the worktree). Gate B0 with the
instrumentation present: **PASS, all rel errors 0.000e+00, step counts identical** (`results/B0_output.txt`).
The reviewer's table reproduced exactly (head 0.01 lb 21.66 MPa @ 9.6 ms, first ign 0.68 ms; 0.001 lb 9.91 @ 63.7 ms,
first ign 28.5 ms; propellant-cannot-ignite 2.03; `saint_robert` 7.27; plenum 0.01 lb 15.97 @ 13.7 ms).

## 0. Summary

1. **The spike is not an ignition-speed problem.** At 4.5 g, slowing first light from 0.7 ms to 2 ms (convective h x0.2,
   Keller-matched), or to 7 ms (plus a 10 ms pyrogen ramp), or removing the DeMar flux, leaves the peak at 20.5-21.7 MPa.
   Cutting the pyrogen feed at 2, 3 or 5 ms (after the whole grain is lit at 1.8 ms) drops the peak to 9.9-10.1 MPa. The peak
   always coincides with pyrogen burn-out. Spike = post-ignition pyrogen delivery (0.2-0.6 kg/s) into a grain that is already
   fully lit, amplified by the erosive response (`saint_robert` 7.3, `mukunda_paul` 14.9 MPa).
2. **Symptom 1 (ignites with ~10x too little pyrogen) is H1 plus the DeMar flux.** Minimum ignition mass in the sim:
   head_basket 0.36 g, forward_plenum 0.029 g. The sim's convective h is 4.8x Keller's at matched G (was 7x). With
   Keller-matched h the head_basket minimum rises to ~1.4 g (4x). The plenum minimum is set by the time-averaged DeMar flux
   (69.4 cal/cm2/s applied as a step on the head cell); without DeMar it is 0.14-0.45 g; with DeMar off and Keller h 1.4-2.3 g.
3. Once any cell lights, the whole 100-cell grain lights in 1.1-2.4 ms in every variant (the propellant's own gas, 400-600 m/s
   front). That near-simultaneity turns the pyrogen tail into a spike. Whether it is physical needs a flame-spread-in-bore
   reference (pending R4).

## 1. Ignition-threshold sweep (Measurement 1)

t_max 0.2 s. "ignites" = fraction of the 100 grain cells ignited. Full tables incl. fine points and the Keller-h variants:
`results/sweep_tables.md`; raw `results/sweep*.jsonl`.

| topology | lb (g) | ignites | t_first ms | t90 ms | peak MPa | t_peak ms | pyro burn ms |
|---|---|---|---|---|---|---|---|
| head_basket | 0.0007 (0.32) | no | - | - | 0.15 | - | 83 (at ~1 atm) |
| head_basket | 0.0008 (0.36) | yes | 64.9 | 66.1 | 9.88 | 101.6 | 69.8 |
| head_basket | 0.001 (0.45) | yes | 28.5 | 30.4 | 9.91 | 63.7 | 38.8 |
| head_basket | 0.003 (1.36) | yes | 2.60 | 4.07 | 10.43 | 20.2 | 14.6 |
| head_basket | 0.00646 (2.93) | yes | 1.27 | 2.46 | 15.34 | 11.5 | 11.5 |
| head_basket | 0.01 (4.54) | yes | 0.68 | 1.81 | 21.66 | 9.6 | 9.6 |
| head_basket | 0.0139 (6.3) | yes | 0.40 | 1.55 | 28.45 | 8.5 | 8.5 |
| head_basket | 0.03 (13.6) | yes | 0.12 | 1.07 | 54.24 | 6.2 | 6.1 |
| forward_plenum | 6e-5 (0.027) | no | - | - | 0.10 | - | 182 |
| forward_plenum | 6.46e-5 (0.029) | yes | 12.9 | 15.5 | 9.89 | 49.8 | 36.4 |
| forward_plenum | 0.00139 (0.63) | yes | 2.21 | 3.17 | 9.98 | 31.7 | 27.3 |
| forward_plenum | 0.003 (1.36) | yes | 0.86 | 2.11 | 10.47 | 22.0 | 25.6 |
| forward_plenum | 0.01 (4.54) | yes | 0.64 | 1.48 | 15.97 | 13.7 | 20.9 |
| forward_plenum | 0.03 (13.6) | yes | 0.29 | 1.06 | 34.44 | 9.4 | 23.9 |

**Minimum ignition mass in the sim: head_basket 0.36 g (0.0008 lb; 0.0007 fails); forward_plenum 0.029 g (6.46e-5 lb; 6e-5
fails)**, versus the real 4.5 g (12x and 155x). Non-ignition cases burn the pyrogen at ~1 atm for 84 ms (r = a*P^n, 3.2 mm
pellets, constant A_burn) without lighting the grain. Peak vs mass: ~1.2x plateau up to ~1.4 g, then monotonically rising
(1.9x at 2.9 g head, 2.6x at 4.5 g).

Threshold sensitivity (diagnostic probes; `results/sweep_hc0.2.jsonl`, `sweep_nodemar.jsonl`):

| variant | head_basket min | forward_plenum min |
|---|---|---|
| baseline | 0.36 g | 0.029 g |
| convective h x0.2 (Keller-matched) | between 0.91 and 1.36 g | unchanged (0.045 g; DeMar ignites it) |
| DeMar flux off | unchanged (cell 4 lights convectively first) | between 0.14 and 0.45 g |
| h x0.2 and DeMar off | between 0.91 and 1.36 g | between 1.36 and 2.27 g |

## 2. Time histories, head_basket (Measurement 2)

Raw: `results/timehist_0.01lb.txt`, `timehist_0.001lb.txt` (analysis `probes/ig1_timehist.py`; snapshots every 20 us / 50 us).
Basket = cells 0-3 (26 mm) at 4.5 g, cell 0 only at 0.45 g.

| | 0.01 lb (4.54 g) | 0.001 lb (0.45 g) |
|---|---|---|
| peak pyrogen gas rate / mean | 0.62 kg/s at 9.6 ms / 0.45 kg/s | 0.034 kg/s at 38.8 ms / 0.011 kg/s |
| consumption t10/t50/t90/t99 | 2.4 / 6.4 / 9.0 / 9.6 ms | 7.2 / 32.7 / 37.7 / 38.7 ms |
| bore P in basket at ignition | ~1 MPa | 0.10-0.18 MPa (~ambient) |
| gas at first-lit cell (cell 4) | T 1890 K, u 425 m/s, G 524 kg/m2s | T 2530 K, u 150 m/s, G 34 kg/m2s |
| first-lit cell flux (convective / DeMar) | h_c 10.8k W/m2K, q 12.3 MW/m2 / 0 | h_c 1.07k, q 1.89 MW/m2 / 0 |
| T_surf 293 -> 756 K | ~0.6 ms | ~28 ms |
| first / 50% / 90% / last ignition | 0.68 / 1.34 / 1.81 / 1.91 ms | 28.5 / 29.7 / 30.4 / 30.6 ms |
| hot-gas front (T>350 K) | 37 mm @0.11 ms, 313 mm @1 ms, 850 mm @2 ms (grain end) | 12 mm @1 ms, 183 mm @3 ms, 500 mm @10 ms |
| P_head peak / time | 21.64 MPa at 9.6 ms (end of pyrogen) | 9.91 MPa at 63.7 ms |
| aft cell (100) G / erosive fraction at peak | 6664 kg/m2s / 0.73 | 2800 / 0.57 |
| aft cell at plateau | G 2885 / 0.57 | G 2466 / 0.53 |

Mass-flow budget at 4.5 g (kg/s; grain rho_p taken as 1800 in the quick script, true 1700, scale 0.94): grain 0.09/0.70/1.64/2.74/3.35
at 1/2/5/8/9.6 ms, pyrogen 0.13/0.21/0.40/0.55/0.62; grain erosive share of its rate 0.37 -> 0.68. The pyrogen is 10-20% of the mass
flow, but the grain flow climbs 2.2x with P and erosive G. The 0.45 g case ignites in an ambient-pressure bore after ~28 ms of
heating at q ~1.8 MW/m2 (inside Keller's 0.84-6.7 band): a plausible, Keller-regime ignition with a mild 1.19x peak. The 4.5 g
case is in a regime (q 11-12 MW/m2, 1.7x Keller's ceiling) that Keller never validated.

## 3. H1 vs Keller-Baer-Ryan 1966 (Measurement 3)

Keller Eq. 11 (read from the page image of `docs/references/keller1966.pdf` p.3, AIAA J. 4(8)):
h = 1.435e-4 * Tg^0.3 * G^0.905 cal/(cm2 s K) (Tg K, G g/cm2 s), 6% scatter. Script `probes/ig1_keller.py`,
output `results/keller_compare.txt`.

| case, first 8 lit cells | G [kg/m2s] | Tg [K] | h_sim [W/m2K] | h_Keller | ratio | q_sim [MW/m2] | q with Keller h |
|---|---|---|---|---|---|---|---|
| 0.01 lb | 343-632 | 1730-1940 | 7.8-10.9k | 1.4-2.4k | **4.5-5.5, median 4.8** | median 11.5 | ~2.4 |
| 0.001 lb | 27-39 | 2340-3140 | 0.65-1.09k | 0.16-0.21k | 4.0-5.5, median 4.5 | median 1.75 | ~0.35 |

All grain cells at their own peak pre-ignition q: ratio median 3.7 (3.3-5.5) at 4.5 g, 3.4 at 0.45 g.
**Update to REOPENED 9.2: the h ratio is ~4.8x for the first-lit cells (was ~7x at G = 72 g/cm2 s in the old tree); sim q is
11.5 MW/m2 vs Keller's band ceiling 6.7 (1.7x; the old 31 MW/m2 / ~10x was pre-T11).** Caveats: the 0.001 lb case
(G = 2.7-3.9 g/cm2 s) is extrapolated below Keller's data; the sim uses static gas T, Keller stagnation Tg; Keller's fitted G/Tg range
is needed (pending R6). Effect of H1 (probe `hc_scale`, scaling only the ignition Goodman h, not erosion): section 6. It postpones first
light, raises the head_basket threshold ~4x, and does **not** reduce the spike.

## 4. H3 checks (Measurement 4)

(a) Depletion vs a bed-spread-limited estimate. Sim: 4.5 g consumed in 9.6 ms (t90 9.0 ms), peak 0.62 kg/s. First-order estimate
(mine, not sourced): bed spread 5-10 ms (repo memory-level figure with no cited source: pending R5) plus a 3.2 mm pellet's 1.6 mm web at
r(7 MPa) = 52 mm/s to r(20 MPa) = 89 mm/s (sim law) = 18-31 ms, so >= 25-40 ms total. The sim is ~3-4x faster, because all pellet area
burns from t = 0 at constant area and r ~ P^0.5 runs away (the 84 ms burn at 1 atm collapses to ~10 ms as the bore pressurizes). But
stretching it (ramp 5/10/20 ms, area shrink ~(m/m0)^(2/3), pressure cap) does not remove the spike (section 6).

(b) Burn-rate law vs DeMar's closed-vessel data (in-repo `docs/references/STARTING YOUR UP-GOERS_pub (1).pdf` p.7, 12): 5.75 in x 2.75 in ID
vessel (0.560 L); BKNO3-Viton peaks at 105 psi at 67 ms, impetus 5000 psi in3/g, so ~0.72 g. The sim law (a 2e-5, n 0.5, 3.2 mm pellets,
ProPep gas) in that vessel (`probes/demar_closed_vessel.py`): burn-out 49.6 ms and 68 psi (impetus 3245 vs 5000) vs measured 67 ms / 105 psi.
The pellet law is within 26% of DeMar's time-to-peak (DeMar's igniters are dipped heads of unknown area; 1 mm grains would give 15 ms, 0.3 mm
powder 4.7 ms). **No evidence the pyrogen burn law is too fast in the sub-MPa regime.**

(c) High-pressure rate: the sim law gives 52 mm/s at 7 MPa and 89 mm/s at 20 MPa. Whether BKNO3 is flat above ~7 MPa and the MIL-P-46994B 31.75 mm/s
anchor are **pending references R1-R3** (only abstracts / secondary mentions seen, so nothing is concluded from them). The sensitivity is
measured regardless: capping the burn-law pressure at 7 MPa (n = 0 above) changes the 4.5 g peak only 21.7 -> 20.6 MPa; the YAML's literature
candidate (3.5e-5, 0.589; abstract-only, **pending reference R2**) makes it worse (28.4 MPa, burn-out 1.2 ms).

## 5. The real Chunc trace (caveat for "no spike")

`static_fire_data/thomas_chunc_firing.csv` starts at t = 0.010 s already at 8.47 MPa (97% of the 8.65 MPa plateau), sample spacing 7 ms. The
ignition rise is **not recorded** and the time origin is not tied to the electrical fire signal. "No spike" is established only for what
follows the first sample. The sim at 4.5 g is at 15 MPa at 20 ms, 12.0 at 30 ms and 10.7 MPa at 40 ms (1.2x plateau, still decaying). The
BALLSstick 2-in subscale shows a real 1.35-1.5x startup over-pressure (static_fire_data README). Whether the real Chunc igniter transient is a
spike that was cropped off cannot be answered from this file (pending R7: the raw trace from the fire signal).

## 6. Isolating sensitivities (Measurement 5; 0.01 lb head_basket, t_max 0.1-0.15)

Full 49-row table: `results/sensitivity_table.md`. Excerpt:

| factor | first ign ms | t90 ms | peak MPa | t_peak ms |
|---|---|---|---|---|
| baseline | 0.68 | 1.81 | 21.66 | 9.6 |
| (i) convective ignition h x0.5 | 1.04 | 2.08 | 21.47 | 9.8 |
| (i) x0.25 | 1.57 | 2.50 | 21.09 | 10.1 |
| (i) x0.2, DeMar flux off (Keller-matched) | 2.02 | 2.86 | 20.70 | 10.4 |
| (i) x0.125, DeMar off | 3.00 | 4.27 | 19.57 | 11.2 |
| (i) x0.05, DeMar off (far outside Keller) | 9.74 | 12.14 | 12.07 | 15.7 |
| (ii) A_burn ramp 5 ms | 2.44 | 3.78 | 21.70 | 12.2 |
| (ii) ramp 10 ms | 3.42 | 4.95 | 21.98 | 14.5 |
| (ii) ramp 20 ms | 4.72 | 6.63 | 22.01 | 19.3 |
| (ii) h x0.2, DeMar off, ramp 10 ms | 6.84 | 7.96 | 20.49 | 16.4 |
| (iii) erosive `mukunda_paul` | 0.68 | 1.86 | 14.88 | 10.6 |
| (iii) MP + h x0.2, DeMar off | 2.02 | 3.01 | 14.32 | 11.4 |
| erosion off (`saint_robert`) | 0.68 | 1.96 | 7.27 | 12.8 |
| DeMar flux off only | 0.68 | 1.81 | 21.65 | 9.6 |
| pyrogen cut at 2 / 3 / 5 ms | 0.68 | 1.81 | 9.89 / 9.88 / 10.10 | 35 / 35 / 23 |
| burn-law pressure cap 7 MPa | 0.68 | 1.81 | 20.59 | 11.8 |
| A_burn shrinks (m/m0)^(2/3) | 0.69 | 1.82 | 18.02 | 14.7 |
| shrink + pressure cap 7 MPa | 0.69 | 1.82 | 16.69 | 16.7 |
| 5 mm pellets | 1.12 | 2.33 | 20.86 | 15.0 |

Other rows: plenum 0.01 lb baseline 15.97 / `saint_robert` 5.73 / MP 10.76 MPa; head 0.45 g `saint_robert` 4.30, MP 6.80; pyrogen alone
(propellant cannot ignite) 2.03 MPa (head, 4.5 g), 0.18 (head, 0.45 g), 2.48 (plenum, 4.5 g). All probe sizes are diagnostics, not proposals.

## 7. Ranked causes

**Symptom 2 (realistic charges spike)**

1. **Post-ignition pyrogen delivery feeding the erosive/pressure loop with the whole grain already lit (H3 tail + H5).** Evidence: cut tests
   (21.66 -> 9.9), peak always at pyrogen burn-out across the variants, peak growing monotonically with charge above ~1.4 g, `saint_robert` 7.3,
   MP 14.9. At 4.5 g the pyrogen delivers ~3.4 of 4.5 g after the grain is lit. Insensitive to ignition timing (section 6).
2. **Near-simultaneous whole-grain lighting after first light (H4 cascade).** 100 cells in 1.1-2.4 ms in every variant up to Keller-matched h; only
   an h reduction of 20x (outside Keller) stretches first-to-t90 and drops the peak to 12 MPa. Driven by the propellant's own gas and the pyrogen
   burst (front 850 mm in 2 ms). Physical bore flame-spread rate: pending R4.
3. **Erosive amplification at transient G (H5)**: aft G 6664 vs 2885 kg/m2s at plateau, erosive share 0.73 vs 0.57. Known and held (MP is the
   user's decision).
4. **Pyrogen burn law / output rate (H3)**: law consistent with DeMar's sealed-vessel data within 26% (4b); cap/area levers move the peak <= 23%;
   ramp does nothing. Constant-area, all-pellets-from-t=0 burning is ~3-4x faster than a first-order bed estimate, but fixing it alone leaves the
   spike. High-pressure rate: pending R1-R3.
5. **H2 ignition criterion**: not separable here; Keller validates the surface-T gate (REOPENED 9.2). No endothermic / finite-energy criterion
   was implemented, so it was not tested.

**Symptom 1 (ignites with 10x too little)**

1. **H1: convective h 4.8x Keller** sets the head_basket minimum (0.36 g -> ~1.4 g with matched h).
2. **DeMar flux (a time-average of 69.4 cal/cm2/s used as a step on the head/basket cells)** sets the forward_plenum minimum (0.029 g; 0.14-0.45 g
   without it).
3. **Self-propagating cascade**: once one cell lights the propellant gas lights the rest, so the ignition mass is just the first-cell mass.
Combined, a Keller-consistent sim needs ~1.4-2.3 g, within a factor of 2-3 of the real 4.5 g.

## 8. Fix options (physically grounded; none implemented; no tuned factors)

**A. Keller-consistent ignition heat transfer (H1) and a non-step DeMar flux.** Change: use, for the *ignition* Goodman heating of unignited cells,
a correlation with the measured Keller form (h = 1.435e-4 Tg^0.3 G^0.905 cal/cm2 s K, or Peretz/Kuo bore data, pending R4/R6); make the DeMar flux
time-dependent (it is a time-average over the action time, not a step). Erosion keeps its own Ma closure. Expected: min ignition mass 0.36 g -> ~1.4 g
(head), 0.029 -> 0.14-2 g (plenum). It will NOT cut the 4.5 g spike (20.7 vs 21.7). Validate: Chunc ignition delay if R7 is obtained; BALLSstick and
Hasegawa A first-light time and peak; the threshold sweep should land in the amateur 0.5-5 g range. Needs Keller's G/Tg range (R6).

**B. Physical pellet-bed model for head_basket (H3).** Change: pellet web consumption (area shrinks; for cylinders/spheres (m/m0)^(2/3), measured
here 21.7 -> 18.0 MPa), bed flame spread along the basket at a sourced speed (pending R5; replaces the 5/10/20 ms probe values) and a pressure law from
BKNO3 data (pending R1-R3; flat above ~7 MPa only if confirmed: 18.0 -> 16.7). Parameters from literature (pellet diameter/density already in the
YAML). Expected: burn-out stretched 9.6 -> ~25-30 ms, spike down ~20-25%, not eliminated. Validate: DeMar's sealed-vessel times (probe included);
for the motor, the Chunc basket drawing and igniter-only data (R7).

**C. Flame-spread-limited whole-grain ignition (H2/H4).** The cascade (100 cells in 2 ms, 400-600 m/s front) is the amplifier. Change: add a
literature-supported finite ignition delay (endothermic AP decomposition / integrated-energy criterion; Linan-Williams plateau estimated at ~35% in
REOPENED 9.2) or a physical bore flame-spread speed (Peretz et al., pending R4). Tested before only as the v0.7.4 flame-front opt-in on the old solver;
here h x0.05 shows what a 5-10x slower spread buys (peak 12 MPa). Validate on Chunc, BALLSstick (must keep its real 1.35-1.5x over-pressure) and
Hasegawa A (must stay at 6.1 MPa).

**D. Erosive closure at transient G (H5): the registry option `mukunda_paul`**, 14.9 MPa (1.8x) at 4.5 g alone. The user's decision (held). Combined with
A-C it is the only lever in this report that moves the peak by more than 25% alone. Validate on all three motors with the existing MP gate.

Suggested order: A (cheap, fixes symptom 1, needs R6) -> B -> C, then re-evaluate D. Do not expect A or B alone to remove the 4.5 g spike. Obtain R7 first:
the committed trace cannot rule a real startup spike in or out.

## 9. Files, instrumentation, open items

Instrumentation (`srm_1d/simulation.py`, probe-only, default-off, B0 PASS bit-identical): kernel arg `probe_hc_scale` + `probe_hq` (records bare h_c and
convective q to snapshot channels 19/20, new snapshot keys `probe_h_c`, `probe_q_conv`); `run_simulation` kwargs `probe_hc_scale`, `probe_ramp_tau`,
`probe_cut_t`, `probe_pcap`, `probe_area23`. They add signature args to `_goodman_ignition_sources_and_mass` and `_run_time_loop`, so 2 tests in
`tests/test_simulation_phase3.py` (direct kernel calls) fail on this branch (19 of 21 pass in that file plus `test_uncontained_pyrogen.py`); this branch is
diagnostic and not intended to be merged as-is.

Probes: `probes/ig1_case.py` (runner, resumable JSONL), `ig1_timehist.py`, `ig1_keller.py`, `demar_closed_vessel.py`, `make_tables.py`,
`make_sens_table.py`. Warm 0.2 s runs take ~4 s each (first run compiles); ~110 runs total. Raw time-history npz (27 MB) are in `C:/tmp/sv_ig1/npz`, not committed.

Open items: paywalled references R1-R7 (`PAYWALLED_REFERENCES.md`); Keller G/Tg fit range; endothermic ignition criterion not tested; the 0.001 lb case lies
below Keller's G range; plenum time histories not taken (only head_basket, per spec); aft_basket not examined.
