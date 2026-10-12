# OB1 report -- a literature-grounded pyrogen burn model (head basket first)

Branch `ob1/pyrogen-bed`, worktree `C:/tmp/sv/ob1`, base `2c08c46`. Committed `motors/machbusterNew.ric` (Chunc), `bpnv`/`mtv`
library pyrogens, roughness 32 um, kappa 0.44, T_ign 756 K, 106 cells, `ma2020`, T11 solver. "Forward combo" =
`ignition_heat_transfer='peretz'` + pyrogen `heat_delivery_mode='none'` (IG2 recommendation, opt-in here);
"current" = shipped defaults. Real igniter (user, R7): 5-6 g of ~3 mm cubes of BKNO3 or MTV in a head basket, no startup
spike in the recorded trace (which starts at 10 ms, so the rise itself is unrecorded). All new behaviour is opt-in; **defaults
are bit-identical** (B0 PASS, section 8).

Preflight: `import srm_1d` -> `C:\tmp\sv\ob1\srm_1d\__init__.py` (inside the worktree).

## 0. Summary

1. **What the three sources say (section 1).** Zhu 2008 (closed bomb, B/KNO3, 2.8-50 MPa): the rate is **not** a power law;
   `r = r0 P` up to 7 MPa and **flat** `r = 7 r0` from 7 to 50 MPa, `r0 = 5.27e-3 .. 6.85e-3 m/s/MPa` (plateau 36.9-47.9 mm/s,
   mean 42.5). MIL-P-46994B: 1.25-1.75 in/s (31.8-44.5 mm/s) closed-bomb *average* with the pressure unstated; standard
   BKNO3 is KNO3 70.7 / B 23.7 / polyester binder 5.6 %. Nakano & Kakami 2016 (the "Mizushima 2016" of the repo YAML):
   `r = 71.1 P^0.589` mm/s is valid **5-200 kPa only** (28:70:2 B:KNO3:binder, 3.2 mm dia), not to 7 MPa. The sim's `2e-5 P^0.5`
   is 1.1-1.4x Zhu at 7 MPa (1.24x the mean), 2.1x the mean plateau at 20 MPa, 2.0x Zhu at 2.8 MPa, and 0.34x Nakano at 0.1 MPa (table T1).
2. **Audit (section 2).** The head basket burns every pellet at **constant area** from t = 0 until the pooled charge is gone. For a cube
   that exhausts the mass after a regression of s/6, while the real web is s/2: **3x too short at any given rate** (this, not the rate
   law, is why 4.5 g burns in 9 ms). The forward plenum already has the exact `(m/m0)^(2/3)` area; the basket does not. The YAML
   citation chain for the BKNO3 law is wrong (author, range, SI conversion 3.50e-5 vs the correct 2.08e-5).
3. **Ranking (section 3, T2).** On the 4.5 g BPNV row (peak 21.6 MPa, 2.57x plateau): (b) shrinking area -16 % -> 18.1; (a) rate
   plateau alone -5..-7 %; (a)+(b) **-25 % -> 16.1-16.3 MPa (1.9x plateau)**; (c) a sequential bed lighting at 3/10/30 m/s changes
   nothing (-16.0 / -17.5 / -16.6 % vs -16.3 % for area alone). With the grounded options the pyrogen burn lasts 36-53 ms (was 9) and
   the peak still reaches 1.9x. **No grounded pyrogen option removes the over-spike.**
4. **Why (T8).** In the best grounded model the pyrogen gas alone gives 1.8 MPa; with the erosive term off the same run peaks at 5.8 MPa
   (1.28x plateau); the 16.1 MPa peak is the Ma quasi-steady erosive response (Mukunda-Paul 10.7). At the Zhu-law peak the pyrogen supplies
   0.09 of 2.1 kg/s. The remaining spike is not a pyrogen-model problem.
5. **Low-pressure branch is unsourced for the sim's regime.** The printed Zhu law (`r = r0 P`, valid >= 2.8 MPa) gives 0.6 mm/s at 1 atm:
   in the head basket first light moves from 0.66 / 0.40 ms to 5.8 / 3.8 ms (4.5 / 6 g), and in the **forward combo the charge does not light at 4.5 or 6 g** (threshold 10.6 g).
   Nakano measured 18 mm/s at 0.1 MPa; the 0.2-2.8 MPa gap between the two papers is unbridged (their compositions differ), so I did
   not stitch a three-branch law. A neutral variant (library law below the cap, capped at Zhu's plateau *rate*, `capZ`) is shipped as a
   configuration of the same opt-in field. The peak is insensitive to which low branch is used (16.1 vs 16.3 MPa).
6. **Ignition threshold (T5).** Head basket, forward combo: 1.38 g (current pyrogen model; IG2 found 1.43 g, same bisection tolerance class) -> **2.7 g** with the area law (and the lit
   boundary is ragged: 2.93 g lights, 3.6 g does too, but only marginally, 99/100 cells at 6-11 ms). Real 4.5-6 g is 1.7-2.2x above it.
7. **Canonical motors (T7).** With the usual Sutton forward-plenum igniters the opt-in plateau changes Chunc / BALLSstick /
   Hasegawa A peaks by 0.0-0.8 % (`capZ`) or +2.8 / +4.3 / 0.0 % (full Zhu law, which delays first light by 10-12 ms). Defaults: unchanged.
8. **Recommendation (section 9):** adopt the shrinking area law (grounded, exact for cubes, already what the plenum does) as the
   head/aft-basket default in a follow-up (IG3-class decision, the user's call) and keep the rate plateau opt-in; do **not** expect either to
   fix the Chunc spike. The remaining lever is the erosive response at transient G (held MP decision) and the whole-grain ignition cascade
   (candidate C of IG1), not the pyrogen. Needs: R5 (bed flame spread), R13-R16 (section 10).

## 1. Part 1 -- what each source says

All three PDFs have unusable text layers for the content (R1 Chinese Founder CID fonts; R2 Japanese, only the English abstract extracts;
R3 an OCR layer too poor to use). Pages were rendered at 100-110 dpi and read as images. Nothing below comes from an abstract or a
secondary mention; where a source is silent I say so.

### 1.1 Zhu, He, Jiang, Zhu, "Investigation on burning rate of B/KNO3 pellet in closed bomb test", Chinese J. Energetic Materials 16(2):204-206 (2008) -- `docs/references/zhu2008.pdf` (R1)

- Method (p.1-2, Fig. 1): 50 mL closed bomb, four runs, loading density ~0.2 g/cm3, **rate from the differential pressure method**
  `p = k W phi R T / (...)` (Eq. 1-2) with the geometric burning law (form function `phi = chi Z (1 + lambda Z)`).
  Pellets are right cylinders (Table 1): No. 1/2: rho 1.750/1.751 g/cm3, D 10.06 mm, H 9.99 mm (charge 9.73 g); No. 3/4: rho
  1.770/1.768, D 6.30, H 6.01/6.02 mm (9.94 g). Composition is **not stated** (the paper cites MIL-STD-1901A; Zhu 2005 and Jiang 2002, refs [2],[3], are the
  likely sources, R15). Peak pressure ~60 MPa (Fig. 2).
- Result (Fig. 3, Section 3.2, Eq. 3, Table 2): linear rate scattered (random) below 2.8 MPa; "between 2.8 and 7 MPa the rate is
  essentially linear in pressure; between 7 and 50 MPa essentially independent of pressure" (the paper attributes the plateau to solid
  products covering the unburnt particles). Fit: `r = r0 p (2.8 <= p <= 7 MPa)`, `r = 7 r0 (7 < p <= 50 MPa)`, with
  `r0 [m/(s MPa)] = 6.848e-3, 6.787e-3, 5.389e-3, 5.271e-3` (runs 1-4; mean 6.074e-3). Plateau rate 47.9, 47.5, 37.7, 36.9 mm/s
  (Fig. 3 shows 0.04-0.05 m/s). Lower density -> faster (runs 3/4 are the denser pellets and slower).
- Conclusion (abstract, Section 4): the rate "does not follow a simple exponential law"; pressure exponent 1 below 7 MPa and 0 above.
  Pellet size affects the product *generation rate* (area), not the linear rate.
- Also: solid fraction of the products at ambient after cooling is "over 85 %"; the authors state that during combustion the products
  are probably still mostly gas, so this does not contradict `gas_mass_fraction`.
- **Limits:** composition unstated (BPNV is B/KNO3/Viton 25/60/15), one bomb volume and loading, pressure range 2.8-50 MPa (nothing below
  2.8 MPa; the printed law through the origin is an extrapolation there), no pellet-size dependence of r.

### 1.2 Nakano & Kakami, "Burning rate measurements of 3.2-mm-dia. BKNO3 pellets in low pressures", J. Japan Soc. Aeronaut. Space Sci. 64(2):139-143 (2016) -- `docs/references/nakano2016.pdf` (R2)

- Method (Section 2): end-face-lit 3.2 mm diameter cylindrical pellets (lengths 2.0 and 4.8 mm, also glued stacks 9.6 and 10 mm),
  pellets glued to a steel rod and to each other with NC lacquer (joint effect on the rate <= 5 % per joint, +/- 12-25 % scatter for the 9.6 mm stack), rate from burn time by high-speed camera; electric igniter in a 0.534 m3 vacuum chamber (Ar 99.99 %, 0.2 kPa up) and a laser
  igniter in a 732 cm3 pressure vessel (N2) for >= atmospheric. B : KNO3 : binder = **28 : 70 : 2** by mass (B ~1 um), supplier's binder
  **undisclosed** (the repo YAML says "NC"; the paper says the binder identity is not disclosed). Pellet temperature 7.2-21.7 C.
- Result (Fig. 4, Section 3, abstract): below 5 kPa the rate is pressure independent, **3.0 +/- 0.5 mm/s**; above 5 kPa it rises, and
  fitting Vieille `r = a p^n` ([mm/s], [MPa]) gives **a = 71.1, n = 0.589 for 5-200 kPa** (the Section 4 summary: "between 5 kPa and
  200 kPa the pressure exponent is 0.589"). Maximum pressure tested: 200 kPa. Fig. 4 also overlays Yano's 1980s data (B mass fraction 0.1-0.8).
- **Limits:** top of the tested range is **0.2 MPa**; at 7 MPa the fit would give 224 mm/s (5x Zhu). The repo YAML says "valid 5 kPa - 7 MPa"
  (wrong, 35x the tested top) and "SI a = 3.50e-5" (the conversion is 71.1e-3 * (1e-6)^0.589 = **2.08e-5**, so the YAML was 1.68x too high).
  Pellet lengths 2.0-10 mm: the "3.2 x 3.2 mm L/D = 1 industry pellet" in the YAML is not from this paper. The authors' own caveat: the rate is
  that of a 3.2 mm cylinder, not a 1-D rate.

### 1.3 MIL-P-46994B, pellets/granules, boron/potassium nitrate -- `docs/references/MIL-P-46994B.pdf` (R3)

- Table I (p.4): standard formulation KNO3 **70.7 +/- 2.0 %**, amorphous boron **23.7 +/- 2.0 %**, binder **5.6 +/- 0.5 %** (Table II: 98 % polyester resin,
  1.5 % MEK peroxide, 0.5 % cobalt naphthenate); Type IV: 72.2 / 22.2 / 5.6. KNO3 15 um max, boron 1.5 um max, 90-92 % pure.
- Table IV (p.7 of the scan): burning rate **1.25-1.75 in/s (31.75-44.45 mm/s)** for Types I, II/III, VI; N/A for IV, V, VII; heat of reaction >= 1500 cal/g
  (1600 Type IV); Type IV density >= 1.60 g/cm3, time to maximum pressure 45 +/- 15 ms.
- Para. 4.6.2 (p.15): the burning rate is a **closed-bomb** measure (10 g sample, MIL-STD-286 method 801.1, six-inch nylon bag, M1A1 squib)
  `R = h/(2 t)`, `h` = average pellet height, `t` = time from first pressure indication to the stabilised maximum. **No pressure is stated**: it is a
  pressure-history *average*, not a rate at a pressure (R13 asks for method 801.1). It numerically brackets Zhu's plateau (36.9-47.9 mm/s) but I use it only as "consistent".
- BPNV (25/60/15 with Viton) is **not** this composition (a binder-heavy 15 % Viton, 25 % B); BKNO3-NC in DeMar's survey is closer.

### 1.4 Ma et al., "A model for igniter mass flow rate history evaluation for solid rocket motors", Int. J. Aerosp. Eng. 2019, 2593602 (the "2020" file)

An **inverse zero-dimensional model**: given a measured igniter-only pressure history it back-computes the gas mass flow rate (with species
properties and heat loss). It has no burn law, no form function and no flame spread: the only charge is a *confined* 1.35 g black-powder
igniter (KNO3 75 / S 15 / C 10, 9 orifices, aluminium case, CEA condensed fraction 43.7 %, gas 0.76 g, plug rupture 0.8 MPa), pressure peak 3.16 MPa at ~5 ms,
decaying to ambient by 26.5 ms, gas-phase mass error 2.6-4.0 %, validated by a 2-D Fluent run at 3.1 % peak error. It does not address pellet baskets or progressive ignition
("different forms of igniter charge ... make it difficult to predict the igniter MFR by estimating the generation rate; neither the burning rate nor the burning area is easy to determine" is
its stated motivation). Useful only as a data point: a small confined igniter delivers its gas over ~20-26 ms. **It supplies nothing for options (b) or (c).**

### 1.5 Comparison with the library

Units: the sim law is `r = a P^n` with `a` in m/s/Pa^n and `P` in Pa (`bpnv.yaml`: a 2e-5, n 0.5; `mtv.yaml`: a 4.4e-5, n 0.35). The Zhu law is
`a = r0/1e6 [m/s/Pa], n = 1, P_plateau = 7e6 Pa`.

| P [MPa] | sim law 2e-5 P^0.5 | Nakano 71.1 P^0.589 (valid 0.005-0.2 MPa) | YAML "candidate" 3.5e-5 P^0.589 | Zhu mean r0 (valid 2.8-50 MPa) | MIL-P-46994B (closed-bomb average) |
|---|---|---|---|---|---|
| 0.1 | 6.3 | 18.3 | 30.8 | (below range) |  |
| 0.2 | 8.9 | 27.6 | 46.4 | (below range) |  |
| 1 | 20.0 | (71, out of range) | 119.7 | (below range) |  |
| 2.8 | 33.5 | (130, out of range) | 219.5 | 17.0 | 31.8-44.5 |
| 7 | 52.9 | (224, out of range) | 376.6 | 42.5 | 31.8-44.5 |
| 20 | 89.4 | (415, out of range) | 698.8 | 42.5 | 31.8-44.5 |

The sim law is closest to Zhu near 7 MPa (1.1-1.4x). Below it (2.8 MPa) it is 2.0x fast (vs the mean r0), above it it keeps rising to 89 mm/s at 20 MPa while the bomb data are flat at 37-48
(1.9-2.4x). At ambient it is 0.34x Nakano. BPNV is not standard BKNO3 (Viton vs polyester/NC binder, 25 vs 23.7-28 % boron); no source covers it.
MTV: only Kubota 1987 (0.3-7 MPa, n ~0.35, in repo); nothing above 7 MPa (R16), so no plateau law is offered for MTV.

## 2. Part 2 -- audit of the sim's pyrogen model

| item | forward_plenum | head_basket / aft_basket |
|---|---|---|
| burn law | `mdot = rho A a max(P_ig,0)^n`, SI (`_plenum_rates`) | `mdot[i] = rho a P[i]^n A_cell` (`_compute_uncontained_pyrogen_mdot`) |
| rate pressure | plenum pressure `P_ig` (solved ODE, orifice-vented) | **local bore pressure of the host cell** (no plenum, no orifice) |
| area evolution | `burn_law '0d'`: `A0 (m/m0)^(2/3)` (exact for cubes, spheres, L/D = 1); `'end_burning'`: constant | **constant** `A_burn_initial / n_cells` until the *pooled* remaining mass hits zero (the `burn_law` field is ignored) |
| depletion | one pool | **one pool shared by all cells**: a cell cannot run out on its own; the last step is scaled uniformly |
| lighting | all burn from t = 0 | **all pellets lit at t = 0** (no bed flame spread) |
| initial area | `A_total = 6 m/(rho d)` (L/D <= 1) or `m (4 lam + 2)/(rho lam d)` | same `A_burn_initial` split evenly over `ceil(L_cart/dx)` cells |
| condensed fraction | solid depletes at the full rate; gas `mdot * gas_mass_fraction` enters the plenum, enthalpy `Cp_pyro T_flame` per kg of **gas** | same: `mdot_arr *= gas_mass_fraction`; mass `(1-gmf)` (21.9 % for BPNV) and its enthalpy (B liquid at 2719 K) are **discarded from the system** (not injected, not radiated, not deposited on the grain) |

Geometric correctness for cubes: all faces regress at `r`, `A = 6 s^2`, `m = rho s^3`, so `A ~ m^(2/3)` exactly (a cylinder with L = D, or a
sphere, gives the same law). The constant-area model consumes the mass `m0 = rho A0 s0/6` after a regression of `s0/6` while the physical exhaustion
is at `s0/2`: **at any given rate it is 3x too fast** for 3 mm cubes. This is the dominant defect I found, and it is the same
`(m/m0)^(2/3)` the plenum already applies. For L/D != 1 the law is an approximation. A second inconsistency: with one pool, the cells at higher pressure burn at a higher rate
but keep their area until the *total* is gone, so cell-to-cell depletion differences are not represented; the per-cell option below fixes that.

Other audit findings:
- The YAML provenance chain for the BKNO3 law is wrong (section 1.2): author "Mizushima" -> Nakano & Kakami; "valid 5 kPa-7 MPa" -> 5-200 kPa; "a = 3.50e-5" -> 2.08e-5; "NC binder" -> undisclosed;
  "3.2 mm L/D = 1" -> not in the paper. IG1's "(3.5e-5, 0.589) makes it worse (28.4 MPa)" probe therefore used a law 1.68x too high *and* far out of range; it should not be read as evidence about the lit law.
  I corrected the comment in `pyrogens/bpnv.yaml` (no value change) and recorded it in DEVNOTES; `docs/v0_7_3/references/pyrogen_physical_params_lit_dive.md` carries the same error (not edited; historical).
- All pellets lit at t = 0: the 84 ms ambient burn in IG1 shows the area is lit from the start; with the real bed, Sutton (RPE 15.3, Fig. 15-14: pellets "ignite in stages") only qualitatively supports staged ignition; no speed is in the repo.
- The first-ms burn at ~1 atm is set by the library law (6 mm/s at 0.1 MPa), 0.35x Nakano's measured rate for a 28 % B pellet; this controls first-light time and the low-charge threshold, not the peak (section 3).

## 3. Part 3 -- options, probes and ranking

Candidates (a)-(d) from the task, each confirmed or rejected:

| | option | grounded? | built as | verdict |
|---|---|---|---|---|
| (a) | rate plateau above 7 MPa | **yes**, Zhu 2008 (B/KNO3, 7-50 MPa, 4 runs) -- but composition differs from BPNV; none for MTV | `Pyrogen.P_plateau` (rate = `a min(P, P_plateau)^n`) | opt-in; small effect |
| (a2) | Zhu piecewise through the origin below 7 MPa | the printed law is valid for >= 2.8 MPa only; below it unsourced and in conflict with Nakano | (same field, `a = r0/1e6, n = 1`) | measured, rejected as a default: no ambient burn |
| (b) | geometric form function `A ~ m^(2/3)` | **yes**, first principles, exact for cubes / spheres / L/D = 1 (and what the plenum already does) | `PyrogenChamber.basket_area_law='shrinking'` (per-cell mass) | opt-in; **largest single effect** |
| (c) | progressive bed ignition / flame spread | **no source in the repo** (Ma 2019 has none; Peretz is bore flame spread; Sutton only qualitative) | `PyrogenChamber.bed_flame_speed_m_s`, DIAGNOSTIC ONLY (default 0) | no effect on the peak; R5 for a sourced value |
| (d) | anything else the sources show | Nakano's ambient rate (3x the sim at 0.1 MPa); the unbridged 0.2-2.8 MPa gap; MIL average rate; Zhu's solids remark | none | listed; no unsourced stitching |

Ranking probes (BPNV 4.5 g, 3 mm cubes, head basket; raw `results/ob1_rank.jsonl`; `capZ` = the library `a, n` capped at Zhu's mean plateau *rate* 42.5 mm/s, i.e. `P_plateau = (0.0425/2e-5)^2 = 4.5 MPa`;
`cap7` = IG1's 7 MPa cap):

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

Ranking by effect on the 4.5 g peak (current combo): (a)+(b) -25 % > (b) -16 % > (a) -5..-7 % > (c) 0 (diagnostic). The gap between the Zhu law and `capZ` is nil (16.27 vs 16.13 MPa), and the measured
spread of Zhu's own r0 (5.27e-3 .. 6.85e-3) moves it only 15.9-16.6 MPa (T8), so the answer does not depend on the unsourced low-pressure branch.

## 4. Part 4 -- implementation (opt-in; defaults unchanged)

| API | default | effect |
|---|---|---|
| `Pyrogen.P_plateau` [Pa] (YAML / `.ric` igniter block / GUI metadata) | 0.0 (off) | rate evaluated at `min(P, P_plateau)` in the plenum and the baskets; Zhu law = `a = r0/1e6, n = 1, P_plateau = 7e6` |
| `PyrogenChamber.basket_area_law` (`run_from_ric(basket_area_law=)`, `data.igniter.basket_area_law`, GUI metadata) | `'constant'` | `'shrinking'`: per-cell mass; `A_cell = A_cell0 (m_cell/m_cell0)^(2/3)`; each cell depletes on its own host pressure; last step exact |
| `PyrogenChamber.bed_flame_speed_m_s` (`run_from_ric(bed_flame_speed_m_s=)` only) | 0.0 | DIAGNOSTIC sequential lighting from the igniter end, needs `'shrinking'`; **no sourced value** |

Files: `srm_1d/propellant.py`, `igniter_plenum.py`, `simulation.py` (`_compute_uncontained_pyrogen_mdot_cells`, kernel branch keyed on `chamber_params_arr[5:7]`, so no kernel signature
change), `openmotor_adapter.py`, `srm1d_plugin.py` (schema metadata only), `pyrogens/bpnv.yaml` (comment correction only). Internal array-length and optional-argument changes are logged in DEVNOTES. No library
YAML for the Zhu law is shipped (composition mismatch; no MTV analogue); a user builds it with `dataclasses.replace(load_pyrogen('bpnv'), a=6.07375e-9, n=1.0, P_plateau=7e6)`.

Tests (`tests/test_ob1_pyrogen_bed.py`, 11): defaults off (fields, array entries, library); plateau off = bit-identical and on = hand-computed `rho A a min(P,P_c)^n`; the Zhu law
(r = r0 P at 2.8 MPa, 42.5 mm/s plateau); plenum `_plenum_rates` at P_ig 30 MPa; `A ~ m^(2/3)` at m/m0 = 1, 0.5, 0.125; per-cell depletion order and exact mass conservation including the gas
fraction (gas delivered = gmf x solid, 1e-12); a single cube's burn time equals the analytic web time `(s/2)/r` to 0.2 %; bed lighting delay; a Chunc run: default = explicit `'constant'` (identical P history and step count)
and `'shrinking'` burns exactly 4.5 g with monotone depletion and a longer burn.

## 5. Part 5 -- measurements

### 5.1 Realistic rows (Chunc head basket, 3 mm cubes L/D 1, 4.5 / 6 g, BPNV / MTV, t_max 0.2 s)

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

t_max 1.0 s (plateau = median 0.4-0.9 s; experiment peak 8.89, plateau 8.66, **ratio 1.027**):

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

MTV has no sourced plateau law, so only the area law is shown. Plateau 8.15-8.19 in every row (exp 8.66): the pyrogen model does not touch the plateau.

### 5.2 Head-basket ignition threshold with the new model (as IG2 4.1)

| option | current combo [g] | forward combo [g] |
|---|---|---|
| base | 0.332 | 1.38 |
| area23 | 0.52 | 2.71 |
| capZ+area23 | 0.52 | 2.71 |
| zhu+area23 | 1.62 | 10.6 |

Ignition boundary is ragged near threshold (the lit criterion is 99 of 100 cells inside 0.2 s): forward combo + area law lights at 2.71 and 2.93 g but not at 2.27, and at threshold first light is 6-11 ms (marginal).
Peak vs charge above the threshold:

| combo | option | 0.6 g | 0.9 g | 1.36 g | 1.8 g | 2.27 g | 2.93 g | 3.6 g | 4.54 g | 6.3 g | 9.0 g |
|---|---|---|---|---|---|---|---|---|---|---|---|
| current | base | 9.9 (1.15) | 10.0 (1.17) | 10.3 (1.21) | 11.3 (1.33) | 12.8 (1.51) | 15.3 (1.81) | 17.9 (2.12) | 21.8 (2.59) | 28.6 (3.42) | 38.5 (4.63) |
| current | area23 | 9.9 (1.15) | 10.0 (1.17) | 10.5 (1.23) | 11.3 (1.33) | 12.2 (1.44) | 13.7 (1.62) | 15.5 (1.83) | 18.2 (2.16) | 23.2 (2.77) | 30.6 (3.68) |
| current | capZ+area23 | 10.0 (1.15) | 10.0 (1.17) | 10.5 (1.23) | 11.2 (1.31) | 11.9 (1.40) | 13.0 (1.54) | 14.3 (1.69) | 16.2 (1.93) | 19.6 (2.34) | 24.3 (2.92) |
| forward | base | NL | NL | NL | 9.9 (1.14) | 10.2 (1.20) | 13.5 (1.59) | 16.4 (1.94) | 20.3 (2.41) | 27.5 (3.28) | 37.7 (4.52) |
| forward | area23 | NL | NL | NL | NL | NL | 12.1 (1.46) | 14.4 (1.70) | 17.1 (2.03) | 22.3 (2.66) | 29.9 (3.59) |
| forward | capZ+area23 | NL | NL | NL | NL | NL | 11.7 (1.41) | 13.2 (1.59) | 15.4 (1.82) | 19.0 (2.26) | 23.8 (2.86) |

The spike exists for every charge above ~1.4 g in every model: peak/plateau 1.4 at 2.27 g -> 1.9 at 4.54 g (best model, current combo). The real 4.5-6 g charge sits 1.7-2.2x above the
forward-combo area-law threshold (2.7 g) but, in this model, in the part of the curve where the sim spikes 1.8-2.3x.

### 5.3 Canonical motors with their usual igniters (library BPNV, Sutton-mass forward plenum)

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

The area law does not apply (plenum). Defaults: B0 identical (Chunc 10.0827, BALLSstick 7.2577, Hasegawa A 6.1073). The opt-in plateau moves peaks by <= +0.8 % (`capZ`), the full Zhu law by +2.8 %
(Chunc), +4.3 % (BALLSstick, pk/trough 1.67 -> 1.74, real ~1.50) and 0.0 % (Hasegawa A); it delays first light by 10-12 ms because the Zhu law is slow at ambient pressure (`capZ`: 0.3-0.5 ms).

### 5.4 Attribution of the residual spike

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

With the best grounded pyrogen model the gas of the charge alone reaches 1.8 MPa; the 16.1 MPa peak is mostly the Ma erosive response (5.8 MPa without it; 10.7 MPa with Mukunda-Paul). A time history of the
Zhu + area-law run shows the pyrogen supplying 0.09 kg/s of the 2.1 kg/s nozzle flow at its 16.3 MPa peak (33.9 ms), after the 0.15-0.24 kg/s delivery at 18-27 ms has primed the erosive/pressure feedback.

## 6. What the probes do and do not establish

- Established: the constant-area model is 3x too fast per unit regression; correcting it (and capping the rate at the bomb plateau) stretches the 4.5 g burn from 9 ms to 36-53 ms (to 116-125 ms for MTV)
  and lowers the peak 25 % (BPNV) / 17-21 % (MTV, area law only). The remaining ~1.9x spike is not caused by how fast or how long the pyrogen burns.
- Not established (pending): whether BPNV follows Zhu's B/KNO3 plateau (composition, R15); the 0.2-2.8 MPa rate; any flame-spread speed (R5); the real Chunc igniter's pressure-time history (R7); MTV above 7 MPa (R16).
  These conclusions are marked pending.
- The Zhu law applied literally is not usable as a head-basket default: it does not light the forward combo at the real charge.

## 7. Files

Probes: `probes/ob1_case.py` (runner), `ob1_sweep.py` (rank / realistic / long / extra / mass / thresh), `ob1_valid.py`, `ob1_hist.py`, `ob1_table.py`, `ob1_tables.py`. Raw: `results/ob1_*.jsonl`, `results/ob1_tables.md`,
`results/ob1_B0_output.txt`. Run from the worktree root with `PYTHONPATH=$PWD`.

## 8. Gates

- **B0 (default tolerances): PASS from a clean compile** (`results/ob1_B0_output.txt`, repeated in `C:/tmp` scratch): Chunc peak rel 3.8e-7, t_peak 6.7e-9, plateau 5.9e-9, steps identical; BALLSstick 7e-16 / 3e-9 / 4e-8, steps identical;
  Hasegawa A 5.5e-10, steps 3.7e-5. The base commit `2c08c46` itself gives 2.4e-7 from a clean compile and *exactly* 0 on a second (warm) run. My branch's clean compile is 3.8e-7 (kernel source changed -> different fastmath codegen).
  **OPEN ITEM:** one warm run on my branch taken AFTER pytest/probe runs had re-cached other specialisations FAILED the Chunc t_peak check (peak 1.5e-6, t_peak 5.97e-5 > 1e-5; the flat peak top makes argmax noise-sensitive). I did not finish establishing whether that is compile-state noise
  or a real default perturbation; the evidence for noise: base is exact only in a pure-B0 warm state, the same-source clean compile passes twice (3.8e-7, repeatable), and default-path code is arithmetically unchanged (the new branches are off; `test_run_default_identity` is bit-identical).
  To close it: run B0 clean then B0 warm again on the branch (a warm run after a clean B0, with no other jobs between) and require exact reproduction like the base.
- **Full suite:** first run on the branch: 750 passed, 1 failed (`test_ric_igniter::test_pyrogen_block_round_trip_every_field`: `_pyrogen_from_mapping` dropped `P_plateau`). Fixed; the 3 igniter test files + the OB1 tests (12) pass (63+). **The full suite has not been re-run after the fix** (pending).


## 9. Recommendation (the user decides defaults)

1. Treat the constant-area head/aft-basket burn as a defect, not a modelling choice: the shrinking `(m/m0)^(2/3)` area (and per-cell depletion) is first-principles for the stated pellet shape and is what the forward plenum already does.
   Candidate for the IG3 default flip together with Peretz h and DeMar off. Cost: a higher ignition threshold (0.33 -> 0.52 g current, 1.38 -> 2.7 g forward), which moves toward the real 4.5-6 g.
2. Keep `P_plateau` opt-in (BPNV composition mismatch, unbridged low-pressure range, no MTV source, forward combo does not light with the literal Zhu law). If a plateau is wanted, the `capZ` form (library `a, n`, cap at the bomb plateau rate) is the conservative reading.
3. Do not expect a pyrogen change to remove the 1.8-2.3x peak at 4.5-6 g: that residual is the erosive response at transient mass flux plus the whole-grain ignition cascade (IG1 candidates C/D), outside this task.
4. Fix the BKNO3 provenance in any docs that still quote "Mizushima / 3.5e-5 / 7 MPa".

## 10. New paywalled / missing references (added to `PAYWALLED_REFERENCES.md`)

R5 pellet-bed flame spread (candidate citations, all unverified DOIs); R13 MIL-STD-286 method 801.1; R14 Mizushima et al. 2024 (J-STAGE, open access, not read); R15 Zhu 2005 / Jiang 2002 (composition behind Zhu 2008);
R16 MTV rate above 7 MPa.
