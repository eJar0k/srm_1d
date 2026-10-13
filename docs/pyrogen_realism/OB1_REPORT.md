# OB1 report -- a literature-grounded pyrogen burn model (head basket first)

Branch `ob1/pyrogen-bed`, rebased onto `openmotor-frontend` (IG3 stage 1 + the B0 t_peak tolerance change `de69787`). Committed `motors/machbusterNew.ric` (Chunc),
library `bpnv`/`mtv` pyrogens, roughness 32 um, kappa 0.44, T_ign 756 K, 106 cells, `ma2020`, T11 solver. **Main's defaults are the "current" combo: Peretz ignition h
and the library DeMar flux ("DeMar on"); "DeMar off" = pyrogen `heat_delivery_mode='none'`.** Real igniter (user, R7): 5-6 g of ~3 mm cubes of BKNO3 or MTV in a head basket,
no startup spike in the recorded trace (which starts at 10 ms, so the rise itself is unrecorded).

**User decisions applied (2026-10-12):** (1) the shrinking pellet area is the head-basket default; `P_plateau` stays opt-in. (2) The library pellet pyrogens use the
MIL-P-46994B Type I-A cylindrical pellet, **0.125 +/- 0.002 in dia x 0.188 +/- 0.004 in = 3.175 x 4.775 mm, L/D = 1.504** (Table III; minimum weight 0.067 g, i.e. 1.77 g/cm3 - a consistency
check against Zhu's 1.75-1.77 g/cm3). For MTV there is no MIL spec: it is the user's stated assumption. (3) Exact right-cylinder form function for the basket default AND the forward plenum's `'0d'`
law. (4) Pre-change library entries are in `SUPERSEDED_PYROGENS`, so committed `.ric` blocks upgrade on load. (5) B0 re-recorded.

Preflight: `import srm_1d` -> `C:\tmp\sv\ob1\srm_1d\__init__.py` (inside the worktree).

## 0. Summary

1. **Sources (section 1).** Zhu 2008 (closed bomb, B/KNO3, 2.8-50 MPa): `r = r0 P` up to 7 MPa and **flat** `7 r0` from 7 to 50 MPa, `r0 = 5.27e-3 .. 6.85e-3 m/s/MPa` (plateau 36.9-47.9 mm/s). MIL-P-46994B: standard
   BKNO3 is KNO3 70.7 / B 23.7 / polyester 5.6 %; burn rate 1.25-1.75 in/s as a closed-bomb average with no pressure stated; Table III gives the pellet geometries. Nakano & Kakami 2016 ("Mizushima" in the old YAML): 71.1 P^0.589 mm/s
   valid **5-200 kPa only**. Ma 2019 igniter paper: an inverse model, no burn law / form function / flame spread.
2. **The dominant defect (section 2).** The head basket burned every pellet at **constant area** until the pooled charge was gone: for a cube/cylinder that exhausts the mass after a regression of ~s/6 instead of s/2 (3x too fast at any rate). Fixed
   by the exact cylinder form function (default).
3. **Effect of the new defaults (current combo, DeMar on, BPNV 4.5 g, I-A pellet).** Peak **20.26 -> 16.72 MPa (-17.5 %)**, peak/plateau 2.40 -> **1.98**, pyrogen burn 10.7 -> 25.4 ms; 6 g: 26.29 -> 20.89 (3.13 -> 2.49); MTV 4.5 g:
   17.92 -> 14.33 (2.13 -> 1.69, pyrogen burn 45 -> 134 ms). DeMar off is within 1.2 % of DeMar on in every row. The plateau is unchanged (8.4 MPa vs exp 8.66). **The over-spike is reduced, not removed** (experiment 1.027).
4. **Why it remains (T8).** New default, 4.5 g: pyrogen gas alone 1.4 MPa; Ma erosive term off -> 5.6 MPa (1.25x plateau); Mukunda-Paul 10.9 MPa (1.62x). The residual is the Ma quasi-steady erosive response, not the pyrogen. Opt-in
   levers on top (rate capped at Zhu's plateau: 15.05 MPa, 1.78x; bed flame speed 3/10/30 m/s: 16.6-17.0, no effect) do not change that. The literal Zhu law (r = r0 P) does **not light** the head basket at 4.5 g under Peretz ignition.
5. **Ignition thresholds (T5).** Head basket: 1.43 g (pre-OB1) -> **3.0 g** (new default; 3 mm cubes 2.7 g), both DeMar on and off. Forward plenum: 0.027 -> 0.029 g with DeMar on, 2.21 -> 2.27 g with DeMar off. The real 4.5-6 g charge is 1.5-2x
   the new head-basket threshold.
6. **Canonical motors (Sutton forward plenum, DeMar on; T7).** Chunc peak 9.951 -> 9.977 MPa (+0.27 %), t_peak 36.7 -> 35.2 ms, plateau unchanged; BALLSstick 7.020 -> 7.061 (+0.58 %), pk/trough 1.62 -> 1.63; Hasegawa A unchanged (6.1073).
   With DeMar off none of the three lights at its Sutton-sized charge (before and after).
7. **Gates.** B0 re-recorded (forward-plenum area changed): clean-compile PASS (max rel error 2.0e-6), **warm run reproduces exactly (0.0)**; full suite 775 passed (section 8). The earlier warm-state t_peak miss is the documented compile-state noise (the pristine
   base fails/passes the same way); main's gate now compares t_peak with the step tolerance.

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

| P [MPa] | sim law 2e-5 P^0.5 | Nakano 71.1 P^0.589 (valid 0.005-0.2 MPa) | Zhu mean r0 (valid 2.8-50 MPa) | MIL-P-46994B (closed-bomb average) |
|---|---|---|---|---|
| 0.1 | 6.3 | 18.3 | (below range) |  |
| 0.2 | 8.9 | 27.6 | (below range) |  |
| 1 | 20.0 | (71, out of range) | (below range) |  |
| 2.8 | 33.5 | (130, out of range) | 17.0 | 31.8-44.5 |
| 7 | 52.9 | (224, out of range) | 42.5 | 31.8-44.5 |
| 20 | 89.4 | (415, out of range) | 42.5 | 31.8-44.5 |

The sim law is closest to Zhu near 7 MPa (1.1-1.4x). Below it (2.8 MPa) it is 2.0x fast (vs the mean r0), above it it keeps rising to 89 mm/s at 20 MPa while the bomb data are flat at 37-48
(1.9-2.4x). At ambient it is 0.34x Nakano. BPNV is not standard BKNO3 (Viton vs polyester/NC binder, 25 vs 23.7-28 % boron); no source covers it.
MTV: only Kubota 1987 (0.3-7 MPa, n ~0.35, in repo); nothing above 7 MPa (R16), so no plateau law is offered for MTV.

## 2. Part 2 -- audit of the sim's pyrogen model (as found, before OB1)

> **Status after the 2026-10-12 user decisions.** The defects listed in this audit have been addressed as follows: the basket area law is now the
> exact cylinder form function by default (`basket_area_law='shrinking'`; `'constant'` stays selectable), the forward plenum's `'0d'` law is the
> same exact function (it was `A0 (m/m0)^(2/3)`, identical only for L/D = 1), the initial-area formula is the single cylinder formula for any
> L/D (the `ld <= 1` sphere branch was wrong for L/D < 1), and the library pellets are the MIL-P-46994B Type I-A cylinder (section 4). Per-cell depletion,
> pooled-mass bookkeeping, lit-at-t=0 and the discarded condensed fraction are as described in the table (the first is fixed by the per-cell path).

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

| | option | grounded? | status |
|---|---|---|---|
| (a) | rate plateau above 7 MPa | **yes**, Zhu 2008 (B/KNO3, 4 runs) -- composition differs from BPNV; none for MTV | opt-in `Pyrogen.P_plateau`; -6.5..-10 % on the peak (T8: 16.72 -> 15.05 with the Zhu-cap form; 15.64 with a flat 7 MPa cap) |
| (a2) | Zhu piecewise through the origin below 7 MPa | valid >= 2.8 MPa only; in conflict with Nakano below | measured: the head basket does **not light** at 4.5 g under Peretz ignition; rejected as a default |
| (b) | exact geometric form function | **yes**, first principles | **default** (this report); largest single effect (-17.5 %) |
| (c) | progressive bed ignition / flame spread | **no source in the repo** | diagnostic knob only; no effect (T8 bed rows) |
| (d) | other | Nakano's ambient rate (3x the sim at 0.1 MPa); the unbridged 0.2-2.8 MPa gap; MIL average rate | listed, not stitched |

Pre-rebase probes (before the IG3 default flip and the geometry change) are archived in `results/pre_rebase/` and `probes/pre_rebase/`; their ranking (area -16 %, plateau -5..-7 %, both -25 %, bed speed 0) is
reproduced by the post-rebase rows of T3/T8 with the I-A pellet (area law -17.5 %, plus the Zhu cap -25 %).

## 4. Part 4 -- implementation

| API | default | effect |
|---|---|---|
| `PyrogenChamber.basket_area_law` (`run_from_ric(basket_area_law=)`, `data.igniter.basket_area_law`, GUI metadata) | **`'shrinking'`** (was absent) | exact right-cylinder form function with `L/D = Pyrogen.particle_LD_ratio`: `D = D0-2s, H = H0-2s, A = pi D H + pi D^2/2, V = pi D^2 H/4`, burn-out at `s = min(D0,H0)/2`; per-cell mass bookkeeping (each cell on its own host pressure); `'constant'` = the pre-OB1 pooled constant-area law (opt-in) |
| forward plenum `burn_law='0d'` | exact cylinder form function (was `A0 (m/m0)^(2/3)`) | identical for L/D = 1; for I-A L/D = 1.504 the area differs; `'end_burning'` unchanged |
| `build_pyrogen_chamber` initial area | `A = m (4 L/D + 2)/(rho (L/D) d)` for any L/D | the former `ld <= 1` sphere branch (`6m/(rho d)`) was wrong for L/D < 1; identical at L/D = 1 |
| `Pyrogen.P_plateau` [Pa] | 0.0 (off) | opt-in rate plateau (YAML / `.ric` block / GUI metadata) |
| `PyrogenChamber.bed_flame_speed_m_s` | 0.0 | DIAGNOSTIC only, no sourced value (R5) |
| library `bpnv.yaml`, `mtv.yaml` | `particle_diameter_m 3.175e-3`, `particle_LD_ratio 1.504` | MIL-P-46994B Table III Type I-A (BPNV: cited; MTV: the user's assumption, stated in the YAML) |
| `SUPERSEDED_PYROGENS` | BPNV and MTV pre-OB1 entries added (each with and without `P_plateau: 0.0`) | every committed `motors/*.ric` embeds the pre-OB1 BPNV and upgrades on load (test); a customised pellet never matches |

The exact function (`igniter_plenum._cyl_area_fraction`) inverts mass fraction to regression depth by bisection and returns the closed form `f^(2/3)` for L/D == 1; `_cyl_area_fraction_general` is the
bisection version so the L/D = 1 agreement is tested. Internal array-length / optional-argument changes are in DEVNOTES. `thermite_cual` and `mtv_fast` are untouched (L/D = 1, so the code change is a no-op for them).

Tests (`tests/test_ob1_pyrogen_bed.py`, 28): hand-checked cylinder geometry for L/D = 0.4, 1, 1.504, 3; endpoints and burn-out (needle / disc remnants); **L/D = 1 equals f^(2/3)** (general routine and
plenum `_burn_area`); plenum area at L/D 1.504; per-cell law at L/D 1.504; analytic web time `(D0/2)/r` for a single cylinder (0.5 %) and a cube (0.2 %); initial-area formula for L/D = 0.5, 1, 1.504, 3; library geometry;
plateau formulas; per-cell mass conservation (gas = gmf x solid); superseded snapshots upgrade (BPNV and MTV, with/without `P_plateau`, customised pellet blocked) and a committed `.ric` block upgrades on load; run-level default = shrinking,
constant-area opt-in conserves mass, default longer burn.

## 5. Part 5 -- measurements (post-rebase; raw `results/ob1r_*.jsonl`, tables `results/ob1r_tables.md`)

Models: **pre-OB1** = 3.2 mm L/D 1 pellet + constant-area basket law (the state before this work; the plenum law is identical there); **new default** = I-A pellet + exact area law; opt-in variants and
3 mm-cube rows are shown for comparison.

### 5.1 Realistic rows (Chunc head basket; 4.5 / 6 g x BPNV / MTV; DeMar on and off; t_max 0.2 s)

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

t_max 1.0 s (plateau = median 0.4-0.9 s; experiment peak 8.89, plateau 8.66, **ratio 1.027**):

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

### 5.2 Ignition thresholds

| topology | model | DeMar on [g] | DeMar off [g] |
|---|---|---|---|
| head_basket | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 1.43 | 1.43 |
| head_basket | I-A pellet, constant-area basket law (opt-in) | 1.51 | 1.51 |
| head_basket | **new default**: I-A pellet, exact cylinder area law | 3.01 | 3.01 |
| head_basket | 3 mm cube, exact area law (L/D 1) | 2.71 | 2.71 |
| forward_plenum | pre-OB1: 3.2 mm L/D 1 pellet, constant-area basket law | 0.0268 | 2.21 |
| forward_plenum | **new default**: I-A pellet, exact cylinder area law | 0.0289 | 2.27 |
| forward_plenum | 3 mm cube, exact area law (L/D 1) | 0.0258 | 2.17 |

Peak vs charge (head basket, BPNV):

| DeMar | model | 0.6 g | 0.9 g | 1.36 g | 1.8 g | 2.27 g | 2.93 g | 3.6 g | 4.54 g | 6.3 g | 9.0 g |
|---|---|---|---|---|---|---|---|---|---|---|---|
| on | pre_ob1 | NL | NL | NL | 9.9 (1.13) | 10.2 (1.19) | 13.4 (1.57) | 16.3 (1.93) | 20.4 (2.42) | 27.5 (3.28) | 37.5 (4.51) |
| on | lib | NL | NL | NL | NL | NL | NL | 14.0 (1.65) | 16.8 (2.00) | 21.8 (2.60) | 29.0 (3.48) |
| on | cube3 | NL | NL | NL | NL | NL | 12.1 (1.46) | 14.4 (1.70) | 17.3 (2.05) | 22.5 (2.69) | 30.1 (3.61) |
| off | pre_ob1 | NL | NL | NL | 9.9 (1.13) | 10.2 (1.19) | 13.4 (1.57) | 16.3 (1.93) | 20.2 (2.39) | 27.3 (3.25) | 37.3 (4.48) |
| off | lib | NL | NL | NL | NL | NL | NL | 14.0 (1.65) | 16.6 (1.97) | 21.5 (2.56) | 28.7 (3.45) |
| off | cube3 | NL | NL | NL | NL | NL | 12.1 (1.46) | 14.4 (1.70) | 17.1 (2.03) | 22.3 (2.66) | 29.9 (3.59) |

### 5.3 Canonical motors with their usual igniters

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

### 5.4 Attribution and ablation of the residual spike

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

## 6. What the probes do and do not establish

- Established: the constant-area basket burn was ~3x too fast per unit regression; with the exact cylinder form function the 4.5 g BPNV burn stretches from 10.7 to 25.4 ms (MTV 45 -> 134 ms) and the peak falls 17-22 %.
  The remaining 1.8-2.0x is the Ma erosive response (erosion off: 1.25x) and does not come from the pyrogen burn model.
- Not established (pending): whether BPNV follows Zhu's B/KNO3 plateau (composition, R15); the 0.2-2.8 MPa rate; a sourced flame-spread speed (R5); the real Chunc igniter's pressure history (R7); MTV above 7 MPa (R16); the MTV pellet
  dimensions (user assumption).
- The thresholds are for 99 of 100 cells lit within 0.2 s; points scanned above the first lit charge stayed lit (no ragged boundary in this set, unlike the pre-IG3 gnielinski runs).

## 7. Files

Probes: `probes/ob1_case.py`, `ob1_sweep.py` (realistic / long / thresh / thresh_plenum_low / mass / extra), `ob1_valid.py`, `ob1_hist.py`, `ob1_table.py`, `ob1_tables.py`. Raw: `results/ob1r_*.jsonl`, `results/ob1r_tables.md`; B0 evidence:
`results/ob1_B0_*.txt`. Pre-rebase material: `results/pre_rebase/`, `probes/pre_rebase/`. Run from the worktree root with `PYTHONPATH=$PWD`.

## 8. Gates

- **Rebase:** onto `openmotor-frontend` (IG3 stage 1, then `de69787`); conflicts were documentation only (DEVNOTES, PAYWALLED_REFERENCES).
- **B0 on the rebased tree BEFORE any default change (item 2).** Clean compile: PASS (Chunc peak 1.2e-7, t_peak 1.0e-9; BALLSstick peak 4.9e-7; Hasegawa A 1.2e-10 / steps 4.3e-7; `results/ob1_B0_rebased_clean.txt`).
  The following warm run FAILED the (then peak-style) t_peak check: Chunc t_peak 0.03673671 vs 0.03673504 (4.55e-5; one time step), peak 2.5e-7, BALLSstick peak 2.5e-6 (`ob1_B0_rebased_warm.txt`).
  To decide whether that was my code, I ran the pristine `bb2f5ab` (archive, own repo): **clean compile FAILS the same way (Chunc t_peak 4.55e-5, the identical value 0.0367367137), warm run reproduces the baseline exactly (0.0)**
  (`ob1_B0_base_bb2f5ab_clean.txt` / `_warm.txt`). So the t_peak miss is the flat-peak argmax flipping one step with the compile state, present in the unmodified base; it is not caused by OB1. `de69787` (main) now compares t_peak with the step tolerance.
- **B0 re-recorded (item 7)** after the defaults changed (the forward-plenum burn area and the library pellet changed): recorded from a WARM-state run, previous IG3 values kept in the `BASELINE` comment block:

  | row | old (IG3) | new |
  |---|---|---|
  | Chunc peak / t_peak / plateau / steps | 9.950855636512404 / 0.036735041862973064 / 8.141968067255977 / 532381 | **9.977396119842126 / 0.03519436869536591 / 8.141928366040833 / 532348** |
  | BALLSstick peak / trough / mid / steps | 7.020318627616664 / 4.343390826697384 / 4.855316005870968 / 843726 | **7.060894562817688 / 4.343382326291922 / 4.855564928020005 / 843680** |
  | Hasegawa A peak / t_peak / steps | 6.107268310023869 / 2.411042627819566 / 2336271 | **6.107263748109281 / 2.4109782579175234 / 2336318** |

  Verification, clean then warm with nothing in between (`results/ob1_B0_final_clean.txt`, `ob1_B0_final_warm.txt`): **clean PASS** (largest relative error 2.0e-6, BALLSstick peak; Chunc peak 4.1e-7, t_peak 1.1e-8, steps identical except Hasegawa A 8.6e-7);
  **warm run reproduces every value exactly: relative error 0.000e+00 on all 11 numbers.** Two further warm runs were bit-identical to each other.
- **Full suite:** `python -m pytest tests/ -q -p no:cacheprovider`: **775 passed, 0 failed** (654 s), run on the committed tree with the final code (later commits touch only probes and documentation). Earlier on this branch a run had 1 failure
  (`test_pyrogen_block_round_trip_every_field`: `_pyrogen_from_mapping` dropped `P_plateau`), fixed with a test.


## 9. Recommendation

1. The exact area law at the basket default is justified and now in place; it is the largest grounded pyrogen-side effect and moves the head-basket threshold to 3.0 g, within a factor 1.5-2 of the real 4.5-6 g.
2. Keep `P_plateau` opt-in (composition mismatch, low-pressure branch unsourced, no MTV source, literal Zhu law does not light the basket).
3. Do not expect further pyrogen-model work to remove the 1.8-2.0x peak at 4.5-6 g. The remaining lever is the erosive response at transient mass flux (held Mukunda-Paul decision) and the whole-grain ignition cascade (IG1 candidate C).
4. Fix the BKNO3 provenance in any doc still quoting "Mizushima / 3.5e-5 / 7 MPa" (`docs/v0_7_3/references/pyrogen_physical_params_lit_dive.md`).

## 10. New paywalled / missing references (added to `PAYWALLED_REFERENCES.md`)

R5 pellet-bed flame spread (candidate citations, unverified DOIs); R13 MIL-STD-286 method 801.1; R14 Mizushima et al. 2024; R15 Zhu 2005 / Jiang 2002; R16 MTV rate above 7 MPa. MIL-P-46994B Table III was read for OB1 (Type I-A, verified from the page image).
