# IG1 -- references needed but not obtainable (abstract-only / paywalled)

Rule (user, 2026-10-11): nothing below is paraphrased from secondary sources in
`IG1_REPORT.md`; every conclusion that would need these is marked
"pending reference". DOIs are given only where verified; "DOI not verified" means
the identifier was NOT guessed -- the user/maintainer should confirm it.
Please drop PDFs into `docs/references/`.

| # | Citation | DOI / URL | Why IG1 needs it | Quantity / figure to extract |
|---|---|---|---|---|
| R1 | Zhu M.-S. et al., "Investigation on Burning Rate of B/KNO3 Pellet in Closed Bomb Test", Chinese Journal of Energetic Materials (Hanneng Cailiao) 16(2):204-206, 2008 (only the abstract was seen: step-wise rate, ~proportional to P below 7 MPa, ~flat 7-50 MPa) | DOI not verified; https://energetic-materials.org.cn/hnclen/article/abstract/20080222 | H3 / fix option B: is the BPNV pellet pressure exponent n=0.5 (a=2e-5) valid at the 7-25 MPa the head-basket pellets see? | Burning rate vs pressure curve (1-50 MPa): numeric r at 1, 3, 7, 10, 20 MPa; pellet diameter, density, composition (B:KNO3:binder) |
| R2 | Mizushima et al. (2016), "Burning Rate Measurements of 3.2-mm-dia. BKNO3 Pellets in Low Pressures", J. Japan Soc. Aeronaut. Space Sci. 64(2) (abstract-only; repo YAML cites r = 71.1 P^0.589 mm/s, 5 kPa-7 MPa, B:KNO3:NC 28:70:2) | DOI not verified; https://www.jstage.jst.go.jp/browse/jjsass/64/2/_contents | H3: the YAML comment's a, n; whether it holds up to 7 MPa and at what binder content | Fitted a, n and their stated validity range; r at 0.1 and 7 MPa; effect of binder (Viton 15% in BPNV vs 2%) |
| R3 | MIL-P-46994B, "Pellets, Boron/Potassium Nitrate" (BPN), min burn rate 31.75 mm/s (seen only as a secondary mention in a thesis abstract, https://tezara.org/theses/587705) | DOI n/a (military spec) | H3 / fix B: absolute rate anchor for BKNO3 pellets | Test pressure and conditions of the 31.75 mm/s minimum (1000 psi?), pellet geometry |
| R4 | Peretz, Kuo, Caveny, Summerfield, "Starting transient of solid-propellant rocket motors with high internal gas velocities", AIAA J. 11(12), 1973 (flame spreading + convective ignition in high-velocity bores; the foundation for fix option A) | DOI not verified | H1/H2/H4: measured convective heat-transfer coefficient and ignition delay in high-G bores; flame-spread speed vs port velocity | h (or Nu) correlation vs G and Tg in bores; ignition delay vs flux; flame-spread rate |
| R5 | A published measurement or model of flame spread through a loose BKNO3 pellet bed / basket (tau_bed = 5-10 ms is a repo memory-level number with no cited source) | not identified | H3(a): bed-spread-limited pyrogen output estimate | Bed flame-spread speed vs pellet size / packing / pressure |
| R6 | Keller, Baer, Ryan, AIAA J. 4(8):1358-1364, 1966 -- ALREADY in `docs/references/keller1966.pdf` (Eq. 11 read from the page image). Only needed in full: the G and Tg range over which Eq. 11 was fitted | DOI not verified for this file | H1: our G = 3-60 g/cm2 s lies at/below the edge of Keller's data; need the fitted range to bound the extrapolation | Range of G [g/cm2 s] and Tg [K] covered by the Eq. 11 data |
| R7 | The user's own Chunc igniter details/pressure data (BKNO3 4.5 g, head basket): pellet size, basket geometry, and the first ~100 ms of the pressure trace (the committed CSV starts at plateau, see report 5) | user data | validation of every fix option | Raw trace from electrical fire signal; basket drawing |

## R6 — RESOLVED by the reviewer (2026-10-11) from `docs/references/keller1966.pdf`

- DOI **10.2514/3.3678** (printed on every page of the PDF) — verified.
- Test ranges (p. 1361): heat flux ~20–160 cal/cm²·s; gas velocity 50–800 m/s
  (Mach 0.07–1.0); gas stagnation temperature 1000–2600 K; pressure 14–25 atm.
  G is not tabulated; from those ranges (N₂, ρ ≈ 1.8–8.5 kg/m³) G spans roughly
  9–680 g/cm²·s, so the sim's low-G cells (~3 g/cm²·s) are an extrapolation.
- Eq. 11 (air/O₂/N₂): h = 1.435e-4 · Tg^0.3 · G^0.905 cal/cm²·s·K (Tg in K,
  G in g/cm²·s), average deviation 6%.
- **Caveat for fix option A:** Eq. 11 was fitted for Keller's own test
  section — a 0.500 × 0.125 in channel (D_h ≈ 5.1 mm) — and the paper
  attributes its strong G-dependence to early laminar→turbulent transition
  from the shock tube's freestream turbulence. It is not a geometry-general
  correlation: transferring it to a ~25 mm motor bore needs a dimensionless
  (Nu–Re) form and should be treated as an extrapolation.
- Also in the PDF: Eq. 17, thermal-ignition temperature
  T_si = 300 K + 286.1·F^0.08 (F in cal/cm²·s), the basis of the 756 K
  check in REOPENED §9.2.

## IG2 additions (2026-10-11) -- needed, not obtainable here

Resolved/available for IG2 (no action): R4 Peretz et al. is in the repo as
`docs/references/peretz1974.pdf` (NTRS 19740005393; NASA CR / Princeton AMS Rep. 1100, Apr 1973,
= Peretz's Ph.D. thesis; the AIAA J. 11(12) 1973 paper cited as R4 is the journal version, which
was not read for IG2).
Eq. III-22..III-31 were read from page images of the PDF (the text layer drops the equations).

| # | Citation | DOI / URL | Why IG2 needs it | Quantity to extract | Conclusions pending |
|---|---|---|---|---|---|
| R8 | Gnielinski V., "On heat transfer in tubes", Int. J. Heat Mass Transf. 63:134-140 (2013) (the correlation Ma 2020 and srm_1d use) | DOI not verified | Part 1(a): the validity range of the property-ratio correction (T_gas/T_wall)^kappa for gases; Ma 2020 only says "kappa = 0.45 for air" with no range | Stated range of T_gas/T_wall (and Re, Pr, L/D) for the K factor and for the entrance term `1+(D/L)^(2/3)` | "ratio 6.5-10.7 at the cold wall is outside Gnielinski's validated range" -- pending. Only the in-repo comparison (Ma's own validation: T_flame 1667-3041 K vs the 1000 K surface, ratio ~1.7-3) is stated |
| R9 | Cohen N.S., Derr R.L., Price C.F., (Peretz ref. 158) correlation of the Rocketdyne igniter-duct data; Carlson & Seader AIAA J. 5(7):1272-1279 (1967); Wrubel & Carlson AFRPL-TR-67-267 (1967); Carlson & Seader AFRPL-TR-65-158 (1965) | Peretz refs 50-53, 158 (DOI not verified) | Peretz fixed a6 = 2.56e-2 and the (A_p/A_t)^0.4 factor from these duct data (A_p/A_t 1.1-3, x/d <= 10, nitrogen and pyrogen gas, 570-620 psig) | Their measured h range (G, Tg, pressure) and the exact A_p/A_t dependence | The Ap/At exponent and a6 validity range (Chunc sim A_p/A_t = 0.56 is below the 1.06-3 tested range) -- pending |
| R10 | Peretz et al., AIAA J. 11(12), 1973 (journal version of R4) | DOI not verified | Whether the journal version states the correlation's tested range (G, Re) more precisely than the thesis | tested Re / G / pressure range for Eq. III-31 | The "validated range" of the Peretz closure is taken from the thesis text only |
| R11 | DeMar J.S., full TRATECH 2021 paper / data (not the slides) | https://www.linkedin.com/in/jsdemar/ (author contact) | Part 3: the slides give q = kE W M / tA / (2 pi R L) but not R, L (the "effective igniter length" and thermocouple distance) or any time record | R, L per igniter; the raw pressure and thermocouple traces | "No time-resolved flux exists" is established for the slides only |
| R12 | Rough-wall heat-transfer: Dipprey & Sabersky, Int. J. Heat Mass Transf. 6:329-353 (1963); Owen & Thomson, J. Fluid Mech. 15:321-334 (1963); Gowen & Smith, Int. J. Heat Mass Transf. 11 (1968) | DOIs not verified | Part 1(c): quantify how much of the Haaland-friction roughness factor (x1.4-1.5 at 32 um, Re 5e4-9e4) survives when the Reynolds analogy is invalid for rough walls (Peretz Ch. III-D states it is invalid and worsens with Re) | St_rough / St_smooth vs k_s+ and Re | The size of the roughness over-statement is not quantified |

## Status after the user's fetch (2026-10-11, later)

Now in `docs/references/` (read them there; conclusions above marked "pending" can be revisited):

| # | File | Note |
|---|---|---|
| R1 | `zhu2008.pdf` | Chinese |
| R2 | `nakano2016.pdf` | Japanese |
| R3 | `MIL-P-46994B.pdf` | military spec |
| R4 | `peretz1974.pdf` | Princeton AMS Rept 1100 (was `19740005393.pdf`) |
| R8 | `gnielinski2013.pdf` | |
| R10 | `peretz1972.pdf` | AIAA Paper 72-1119 (the conference paper that preceded AIAA J. 11(12) 1973; DOI 10.2514/3.50676 is cited for the journal version in `docs/core_loop_opt/LITERATURE_ALLSPEED_PISO.md`, unverified). Its Eq. 14 prints a constant of 1.78e-3 where the thesis Eq. III-31 has 1.56e-3 -- CHECKED in IG3 stage 1 from the page images (paper p. 6, thesis p. 40): the SAME law (same variables, exponents, CGS units: h = a Pr^-0.6 cp (pu/R)^0.8 W^0.1 T_af^-0.67 (x_r d_h)^-0.1 (A_p/A_t)^0.4), with a different numerical constant -- 1.78e-3 (paper; no derivation given, equivalent to a6 = 2.92e-2) vs 1.56e-3 (thesis, a6 = 2.56e-2 set to match the Rocketdyne duct data to x/d = 10). The thesis later varies the constant (1.56e-3 -> 1.704e-3, +9.2 %) as a sensitivity study. The code carries the thesis constant, unchanged. The paper's constant gives h 1.141x higher. Pending: R9 / R10 (journal version) to see which calibration the authors finally endorsed |
| R11 | `STARTING YOUR UP-GOERS_pub (1).pdf` | The slides are the whole DeMar source: no R, L or time record exists beyond them |
| R12 | `dipprey1963.pdf`, `owen1963.pdf`, `gowen1968.pdf` | |

Still open: **R5** (pellet-bed flame spread; not found), **R7** (raw Chunc ignition trace; unavailable),
**R9** (Rocketdyne igniter-duct data; see `R9_SOURCES.md`: AIAA J. 5(7) 1967 DOI 10.2514/3.4183 is
verified and fetchable; the two AFRPL reports and Cohen-Derr-Price are not online).

## OB1 additions (2026-10-11)

R1, R2, R3 were READ for OB1 (page images; the text layers are unusable for R1/R2). What they do and do not
support is in `OB1_REPORT.md` section 1. Remaining / new needs:

| # | Citation | DOI / URL | Why OB1 needs it | Quantity to extract | Conclusions pending |
|---|---|---|---|---|---|
| R5 | Pellet-bed / porous-charge flame spread. Search candidates (NONE read, NONE verified, none in repo): Kuo K.K., Vichnevetsky R., Summerfield M., "Theory of flame front propagation in porous propellant charges under confinement", AIAA J. 11(4):444-451 (1973); Gough P.S. & Zwarts F.J., "Modeling heterogeneous two-phase reacting flow", AIAA J. 17(1):17-25 (1979); Baer M.R. & Nunziato J.W., "A two-phase mixture theory for the deflagration-to-detonation transition in reactive granular materials", Int. J. Multiphase Flow 12(6):861-889 (1986); any measurement of flame spread through a loose BKNO3 / MTV pellet basket | DOIs not verified | Option (c): a sourced bed flame-spread speed (the repo's tau_bed = 5-10 ms is memory-level only). OB1 probes show a 3-30 m/s sequential lighting changes the 4.5 g peak by <2 %, so this is lower priority than it looked | Flame-spread speed vs pellet size / packing / pressure for BKNO3 or MTV beds | Whether a sourced bed model can reduce the spike: pending (diagnostic result only) |
| R13 | MIL-STD-286, Method 801.1 (closed-bomb burning rate) | DOI n/a (military std) | MIL-P-46994B 4.6.2 defines the 1.25-1.75 in/s (31.8-44.5 mm/s) minimum as R = h/(2t) from a 10 g closed-bomb run, t = first pressure indication to stabilised maximum, with NO pressure stated | Bomb volume, loading, and the pressure range the average rate refers to | Whether the MIL rate is evidence for the Zhu plateau (they agree numerically: 36.9-47.9 vs 31.8-44.5 mm/s) is pending; stated only as consistent |
| R14 | Mizushima et al., "Measurements of Burning Rates of BKNO3 Pellets with a Diameter of 3.2 mm under Vacuum Conditions", Trans. JSASS Aerospace Tech. Japan 22 (2024) (cited in `docs/v0_7_3/references/pyrogen_physical_params_lit_dive.md`; J-STAGE open access, https://www.jstage.jst.go.jp/article/tastj/22/0/22_22.25/_pdf) -- NOT read | URL not verified | Possible second low-pressure data set / the binder identity | Rate vs pressure, binder, pellet length | The 0.2-2.8 MPa gap between Nakano (<= 0.2 MPa) and Zhu (>= 2.8 MPa) stays unbridged: pending |
| R15 | Zhu M.-S., Long X.-P., Jiang X.-H., et al., "B/KNO3 burning performance parameter calculation", Acta Armamentarii (Bingong Xuebao) 26(3):320-322 (2005) and Jiang X.-H., He B., Chen L., "Study on burning characteristics of boron potassium", Proc. Symposium on Propellants, Explosives and Insensitive Munitions, Zhuhai (2002) (refs [2], [3] of Zhu 2008) | DOIs not verified | Zhu 2008 does not state the B/KNO3 composition or binder (it cites MIL-STD-1901A); the earlier papers may | Composition, binder, density | Whether Zhu's r0 applies to BPNV (25:60:15 B/KNO3/Viton) or MIL-type (24:71:5): pending |
| R16 | MTV (Mg/Teflon/Viton) burn rate ABOVE 7 MPa | Kubota & Serizawa 1987 (in repo, 0.3-7 MPa only) | OB1's plateau option has a sourced law only for B/KNO3; MTV has none | Strand or closed-bomb MTV rate 7-50 MPa | Any MTV plateau: pending (not implemented) |
