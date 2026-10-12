# IG2 — fix option A: literature-grounded ignition heat transfer + the DeMar question

**User go-ahead 2026-10-11** for IG1's option A, with an explicit condition:
*"Going to need good justification for time variance of DeMar flux, or
dropping the DeMar values altogether, since we're now fitting time-averaged
data to a time-dependent function."* IG2 implements OPT-IN modes (defaults
unchanged), measures them, and recommends — the user decides defaults.

Read first: `IG1_REPORT.md`, `PAYWALLED_REFERENCES.md` (R6 resolution and
the Keller geometry caveat), `docs/pyrogen_realism/TASKS.md`.

## New references available in `docs/references/`

- `peretz1974.pdf` — R4, Peretz/Kuo/Caveny/Summerfield, starting transients
  with high internal gas velocities. **The key source for A.**
- `keller1966.pdf` — DOI 10.2514/3.3678; Eq. 11 h-correlation (channel-
  specific, D_h ≈ 5.1 mm, early BL transition), Eq. 17 T_si; ranges recorded.
- `STARTING YOUR UP-GOERS_pub (1).pdf` — DeMar 2021 (TRATECH), the source of
  the 69.4 cal/cm²·s BPNV flux; repo summary in
  `docs/v0_7_0/references/primary_sources_summary.md` §2.5 (sealed 5.75 ×
  2.75 in vessel, one K-type thermocouple + 1000 psi transducer; ONE flux
  number per pyrogen).
- `zhu2008.pdf` (Chinese), `nakano2016.pdf` (Japanese), `MIL-P-46994B.pdf` —
  R1–R3, pyrogen burn rate (option B; not needed for IG2 — skim only if a
  number bears on A). R5 not found.
- PDF text extraction works with the Read tool (no page images:
  `pdftoppm` is missing). If a PDF yields no text, say so; do not guess.

## Real igniter (user, R7)

The real Chunc was lit with **5–6 g of ~3 mm cubes of BKNO3 or MTV**
(uncertain which). Use 4.5 g and 6 g, BPNV and MTV library materials,
`particle_diameter_m=3e-3`, `particle_LD_ratio=1`, head basket, for the
realistic-charge rows.

## Part 1 — why is the ignition h ~4.8× Keller? (decompose before fixing)

Ignition heat transfer (`_bare_heat_transfer_coeff` and the Goodman path in
`simulation.py`) uses the same Gnielinski chain as Ma's erosive model. At the
first cells to ignite (IG1 cases 0.45 g and 4.5 g, head basket), report the
factor each term contributes relative to a plain developed-flow pipe value:
(a) the property-ratio correction (T_gas/T_surf)^kappa — at a COLD wall the
ratio is ~10, far outside the range Gnielinski validated (find and cite that
range: Gnielinski / VDI or the Ma 2020 paper in `docs/references/`); (b) the
entrance term (D/L)^(2/3) with L = x_from_head; (c) roughness via Haaland
friction (32 µm); (d) the laminar/turbulent blend. Report which terms are
outside their published validity at ignition conditions.

## Part 2 — choose a literature-grounded ignition heat-transfer closure

From Peretz 1974 (primary) and Keller 1966: what correlation does Peretz
use / validate for convective ignition in SRM bores during starting
transients, with what data and range? Candidates: (i) Peretz's bore
correlation; (ii) Keller Eq. 11 in Nu–Re form using Keller's D_h (state the
extrapolation); (iii) Gnielinski restricted to its validated property-ratio
range (only if the source says so). Implement the chosen one(s) as opt-in
`ignition_heat_transfer=` modes (default = today's, bit-identical). The Ma
erosive path is NOT changed — Ma 2020 stays the erosive model of record.

## Part 3 — the DeMar flux: justify a time shape, or drop it

Read DeMar 2021 and record exactly what the 69.4 cal/cm²·s is: how it was
derived (thermocouple heating rate? energy/area over which window?), the
geometry, and whether any time-resolved record exists. Then:
- a time-dependent form is acceptable ONLY if it is derived from DeMar's own
  measurement (e.g. a published time history, or a definitional relation
  that preserves the measured average for the same configuration) — no
  assumed shapes fitted to the average;
- otherwise the recommendation is to DROP the DeMar flux and let ignition
  come from the physically modelled channels (convective heating by the
  igniter gas with the Part-2 h, plus radiation only where already modelled).
  Dropping needs no new code: `Pyrogen.heat_delivery_mode='none'` exists.
State the conclusion with citations.

## Part 4 — measure (opt-in modes only; defaults unchanged)

Combos: ignition h ∈ {current, chosen} × DeMar ∈ {step (current), off,
(time form only if Part 3 justifies it)}.
1. Ignition-threshold sweep (as IG1 §1), head basket and forward plenum.
2. Realistic charges: 4.5 g and 6 g, BPNV and MTV, 3 mm cubes, head basket —
   first ignition, t90, peak, t_peak, plateau.
3. Validation, canonical motors with their usual igniters: Chunc,
   BALLSstick, Hasegawa A — peak / plateau / ratios vs experiment.
4. Gates: B0 PASS with defaults (default tolerances); full suite green.

## Deliverables

Opt-in implementation + tests (incl. default identity) on your branch;
`docs/pyrogen_realism/IG2_REPORT.md` with the decomposition, the chosen
closure and its justification, the DeMar conclusion, all tables, and a
recommendation for defaults. Paywalled-reference rule applies: list new
needs in `PAYWALLED_REFERENCES.md`; mark dependent conclusions "pending".
