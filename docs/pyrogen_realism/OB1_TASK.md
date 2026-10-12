# OB1 — fix option B: a literature-grounded pyrogen burn model (head basket first)

**User go-ahead 2026-10-11 ("Attempt option B").** IG1 and IG2 showed the
realistic-charge over-spike is set by pyrogen gas still feeding a grain that
is already lit (the peak sits at pyrogen burn-out; cutting the pyrogen at 2 ms
gives 9.9 MPa instead of 21.7). Fixing the ignition heat transfer (IG2) did
not change it. OB1 targets the pyrogen side: how fast, how long and at what
rate the igniter charge actually produces gas.

Read first: `IG1_REPORT.md` (sections 0, 4, 6), `IG2_REPORT.md` (sections 0, 4.2),
`PAYWALLED_REFERENCES.md` (status table at the end).

## Decided elsewhere (do NOT implement here)

The user adopted Peretz as the default ignition h and dropping the DeMar flux;
those flips land separately (IG3, together with a new igniter-sizing default).
In OB1 treat them as an opt-in combo: `ignition_heat_transfer='peretz'` +
pyrogen `heat_delivery_mode='none'` ("forward combo"). Measure both the current
defaults and the forward combo.

## Real igniter (user, R7)

Real Chunc: **5–6 g of ~3 mm cubes of BKNO3 or MTV** (uncertain which), head
basket, no startup spike in the real trace (flat ~8.8 MPa; the committed CSV
starts at 10 ms, so the ignition rise itself is unrecorded). Realistic rows:
4.5 and 6 g, `bpnv` and `mtv`, `particle_diameter_m=3e-3`, `particle_LD_ratio=1`,
`injection_topology='head_basket'`.

## Part 1 — literature extraction (token-light)

Sources in `docs/references/`: `zhu2008.pdf` (R1, Chinese), `nakano2016.pdf`
(R2, Japanese), `MIL-P-46994B.pdf` (R3), and the Ma 2020 igniter paper
(`International Journal of Aerospace Engineering - 2020 - Ma - A Model for
Igniter Mass Flow Rate History Evaluation for.pdf`), plus anything else already
in `docs/references/` that models igniter pellet beds or porous-bed flame spread.

Extract: BKNO3 (and MTV if present) burn rate r(P) over ~0.1–30 MPa with the
stated range, composition, density and pellet geometry; any igniter mass-flow
history model (Ma 2020) and its assumptions (progressive pellet ignition? form
function? validated how?). Compare with the library `bpnv.yaml` (a 2e-5, n 0.5)
and `mtv.yaml` (a 4.4e-5, n 0.35): check the units the sim uses for a, and note
that BPNV (B/KNO3/Viton 25/60/15) is not the standard BKNO3 composition.

**Token rules for the translations (user):** extract text locally with
`pdftotext -layout` (`/mingw64/bin/pdftotext`) or PyMuPDF (`fitz`, in the pyenv
python) into files in your own scratch directory; `grep`/`sed` the extracted
text and read only the relevant passages (abstract, method, results, tables,
figure captions); translate those passages yourself. **Never call Read on a
whole PDF.** For a figure or table you need, render that single page with
`fitz` at <= 110 dpi, crop to the region, and Read the PNG. If a PDF has no
text layer, say so and render only the pages you need.

## Part 2 — audit the sim's pyrogen model

For each topology (`forward_plenum`, `head_basket`): the burn law and its units;
how the burning area evolves (the forward plenum uses `A0 (m/m0)^(2/3)` in
`igniter_plenum._burn_area`; the uncontained path in
`simulation._compute_uncontained_pyrogen_mdot` appears to use a constant
`A_burn_per_cell` until depletion — confirm, and say whether that is
geometrically right for cubes); which pressure drives the rate (local bore P);
all pellets lit at t = 0; and where the condensed fraction
(`gas_mass_fraction` 0.781 for BPNV) goes — mass and enthalpy.

## Part 3 — model options (each grounded in a source or first principles)

Candidates, to be confirmed or rejected by Parts 1–2:
(a) burn-rate law from R1–R3 (for example a plateau above ~7 MPa, if zhu2008
supports it); (b) a geometric form function for the head-basket pellets;
(c) progressive bed ignition / flame spread through the basket — only if a
source gives it (Ma 2020 igniter paper; porous-bed flame-spread literature); if
none is in the repo, list candidates for R5 in `PAYWALLED_REFERENCES.md` with
DOIs marked unverified; (d) anything else the sources show matters. Rank them
by their effect on the over-spike using quick probes before building anything.

## Part 4 — implement (opt-in only)

Implement the grounded options as opt-in choices (for example a new pyrogen
burn-law value or library YAML fields). Defaults stay bit-identical. Add tests:
hand-checked formulas, default identity, mass conservation of the pyrogen.

## Part 5 — measure

1. Chunc head basket, realistic rows (4.5 / 6 g × `bpnv` / `mtv`), t_max 0.2 s
   (and 1.0 s for peak/plateau), under {current defaults, forward combo} ×
   {current pyrogen model, each new option, all combined}: first ignition, t90,
   peak, t_peak, plateau, peak/plateau, pyrogen burn duration.
2. Head-basket ignition threshold with the new model (as IG2 section 4.1).
3. Canonical motors with their usual igniters (Chunc, BALLSstick, Hasegawa A;
   Sutton-sized forward plenum): unaffected, or by how much.
4. Gates: B0 PASS with defaults (default tolerances); full suite green.

## Deliverables

`docs/pyrogen_realism/OB1_REPORT.md` (sources and what each says, the audit,
the options and their ranking, all tables, a recommendation — the user
decides defaults); opt-in implementation + tests + probes + raw results on your
branch. Paywalled-reference rule applies: list new needs, mark dependent
conclusions "pending".
