# v0.7.0 References — Index

This directory contains paper extractions and source-document summaries for the v0.7.0 hot-gas plenum igniter model. A fresh agent should read these to understand the literature foundations of the [DESIGN.md](../DESIGN.md) architectural decisions without re-extracting from the original PDFs.

## Documents

### Extractions of post-fetch papers (deep, equation-level)

- **[extraction_ma2019.md](extraction_ma2019.md)** — Ma et al. 2019, IJAE. *0D igniter chamber inverse model.* CRITICAL: Ma 2019 is an **inverse problem** (back-solves mdot_ig from measured P), not a forward predictor. Useful as: NASA-CEA 7-coefficient C_p(T) machinery, choked-throat closure form, multi-species mixing rules. NOT useful as: drop-in replacement for v0.6.0's igniter — v0.7.0 needs a forward 0D plenum (which Ma explicitly chose not to model).

- **[extraction_salita_wang_dagostino_2001.md](extraction_salita_wang_dagostino_2001.md)** — Three 2001 conference papers: Salita's framework introduction (AIAA 2001-3443), Wang's CFD outlook (AIAA 2001-3447), and d'Agostino's validated quasi-1D solver (AIAA 2001-3449). Salita's Baer/Ryan T_ign correlation (B≈5.5, C≈0.92) and recommendations on radiation handling. d'Agostino's sub-1% match on Ariane 4 & 5 with prescribed mdot_ig + Lenoir-Robillard erosive.

- **[extraction_peretz_pardue_cavallini.md](extraction_peretz_pardue_cavallini.md)** — Three foundational documents spanning 1973-2009: Peretz's foundational 1D HVT model, Pardue & Han's two-phase Shuttle extension, Cavallini's SPINBALL/SPIT (Vega program). The 1973→1992→2009 timeline of what stayed fixed and what evolved. **Peretz's Goodman cubic-polynomial integral method is the load-bearing kernel for v0.7.0** — see equations doc.

### Equation derivations

- **[equations_goodman_integral.md](equations_goodman_integral.md)** — Full derivation of the Goodman cubic-polynomial heat-balance integral method for per-cell solid-phase conduction. Reduces a per-cell PDE to a single ODE in penetration depth `δ(t)`. Numba-friendly, ~5% error vs. exact PDE. Derives DESIGN.md Eqs. 7 and 8 from first principles. Includes asymptotic checks against analytical limits.

### Primary-source summaries (textbook + amateur reference)

- **[primary_sources_summary.md](primary_sources_summary.md)** — Sutton 9e §14.2 + §15.3 (the textbook three-phase ignition framework + the empirical sizing equation `m = 0.12·V_F^0.7`) plus DeMar 2021 amateur experimental survey (`M = V·P/W` formula with measured impetus W for 8 pyrogen compositions). The Sutton + DeMar materials are the empirical / amateur-friendly grounding; the academic papers are the modeling depth.

## Source PDFs (in repo root)

All in `c:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/`:

| File | Reference |
|---|---|
| `Rocket Propulsion Elements.pdf` | Sutton 9e |
| `STARTING YOUR UP-GOERS_pub (1).pdf` | DeMar 2021 amateur deck |
| `International Journal of Aerospace Engineering - 2020 - Ma - A Model for Igniter Mass Flow Rate History Evaluation for.pdf` | Ma et al. 2019 |
| `salita2001.pdf` | Salita 2001 (AIAA 2001-3443) |
| `wang2001.pdf` | Wang 2001 (AIAA 2001-3447) |
| `10.2514@6.2001-3449.pdf` | d'Agostino, Biagioni, Lamberti 2001 (AIAA 2001-3449) |
| `19740005393.pdf` | Peretz, Caveny, Kuo, Summerfield 1973 (Princeton AMS-1100) |
| `pardue1992.pdf` | Pardue & Han 1992 (AIAA 92-3277) |
| `74323997.pdf` | Cavallini 2009 (Sapienza PhD thesis) |
| `International Journal of Aerospace Engineering - 2020 - Ma - A New Erosive Burning Model of Solid Propellant Based on Heat.pdf` | Ma 2020 (the erosive paper srm_1d already uses, NOT the igniter paper) |
| `hasegawa2006.pdf` | Hasegawa Motor A validation data |

## Reading order for a fresh agent

1. **Start here**: [../DESIGN.md](../DESIGN.md) — the v0.7.0 architecture decisions, with rationale.
2. **Implementation tasks**: [../TASKS.md](../TASKS.md) — concrete file-level breakdown.
3. **Equation depth**: [equations_goodman_integral.md](equations_goodman_integral.md) — derivation of the per-cell solid-phase conduction kernel. This is the single most equation-heavy doc and the only one you must understand fully to implement Phase 2.
4. **Architectural depth**: [extraction_peretz_pardue_cavallini.md](extraction_peretz_pardue_cavallini.md) — the 1973-2009 timeline of how the 1D ignition transient field evolved. Most relevant for understanding why we chose what we chose.
5. **Empirical grounding**: [primary_sources_summary.md](primary_sources_summary.md) — the Sutton sizing equation and DeMar's measured pyrogen impetus values. Essential for pyrogen YAML defaults.
6. **Modern alternatives considered and rejected**: [extraction_salita_wang_dagostino_2001.md](extraction_salita_wang_dagostino_2001.md) and [extraction_ma2019.md](extraction_ma2019.md). Read if you want to deviate from the chosen architecture; otherwise skim.

## Verifying source claims

If you encounter a claim in any extraction that seems wrong or under-supported, the source PDFs are in repo root. Use:

```bash
pdftotext -layout "<filename>.pdf" /tmp/extracted.txt
grep -n "<keyword>" /tmp/extracted.txt
```

For Cavallini's thesis (49 MB), pdftotext is mandatory — Read tool will fail.

For Peretz 1973, use `pdftotext -raw` (the default layout mode garbles the text due to scan layout issues).

## What's NOT in these references

- **Anything from Salita 2001 Parts 2-4 (AIAA 2001-3444, -3445, -3446)** — these were not available to fetch. The user has Parts 1 and 5 only.
- **JANNAF/SPP-12 IGT module documentation** — proprietary (SEA Inc. commercial).
- **Detailed Shuttle RSRM igniter mdot(t) traces** — sparse in open literature.
- **Multi-species coupling stability analysis** — not directly addressed in any reviewed paper.

If those gaps become load-bearing, surface them as PDF requests to the user (per session convention: don't quietly skip paywalled material).
