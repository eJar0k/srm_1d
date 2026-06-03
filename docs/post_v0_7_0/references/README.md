# Post-v0.7.0 References

Research material supporting the next physical-model iteration. v0.7.0 is tagged and frozen; this folder hosts work toward the spike-taildown residual that motivates v0.7.0.x / v0.7.1.

## What's here

| Document | One-line purpose |
|---|---|
| [spinball_walkthrough.md](spinball_walkthrough.md) | **Start here.** SPINBALL vs srm_1d v0.7.0 mechanism-by-mechanism, with recommendation among Z-N / Al2O3 two-phase / Peretz participation / per-cell variable γ. |
| [extraction_spinball_2009.md](extraction_spinball_2009.md) | Cavallini, Favini, Di Giacinto, Serraglia (AIAA 2009-5512). Canonical SPINBALL equation set (Eq. 1), SPIT-vs-SPINBALL contrast (Eq. 1 vs Eq. 2), Z23 results. Contains §10 correction to the v0.7.0 thesis extraction. |
| [extraction_spit_to_spinball_2008.md](extraction_spit_to_spinball_2008.md) | Favini, Cavallini, Di Giacinto, Serraglia (AIAA 2008-5141). The SPIT→SPINBALL transition; GREG burnback module; Z9 entire-burn validation. |

PDFs are in the parent directory (`srm_1d/docs/post_v0_7_0/`):
- `cavallini2009.pdf` — AIAA 2009-5512 (SPINBALL)
- `digiacinto2008.pdf` — AIAA 2008-5141 (SPIT→SPINBALL + GREG)

## What's NOT here

The SPIT-lineage equations (igniter sub-model, ignition criterion) deferred by both Sapienza papers to:
- Di Giacinto & Serraglia 2001, AIAA 2001-3448
- Serraglia 2003 PhD thesis
- Favini et al. 2005, EUCASS

These are not required for the current spike-taildown research direction — see [spinball_walkthrough.md](spinball_walkthrough.md) §6.

## Relation to v0.7.0 references

The frozen v0.7.0 reference set is in [`../../v0_7_0/references/`](../../v0_7_0/references/). This `post_v0_7_0` folder is the next-iteration parallel. Cross-references with `[[name]]` link both ways.

The v0.7.0 extraction `extraction_peretz_pardue_cavallini.md` contains the 2009 Cavallini thesis summary; the §3.1 description of "6 mass-conservation equations" is superseded by the 2009 conference paper formulation (Eq. 1: ONE mass equation). Correction is documented in [extraction_spinball_2009.md §10](extraction_spinball_2009.md).
