# Post-v0.7.1 References

Research material supporting work beyond the v0.7.1 N-species
release. v0.7.1 introduced per-cell γ/Cp/R, sensible-enthalpy
advection, per-species Cp at source sites, and strict per-cell
T_ceiling. Phase 5 LHS re-calibrated Hasegawa Motor A against the
new physics; this folder hosts the literature that informed that
re-calibration and the candidates queued for v0.7.2+.

## What's here

| Document | One-line purpose |
|---|---|
| [pyrogen_heat_flux_literature.md](pyrogen_heat_flux_literature.md) | DeMar / Sandia LDRD / AIAA pyrogen heat-flux bounds. **The "DeMar 69.4 is time-averaged, not peak" finding** that reframes elevated LHS heat-flux values (e.g. v0.7.1 Phase 5 rank-1 at 232 cal/cm²/s ≈ 9.7 MW/m²) as defensible peak-transient interpretations rather than missing-physics compensation. Cites Sandia LDRD 2022, MDPI 2020, AIAA JSR. |
| [copper_thermite_igniter_literature.md](copper_thermite_igniter_literature.md) | CuO/Al thermite chemistry + heat delivery + ignition mechanism. **Architectural mismatch with v0.7.0+ 0D pyrogen plenum** — would need a surface-flux-only boundary path with no mass injection. Deferred to v0.8. Cites Reese 2015 (AIAA JPP, gold standard), Nishii 2024 (only published 1D thermite + grain model), Nakka, US Patent 4,464,989. |

## What's NOT here

Items the post-v0.7.1 research touched on but did not include in this folder:

- **Z-N dynamic burn rate**: still queued (covered in
  [`../../post_v0_7_0/references/spinball_walkthrough.md`](../../post_v0_7_0/references/spinball_walkthrough.md)
  and recorded in the `project_spinball_research_state` memory).
  The v0.7.1 N-species infrastructure is the prerequisite; Z-N
  becomes a v0.7.2+ candidate once Phase 5 calibration is final.
- **Per-cell transport** (k_thermal, μ_gas, Pr) frozen-vs-effective:
  flagged by the user as the next investigation. Not yet researched.
- **Time-varying pyrogen heat-flux profile**: surfaced as an open
  improvement in the BPNV heat-flux digest; deferred to v0.7.2+ if
  Phase 5 calibration ergonomics suggest the scalar-flux limitation
  is binding.

## Provenance

Both digests in this folder originate from haiku-model general-purpose
subagent runs during the v0.7.1 Phase 5 session on 2026-05-22. The
agents performed open WebSearch / WebFetch literature dives; their
raw JSONL transcripts were transient and are not preserved. The
digests here include URLs for all primary citations so future
sessions can re-fetch the sources directly if needed.

## Relation to other reference folders

- [`../../post_v0_7_0/references/`](../../post_v0_7_0/references/) —
  the SPINBALL research that motivated v0.7.1 itself (N-species
  + Z-N candidacy). Cross-references with `[[name]]` link both ways.
- [`../../v0_7_1/`](../../v0_7_1/) — design and tasks docs for the
  v0.7.1 work itself.

## Related memories

- `[[srm-1d-pyrogen-heat-flux-literature-bounds]]` — quick-lookup
  memory equivalent of pyrogen_heat_flux_literature.md.
- `[[srm-1d-copper-thermite-igniter-research]]` — quick-lookup
  memory equivalent of copper_thermite_igniter_literature.md.
- `[[srm-1d-v0-7-1-progress-state]]` — v0.7.1 implementation arc.
