# srm_1d v0.7.2 — Ignition Model Rework Design Package

**Status**: design phase (2026-05-23). v0.7.1.1 just shipped; v0.7.2
implementation has not started. This folder collects the design docs
for the four candidate fixes to the structural ignition-kernel
artifact identified in the v0.7.1 Phase 5 close-out + the v0.7.1.1
cross-motor cleanup.

## Problem statement (the v0.7.2 target)

At default knobs the model produces an unphysical ignition-spike
artifact for every fired motor in the library:

| Motor       | Frozen P_peak | Effective P_peak | Experimental peak |
|-------------|---------------|------------------|---------------------|
| Hasegawa A  | 5.84 MPa      | 8.27 MPa         | 6.5 MPa @ t=1.1 s   |
| Zerox       | 7.85 MPa      | 10.20 MPa        | ~4.0 MPa @ t=0.2 s  |
| BALLSstick  | 9.33 MPa      | 14.48 MPa        | n/a in repo         |
| Chunc       | 13.14 MPa     | 20.27 MPa        | n/a in repo         |

Effective transport AMPLIFIES the artifact by +30-55% vs frozen across
all 4 motors. Gas-transport choice is not the lever; the **structural
cause is hypothesized to be the ignition kernel itself**: each bore
cell solves its Goodman cubic-polynomial pre-ignition heat-balance
integral independently, and at typical parameter settings all cells
cross `T_surf > T_ignition` within a few ms of each other →
simultaneous bore ignition → sharp pressure transient that the
experimental traces don't show.

Per the LHS path B finding (`hasegawa_a_lhs_effective.py` rank-1
fitness 0.1933 with k_solid at literature center 0.331 W/(m·K)), no
combination of knobs within physical-realism bounds simultaneously
matches the spike, plateau, and tail-off. The fix has to be
structural.

## Four candidate approaches under design-doc consideration

| # | Candidate | One-line summary | Implementation footprint |
|---|-----------|------------------|--------------------------|
| 1 | [Z-N dynamic burn rate](candidates/01_z_n_dynamic_burn_rate.md) | Relaxation ODE on steady-state Saint-Robert burn rate per cell; lag smooths pressure transients | Smallest — one new per-cell state + ODE step |
| 2 | [Spatial ignition-front coupling](candidates/02_spatial_ignition_front_coupling.md) | Couple adjacent unignited cells via mass-flux-driven heat transfer / lateral conduction / radiation so cells don't all ignite simultaneously | Medium — modifies Goodman call site, adds cell-to-cell coupling term |
| 3 | [Pyrogen spatial distribution refactor](candidates/03_pyrogen_spatial_distribution.md) | Replace cell-0-only pyrogen injection with an axial profile (decay length / jet-penetration model); spreads early energy over the plume's actual reach | Medium — adds axial profile to `igniter_plenum.py`; per-pyrogen parameters |
| 4 | [Submerged pyrogen modes](candidates/04_submerged_pyrogen_modes.md) | Support alternative igniter topologies: head-end basket inside grain core (a); nozzle-inserted cartridge like Super Loki (b). Enables ISP Super Loki validation | Largest — new boundary conditions / source-term placements, new YAML schema for igniter topology |

These are **not mutually exclusive**. Z-N (1) is fundamentally a
burn-rate physics fix; (2) and (3) and (4) are ignition-stage fixes.
A v0.7.2 release could ship any subset; the natural pairing is (1) +
one of (2)/(3)/(4) based on which best matches the cross-motor
diagnostic data.

## Decision criteria

Per the data we have, the right candidate should:

1. **Reduce the at-default-knobs spike-to-plateau ratio for all 4
   fired motors** toward physically reasonable values (~1.0-1.5
   based on experimental traces of fired SRMs).
2. **Preserve the v0.7.1 effective-transport plateau / erosive-peak
   shape match** for Hasegawa A — that part is good and should not
   regress.
3. **Have published validation precedent** in 1D SRM internal
   ballistics codes (not invented from scratch).
4. **Stay within srm_1d's pure-numerics solver / data-defined motor
   philosophy** (see CLAUDE.md). Knobs should map to physical
   quantities with literature ranges, not free fit parameters.
5. **Reasonable implementation cost** — v0.7.2 should ship within a
   few sessions, not become a multi-month rewrite.

## Decision process

After all 4 design docs land + literature digests are filed:
1. Score each candidate against the 5 criteria above
2. Pick 1-2 to implement for v0.7.2
3. Defer the others to v0.7.3+ or document as alternatives
4. Begin implementation on the chosen candidate(s); create
   `srm_1d/docs/v0_7_2/TASKS.md` with phase breakdown

## Relevant repo context

- `srm_1d/solid_thermal.py` — current per-cell Goodman solver
  (independent cells; the model that all 4 candidates modify)
- `srm_1d/igniter_plenum.py` — current forward-plenum pyrogen
  injection into cell 0 (the model candidates 3 and 4 modify)
- `srm_1d/simulation.py` — `_run_time_loop`; the integration point
  for any per-cell ODE state (candidate 1) or coupling term (2)
- `srm_1d/burn_rate.py` — Ma 2020 erosive-burning chain; Z-N (1)
  would wrap the steady-state r₀(P) lookup here
- `srm_1d/docs/post_v0_7_0/references/spinball_walkthrough.md` —
  prior literature decision doc; flagged Z-N as the spike-taildown
  candidate
- `srm_1d/docs/post_v0_7_1/references/` — pyrogen heat flux,
  AP/HTPB k_solid, copper thermite literature (extracted during
  Phase 5)

## Related memories

- `[[project_v0_7_1_phase5_task1_task2_findings]]` — the diagnostic
  data this design package responds to.
- `[[project_hasegawa_calibration_state]]` — Hasegawa A calibration
  state to preserve.
- `[[project_fired_motor_set]]` — validation set, now includes ISP
  Super Loki for the candidate 4 design.
- `[[project_spinball_research_state]]` — earlier SPINBALL extract
  that motivated Z-N as the primary candidate.
- `[[roughness-kappa-physical-bounds]]` — physical-realism bounds to
  preserve on any new knobs introduced by v0.7.2.
