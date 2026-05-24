# Submerged Pyrogen Modes — Extended Literature Digest

*Subagent literature dive, 2026-05-23. Extends what srm_1d already had
documented (forward-plenum-only `igniter_plenum.py`,
`docs/v0_7_0/DESIGN.md`, and the post-v0.7.1 igniter literature digests).*

## What srm_1d already has documented

v0.7.0+ ships **one** igniter topology: a forward 0D pyrogen plenum
(`igniter_plenum.py`) with a Saint-Robert burn law and a choked/
subsonic orifice that injects mass, sensible enthalpy, and (optionally)
axial momentum into bore cell 0. Pyrogen formulations covered are
BPNV / BKNO3-V / MTV-class via DeMar 1995 heat-flux calibration
(`pyrogen_heat_flux_literature.md`). Surface-flux-only thermites are
scoped for v0.8 in a separate boundary-condition path
(`copper_thermite_igniter_literature.md`). The DESIGN doc explicitly
defers "squib stage" and "non-forward igniter geometries" — no model
exists for submerged baskets, nozzle-inserted cartridges, or aft-end
igniters.

## New findings from this dive

**Industry/military topology taxonomy.** Open-literature SRM igniter
topologies fall into four buckets that map cleanly to distinct 1D
boundary treatments:
1. *Forward pyrogen* — what srm_1d already models (Sutton 9e §15,
   NASA SP-8051).
2. *Head-end basket of BKNO3/BPN pellets* — the Apogee/ProCast amateur
   path and the canonical "pre-pressurize fast" military approach.
   Topologically equivalent to a forward plenum but with effectively
   zero residence time and very high A_burn/V_plenum, so the existing
   module covers it with parameter changes only.
3. *Spherical/internal pyrogen suspended in the head-end bore*
   (Thiokol US5150654A) — a small SRM that lives *inside* the main
   grain bore, with discharge ports oriented forward-out. This is the
   canonical "submerged head-end pyrogen" — burns inside the main bore
   volume.
4. *Aft-end / nozzle-cavity igniter* (US4023497A, US4503773A) —
   neutral-thrust pressure vessel removably bonded into an aft-grain
   cavity, with discharge ports oriented forward into the bore. The
   cylindrical-grain "aft pyrogen" is *not* a cartridge that ejects;
   it's a consumable filament-wound pressure vessel that burns out in
   place. ASRM (NASA NTRS 19950017219) used a 12-port multi-port
   igniter that vented hot gas radially into circumferential
   propellant slots.

**Super Loki / amateur "cartridge-in-nozzle".** The MIT/RRS Super Loki
Dart report and Loki Research reload docs (38mm Red, I110) confirm
the amateur/sounding-rocket convention is *not* a cartridge that sits
inside the nozzle and ejects — it's a thermalite or e-match-plus-
pyrogen string fed *through* the nozzle from outside the pad, with
the active pyrogen sitting up at the head end of the bore. The
"ejection" is just the unburned igniter housing falling out the
nozzle once flow starts; throat area is not significantly modulated.
John Coker's *Forward Igniter* page (jcrocket.com) documents the
inverse case: a head-end initiator screwed into the forward closure,
with a secondary pyrogen charge placed *inside* the top grain bore —
explicitly a submerged-pyrogen amateur topology.

**Boundary-condition formulations in published 1D codes.** SPINBALL
(Cavallini/Favini/Di Giacinto, AIAA 2009-5512 and Cavallini PhD
thesis 2014) treats the igniter as a **distributed Q1D source term**:
mass, energy, and (optionally) axial momentum are added across the
cells the igniter physically occupies, not as a head-end boundary
condition. The recent ScienceDirect 2025 parametric study
(S1290072925006453) systematically explored axial-port vs canted
multi-port igniter discharge with full 3D CFD and confirms that
**jet Mach number, stagnation T/P, and injection angle dominate
induction-interval physics** — momentum direction is first-order, not
a tuning detail.

**Cartridge-in-nozzle dynamics in real designs.** The "true"
cartridge-in-throat case (US4751881, US5007236, US5062206 "removable
rocket motor igniter") *is* used on military liquid engines and some
pyrotechnic motors. These igniters have rupture-pressure positioning
rings; the throat is partially blocked during the first 5-50 ms and
the cartridge then ejects. There is no published open closed-form
transient-throat model, only patent prose. For amateur SRMs this case
is rare — the cartridge mass is small relative to the throat area.

## Implementation guidance for srm_1d

**Topology (3) — submerged head-end pyrogen basket.** The lowest-cost
extension. Reuse the existing `PyrogenChamber` ODE wholesale, but
allow `V_plenum` to refer to a *bore-internal* volume and add an
`injection_cell_range = (i_start, i_end)` parameter. Replace the
single-cell mass/enthalpy/momentum sink at cell 0 with a
partition-of-unity hat function (already used for end-face injection)
distributed across `injection_cell_range`. The choked orifice still
couples plenum P to a representative bore P (volume-weighted average
of the spanned cells). User-supplied inputs: `axial_extent_m`,
`i_start_fraction`, `port_orientation` (forward / radial / aft). **One
day of work.**

**Topology (4) — aft-end / nozzle-cavity igniter.** Use the same
distributed-source machinery as (3), but locate `injection_cell_range`
near `cell N-1` and *reverse the momentum sign* (forward-firing jet
means negative axial momentum into the bore). The aft pyrogen
pressure-feedback must close on cell N-1 bore pressure, not cell 0.
This also changes the nozzle BC handling: until the cartridge burns
out, the effective throat area is
`A_throat - A_cartridge_blockage(t)`. Add a `throat_blockage_history(t)`
callable. **Two to three days of work.**

**Topology (2) — pre-pressurization basket.** No code changes needed.
Document a recipe: set `V_plenum` very small, `A_burn_initial` very
large, `burn_law = "0d"`. Falls out of existing model.

## Known limitations

- No open-literature instantaneous heat-flux data exists for ISP
  Super Loki — only chamber-averaged thrust and burn-time. Cannot
  validate ignition-spike shape against Super Loki without
  proprietary Sandia/ISP test data.
- Cartridge ejection dynamics (US4751881 class) require a 2-phase
  model (constrained cartridge → free-flight in nozzle exit) that
  1D PISO cannot capture cleanly; treat as a step-change in
  `A_throat` at a user-specified ejection pressure.
- ASRM-style canted multi-port igniters inject momentum at angles to
  the bore axis. A 1D code can only model the axial component;
  circumferential flame-spread benefits are lost. Document as a
  known fidelity limit.
- Thermite-style igniters (Reese 2015) need a separate surface-flux-
  only path (already scoped for v0.8); the distributed-source
  machinery here doesn't apply.
- The "no published instantaneous transient-throat-area model" gap
  means cartridge-in-nozzle support will be patent-prose-derived,
  not literature-calibrated. Recommend deferring (4) until a user
  supplies experimental pressure trace data.

## Primary citations (with URLs)

1. **US Patent 5,150,654A** (Thiokol, 1992) — Spherical igniter for
   full head-end web rocket motors. The canonical "submerged
   head-end pyrogen" patent. https://patents.google.com/patent/US5150654A/en
2. **US Patent 4,503,773A** (1985) — Aft-end igniter for full
   head-end web SRMs. Consumable filament-wound pressure vessel
   bonded into aft-grain cavity.
   https://patents.google.com/patent/US4503773A/en
3. **US Patent 4,023,497A** (1977) — Aft-end ignition system for
   rocket motor. Earlier aft-end design.
   https://patents.google.com/patent/US4023497A/en
4. **US Patent 4,751,881** — Igniter capable of being fitted in the
   nozzle of a propulsion unit. Cartridge-in-nozzle baseline patent.
   https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/4751881
5. **NASA SP-8051** — Solid Rocket Motor Igniters (NASA monograph).
   Canonical taxonomy.
   http://www.vibrationdata.com/tutorials_alt/NASA_SP8051.pdf
6. **NTRS 19950017219** — ASRM Multi-Port Igniter Flow Field
   Analysis. 12-port igniter venting into circumferential propellant
   slots. https://ntrs.nasa.gov/citations/19950017219
7. **NASA CR-61238 (1968)** — Super Loki Dart Meteorological Rocket
   System. The Super Loki motor + igniter system description.
   https://www.rrs.org/wp-content/uploads/2014/01/Super-Loki-Dart-Meteorological-Rocket-System-1968.pdf
8. **MIT Super Loki Report** — companion technical description with
   motor cross-section and igniter installation.
   https://wikis.mit.edu/confluence/download/attachments/122633025/SuperLokiReport.pdf
9. **Cavallini et al., AIAA 2009-5512** — SRM Internal Ballistic
   Numerical Simulation by SPINBALL Model. Distributed igniter
   source-term formulation.
   https://arc.aiaa.org/doi/10.2514/6.2009-5512
10. **ScienceDirect S1290072925006453 (2025)** — Parametric study of
    igniter design on ignition transient performance. Axial-port vs
    canted multi-port; jet Mach/T/P/angle sensitivity.
    https://www.sciencedirect.com/science/article/abs/pii/S1290072925006453
11. **John Coker — Forward Igniter** (amateur reference). Head-end
    initiator with secondary charge submerged inside top grain bore.
    http://jcrocket.com/forward-igniter.shtml
12. **Richard Nakka — Igniter Systems**. Authoritative amateur design
    reference for pyrogen vs basket igniters.
    https://www.nakka-rocketry.net/igniter.html
