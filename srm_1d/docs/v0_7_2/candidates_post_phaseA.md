# srm_1d v0.7.3 Candidate Breakdown — post-Phase-B analysis

**Context**: v0.7.2-phaseA shipped (commit `e507c09`, 2026-05-24)
with Phase A (pyrogen axial distribution) as the load-bearing
contribution and Phase B (cumulative-G v1 + flame-front v2)
disabled-by-default after both empirically amplified rather than
smoothed the simultaneous-ignition spike. The Phase B negative
finding clarifies what the artifact ISN'T (a missing-augmentation
problem like the literature describes), which narrows the candidate
space for v0.7.3.

This doc analyzes each candidate the user flagged as a plausible
next-direction angle. Not a scope decision — that's pending the
user's read of this breakdown. Each section follows the same
structure: physical motivation, implementation sketch, what it
addresses, what it WON'T address, implementation cost, dependencies.

---

## 1. Z-N dynamic burn rate (original candidate 1)

### Physical motivation

Solid propellant pyrolysis has a transient lag (~ms scale) between
pressure changes and the burn rate equilibrating to the new
steady-state value. The current Saint-Robert + Ma erosive chain
treats `r_b(P, Re)` as instantaneous — any pressure spike that the
PISO solver computes immediately drives a matching burn-rate spike.
A Z-N relaxation ODE per cell with `τ_ZN = κ · α / r²` adds the
physical lag, blunting transient pressure spikes by preventing
instantaneous burn-rate response.

### Implementation sketch

Per cell `i`, advance a dynamic burn rate state:
```
dr_dyn[i]/dt = (r_steady[i] - r_dyn[i]) / τ_ZN[i]
τ_ZN[i] = κ · α_solid / max(r_dyn[i], r_floor)²
```
where `r_steady[i]` is the existing Ma 2020 erosive-augmented burn
rate and `α_solid = k_solid / (ρ_propellant · Cps)`. The mass-source
and regression-rate paths use `r_dyn[i]` instead of `r_steady[i]`.

Cleanest insertion: new Numba kernel `_advance_zn_burn_rate(...)`
called once per step inside `_run_time_loop` after the existing
`compute_burn_rates(...)` block; replaces direct use of `r_total`
in mass / regression paths with `r_dyn`. New per-cell state array
`r_dyn[N]` allocated at sim init. New knob `kappa_zn: float = 1.0`
on Propellant (no fit constant if κ ≈ 1 per Greatrix 2008).

### What it addresses

- **Burn-rate ramp-up rate during ignition transient**: each cell's
  r_b lags r_steady by τ_ZN, smoothing the leading edge of the
  pressure spike even when ignition timing is unchanged.
- **Pressure dynamics during throat-erosion recovery and end-of-burn
  knee**: same lag mechanism, different transient.

### What it WON'T address

- **The simultaneous-ignition event itself** — Z-N affects post-
  ignition burn rate, not when cells cross `T_surf > T_ignition`.
  If cells all ignite within a few ms (current artifact), Z-N
  smooths each cell's rate ramp but the cells still all start
  ramping near-simultaneously. Net effect: spike is shorter and
  rounder, but its existence and rough magnitude carry over.

### Implementation cost

Smallest of the candidates. ~150-250 LOC across simulation.py +
propellant.py + ~5-8 new tests. 1-2 sessions including
re-calibration against Hasegawa A baseline.

### Dependencies / stacking

Independent — works alongside Phase A and any igniter-topology
candidate. **Best-case stack**: Z-N + 4b (aft-inserted) would
attack the artifact from both ends: 4b breaks up simultaneous
ignition via topology; Z-N smooths whatever spike remains.

---

## 2. Submerged pyrogen — head-end basket (candidate 4a)

### Physical motivation

Industry-standard topology for medium and large SRMs (Thiokol
US5150654A; ASRM NTRS 19950017219). Pyrogen burns INSIDE the main
bore volume rather than upstream of it. Mass and enthalpy enter
through ports oriented forward and/or radially, distributed across
the first 5-15% of bore length. Provides faster pressurization
than a forward plenum and (potentially) more uniform initial bore
heating.

### Implementation sketch

Extends existing `PyrogenChamber` ODE machinery rather than
replacing it. Add `injection_topology: str = 'forward'` field;
new value `'submerged_head'` plus `axial_extent_m` and
`i_start_fraction` parameters. Pyrogen-chamber pressure-coupling
closes against volume-weighted average bore pressure across the
injection cell range. Mass / enthalpy injection reuses Phase A's
`_compute_pyrogen_axial_weights` machinery (top-hat instead of
exponential decay, but same kernel pattern).

### What it addresses

- **Initial bore-heating distribution**: energy spreads naturally
  across multiple cells from t=0, not just cell 0. Cells receive
  pyrogen mass simultaneously across the basket extent rather
  than relying on Phase A's axial-decay model.
- **Plenum-as-option foundation**: lays groundwork for the unified
  igniter-architecture API (candidate 6).

### What it WON'T address

- **The Goodman per-cell solver's tendency to push all cells past
  T_ignition near-simultaneously** — if the artifact lives in the
  ignition kernel itself (rather than the energy-distribution
  model), spreading the energy across more cells via 4a may
  actually MAKE IT WORSE by giving more cells a fast ramp.
- **The fundamental simultaneous-ignition artifact**: 4a changes
  the spatial energy profile but the cells in the basket extent
  still see similar local h_c and ignite in a tight time window.

### Implementation cost

Medium. ~200-300 LOC across `igniter_plenum.py` +
`simulation.py` + new Propellant/Pyrogen YAML schema + ~6-8 tests.
1-2 sessions.

### Dependencies / stacking

Depends on Phase A's axial-weight kernel infrastructure (already
shipped). Stacks with Z-N (the lag still applies regardless of
where energy enters). Stacks with 4b — same architecture, different
injection range.

---

## 3. Submerged pyrogen — aft-inserted impinging cartridge (candidate 4b)

### Physical motivation

Super Loki / ISP Corporation class. A short propellant cartridge
sits in the aft of the grain bore (or partially in the nozzle
throat) and fires FORWARD. The pyrogen plume impinges on grain
surfaces along the cartridge's axial extent — the burn propagates
back→front rather than front→back. Common in amateur high-power
designs (e-match + thermalite-on-PVC strings sit in the nozzle).
**User-flagged as a clean test of whether the simultaneous-ignition
artifact is caused by current concentrated head-end mass injection**:
if the artifact persists with mass-injection topology completely
reversed, the artifact lives in the Goodman/h_c kernel rather than
in the pyrogen model.

### Implementation sketch

Reuses the same `PyrogenChamber` + axial-weight machinery as 4a
but with:
- `injection_topology = 'aft_cavity'`
- Injection cell range = `(i_aft_start, N-1)`. **Cartridge length
  can be either user-specified (`axial_extent_m`) or pyrogen-mass-
  defined**: derive extent from pyrogen mass and cross-section
  geometry (`extent = m_pyrogen / (rho_p · A_port[aft])` for a
  full-cross-section plug).
- Momentum sign reverses: `momentum_sign = -1.0` (forward-firing
  jet from aft cartridge contributes negative axial momentum
  into the bore).
- Pressure-feedback closes on `P[i_aft_start..N-1]` volume-weighted
  average bore pressure (cell N-1 region), not P[0].
- **Optional transient throat blockage**: until cartridge ejects
  (user-specified ejection pressure or time), effective throat
  area = `A_throat - A_cartridge_blockage`. Step-change at
  ejection.

### What it addresses

- **The simultaneous-ignition diagnostic question**: if the spike
  artifact persists under reversed topology, the artifact is in
  the per-cell ignition kernel. If the spike disappears or
  qualitatively changes, the artifact was driven by head-end
  energy concentration. **Either outcome is informative**.
- **Super Loki / amateur SRM validation set expansion**: ISP
  Super Loki .ric is in the repo with experimental data already
  embedded (commented-out in `examples/ISP_Super_Loki.py`).
  Activating that test case unlocks a 5th fired-motor data point
  with a fundamentally different igniter topology.

### What it WON'T address

- **The same caveat as 4a**: if the Goodman per-cell solver
  produces near-simultaneous ignition under any reasonable
  energy distribution (uniform-h_c, similar surface temperatures
  pre-ignition), the artifact remains.
- **The Hasegawa A example trace specifically**: Hasegawa A has
  a forward-plenum igniter per the .ric, not an aft cartridge.
  4b would only validate against Super Loki (and any other
  motors the user wants to add aft cartridges to).

### Implementation cost

Largest of the post-A candidates if cartridge-ejection dynamics
are included. ~300-450 LOC across `igniter_plenum.py` +
`simulation.py` + `nozzle.py` (transient throat blockage) +
YAML schema + ~8-10 tests. **2-3 sessions** for the topology +
ISP Super Loki experimental validation.

### Dependencies / stacking

Reuses Phase A's axial-weight machinery and 4a's
`injection_topology` enum (build 4a first → 4b becomes much
cheaper). Stacks with Z-N. **Conceptually orthogonal to per-cell
coupling alternatives** — 4b tests a topology hypothesis; if it
fails to fix the artifact, that's evidence the per-cell coupling
work is the right next step.

---

## 4. Per-cell coupling alternatives (post-Phase-B redesign)

### Physical motivation

Phase B-v1 (cumulative-G) and Phase B-v2 (flame-front gating)
both showed that POSITIVE h_c augmentation amplifies the cascade
because PISO already captures upstream contributions to local h_c.
The complementary direction is **NEGATIVE augmentation** (damping):
reduce h_c at cells far from any recently-ignited cell, leaving
only cells immediately adjacent to the burning front at full
local h_c. Forces sequential propagation by SLOWING distant cells
rather than ACCELERATING adjacent ones.

Other coupling mechanisms worth considering:
- **Solid-phase axial conduction**: heat propagates along the
  propellant surface from burning cells to unignited neighbors
  via the solid itself (current per-cell Goodman treats cells as
  axially isolated).
- **Goodman per-cell shared boundary layer**: the boundary layer
  thickness `δ` at adjacent cells should be related; current
  per-cell Goodman lets each cell's δ evolve independently.

### Implementation sketch

**Sub-option 4A (reverse-polarity damping)**: minimal change to
existing Phase B infrastructure. Replace `flame_spread_boost = 3.0`
with `flame_spread_damp = 0.3` (or similar < 1.0) and **invert the
gating logic**: damp h_c at cells NOT downstream of a recently-
ignited cell. The cell that IS the immediate neighbor of the
burning front keeps full local h_c (factor = 1.0); all other
unignited cells get factor < 1.0. Same kernel
(`_compute_flame_front_augment`) with different output values.

**Sub-option 4B (solid-phase axial conduction)**: add a per-cell
axial conduction term to the Goodman ODE:
```
∂T_surf/∂t = (existing Goodman RHS) + k_solid / (ρ_p · Cps) · ∂²T_surf/∂x²
```
where the second-derivative is approximated by central differences
on the cell-centered T_surf array. Adds modest computational cost
(O(N) per step) and one boundary-condition decision at the ends.

**Sub-option 4C (boundary-layer shared)**: extend Goodman's δ
state to couple `δ[i]` to `δ[i±1]` via a small lateral-conduction
term, analogous to 4B but on δ instead of T_surf. More complex,
less immediately defensible.

### What it addresses

- **The simultaneous-ignition artifact directly** (reverse-polarity
  damping in particular): slows distant cells, forces front
  propagation timescale to be set by the front itself rather than
  by global pyrogen-driven uniform heating.

### What it WON'T address

- **Hasegawa A pyrogen-driven spike specifically**: if the spike
  is dominated by the cell-0-pyrogen-flux mechanism (Goodman
  surface heating), per-cell coupling between unignited cells
  downstream may not change the cell-0 ignition timing.
- **Burn-rate dynamics post-ignition** (Z-N's territory).

### Implementation cost

- Sub-option 4A (reverse-polarity damping): **smallest** — one
  numeric knob flip on existing Phase B machinery, ~30 LOC + 3-5
  tests. ~0.5 session.
- Sub-option 4B (solid-phase axial conduction): medium — adds
  one new term to Goodman + per-cell discretization + BC choice.
  ~150-200 LOC. 1-2 sessions.
- Sub-option 4C (boundary-layer shared): largest — modifies
  Goodman's internal state coupling. ~250 LOC + careful
  validation. 2+ sessions. Likely diminishing returns vs 4B.

### Dependencies / stacking

Sub-option 4A reuses Phase B infrastructure directly (just flip
the polarity). 4B/4C are independent of Phase B and stack with
Z-N. **Promising stack**: 4A + 4b (aft-inserted) — different
igniter topology + reverse-polarity coupling to slow distant
cells.

---

## 5. Different heating modes

### Physical motivation

The current model uses three heat-transfer mechanisms to drive
surface temperature: Bartz convective h_c (local, dominant);
DeMar pyrogen surface heat flux (cell 0 only, transient); and
optional `radiation_emissivity` (adjacent-cell radiation, already
wired). Additional mechanisms from the literature:

- **Pardue 1992 Al2O3 two-phase condensation**: in aluminized
  propellants, molten Al2O3 droplets condense on cooler bore
  surfaces, releasing latent heat that preferentially heats
  UNIGNITED downstream cells. This is the physical mechanism
  flagged in the SPINBALL walkthrough as the secondary spike-
  taildown candidate (after Z-N).
- **Enhanced radiation at distance**: current radiation path
  treats only immediate neighbors (`radiation_emitter[i±1]`).
  Real SRM bore radiation has line-of-sight to many cells; a
  multi-cell radiation matrix would couple cells over longer
  distances.

### Implementation sketch

**Pardue 1992 Al2O3**: new propellant property `aluminum_fraction:
float = 0.0` (Hasegawa A has 14%). Compute Al2O3 mass production
rate per burning cell; transport along bore via PISO (need new
species index for Al2O3 droplets); deposit on cells where local
T_surf < T_Al2O3_melt (≈2327 K). Latent heat of fusion adds to
surface heating. Substantial implementation: ~400 LOC including
N-species infrastructure extension to 4 species (igniter + grain +
ambient + Al2O3 condensable).

**Multi-cell radiation matrix**: extend existing emitter/receiver
infrastructure. For each unignited cell j, sum radiation
contributions from all `radiation_emitter[i]` cells with line-of-
sight (in 1D, all cells in the bore). View-factor approximation:
F_ij = 1 / (4·(x_j - x_i)² + 1) — empirical. ~150 LOC + tests.

### What it addresses

- **Pardue Al2O3** specifically the Hasegawa A late-tail-off shape
  mismatch (the SPINBALL walkthrough's documented motivation).
  Also potentially the simultaneous-ignition artifact IF the
  condensation heating mechanism preferentially heats downstream
  cells while upstream cells are too hot for condensation.
- **Multi-cell radiation** — propagates heat further axially than
  current 2-cell-radius radiation model. May help propagate
  ignition front.

### What it WON'T address

- Spike-magnitude mismatch from the pyrogen-driven ignition
  transient (which is faster than radiation/condensation
  timescales).

### Implementation cost

Substantial. Pardue alone is comparable to v0.7.1 N-species
infrastructure scope (~2-3 sessions including validation).
Multi-cell radiation is more contained (~1 session). Both
require empirical/literature parameter calibration that's harder
to physically defend than Z-N's `κ ≈ 1`.

### Dependencies / stacking

Pardue 1992 ALMOST requires the N-species infrastructure to be
extended (4th species for Al2O3 droplets). Stacks with Z-N and
all topology candidates.

---

## 6. Plenum-as-option refactor (architecture unification)

### Physical motivation

User's stated long-term goal: forward-plenum, head-end basket, and
aft-inserted cartridge should all be SELECTABLE per motor via
YAML configuration rather than baked into separate code paths.
Foundation for general igniter architecture covering current and
future topologies (e.g., distributed-port igniters, multi-port
canted designs).

### Implementation sketch

Unify under a single `Igniter` interface:
```python
@dataclass
class Igniter:
    topology: str         # 'forward' | 'submerged_head' | 'aft_cavity'
    pyrogen: Pyrogen      # composition / kinetics
    m_pyrogen: float      # total mass
    chamber_geometry: dict  # topology-specific
    injection_distribution: 'AxialDistribution'  # weight kernel + range
    momentum_sign: float = 1.0
```
Each topology supplies its own pressure-coupling rule
(forward → P[0]; submerged_head → volume-avg P[range];
aft_cavity → P[N-1]) and momentum delivery face (forward →
face 1; submerged_head → all faces in range; aft_cavity →
last face).

Motor YAML grows an `igniter:` block with `topology` and
topology-specific parameters.

### What it addresses

- **Long-term maintainability**: replaces a growing branch tree
  in `_run_time_loop` (one for each topology) with a polymorphic
  dispatch table.
- **Cross-motor validation framework expansion**: enables
  per-motor igniter topology selection, unlocking systematic
  cross-topology studies (e.g., "how does Hasegawa A respond if
  we model it as an aft-cartridge?").

### What it WON'T address

- **The structural ignition artifact directly** — pure refactor,
  no new physics.

### Implementation cost

Substantial but mostly mechanical. ~400-500 LOC reorganizing
`igniter_plenum.py` + `simulation.py` + new `Igniter` dataclass
+ YAML schema + migration of existing motor YAMLs + tests.
2-3 sessions.

### Dependencies / stacking

**Best done AFTER 4a (head-end basket) and 4b (aft-inserted)
ship** as topology-specific implementations; the refactor then
extracts the common pattern. Doing it first risks designing for
hypothetical future use cases without empirical guidance.

---

## Stacking matrix (which can run together)

|     | Z-N | 4a head-end | 4b aft-insert | 4-coupling | 5-heating | 6-refactor |
|-----|-----|-------------|---------------|------------|-----------|------------|
| **Z-N**       | — | ✓ | ✓ | ✓ | ✓ | ✓ |
| **4a basket** | ✓ | — | ✓ (separate motors) | ✓ | ✓ | (build 4a first, then refactor) |
| **4b aft-insert** | ✓ | ✓ | — | ✓ | ✓ | (build 4b first, then refactor) |
| **4 coupling**    | ✓ | ✓ | ✓ | — | ✓ | ✓ |
| **5 heating**     | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| **6 refactor**    | ✓ | (after) | (after) | ✓ | ✓ | — |

All candidates stack cleanly with Z-N; (6) refactor is best done
after at least one new topology lands.

## Recommended attack ordering

If user wants to pursue MULTIPLE candidates, my read of the
trade-offs:

1. **First**: candidate **4b (aft-inserted impinging cartridge)**.
   The user flagged this as the cleanest diagnostic test —
   reversing the mass-injection topology will produce either a
   qualitatively-different trace (artifact was in pyrogen-source
   model) or essentially the same artifact (artifact is in
   per-cell Goodman). Either outcome informs the next step
   decisively. Bonus: unlocks ISP Super Loki validation.

2. **Second**: depending on 4b's outcome —
   - If 4b's trace is meaningfully different (topology DID drive
     the artifact): proceed to **4a (head-end basket)** for
     architectural completeness + maybe **6 (refactor)** to
     unify the three topologies.
   - If 4b's trace shows the same artifact: the per-cell Goodman
     coupling is the load-bearing fix. Proceed to **4 sub-option
     A (reverse-polarity damping)** — cheapest experiment using
     existing Phase B machinery. Then **Z-N** to smooth whatever
     spike remains.

3. **Third (parallel or after)**: **Z-N** if not already done —
   it's smallest scope and stacks with everything, and a clean
   win regardless of which other path proves dominant.

Items 5 (Al2O3 condensation) and the refactor stay deferred until
the primary artifact is resolved.

## Pending user scope decision

Which subset to implement in v0.7.3? My recommendation: **4b alone**
as the v0.7.3 diagnostic milestone (informs everything else), then
follow-up scope after seeing the result. But the user may prefer:
- Stack 4b + Z-N for a bigger v0.7.3 ship
- Skip topology testing, go straight to 4 sub-option A (cheapest
  experiment) + Z-N
- Build 4a + 4b + refactor as a complete igniter-architecture v0.8
  milestone, deferring the spike artifact to v0.9
