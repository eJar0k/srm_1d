# Candidate 4: Submerged Pyrogen Modes

**One-line**: Extend `igniter_plenum.py` to support alternative
igniter topologies — head-end basket physically inside the upper
grain core, and nozzle-cavity / aft-end igniters like the ISP Super
Loki cartridge — by parameterizing the injection cell range and
momentum-direction sign.

## Motivation

srm_1d currently models ONE igniter topology: a forward 0D pyrogen
plenum that injects into bore cell 0. Real SRMs use multiple
topologies (Sutton 9e §15, NASA SP-8051) which need different
boundary-condition treatments:

1. **Forward pyrogen** (current model) — plenum upstream of cell 0
2. **Head-end basket of BKNO3/BPN pellets** — topologically
   equivalent to (1); parameter tweaks only, no code change
3. **Spherical/submerged head-end pyrogen** (Thiokol US5150654A) —
   pyrogen burns *inside* the main bore volume, distributed across
   the first 5-15% of the bore
4. **Aft-end / nozzle-cavity igniter** (US4023497A, US4503773A;
   Super Loki ISP Corporation) — pyrogen burns inside the bore
   near the aft end, forward-firing jet

The ISP Super Loki motor is already in the repo
(`srm_1d/motors/ISP_Super_Loki.ric`) with embedded (commented-out)
experimental pressure trace; it's the only fired motor in the
library with a non-forward igniter topology. Building topology
support unlocks Super Loki as a validation point and enables future
amateur / sounding rocket use cases.

This candidate is **largest in implementation scope** but is also
the cleanest extension of candidate 3's machinery: once the
axial-distribution kernel exists (candidate 3), supporting
alternative cell ranges is a straightforward parameterization.

## Physical model

Generalize the pyrogen injection model from "weighted by exponential
decay starting at cell 0" (candidate 3) to "weighted by an axial
profile across a user-specified cell range":

```
For topology (3) submerged head-end basket:
  injection_range = (0, i_basket_end)
  w_i = (1 / N_inject) for i in injection_range, else 0   # uniform top-hat

For topology (4) aft / nozzle-cavity:
  injection_range = (i_aft_start, N-1)
  w_i = uniform top-hat across injection_range
  momentum_sign = -1.0     # forward-firing jet → negative axial momentum
```

Pyrogen chamber pressure-coupling closes against the **volume-weighted
average bore pressure across `injection_range`**, not against cell 0:
```
P_bore_ref = sum_{i in range} (V_cell_i * P[i]) / sum_{i in range} V_cell_i
```

For topology (4), the nozzle boundary condition also needs to track
**transient throat-area blockage** if the cartridge ejects (e.g., the
Super Loki cartridge sits in the nozzle until pressure ejects it):
```
A_throat_effective(t) = A_throat_nominal - A_cartridge_blockage(t)
```
This is a step-change at the ejection event (typically 5-50 ms into
the burn), modeled as a user-supplied threshold pressure or time.

## Implementation interface

Affects:
- `srm_1d/igniter_plenum.py` — `PyrogenChamber` gains
  `injection_topology` parameter with three values: `'forward'` (current
  default), `'submerged_head'`, `'aft_cavity'`
- `srm_1d/igniter_plenum.py` — new method
  `_compute_injection_weights(self, x_centers, dx, N)` returning the
  per-cell weight vector for the configured topology
- `srm_1d/simulation.py:_run_time_loop` — pyrogen injection block uses
  the chamber's weight vector + topology's momentum sign + topology's
  P_bore_ref
- `srm_1d/nozzle.py` — `Nozzle` gains optional
  `transient_throat_blockage` callable for topology (4)
- `srm_1d/motors/pyrogens/*.yaml` — pyrogen YAMLs grow optional
  `injection_topology`, `axial_extent_m`, `i_start_fraction`,
  `momentum_sign` keys with safe defaults

Proposed config schema additions:

```yaml
# srm_1d/motors/pyrogens/<name>.yaml
name: my_pyrogen
# ... existing fields ...

# v0.7.2 additions for non-forward topologies:
injection_topology: 'forward'     # 'forward' | 'submerged_head' | 'aft_cavity'
axial_extent_m: 0.020             # length of injection region (default: cell 0 only)
i_start_fraction: 0.0             # bore fraction where injection begins (0 = head end)
momentum_sign: 1.0                # +1 = downstream-firing; -1 = upstream-firing
```

For topology (4)'s transient throat blockage, the `Nozzle` YAML grows:

```yaml
# srm_1d/motors/<motor>.transport.yaml or directly in <motor>.ric
nozzle:
  # ... existing fields ...
  transient_throat_blockage:
    initial_area_fraction: 0.7    # cartridge blocks 30% of throat at t=0
    ejection_pressure_Pa: 5.0e6   # cartridge ejects when bore P crosses this
```

ISP_Super_Loki.ric / transport.yaml would carry these new fields to
demonstrate the topology.

## Validation strategy

1. **Unit tests**:
   - `injection_topology='forward'` with `axial_extent_m=0.0` recovers
     v0.7.1 behavior byte-for-byte (regression gate)
   - Conservation: sum of injected mass / enthalpy / momentum across
     cells equals plenum vent rate × dt for all topologies
   - For topology (4): momentum sign reverses correctly and PISO
     pressure predictor responds appropriately
2. **ISP Super Loki**: re-run with the new topology (cartridge-in-
   nozzle) and compare against the embedded experimental trace
   (must first rename + activate the experimental array in
   `examples/ISP_Super_Loki.py`).
3. **Cross-motor regression**: re-run
   `cross_motor_frozen_vs_effective.py`; the 4 forward-pyrogen motors
   (Hasegawa A / Zerox / BALLSstick / Chunc) should be unaffected
   when topology defaults to `'forward'`.
4. **Submerged-head trace shape**: synthetic test with a Hasegawa-A-like
   motor + `injection_topology='submerged_head'` should show different
   axial pressure development than the forward case (peak appears
   later along the bore, not at cell 0).

## Risks / open questions

- **No open instantaneous flux data for ISP Super Loki** — only
  chamber-averaged thrust + burn time exist publicly. Cartridge-in-
  nozzle ignition-spike validation requires user-supplied data
  beyond what's already commented in `examples/ISP_Super_Loki.py`.
- **Cartridge ejection dynamics** (US4751881 class): no published
  open closed-form transient-throat model. The step-change-A_throat
  approximation is patent-prose-derived; sensitivity to ejection
  pressure may not be well-constrained.
- **ASRM-style canted multi-port igniters**: 12-port radial venting
  with axial components (NTRS 19950017219). A 1D code can model the
  axial momentum component but loses circumferential flame-spread
  benefits. Document as a 1D-fidelity limit; users wanting full
  fidelity need 3D CFD.
- **Aft-cavity pyrogen pressure feedback**: when the pyrogen is at
  the aft end, P_plenum couples to cell N-1 bore pressure which is
  itself coupled to the nozzle. This creates a tighter feedback loop
  than the forward case where cell 0 is far from the nozzle. May
  show new stability issues during ignition transient — needs CFL
  audit.
- **Backwards compat for current motor YAMLs**: every existing motor
  YAML lacks the new `injection_topology` etc. fields. Need explicit
  default-to-forward handling in the YAML loader so v0.7.0 / v0.7.1
  motor configs continue to work unchanged.

## Literature

See [../references/04_submerged_pyrogen_modes.md](../references/04_submerged_pyrogen_modes.md)
(extended lit dive — populated by subagent 2026-05-23). Key
references:
- US Patent 5,150,654A (Thiokol) — canonical submerged head-end
  pyrogen design.
- US Patents 4,503,773A / 4,023,497A — aft-end igniter designs.
- US Patent 4,751,881 — cartridge-in-nozzle baseline.
- NASA SP-8051 — Solid Rocket Motor Igniters monograph (canonical
  taxonomy).
- NASA CR-61238 (1968) + MIT Super Loki Report — Super Loki motor +
  igniter system documentation.
- Cavallini AIAA 2009-5512 — SPINBALL distributed igniter source-term
  formulation (covers head-end basket case).
- John Coker / Richard Nakka amateur references — submerged-head
  igniter design conventions for high-power amateur SRMs.

## Estimated implementation cost

**Largest of the 4 candidates**. ~250-400 LOC across `igniter_plenum.py`
+ `simulation.py` + `nozzle.py` + `propellant.py` + YAMLs + ~8-10
new tests. **2-3 sessions for topology (3) submerged head; another
2-3 sessions for topology (4) with cartridge-in-nozzle** including
ISP Super Loki validation runs.

**Recommendation**: build on top of candidate 3 (pyrogen spatial
distribution). Candidate 3 delivers the axial-weight kernel; this
candidate just parameterizes the injection cell range and adds
topology-specific momentum / pressure-coupling logic. Doing 3 first
makes 4 substantially cheaper.

## Related candidates

- **Strong dependency on candidate 3**: this candidate is essentially
  "candidate 3 + parameterized injection range + topology-specific
  momentum / nozzle treatments." Building 3 first is highly
  recommended.
- **Independent of Z-N (candidate 1)** and **spatial coupling
  (candidate 2)** — works alongside both without conflicts.
- **Enables ISP Super Loki validation** as a fifth fired-motor data
  point for the cross-motor library. This is the only candidate that
  expands the validation set (the others only improve fit on the
  existing 4 motors).
