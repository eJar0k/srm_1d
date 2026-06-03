# Candidate 3: Pyrogen-Injection Spatial Distribution

**One-line**: Replace cell-0-only pyrogen injection with an axial
profile (exponential decay length `L_jet = κ_jet · d_throat`) so
pyrogen mass + enthalpy + momentum deposit along the bore over the
plume's actual physical reach, not as a delta function at the head
end.

## Motivation

The current `igniter_plenum.py` pumps 100% of pyrogen mdot, sensible
enthalpy, and axial momentum into cell 0 (the first head-end bore
cell). Real pyrogen jets penetrate some distance along the bore
before equilibrating with the bore flow — for an axial-vent sonic
pyrogen, jet potential-core length is ~6-7 throat diameters
(Witze 1974); for an SRM bore of length 0.3 m and throat diameter
7 mm, that's ~5 cm or roughly the first 10-25% of the bore.

Concentrating all pyrogen energy into cell 0 creates an unrealistic
localized hot spot that:
1. Drives cell 0 to ignite far before downstream cells, contributing
   to the simultaneous-bore-ignition artifact when downstream cells
   catch up
2. Generates an unphysically large head-end pressure pulse during
   the first ms of pyrogen burn (visible as the early spike in all
   4 fired motors at default knobs)

Distributing pyrogen energy axially fixes both: the spike is blunted
because the energy density per cell is lower, and cells in the
impingement region heat up together while still being faster than
cells outside the impingement region (preserving the propagation
mechanism candidate 2 also addresses).

## Physical model

For each bore cell `i` with axial center `x_i`:
```
w_i = exp(-x_i / L_jet) * dx[i] / sum_j[exp(-x_j / L_jet) * dx[j]]
mdot_pyrogen_i = w_i * mdot_plenum
hdot_pyrogen_i = w_i * hdot_plenum
momdot_pyrogen_i = w_i * momdot_plenum

L_jet = kappa_jet * d_throat_pyrogen
```

where:
- `L_jet` is the characteristic exponential decay length of the
  pyrogen jet's energy deposition
- `kappa_jet` is the dimensionless decay length, one new physical
  knob per pyrogen
- `d_throat_pyrogen` is the existing `PyrogenChamber.D_throat`

**Defensible defaults** from the jet-penetration literature:
- `kappa_jet ≈ 6-10` for choked sonic axial-vent pyrogens
  (Witze coaxial-jet theory + Hersch/Rieser sonic injection data)
- `kappa_jet ≈ 2-4` for predominantly radial-vent pyrogens (JICF
  with q < 5 has modest lateral penetration)

**Conservation guarantees**:
- `sum_i w_i = 1` by construction → exact mass / enthalpy / momentum
  conservation
- In the limit `L_jet → 0` the model reduces to current cell-0-only
  behavior (the current state becomes a special case)

## Implementation interface

Affects:
- `srm_1d/igniter_plenum.py` — add `L_jet` parameter to
  `PyrogenChamber` (or to the injection helper in `simulation.py`)
- `srm_1d/simulation.py:_run_time_loop` — pyrogen injection block
  changes from cell-0-only writes to a weighted loop over the
  affected cell range
- `srm_1d/propellant.py:Pyrogen` — new field `kappa_jet: float = 8.0`
  (default for axial-vent BPNV/MTV pyrogens)
- `srm_1d/motors/pyrogens/*.yaml` — each pyrogen YAML grows a
  `kappa_jet` entry (with safe default if missing)

Proposed kernel (Numba):

```python
@njit(cache=True)
def _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N):
    """Compute exponential-decay weights for pyrogen distribution.

    w_i = exp(-x_i / L_jet) * dx[i] / sum_j[exp(-x_j / L_jet) * dx[j]]

    L_jet=0 special case: all weight in cell 0 (returns current
    cell-0-only behavior).
    """
    w = np.zeros(N)
    if L_jet <= 1e-12:
        w[0] = 1.0
        return w
    total = 0.0
    for i in range(N):
        w[i] = np.exp(-x_centers[i] / L_jet) * dx[i]
        total += w[i]
    if total <= 0.0:
        w[0] = 1.0
        return w
    for i in range(N):
        w[i] /= total
    return w
```

Existing pyrogen-injection call sites in `_run_time_loop` change
from:
```python
mass_source[0] += mdot_pyrogen / dx[0]
mass_source_by_species[0, _SPECIES_IGNITER] += mdot_pyrogen / dx[0]
thermal_source[0] += mdot_pyrogen * Cp_pyrogen * T_pyrogen / dx[0]
# (and similarly for momentum)
```
to:
```python
w = _compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)
for i in range(N):
    if w[i] <= 0.0:
        continue
    mass_source[i] += w[i] * mdot_pyrogen / dx[i]
    mass_source_by_species[i, _SPECIES_IGNITER] += w[i] * mdot_pyrogen / dx[i]
    thermal_source[i] += w[i] * mdot_pyrogen * Cp_pyrogen * T_pyrogen / dx[i]
    # momentum: face-centered, slightly trickier — weight face injection by an
    # adjacent-cell pair-of-weights average. See implementation notes.
```

The weights `w` are constant per simulation (assuming bore geometry
doesn't change during pyrogen burn — true for the first ~100 ms of
typical operation), so the weight array can be computed once at
simulation start and reused every step.

## Validation strategy

1. **Unit tests**:
   - `L_jet=0` recovers cell-0-only behavior byte-for-byte (regression
     gate for the current model)
   - `sum(w) == 1` for all `L_jet > 0`
   - Hand-calculated weight distribution matches kernel output for a
     fixed 10-cell test geometry
2. **Hasegawa A**: re-run `hasegawa_motor_a.py`; target: head-end
   pressure spike drops toward experimental, axial pressure
   snapshot at t=10 ms shows uniform pressure across the impingement
   region rather than a sharp head-end spike.
3. **Cross-motor**: re-run `cross_motor_frozen_vs_effective.py`;
   target: all 4 motors show meaningful spike reduction (15-30%) at
   `kappa_jet = 8`.
4. **Mass / energy conservation**: pytest gate that overall mass
   balance error stays under 1% with the new distribution (should
   be invariant since weights sum to 1).
5. **Limit-case tests**: simulations at `kappa_jet = 0.5` (very
   localized) should reproduce current behavior; `kappa_jet = 20`
   (very spread) should show essentially uniform pyrogen deposition
   across the first half of the bore.

## Risks / open questions

- **κ_jet calibration**: literature gives a range (2-10) but the
  right value for a specific pyrogen depends on its throat geometry,
  port count, and orientation (axial vs canted vs radial). Per the
  ScienceDirect 2025 parametric study, jet Mach number and angle are
  first-order. Initial v0.7.2 ships a single scalar per pyrogen;
  more sophisticated models (per-port distribution) are v0.7.3+.
- **Momentum injection** is face-centered in PISO, not cell-centered.
  The axial profile needs careful translation: face k between cells
  i-1 and i receives a weight derived from `0.5*(w[i-1] + w[i])`
  (linear interpolation between adjacent cell weights). Need a smoke
  test that pyrogen momentum still balances against PISO's pressure
  predictor as the v0.7.0 audit data confirmed.
- **Interaction with v0.7.1 N-species**: the pyrogen species mass
  source now gets distributed too, which means cells far from the
  head end get a small pyrogen fraction earlier than in the current
  model. Phase 3.5's per-species Cp lookup handles this correctly
  (each cell's mixture Cp updates per advected Y), so no kernel
  change is needed.
- **Interaction with the relaxed T_ceiling**: candidate 2's strict
  T_ceiling formula correctly handles the case where a cell has
  significant pyrogen mass fraction; the distributed-pyrogen change
  may make this case more common (cells in the impingement band see
  Y[i, pyrogen] > 0.05 earlier). Verify the strict ceiling doesn't
  artificially clip the trace.

## Literature

See [../references/03_pyrogen_spatial_distribution.md](../references/03_pyrogen_spatial_distribution.md)
(extended lit dive — populated by subagent 2026-05-23). Key
references:
- Peretz, Kuo, Caveny, Summerfield 1973 (NASA TN, NTRS 19740005393)
  — foundational head-end-boundary + flame-spreading framework.
- Cavallini SPINBALL impingement-region — the canonical published
  axial-distribution model.
- Witze 1974 — coaxial-jet potential-core length (sets `kappa_jet`
  default for axial-vent pyrogens).
- Hersch & Rieser 1971 (NTRS 19710018794) — empirical L_pen/d_j vs
  pressure ratio and Mach number.
- Yıldız et al. 2025 (ScienceDirect S1290072925006453) — recent
  parametric study confirming jet angle is first-order.

## Estimated implementation cost

Medium. ~80-150 LOC across `igniter_plenum.py` +
`simulation.py` + `propellant.py` + pyrogen YAMLs + ~5 tests. 1-2
sessions. Smaller than candidate 4 (submerged pyrogen) which adds
new topologies + boundary conditions.

## Related candidates

- This candidate is the structural ignition-stage fix attacking the
  same artifact as candidate 2, but from the pyrogen-source side
  rather than the gas-side h_c side. **They are complementary**:
  distributing pyrogen energy axially AND coupling cell-to-cell
  h_c are the two halves of the SPINBALL / Han 2017 architecture.
- **Z-N (candidate 1)** is independent — Z-N affects burn-rate
  dynamics post-ignition, this affects ignition timing. Stacks
  cleanly.
- **Candidate 4 (submerged pyrogen)** uses essentially the same
  axial-distribution machinery this candidate builds, just with the
  injection cell range located mid-bore or aft instead of always
  at the head end. Candidate 4 should be built ON TOP OF candidate
  3's machinery — pick this one first and 4 becomes much cheaper.
