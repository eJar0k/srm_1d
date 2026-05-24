# srm_1d v0.7.2 — Implementation Plan

**Scope (user-confirmed 2026-05-23)**: candidates 2 (spatial
ignition-front coupling) + 3 (pyrogen axial distribution) as a
paired ship. Z-N (candidate 1) and submerged pyrogen modes (4)
deferred to v0.7.3+.

**Rationale**: 2+3 together replicate the published SPP / SPINBALL /
Han 2017 architecture. Distributing pyrogen energy axially (3)
eliminates the head-end hot spot; coupling cell-to-cell h_c via
cumulative upstream mass flux (2) eliminates the simultaneous-bore-
ignition artifact. Both attack the v0.7.1.1 cross-motor cleanup
finding from complementary physical levers and introduce no new
fitted constants (defaults from Witze coaxial-jet theory + Dittus-
Boelter Re^0.8 scaling).

## Phase A — Pyrogen axial distribution (candidate 3) [FIRST]

Smaller scope; purely additive (current cell-0-only is `L_jet → 0`
limit); establishes the per-cell weight kernel that improves
candidate 2's `G_igniter` input quality.

### Phase A.1 — kernel + plenum extension
- [ ] Add `_compute_pyrogen_axial_weights(x_centers, dx, L_jet, N)`
      Numba kernel in `srm_1d/igniter_plenum.py` (or
      `srm_1d/simulation.py` — TBD which module is the cleaner home).
- [ ] `Pyrogen` dataclass in `srm_1d/propellant.py` gains
      `kappa_jet: float = 8.0` (axial-vent default per Witze 1974 +
      Hersch/Rieser 1971).
- [ ] `PyrogenChamber` exposes a method to derive `L_jet = kappa_jet *
      d_throat_pyrogen` for the current chamber configuration.
- [ ] Pyrogen YAMLs (`srm_1d/motors/pyrogens/*.yaml`) gain optional
      `kappa_jet` key with safe default.

### Phase A.2 — wire weights into `_run_time_loop`
- [ ] Compute weight array once at simulation init (geometry doesn't
      change during pyrogen burn — safe).
- [ ] Replace cell-0-only pyrogen mass / enthalpy writes with the
      weighted-loop pattern.
- [ ] Momentum injection (face-centered): translate cell-centered
      weights to face weights via linear interpolation
      (`0.5*(w[i-1] + w[i])` for face k between cells i-1 and i).
- [ ] `mass_source_by_species` updates: distribute pyrogen species
      mass per weight (v0.7.1 N-species infrastructure handles the
      mixture Cp lookup automatically per cell).
- [ ] Clear `srm_1d/__pycache__/` + `.nbi`/`.nbc` after each kernel
      change (CLAUDE.md gotcha #1).

### Phase A.3 — gates
- [ ] **Regression**: `kappa_jet=0` (or `L_jet=0` edge case)
      recovers v0.7.1.1 byte-for-byte for Hasegawa A. Pytest gate.
- [ ] **Conservation**: pytest gate that sum-of-weights = 1 ±1e-12
      for typical and edge-case geometries.
- [ ] **Hand-calc test**: 10-cell uniform-dx test geometry; verify
      exponential decay matches by-hand calculation.
- [ ] **Momentum balance**: pyrogen-momentum ledger (already in
      v0.7.0) closes to <1% after distribution.

### Phase A.4 — validation
- [ ] Re-run `hasegawa_motor_a.py`; capture trace + flow snapshot
      at t=10 ms (should show distributed pressure across head-end
      cells rather than localized spike at cell 0).
- [ ] Re-run `cross_motor_frozen_vs_effective.py`; tabulate spike
      reduction per motor.
- [ ] Save artifacts under `artifacts/v0_7_2_phaseA_*` for diffing
      vs v0.7.1.1 baseline.

## Phase B — Spatial ignition-front coupling (candidate 2)

After Phase A lands so Phase B's `G_igniter` is informed by the
distributed pyrogen mdot.

### Phase B.1 — cumulative-G kernel
- [ ] Add `_compute_cumulative_mass_flux(G_igniter, rho_p, r_b, P_b,
      A_p, is_burning, dx, N)` Numba kernel.
- [ ] Add `_blowing_augmentation(G_ratio)` Numba kernel (Dittus-
      Boelter Re^0.8 form).
- [ ] Capture `G_ref` constant at simulation init from the pyrogen
      plenum's expected first-vent mass flux.

### Phase B.2 — wire into Goodman call site
- [ ] `_run_time_loop` computes `G_cum[N]` per timestep before the
      Goodman ODE step.
- [ ] **Gate on `is_burning[i] == 0`** when applying augmentation —
      augment h_c only for unignited cells (avoids double-counting
      with Ma 2020 erosive burn rate which also uses local Re).
- [ ] Add `flame_spread_enabled: bool = True` knob on `Propellant`
      for diagnostic A/B comparisons (set False to recover Phase A
      behavior without coupling).

### Phase B.3 — gates
- [ ] **Regression**: `flame_spread_enabled=False` recovers Phase A
      byte-for-byte. Pytest gate.
- [ ] **Kernel test**: `G_cum[i]` is monotonically increasing in i
      when upstream cells are burning.
- [ ] **Kernel test**: `_blowing_augmentation(0) == 1.0` and is
      monotonically increasing in G_ratio.
- [ ] **No-burning baseline**: at t=0 before any cell ignites,
      G_cum is exactly G_igniter / A_port[0] for all cells.

### Phase B.4 — validation
- [ ] Re-run `hasegawa_motor_a.py`; target: spike drops toward
      experimental 6.5 MPa; ignition propagation visible in snapshot
      history (`is_burning` ignites cells sequentially over 50-200
      ms, not simultaneously).
- [ ] Re-run `cross_motor_frozen_vs_effective.py`; target: all 4
      motors show spike-to-plateau ratio < 1.5; Zerox / BALLSstick /
      Chunc over-prediction shrinks meaningfully.
- [ ] Compare Phase A vs Phase B vs (A+B) traces to attribute the
      improvement contributions.

## Phase C — Integration validation + optional re-calibration

After both Phases A and B land.

### Phase C.1 — cross-motor sanity check
- [ ] Full `cross_motor_frozen_vs_effective.py` re-run on the v0.7.2
      build for diff vs v0.7.1.1.
- [ ] Re-run `pytest srm_1d/tests/` — all 213 baseline + new Phase A
      and B tests must pass.

### Phase C.2 — Hasegawa A re-calibration (only if needed)
- [ ] Score `hasegawa_motor_a.py` default trace against experimental;
      if P_peak error > 15% or shape MSE > 0.15 MPa^2, run a small
      LHS sweep over (roughness, kappa, T_ignition, k_solid, κ_jet)
      with the v0.7.2 build to find new rank-1.
- [ ] If LHS rank-1 lies within physical-realism bounds (roughness
      ≥ 15 µm, kappa near 0.45, k_solid in 0.20-0.40), canonize and
      update `hasegawa_motor_a.py`. If outside, retain v0.7.1.1 knobs
      and document.

### Phase C.3 — tag + memory updates
- [ ] DEVNOTES Calibration State updated with v0.7.2 cross-motor
      results.
- [ ] CLAUDE.md gotcha #5 updated to reflect that the structural
      fix shipped.
- [ ] `project_hasegawa_calibration_state` memory updated.
- [ ] `project_v0_7_2_design_package` memory annotated with what
      shipped vs what deferred.
- [ ] Tag `v0.7.2`.

## Tag criteria

`v0.7.2` ships when:

- All 213 baseline pytest tests pass plus the new Phase A+B tests
  (target +8-12 new tests).
- Hasegawa A `hasegawa_motor_a.py` P_peak under-prediction shrinks
  to ≤ 5% (vs current 31% over-prediction at the effective default).
- All 4 fired motors in `cross_motor_frozen_vs_effective.py` show
  spike-to-plateau ratio < 1.5 at default knobs.
- No new fitted constants outside literature-defensible ranges
  (`kappa_jet` ∈ [2, 12]; cumulative-G coupling uses Dittus-Boelter
  Re^0.8 with no free parameter).
- DEVNOTES API-breaking-change log updated (new `kappa_jet` field
  on `Pyrogen`; new `flame_spread_enabled` field on `Propellant`).

## Deferred to v0.7.3+

- **Candidate 1 (Z-N dynamic burn rate)** — if v0.7.2's 2+3 ship
  does not fully close the spike artifact, Z-N is the v0.7.3
  first-pick because it addresses a different phenomenon (burn-rate
  ramp lag) and stacks cleanly.
- **Candidate 4 (submerged pyrogen modes)** — depends on Phase A's
  axial-weight kernel, so v0.7.2 unblocks it. Defer to v0.7.3+
  unless the user wants ISP Super Loki validation prioritized.
- **Cross-motor effective-transport LHS recalibration** for
  Zerox / BALLSstick / Chunc / Super Loki — if the v0.7.2 build
  matches experimental shape qualitatively, per-motor LHS is
  optional polish, not a tag blocker.

## Current state (2026-05-23 design-phase close-out)

- Design package shipped at commit `52bf32c`.
- Implementation has NOT started.
- This TASKS.md is the v0.7.2 implementation plan; pending user
  sign-off on Phase A start.
