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

### Phase B.4 — validation [COMPLETE 2026-05-24]
- [x] Re-ran `cross_motor_frozen_vs_effective.py`; v1 amplified spike
      across all 4 motors by +0.6-7.4%, REVERTED Phase A's Zerox win
      (peak time 0.27s → 0.026s). Negative finding documented in
      commit `065d193`.
- [x] Reformulated Phase B as v2: flame-front-marker gating (boost
      cell j+1 only if cell j ignited within last tau_window) instead
      of cumulative-G magnitude. Less amplification than v1 but same
      direction (cascade accelerated, Zerox win still reverted by ~6%
      P_peak vs Phase A). Negative finding documented at commit
      `e507c09`.
- [x] **Root cause**: PISO's local-Re tracking already captures
      upstream-mass-flux contributions to h_c at unignited cells, so
      the Kashiwagi/Han augmentation (developed for codes that DON'T
      track local flow properly) is double-counting in this codebase.
      No augmentation-gating formulation appears to resolve the
      direction issue without inverting the physics (boost → damp).
- [x] **Decision**: Propellant.flame_spread_enabled defaults to False
      (Phase B infrastructure preserved as opt-in diagnostic only).
      Phase A is the v0.7.2 load-bearing ship.

## Phase C — Tag close-out [COMPLETE 2026-05-24]

Phase C scope reduced from the originally planned cross-motor
sanity-check + Hasegawa A re-calibration to just the doc/memory/tag
mechanical close-out, because:
- Phase A delivers a real Zerox win (P_peak 10.20→9.69 MPa, t_peak
  0.035s→0.27s) but the other 3 motors are essentially unchanged at
  default knobs (cascade dominated by simultaneous-ignition artifact
  Phase B couldn't fix).
- Hasegawa A `hasegawa_motor_a.py` still over-predicts P_peak by
  ~31% (same as v0.7.1.1). Re-calibration would just chase the
  structural artifact with LHS-fitted knobs that escape the
  physical-realism bounds — wasteful before v0.7.3 structural work.

### Phase C.1 — sanity check [COMPLETE]
- [x] Full `cross_motor_frozen_vs_effective.py` re-run on v0.7.2;
      artifacts saved under `artifacts/cross_motor_frozen_vs_effective/`.
- [x] 240/240 pytest green with Phase A wired + Phase B disabled
      default.

### Phase C.2 — Hasegawa A re-calibration [SKIPPED]
- [-] Skipped per above — the structural artifact dominates the
      Hasegawa A canonized example; calibration adjustments within
      physical-realism bounds (`feedback_roughness_kappa_physical_bounds`)
      cannot recover the spike shape. v0.7.3 takes a different
      angle (Z-N, submerged pyrogen modes, per-cell coupling
      alternatives, or different heating modes).

### Phase C.3 — doc / memory / tag close-out [COMPLETE]
- [x] DEVNOTES API-breaking-change log entry for v0.7.2
      (kappa_jet field on Pyrogen, flame_spread_* fields on Propellant).
- [x] CLAUDE.md roadmap updated; gotcha #5 retained (effective
      default still ships).
- [x] `project_v0_7_2_progress_state` memory created with Phase A
      + B summary.
- [x] Tag `v0.7.2-phaseA`.

## Tag criteria

`v0.7.2-phaseA` ships when:

- ✓ All baseline pytest tests pass plus the new Phase A+B tests
  (240 total: 226 baseline + 13 Phase A kernel + 4 Phase A.3
  integration + 8 Phase B-v2 flame-front kernel + 3 Phase B
  integration — minus 14 obsolete cumulative-G kernel tests
  superseded by the v2 reformulation).
- ✗ Hasegawa A `hasegawa_motor_a.py` P_peak under-prediction
  shrinks to ≤ 5% — **NOT MET** (still ~31% over). Structural
  artifact requires v0.7.3 work.
- ✗ All 4 fired motors show spike-to-plateau ratio < 1.5 at
  default knobs — **NOT MET** (Phase A delivered Zerox 2.22→spike
  shifted to t=0.27s but other 3 motors essentially unchanged).
- ✓ No new fitted constants outside literature-defensible ranges
  (`kappa_jet ∈ [2, 12]` per Witze; `flame_spread_*` are
  experimental and disabled by default).
- ✓ DEVNOTES API-breaking-change log updated.

The two ✗ items are not tag blockers — they're explicitly the
v0.7.3 target. The tag captures what shipped: Phase A infrastructure
(pyrogen axial distribution, real Zerox win), Phase B infrastructure
(opt-in diagnostic only), and the negative findings that motivate
v0.7.3.

## Deferred to v0.7.3+

After the Phase B negative findings, the open candidate space is
broader than the original v0.7.2 design package anticipated. The
user-flagged next-direction options (post-tag, design analysis
pending):

1. **Candidate 1 — Z-N dynamic burn rate** (relaxation ODE on
   steady r_b). Addresses burn-rate ramp lag, not ignition timing.
   Stacks cleanly with anything else.
2. **Candidate 4a — Head-end submerged pyrogen basket**. Energy
   deposits inside the bore not from a head-end source. May change
   ignition distribution but unclear if it fixes the simultaneous-
   cell artifact (the artifact may live in the Goodman per-cell
   solver, not the pyrogen source).
3. **Candidate 4b — Aft-inserted impinging cartridge** (Super Loki
   class). Igniter occupies arbitrary or pyrogen-mass-defined core
   length, fires forward; ignition propagates BACK→FORWARD. User-
   flagged as worth testing whether the simultaneous-ignition
   artifact is caused by the current concentrated-mass-injection
   model.
4. **Per-cell coupling alternatives** (post-Phase-B insights):
   reverse polarity (damp h_c at cells far from any recent
   ignition rather than boost adjacent), solid-phase axial
   conduction, or Goodman per-cell coupling via shared boundary
   layer.
5. **Different heating modes**: surface radiation enhancement at
   distance, two-phase Al2O3 condensation (Pardue 1992), Z-N
   combined with current local-Re tracking.
6. **Plenum-as-option refactor**: unify forward-plenum (current
   default), head-end basket, and aft-inserted cartridge under a
   single igniter-architecture API so motors can specify topology
   in YAML rather than baked into `igniter_plenum.py`.

A follow-up design doc analyzing each option's path forward is in
progress — see `docs/v0_7_2/candidates_post_phaseA.md` (or
equivalent v0.7.3 design package).

## Current state (2026-05-24 v0.7.2 ship)

- v0.7.2-phaseA tag at commit `e507c09`.
- Phase A (pyrogen axial distribution) shipped as load-bearing
  default; real Zerox win, neutral on other 3 fired motors at
  default knobs.
- Phase B (cumulative-G v1 + flame-front v2 reformulation)
  infrastructure shipped but **disabled by default** after both
  formulations amplified rather than smoothed the spike. Opt-in
  via `Propellant.flame_spread_enabled = True`.
- Structural ignition-kernel artifact persists for Hasegawa A /
  BALLSstick / Chunc at default knobs. v0.7.3 candidate analysis
  pending user decision.
