# v0.7.5 — cross-motor re-LHS result (the v0.7.4 closing task)

Ran 2026-06-03 on a 16-core box, **N=3000 samples/motor**, ~13.8 h wall
(parallel efficiency η≈0.52 ≈ 8 physical cores + HT). Driver:
`examples/cross_motor_lhs_v075.py`; per-sample CSV checkpointing in
`tools/sensitivity.py` (crash-safe). Artifacts (gitignored) under
`artifacts/cross_motor_lhs_v075/`.

## Setup (locked per the v0.7.4 TASKS task)

- **Motors:** the fired motors with experimental traces — Hasegawa A, Zerox,
  Chunc (machbusterNew).
- **Transport:** FROZEN (each motor's `<motor>.frozen.transport.yaml`) — the
  v0.7.3.3 frozen-wins finding.
- **Spike fix ENABLED + FIXED:** `flame_front_enabled=True`, `zn_enabled=True`,
  `kappa_zn=1.0` (Z-N strength is not an LHS variable).
- **Igniter:** explicit `pyrogen='bpnv'` (v0.7.x has no .ric igniter block).
- **Swept knobs** (shared across all 3 motors), physical bounds enforced:
  `roughness ∈ [15, 100] µm`, `kappa ∈ [0.40, 0.50]`,
  `T_ignition ∈ [750, 950] K`, `k_solid ∈ [0.26, 0.32] W/(m·K)`.
- **Combine:** fixed-seed LHS evaluates the SAME knob sets for every motor;
  each motor's fitness is median-normalized (equal weight) and summed; lowest
  combined = the cross-motor optimum.

## Rank-1 cross-motor optimum (combined 1.633)

| knob | value | sanity |
|------|-------|--------|
| **roughness** | **32.2 µm** | physical (>15; near Hasegawa A's 37) |
| **kappa** | **0.439** | physical (≈ the 0.45 Gnielinski center) |
| **T_ignition** | **756 K** | in-band (low end) |
| **k_solid** | **0.271 W/(m·K)** | in the AP/HTPB+Al band (0.26–0.32) |

Per-motor fit: **hasegawa_a 0.374** (good), **zerox 1.484** (moderate),
**chunc 6.559** (poor).

### Top-5 (tightly clustered → robust optimum)

| Rank | combined | roughness [µm] | kappa | T_ign [K] | k_solid | fit_has | fit_zerox | fit_chunc |
|------|----------|----------------|-------|-----------|---------|---------|-----------|-----------|
| 1 | 1.633 | 32.2 | 0.439 | 756 | 0.271 | 0.374 | 1.484 | 6.559 |
| 2 | 1.646 | 29.0 | 0.467 | 766 | 0.269 | 0.370 | 1.531 | 6.505 |
| 3 | 1.655 | 34.1 | 0.419 | 756 | 0.271 | 0.383 | 1.507 | 6.632 |
| 4 | 1.663 | 36.0 | 0.445 | 764 | 0.276 | 0.377 | 1.614 | 6.347 |
| 5 | 1.664 | 30.4 | 0.488 | 782 | 0.260 | 0.386 | 1.628 | 6.297 |

## Reading

- **A clean, fully-physical cross-motor optimum.** The top-5 sit in a narrow
  box (roughness 29–36 µm, kappa 0.42–0.49, T_ign 756–782 K, k_solid
  0.26–0.28) with no knob pegging a bound — i.e. not a degenerate fit, and
  none of the `feedback_roughness_kappa_physical_bounds` rejection criteria
  trip.
- **Chunc stays the worst fit (≈6).** Consistent with the documented known
  limitation: Chunc is a high-L/D motor whose ignition-transient mass-flux
  regime Ma's quasi-steady erosive closure never benchmarked
  (`IGNITION_SPIKE_CLOSEOUT.md`). The median-normalized combine keeps the
  ranking meaningful despite the scale gap.

## Folding in

These are the calibrated **shared** Ma+Goodman knobs for the fired-motor set
on the v0.7.x base. They are committable as-is (this doc + the driver +
checkpointing). Applying them to the canonical defaults / examples
(`run_simulation` defaults, `hasegawa_motor_a.py`, etc.) is a follow-up
calibration step — recommend updating to **roughness 32 µm, kappa 0.44,
T_ignition 756 K, k_solid 0.271** and re-checking each fired motor's trace.
Per the v0.8.0 tag gate, cut v0.8.0 only from a base containing this v0.7.5.

## Folded in — v0.8.0 base (2026-06-05)

Knobs applied to the canonical defaults (`run_simulation`,
`run_from_ric`, `Propellant.k_solid`) and to the three fired-motor
validation examples (`hasegawa_motor_a.py`, `zerox.py`,
`machbusterNew.py`). Re-ran each on the v0.8.0 flat base — note these
use each motor's **embedded `.ric` transport** and the **default
spike-fix state (F+Z OFF)**, NOT the LHS basis (frozen sidecars + F+Z
on), so they are the as-shipped user-path numbers, not a reproduction
of the LHS fitnesses above.

| motor | exp peak [MPa] | sim P_peak [MPa] | ratio | note |
|-------|----------------|------------------|-------|------|
| Hasegawa A | 6.44 | **6.14** @ 2.38 s | 0.95× | excellent (old effective default over-predicted ~1.31×) |
| Zerox      | 3.99 | 7.07 @ 0.017 s   | 1.77× | ignition-spike overshoot; replaces the old per-motor kappa=0.329 (sub-physical) fit with the shared physical optimum |
| Chunc      | 8.88 | 12.65 @ 0.012 s  | 1.42× | known high-L/D QS-erosive limit; down from ~1.9× pre-recal |

The residual Zerox/Chunc over-prediction is the documented ignition-
transient QS-erosive limitation (`IGNITION_SPIKE_CLOSEOUT.md`), not a
calibration miss. The recal is a net win: Hasegawa A lands on target and
the worst-case Chunc spike drops ~25%, all with fully-physical shared
knobs. 373/373 pytest green (one default-pinning assert in
`test_adapter.py` updated 850→756).
