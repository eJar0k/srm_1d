# Ignition-tuning Cartesian sweep audit (Step 4) — 2026-05-21

## Setup

- Branch: `v0.7.0-phase4` at commit `2e20355` (post `--mode
  ignition-tuning` addition).
- Motor: Hasegawa A, `t_max = 5.5 s`.
- CLI: `ignition_spike_diagnostic.py --mode ignition-tuning --t-max
  5.5`.
- Matrix: 18-variant Cartesian sweep over `T_ignition ∈ {650, 750,
  850} K`, `k_solid ∈ {0.2, 0.3, 0.5} W/(m·K)`, `radiation_emissivity
  ∈ {0.0, 0.45}`.
- Artifacts: `artifacts/ignition_diagnostics/hasegawa_a_ignition_tuning/`.

## Result

| rank | variant            | MSE [MPa²] | P_peak [MPa] | t_peak [s] |
| ---- | ------------------ | ---------- | ------------ | ---------- |
| 1    | T650_k0p2_eps0p45  | **1.3456** | 12.18        | 0.0242     |
| 2    | T650_k0p3_eps0p45  | 1.3461     | 12.17        | 0.0242     |
| 3    | T750_k0p2_eps0p45  | 1.3463     | 12.17        | 0.0243     |
| …    | …                  | …          | …            | …          |
| 18   | T850_k0p5_eps0p00  | 1.3496     | 12.09        | 0.0246     |

**All variants fail.** MSE 1.34–1.35 MPa² is ~9× the
plan's 0.15 MPa² pass threshold.

## Diagnosis

The peak pressure is essentially constant across the sweep
(12.09–12.18 MPa, range 0.09 MPa = 0.7 %). T_ignition shifts the
ignition timing by < 1 ms and k_solid has almost no observable effect
on peak. Radiation emissivity (0.0 vs 0.45) shifts peak by < 0.05 MPa.

This is the **expected** physical outcome: T_ignition and k_solid
control *when* cells ignite, not *how much* gas they produce after
ignition. With instantaneous ignition (full r_total at first burn),
the chamber receives the same total mass flux regardless of when each
cell crossed its ignition criterion. The peak is set by:

- Saint-Robert burn rate `a · P^n`
- Bore Klemmung `Kn = A_burn / A_throat`
- Erosive amplification (`kappa`, `roughness`)

None of which are in the Step 4 matrix.

## Decision

**Escalate to Step 5.** Per the plan, this triggers the
Peretz-aligned thermal-layer establishment (participation fraction)
model. The structural argument is that the 12 MPa instantaneous spike
in the simulation reflects *every* grain cell contributing its full
steady-state burn rate as soon as `T_surf > T_ignition` — but
physically, a freshly-ignited surface has a thermal layer thickness
δ << δ_steady and cannot yet drive its full mass flux. The
participation fraction `φ[i] = min(1, δ_burned / δ_steady)` emerges
from the same Goodman ODE that drives ignition, so it introduces no
new fitted constant.

## Side observation — radiation effect tiny

The sweep included an eps=0 vs eps=0.45 axis to confirm that
radiation is operating correctly. Within each (T_ignition, k_solid)
pair, eps=0.45 produces slightly higher peak (~12.17 vs 12.14 MPa) —
consistent with marginally faster ignition spread leading to slightly
faster pressure rise, but well within the same instant-spike regime.
This confirms that the radiation kernel is wired in correctly and
the residual physical mismatch is in the burn-rate model, not the
ignition spread model.
