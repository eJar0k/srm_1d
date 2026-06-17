# Ignition-Spike Investigation — RE-OPENED (2026-06-16/17)

**This document updates `IGNITION_SPIKE_CLOSEOUT.md`.** The user re-opened the
question with the explicit instruction to *suppress the closeout's prior
assumptions* and test whether the spike is a structural solver artifact rather
than fundamental burn-rate physics. The verdict after a full multi-angle
re-investigation:

> **The spike splits into two distinct phenomena. The dramatic *sharp early
> 1.54× spike* IS a structural solver artifact (and the user was right about
> it). The smaller *residual ~1.25× erosive hump* is largely genuine quasi-steady
> erosive physics — the closeout's core conclusion stands, now on far stronger,
> quantitative, multi-source evidence.**

One real bug was found and fixed along the way (the Phase F flame-front was
silently bypassed for `forward_plenum`). Two physically-grounded fidelity
upgrades were recorded for future implementation. The diagnostic motor
throughout is Chunc / `machbusterNew` (high-L/D, canonical knobs: bpnv pyrogen,
Sutton sizing, roughness 32 µm, kappa 0.44, T_ignition 756 K).

---

## 1. Baseline & decomposition

| quantity | value |
|---|---|
| P_peak (baseline) | 12.65 MPa @ 11.6 ms |
| plateau (0.4–0.9 s) | 8.23 MPa |
| spike ratio | **1.54×** |
| experimental (Chunc) | ~8.5 MPa by 10 ms, flat ~8.8, **no spike** |
| erosive-disabled | early spike **gone** (peak at end) → the spike is 100% erosive |

The plateau and decay already match experiment well; the entire model error is
the 0–25 ms transient.

---

## 2. The sharp early spike is structural (the supersonic-fill ignition gate)

**Ignition is near-simultaneous: all 100 grain cells light in 0.61–1.47 ms**
(0.86 ms spread). A fill-resolved probe shows why — a hot (~3100 K) gas-dynamic
contact front sweeps the bore at **~800 m/s** and the convective Goodman gate
ignites each cell the instant the hot gas arrives.

Tracing one mid-grain cell (x = 435 mm) through its ignition:

```
 t_ms  Mach   u(m/s)     Re   Tgas  Tsurf   q''(MW/m2)
 0.92  0.06     21      5290   293   293       0.0    (cold, quiescent)
 0.97  0.66    226    136610   296   294       0.0    (fast gas, still COLD)
 1.02  1.55   1703    132780  3119   820      39.0    IGNITES, one step
```

The cell ignites at **Mach 1.55, Re ≈ 133,000**. Over the 0.6–1.5 ms ignition
window the bore is transonic-to-supersonic (**max Mach 1.3–3.3, Re up to
330,000**), and the velocities are continuity-consistent (`ρuA` = 0.5–0.8× the
cumulative source — honest gas dynamics, not a PISO blow-up).

**Root cause:** the ignition gate sets surface heat flux by feeding these
supersonic fill velocities into the **Gnielinski/Haaland steady, low-Mach,
developed-pipe-flow correlation** (`_bare_heat_transfer_coeff` → `gnielinski_nusselt`).
At Re ≈ 10⁵ the `Nu ∝ Re^0.8` scaling returns q″ ≈ 40–78 MW/m², overwhelming the
Goodman thermal-inertia lag and lighting the whole grain in ~1 ms. The
correlation is being used 1–3 orders outside its validity. This is structural,
not burn-rate physics — and it is the dominant cause of the *sharp* spike.

### Gate-Mach-clamp experiment (negative for the principled bound)
Bounding the ignition-gate Re to the Gnielinski validity Mach:

| clamp | P_peak | ratio | full ignition |
|---|---|---|---|
| off (baseline) | 12.65 | 1.54× | 1.5 ms |
| Mach 0.30 (validity bound) | 12.56 | 1.53× | 1.6 ms |
| Mach 0.05 | 11.66 | 1.42× | 4.4 ms |
| Mach 0.02 | 10.61 | 1.29× | 11.2 ms |

At the *defensible* Mach 0.3 the clamp is a no-op. Lesson: **hot-gas presence,
not flux magnitude, is the bottleneck** — once 3100 K gas occupies a cell (which
real gas dynamics deliver everywhere in ~1.5 ms), even a validity-clamped flux
ignites the cold surface within ~1 ms. The clamp was reverted (not kept).

---

## 3. The Phase F flame-front bug (found and FIXED)

The principled way to spread ignition is the v0.7.4 Phase F flame-front (withhold
surface heating from cells the front hasn't reached). But `flame_front_velocity`
from 15→100 m/s gave **identical** results (full ignition fixed at ~4.6 ms).

**Bug:** for `forward_plenum`, the DeMar surface-flux target is `head_grain_cell`
= the *first unignited grain cell*, and the gate exempted any cell carrying DeMar
flux. As cells ignite, that target **marches down the grain** — a self-propagating
DeMar cascade that bypassed the front entirely. (This is why the v0.7.4 closeout
found flame-front "ineffective on Chunc" — it was testing a bypassed front.)

**Fix (landed in the working tree, `simulation.py`,
`_goodman_ignition_sources_and_mass`; not yet committed):** removed the DeMar-flux
exemption from the `heat_cell` gate. The seed grain cell is
always in `ignitable`, so it still lights from the head-end plume; DeMar flux is
applied only at ignitable cells, so it cannot ignite ahead of the front.

Verification (default-off, no regression):

| case | P_peak | ratio | full ignition |
|---|---|---|---|
| default (ff off) | 12.65 | 1.54× | 1.5 ms (unchanged) |
| ff v=100 | 11.18 | 1.36× | 8.9 ms (≈ 0.86 m / 100) |
| ff v=50 | 10.47 | 1.27× | 17.0 ms |
| ff v=15 | 10.32 | 1.25× | 54.5 ms |

The front is now authoritative and `flame_front_velocity` controls ignition
spread as documented — **without** the `diagnostic_disable_pyrogen_surface_heating`
crutch. 406/406 non-working-tree tests pass.

---

## 4. The residual ~1.25× erosive hump is mostly physical

Flame-spread (any speed) plus a per-cell burn-establishment ramp both bottom out
at a hard **~1.25× floor** (10.3 MPa); the establishment ramp adds essentially
nothing on top of flame spread. `diagnostic_disable_erosive` removes the hump →
**it is erosive.** Two papers (now in `docs/references/`) settle its nature:

**Mukunda & Paul 1997 (`mukunda1997.pdf`) — the hump is largely real.** Their
universal, **Hasegawa-validated** law `η = 1 + 0.023(g^0.8 − 35^0.8)·H(g−35)`,
`g = (G/ρ_p r₀)(Re₀/1000)^−0.125`, computed per-cell on Chunc:

| | η (sim, Ma) | η (Mukunda-Paul) | Ma/MP | cells g>35 |
|---|---|---|---|---|
| spike (12 ms) | 2.76 | **2.08** | 1.33 | 89/100 (g≈181) |
| plateau (600 ms) | 1.27 | 1.19 | 1.06 | 53/100 |

The threshold (g_th=35) is **not tripped-suppressed** at the transient (g≈181).
The validated universal law itself predicts η≈2.1 at the spike condition — i.e.,
substantial erosive enhancement is genuine QS physics (their Fig. 10 shows the law
reproducing real-motor erosive humps). **Ma over-predicts only ~33% vs the
universal law at high g** (agreeing within 6% at the plateau).

**Liñán & Williams 1971 (`linan1971.pdf`) — the Goodman gate is validated.**
Ignition is a thermal runaway, but *to leading order it occurs when the inert
surface temperature reaches a chemistry-set threshold* — exactly the Goodman
`T_surf > T_ignition` structure. The reactive correction is one order-unity factor
(e^b = 0.65, ~35%). So the Goodman kernel is leading-order correct; an
integrated-energy/AP-melt upgrade is a *fidelity* refinement, not a spike fix.

---

## 5. The last floor-breaking lever (relaminarization) — TESTED AND DEAD

The only remaining non-tuning candidate to suppress the erosive term transiently
was favorable-pressure-gradient relaminarization (turbulence cannot be sustained
under strong flow acceleration; erosive burning is a turbulent effect). Criterion:
`K = (ν/u²)(du/dx) > ~3×10⁻⁶` (Kays-Crawford / Narasimha-Sreenivasan).

Instrumented on Chunc at the erosive-active cells:

| | peak-erosive cell K | median K | % above 3×10⁻⁶ |
|---|---|---|---|
| spike (12 ms) | **−4.6×10⁻⁷** | 7.6×10⁻⁸ | **0%** |
| plateau (600 ms) | −1.6×10⁻⁷ | 1.4×10⁻⁸ | 0% |

K ≈ 10⁻⁷ — two orders below threshold — at both spike and plateau, and **negative**
(adverse/decelerating, du/dx ≈ −6800 s⁻¹) at the peak-erosive aft cell (x = 825 mm).
The spike flow is fully turbulent. **Relaminarization is not active; the lever is
dead.** (The sibling Beddini head-end-laminar spatial gate is dead too: the spike
is at the aft, far past any head-end laminar zone, and it is a steady effect that
would alter the plateau.)

---

## 6. Conclusion & disposition

- **No parameter-free closure removes the residual erosive hump.** Every
  physically-grounded lever — ignition timing (flame spread), burn buildup
  (establishment ramp), ignition criterion (Liñán-Williams validates Goodman),
  erosive threshold (Mukunda-Paul: not tripped), relaminarization (K dead) — has
  been eliminated. The hump is ~⅔ genuine QS erosive physics + ~⅓ Ma high-g
  over-prediction. **The v0.7.4 closeout stands, reinforced.**
- **The sharp early spike is structural** (supersonic-fill convective ignition
  gate + the now-fixed DeMar-cascade flame-front bug). With the flame-front bug
  fixed and the front enabled, the 0–10 ms rise reshapes toward the experimental
  trace; the sharp spike becomes a modest delayed hump.
- **Landed (working tree, uncommitted):** the DeMar-cascade flame-front fix
  (correctness; opt-in Phase F now functional; default-off, no regression).
- **Recorded for future (fidelity, NOT spike fixes):** (1) de-tune Phase Z to the
  parameter-free `τ_cond = ln(10²)·α_s/r² ≈ 4.6 α_s/r²` (Lengellé thermal-wave
  time); (2) integrated-energy / AP-melt ignition criterion with the endothermic
  AP-decomposition debit (Liñán-Williams §5 plateau mechanism). See the
  `project-ignition-fidelity-candidates` memory.

Diagnostic scripts used live under `c:/tmp/` (chunc_spike_probe, chunc_fill_probe,
chunc_gate_audit, chunc_mach_audit, chunc_cell_history, run_gate_clamp_experiment,
run_flamespread_experiment, run_flamespread_combo, chunc_mukunda_compare,
chunc_K_relaminarization); comparison plots under
`artifacts/ignition_mach_clamp/`, `artifacts/ignition_flamespread*/`.
