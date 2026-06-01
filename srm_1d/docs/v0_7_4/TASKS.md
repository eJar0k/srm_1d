# v0.7.4 — Ignition-Transient Spike: Implementation Plan

> **Implementation outcome (2026-05-31) — read this first.**
> Both phases shipped (opt-in, default OFF) and validated on Chunc +
> Hasegawa A. Key corrections vs the plan below, found during validation:
>
> 1. **Phase F `v_flame` derivation was wrong.** The planned per-step
>    `v_flame = q''/(ρ·Cps·ΔT)` is dimensionally a *burn/regression*
>    velocity (~1–2 mm/s), not a lateral *flame-spread* velocity (~1–10
>    m/s) — it starved the grain. Replaced with a **bounded literature
>    constant** `Propellant.flame_front_velocity` (default 3.0 m/s; AP/HTPB
>    spread 1–10 m/s, Peretz-Kuo-Caveny-Summerfield 1973, Kumar & Kuo
>    1984). The derived de Ris / successive-heating forms remain options
>    if a parameter-free spread velocity is wanted later.
> 2. **`head_basket` seed bug.** The front seed (`cart_i_start`) maps to
>    the non-grain motor-head cavity, so the grain (starting a few cells
>    aft) was fully gated → 0 cells ignited → motor never lit (velocity-
>    independent 0.52 MPa). Fixed by **snapping the seed to the first
>    grain cell** in the propagation direction.
> 3. **`head_basket` gate-bypass.** The pyrogen-flux exemption (meant for
>    forward-plenum's single DeMar cell) was exempting *every* cell under
>    `head_basket` radiation (pellets radiate bore-wide) → gate a no-op.
>    Fixed by making the exemption **topology-aware** (forward_plenum only).
>
> **Results (frozen transport, canon knobs, v_flame=3, kappa_zn=1):**
> | Motor | baseline | F | Z | F+Z |
> |---|---|---|---|---|
> | Chunc P_peak [MPa] | 16.93 | 10.66 | 11.50 | 10.62 |
> | Chunc ign spread | 0.8 ms | 205 ms | 0.9 ms | 214 ms |
> | Hasegawa P_peak | 6.24 | 6.41 | 6.19 | 6.41 |
>
> - **Phase F works for both topologies** (Hasegawa clean 481 ms front;
>   Chunc front runs, peak shifts 0.008 s→0.278 s — a gradual mass hump,
>   not an ignition spike). **Phase Z preserves the plateau**, Hasegawa
>   un-regressed.
> - **Both F and Z cut Chunc's spike ~35%** but **do NOT stack** (F+Z ≈ F):
>   once F spreads ignition over ~205 ms, Z-N's ~4 ms lag is irrelevant.
> - **Neither reaches 8.5 MPa.** Spike overshoot drops from ~2× (16.93) to
>   **~1.25× (10.6)** of the plateau — a substantial partial fix, residual
>   ~25% hump remains (candidate causes: erosive-augmentation ramp as the
>   front passes, conical geometry, Al₂O₃ two-phase lag; or push v_flame to
>   the 1 m/s floor / kappa_zn to 2 within physical bounds).
>
> **Ignition-model + energy-balance audit (2026-05-31).** Prompted by the
> partial result, we audited the ignition side from scratch. Findings,
> all by measurement: the Goodman kernel is accurate (2.3% vs exact PDE);
> the `T_ceiling` clip discards only 0.4% of injected enthalpy (0% at the
> spike peak) — not the cause; ignition is slaved to gas arrival (`T_surf`
> jumps 293 K → past `T_ignition` in one 0.5 ms step on front arrival).
> A convective wall-loss sink and a fix for the pyrogen-radiation
> double-count were added (energy-conservation completions, **always on**),
> but both are quantitatively negligible — the wall heat-loss *power* is
> small because `h_c` is laminar-low during the low-Re fill. **Conclusion:
> the bore gas reaching ~flame temperature is physically correct (both
> sources ~3000 K, wall losses tiny), so the spike is NOT an ignition or
> energy artifact — it is the erosive burn-rate over-response (Root B),**
> Ma's quasi-steady erosive firing instantly off the genuine peak-G at the
> smallest-bore condition (erosive fraction 0.13→0.61→0.13). The remaining
> lever is a transient/unsteady erosive closure. See the
> `project_ignition_model_audit` memory for the full trace.
>
> Harness: `examples/chunc_ignition_2x2.py`. Below is the original plan.

---

Companion to [`README.md`](README.md) (research synthesis). Attack order
fixed by user decision (2026-05-30): **Phase F (flame-front, bottom-up)
first → Phase Z (Z-N relaxation, top-down) → Phase FZ (combined +
re-calibration)**, evaluating each independently before combining.

**Diagnostic motor**: Chunc (`machbusterNew.ric`), `run_chunc_frozen.py`.
**Validation criterion (all phases)**: ignition spike drops toward
~8.5 MPa (no overshoot) **while plateau magnitude and taildown stay put**.
Cross-check the fired set (Hasegawa A, Zerox, BALLSstick) for no
regression.

Global reminders: delete `srm_1d/__pycache__/` after every `@njit` edit
(gotcha #1); `f_active` burnout ramp must scale both `C_burn` and
regression (gotcha #2); never hand-edit `.ric` files.

---

## Phase F — Flame-spread front (bottom-up, FIRST)

### Physical basis
Real grains light as a **front** propagating from the igniter outward by
**successive heating-to-ignition** (Peretz-Kuo-Caveny-Summerfield 1973):
each cell ignites only once the *local* convective+radiative flux from
the established flame just upstream has heated its surface to the
autoignition temperature. Front velocity (~1–10 m/s) is set by this
heating sequence and is **decoupled from the acoustic fill** (~300–1000
m/s). Our model currently lets every cell ignite from the bulk hot-gas
fill simultaneously (README §1, Root A).

### Design — front-gated ignition (recommended: F-1)
Add an **ignition-front position** `x_front` (m) and gate
`has_ignited[i]`: a grain cell may cross to ignited **only if the front
has reached it**. The front advances when the frontmost burning cell's
*immediate downstream neighbor* reaches `T_ignition` under heating that
includes the upstream flame (the existing adjacent-cell radiation +
convective `h_c` already provide this localized flux); cells far ahead of
the front are gated out, so they cannot ignite off the bulk fill.

- Front seed: at the cartridge/igniter location.
  - `forward_plenum` / `head_basket`: `x_front` starts at the head grain
    cell, spreads fore→aft.
  - `aft_basket`: seed at the aft cartridge, spreads aft→fore.
- Advance rule (parameter-free): the next unignited neighbor heats via
  its *own* Goodman ODE (already in
  [`_goodman_ignition_sources_and_mass`](../../simulation.py) ~L1000),
  but the **ignition gate is ANDed with `x within front`**. When that
  neighbor ignites, `x_front` steps to it. This yields a successive
  heating-to-ignition velocity emergent from the physics — no prescribed
  v_flame constant.
- Anti-stall floor: if no cell has ignited yet (cold start) the front
  seed cell must remain ignitable from the igniter heat flux directly
  (it already receives pyrogen surface flux / DeMar / radiation via
  [`_compute_pyrogen_heat_flux_arr`](../../simulation.py) L1509). Guard:
  the seed cell is always "within front."

### Code touch-points
- **Reuse the scaffold**: [`_compute_flame_front_augment`](../../simulation.py#L723)
  already tracks "cell immediately downstream of a recently-ignited
  cell." Phase B used it to *boost* `h_c` (failed — amplified the spike).
  Phase F **does not touch `h_c`**; instead it produces a per-cell
  `ignitable[i]` mask (or compares `x_centers[i]` to `x_front`).
- **Gate the ignition event**: in `_goodman_ignition_sources_and_mass`
  at the `_surface_has_ignited` check ([simulation.py:1092](../../simulation.py#L1092)),
  AND the condition with `ignitable[i]`. Keep computing `T_surf[i]`
  (cells ahead may pre-heat) but withhold the `has_ignited[i]=True`
  transition until the front arrives.
- **New state**: `x_front` scalar (or `ignition_front_idx` int) tracked
  across steps in `_run_time_loop`; `ignitable[N]` scratch.
- **Topology-aware seed**: branch on `topology_code` (already plumbed,
  L1230) + `cart_i_start/cart_i_end`.
- **New knob (off by default)**: `Propellant.flame_front_enabled: bool = False`
  so disabled runs reduce **byte-for-byte** to current behavior (mirror
  the `flame_spread_enabled` gating pattern, propellant.py:283).

### Tests
- `tests/test_flame_front_gate.py`: (a) disabled → identical
  `ignition_time[i]` to baseline; (b) enabled → `ignition_time[i]` is
  **monotonic** in distance from the seed (front, not simultaneous);
  (c) front never stalls on a canonical motor (all grain cells ignite);
  (d) aft_basket seeds at aft, spreads forward.
- Add an `ignition_time[i]` spread diagnostic to
  [`tools/ignition_diagnostics.py`](../../tools/ignition_diagnostics.py):
  baseline Chunc should show <5 ms spread (confirming Root A);
  Phase F should widen it to tens of ms.

### Validation
Run `run_chunc_frozen.py`. Expect the spike to drop materially as
mass-addition de-synchronizes. Record the new `ignition_time[i]` spread
and P_peak. **If the spike persists even with tens-of-ms ignition
spread, Root B (rate) dominates → Phase Z is the load-bearing fix.**
Either outcome is decisive (this is also the 4b-style diagnostic the
v0.7.3 candidate doc wanted).

### Risk / scope
~150–250 LOC. Main risk is the anti-stall guard and the aft→fore seed
direction. Does **not** require re-calibration by itself if it only
re-times ignition (mass conservation unchanged).

---

## Phase Z — Zeldovich-Novozhilov dynamic burn rate (top-down, SECOND)

### Physical basis
The burn rate has no memory (README §1, Root B). Z-N supplies the
condensed-phase thermal-inertia lag: `r` relaxes toward the quasi-steady
`r_qs` over the thermal-wave penetration/diffusion time. **Lumped form**
(README §4 — collapses Greatrix's two-layer model to the physical
relaxation; we do NOT solve the per-cell solid conduction PDE):

```
dr_dyn/dt = ( r_qs − r_dyn ) / τ_ZN
τ_ZN      = κ_zn · α_solid / max(r_dyn, r_floor)²
α_solid   = k_solid / (ρ_p · Cps)          (already computed, simulation.py:1299)
```

- `r_qs` = existing Ma-2020-augmented total from `compute_burn_rates`
  (Saint-Robert + Ma erosive). **Ma 2020 is the erosive model of record —
  Z-N lags its output, it does NOT replace it.** Do not substitute
  Lenoir-Robillard / King / Lengellé erosive correlations (inferior to
  Ma). One ODE on the **total** lags both the normal and erosive
  contributions — defensible because the surface thermal layer is the
  shared bottleneck for any surface-heat-flux change (pressure- or
  crossflow-driven). Lengellé 1975 (DOI 10.2514/3.49697) is cited only
  as evidence that the velocity-coupled (erosive) rate also relaxes —
  not as a model to adopt.
- `κ_zn ≈ 1` (parameter-free; Greatrix Eq. 14 lag form, τ set by
  penetration physics Eq. 19, confirmed τ ≈ 4 ms @ 5 mm/s by Lengellé).
- **Numerics**: use the analytic relaxation update for stiffness safety,
  `r_dyn = r_qs + (r_dyn − r_qs)·exp(−dt/τ_ZN)` (dt ≈ µs ≪ τ ≈ ms, so
  explicit is stable too, but the exp form can never overshoot). `r_floor`
  (~0.1 mm/s) avoids divide-by-zero at ignition.
- **Self-attenuation**: τ ∝ 1/r² shrinks as r rises → negligible lag at
  the high-r plateau → steady-state preserved (verify in tests).

### Code touch-points
- **New kernel** `_advance_zn_burn_rate(r_dyn, r_total, alpha_solid, kappa_zn, r_floor, is_burning, dt, N)`
  in [simulation.py](../../simulation.py); call it once per step **after**
  `compute_burn_rates` ([simulation.py:1389](../../simulation.py#L1389))
  and before source assembly.
- **Route `r_dyn`** into the two consumers, replacing `r_total`:
  mass-source assembly ([simulation.py:1114](../../simulation.py#L1114))
  and `advance_bore_regression` ([simulation.py:1362](../../simulation.py#L1362)).
  Keep `r_erosive` proportional (`r_erosive_dyn = r_dyn · r_erosive/r_total`)
  for diagnostics/snapshots.
- **New per-cell state**: `r_dyn[N]`, allocated at sim init in
  `run_simulation`, seeded to `r_qs` (or 0) at t=0; persisted across steps.
- **New knob (off by default)**: `Propellant.zn_enabled: bool = False`,
  `Propellant.kappa_zn: float = 1.0`. Disabled → `r_dyn ≡ r_total`
  byte-for-byte.
- Interaction with `tau_establishment`: Z-N supersedes it physically;
  keep `tau_establishment=0.0` whenever `zn_enabled` (document, don't stack).

### Tests
- `tests/test_zn_burn_rate.py`: (a) disabled → identical trace;
  (b) step change in `r_qs` → `r_dyn` relaxes with the correct τ
  (analytic check); (c) constant `r_qs` for ≫τ → `r_dyn → r_qs`
  (plateau preserved); (d) mass conservation across a transient
  unchanged within tolerance.

### Validation
Run `run_chunc_frozen.py` with Phase F **off**, Z-N **on** (independent
evaluation). Expect a shorter, rounder, lower spike. Compare the spike
reduction to Phase-F-only.

### Risk / scope
~150–250 LOC. Risk: over-/under-damping → one physically-bounded dial
`κ_zn ∈ [0.5, 2]` before escalating. Plateau depression if τ self-
attenuation is wrong (test c guards it).

---

## Phase FZ — Combined + re-calibration (THIRD)

### Evaluate independently, then combined
Produce four Chunc traces against `CHUNC_EXPERIMENTAL`:
1. baseline (both off) — the 16.85 MPa reference
2. Phase F only
3. Phase Z only
4. Phase F + Z

This isolates each mechanism's contribution and tests the README §3
hypothesis (F fixes spatial simultaneity, Z polishes the per-cell rate;
combined should be best). Put the 2×2 in a new
`examples/chunc_ignition_2x2.py` (mirror
`hasegawa_a_lhs_mode_transport_2x2.py`).

### Cross-motor regression + re-LHS
Existing knobs were calibrated to the spiky baseline. After F+Z:
- Re-run the fired set (Hasegawa A, Zerox, Chunc, BALLSstick —
  `project_fired_motor_set`) at default knobs; confirm no new
  divergence.
- If the spike fix shifts the optimum, re-LHS with **frozen** transport
  (per the v0.7.3.3 frozen-wins finding,
  `project_v0_7_3_post_phaseB_state`) and physical bounds enforced
  (roughness ≥ 15 µm, kappa ≈ 0.45 — `feedback_roughness_kappa_physical_bounds`).
  κ_zn stays ≈ 1; do not let the LHS treat it as a free knob.

### Decision points
- If **F alone** closes Chunc → Z is optional polish (still ship for the
  throat-erosion/end-of-burn transients it also improves).
- If **F alone** does not → confirms Root B dominance; Z is load-bearing.
- If **neither nor both** closes it → revisit the diagnostic (e.g.
  Pardue Al₂O₃ condensation, the deferred secondary candidate).

---

## Sequencing summary

| Phase | Root | Mechanism | New knob (default) | LOC | Gate |
|---|---|---|---|---|---|
| **F** | A (spatial) | flame-spread front gate | `flame_front_enabled=False` | 150–250 | ignition_time monotonic; spike drops; no stall |
| **Z** | B (temporal) | Z-N lumped relaxation | `zn_enabled=False`, `kappa_zn=1.0` | 150–250 | relaxes at τ=α/r²; plateau preserved |
| **FZ** | A+B | combine + re-cal | — | glue + LHS | Chunc spike→~8.5 MPa, plateau+tail intact; fired set no-regress |

Both knobs default **off** → tag-time behavior is unchanged until a motor
opts in, consistent with `feedback_defaults_reflect_majority_use`.
