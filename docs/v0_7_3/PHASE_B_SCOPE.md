# v0.7.3 Phase B Scope — heat-flux completeness for uncontained ignition

**Status (2026-05-25)**: scoping draft pending user review. Built on
the Phase A.3 diagnostic D finding (see `TASKS.md`): bore stays
stagnant under uncontained topologies because Bartz `h_c` is
Reynolds-limited and three coupled gaps prevent surface ignition
despite the bore being at T_flame from t=0.

## Diagnostic recap (audit-corrected, 2026-05-25)

The aft_basket heatmap at t≤0.5s shows:
- `T_gas` ≈ 2800 K everywhere from t=0 (v0.7.0/v0.7.1 IC sets
  `T_initial_gas = rep_tab.T_flame`)
- `u` ≈ 0 m/s (no axial flow; brief transients <7 m/s during
  pyrogen burn)
- `T_surf` rises 300→440 K in 437 ms (never reaches T_ignition=850 K)
- `P` stays at 0.10 MPa (no pressure gradient)

**The audit-revised gap diagnosis** (after reading the existing
heat-flux pipeline in `_goodman_ignition_sources_and_mass`):

1. **IC short-circuit**: `T_gas = T_flame` at t=0 across the bore
   means pyrogen mass injection at T_flame_pyrogen does not create
   a temperature gradient. Without a T gradient there is no
   density gradient, no pressure gradient, no flow. Forward_plenum
   hides this because the choked orifice produces a localized high-
   P injection at cell 0 that creates flow regardless of the bore
   IC. (B.0 below)
2. **Cell-to-cell radiation is gated on `is_burning[i]`** (only
   propellant-burning cells emit), so pyrogen-hot cells with
   T_gas=T_flame_pyrogen don't contribute even though
   `Propellant.radiation_emissivity` would let them in principle.
   (B.2 below)
3. **No pellet-to-surface radiation pathway**. Real BKNO3/MTV
   pellets radiate ~50-100 W/cm² to adjacent propellant; we have
   no representation of the pellets AS emitters. Forward_plenum
   covers this via DeMar at cell 0 only; uncontained has nothing.
   (B.4 below — primary structural fix)

**Audit-rejected gaps** (originally proposed, then withdrawn):

- ~~h_c floor missing~~ — already exists.
  `gnielinski_nusselt` returns `max(Nu, 3.66)` (the laminar
  fully-developed limit). At our diagnostic conditions, the floor
  delivers ~15 W/cm² which exactly matches the observed 140K
  T_surf rise over 437ms. The floor is correctly bounding heat
  transfer; the issue is that 15 W/cm² isn't enough flux to
  ignite within the pyrogen burn window. The fix is more flux,
  not a higher floor.
- ~~No gas-to-surface radiation~~ — cell-to-cell already IS
  gas-to-surface (it uses `T[i±1]`, the gas T of neighbors); it's
  just gated wrong. See B.2.

## Phase B scope

### B.0 — IC fix (CRITICAL, must land first)

**Diagnosis**: `T_initial_gas = rep_tab.T_flame` is a v0.7.0
"numerical-stability shortcut" preserved into v0.7.1 (see
`docs/v0_7_1/DESIGN.md` §3, simulation.py:1842-1854). The
docstring explicitly says "physically-correct ambient
initialization is a follow-up improvement". This is now the
v0.7.3 follow-up.

**Change**:
- Default `T_initial_gas = T_ambient` (300 K or P_ambient/(R·T))
  rather than `rep_tab.T_flame`.
- Preserve override via `initial_gas_temperature=` kwarg.
- For backward compatibility on existing motor calibrations
  (Hasegawa A, Zerox), test whether the change breaks
  forward_plenum traces. If it does, document the regression and
  decide whether to bake in `initial_gas_temperature=T_flame` per-
  motor or accept the recalibration.

**Risk**:
- The v0.7.0 numerical-stability rationale wasn't documented in
  detail. Likely was about avoiding cold-bore Mach number /
  acoustic issues at startup; PISO + adaptive CFL should handle
  this now (v0.7.0 had different stability controls). Need to
  test, not just change.
- Will likely change forward_plenum ignition traces because the
  cell-0 pyrogen jet now meets cold gas instead of hot gas. Could
  affect Hasegawa A v0.7.0 calibration. The post-Phase-B re-
  calibration is in scope.

**Test plan**:
- Run `hasegawa_motor_a.py` (forward_plenum) at T_initial_gas=300
  vs 2800. Diff the pressure traces and verify the simulation
  doesn't trip the collapse detector. If it diverges, either
  shrink the dt or document the recalibration requirement.
- Run `hasegawa_motor_a_aft_basket.py` at the new IC and inspect
  the u_cell heatmap for the expected T-gradient-driven flow.

**Files**:
- `simulation.py` L1842-1854 — change default `T_initial_gas`
- Update relevant docstrings + `docs/v0_7_1/DESIGN.md` §3
  follow-up note

**Scope**: ~50 LOC + ~3 LHS sanity-check tests.

### B.1 — h_c conductive floor [DROPPED — ALREADY EXISTS]

Audit (2026-05-25) found that `gnielinski_nusselt` already returns
`max(Nu, 3.66)` — the fully-developed laminar limit. At our
diagnostic conditions (k_gas=0.5 W/m/K, D_hyd=30 mm):
`h_c_floor ≈ 61 W/m²/K`, `q ≈ 15 W/cm²` at ΔT=2500 K, which
matches the observed 140 K T_surf rise over 437 ms exactly. The
model is correctly predicting that floor-Bartz alone isn't enough
to ignite in pyrogen-burn time. **Fix is more flux (B.4/B.5),
not a floor change.**

### B.2 — extend radiation_emitter criterion to pyrogen-hot cells

**Audit-corrected diagnosis**: Cell-to-cell radiation `σ·ε·(T_gas
[i±1]^4 - T_surf[i]^4)` ALREADY EXISTS in
`_goodman_ignition_sources_and_mass` (see simulation.py:826-871).
But `radiation_emitter[i] = is_burning[i]` (simulation.py:773),
which means **only propellant-burning cells emit**. In our
uncontained diagnostic with T_gas[cartridge_cells] ≈ T_flame_pyrogen
≈ 2700 K, those cells should be emitting strongly — but they
aren't, because their `is_burning` flag is False (no propellant
ignited yet).

The fix is a one-line change to the emitter criterion: cells with
substantial pyrogen species mass fraction should also count as
emitters. Pyrogen-burning cells produce gas at T_flame_pyrogen
which radiates exactly like propellant gas at T_flame_propellant.

**Change**:
- Extend `radiation_emitter[i] = is_burning[i]` to
  `radiation_emitter[i] = is_burning[i] OR (Y_species[i, IGNITER]
  > Y_emit_threshold)` where `Y_emit_threshold = 0.5` (cell is
  majority-pyrogen gas).
- This is a single-line addition in the per-step radiation_emitter
  refresh; no new API surface.
- Caveat: requires `Propellant.radiation_emissivity > 0` for any
  effect. Hasegawa A baseline ships with `radiation_emissivity=0`
  per the v0.7.0 numerical-stability rationale
  (`_default_radiation_emissivity` returns 0.0 in
  `openmotor_adapter.py:291`). So this fix is harmless to existing
  calibrations but only activates if the user opts in to
  radiation.

**Files**:
- `simulation.py` — update `radiation_emitter` refresh line
- `tests/test_radiation_emitter.py` (new) — verify pyrogen-hot
  cells emit when Y_pyrogen > threshold

**Scope**: ~30 LOC + ~3 new tests.

### B.2b — pellet-to-surface radiation (per-cartridge-cell)

**Audit finding**: The most under-modeled pathway for amateur HPR
ignition. Real pyrogen pellets (BKNO3, MTV) radiate strongly to
adjacent propellant surfaces; this is the documented secondary
ignition mechanism per Sutton 9e §15.3 ("primarily convective" =
acknowledges secondary radiative role) and the Sandia 2022 LDRD
peak flux of ~1 GW/m² for pyrotechnic emission. Currently the
model has no representation of the pellets AS radiative emitters
— only the gas in the cartridge cells emits, and only if B.2 is
applied.

This is functionally **B.4** in the original scope but reframed:
"per-cell DeMar" is really "per-cell pellet radiative flux" — see
B.4 below.

### B.3 — pyrogen form archetypes

**Diagnosis**: `build_pyrogen_chamber` defaults `A_burn_initial =
4·π·r²` for a single sphere of the pyrogen mass (~50 cm² for 4.8 g
BKNO3). Real amateur pyrogen is typically pellets or slivers with
much higher total surface area (5-20×). Per the user-supplied
archetype memory, there are three archetypes with characteristic
burn timescales.

**Change**:
- New `Pyrogen.form: str = 'pellets'` (default reflects amateur
  HPR majority) with values `'powder' | 'pellets' | 'chunks'`.
- Multipliers on `A_burn_initial` per archetype:
  - `'powder'`: ×20 (fine powder, ~1-10 ms)
  - `'pellets'`: ×5 (standard BKNO3/MTV pellets, ~100 ms)
  - `'chunks'`: ×1 (single-sphere geometric, ~1000 ms — the
    pre-Phase-B default for backward compat documentation)
- Override mechanism: explicit `pyrogen_burn_area=` kwarg
  continues to wins (literal user-specified value).

**Test plan**:
- Verify the three archetypes give monotonically-decreasing burn
  timescales for a fixed mass.
- Verify `'pellets'` default + `pyrogen_burn_area=None` gives the
  new multiplied default.
- Verify `pyrogen_burn_area=X` (explicit) overrides the multiplier.

**Files**:
- `propellant.py` — add `Pyrogen.form` field
- `openmotor_adapter.py` — `build_pyrogen_chamber` consumes
  `pyrogen.form` for default A_burn

**Scope**: ~50 LOC + ~3 new tests.

### B.4 — pyrogen-to-surface heat delivery (mutually exclusive modes)

**Energy-balance audit finding**: DeMar's time-averaged heat flux
(69.4 cal/cm²/s for BPNV) is the lumped TOTAL flux to a
calorimeter surface adjacent to a burning pellet stack —
convective + radiative + particle. Estimating the radiative
component alone: σ·0.7·2800⁴ ≈ 240 W/cm² (~83% of DeMar's 290
W/cm²). If we apply BOTH DeMar AND a separate pellet-radiation
model, we **double-count the radiation**. They must be mutually
exclusive at the implementation level.

**Form-extensibility audit**:

| Form | DeMar | Pellet-radiation |
|---|---|---|
| Pellets (BKNO3, MTV) | ✓ measured 2021 | ✓ σ·ε·T⁴, lit-bounded ε |
| Powder (Cu/Al, Fe/Al thermite) | △ Sandia 2022 has peak/sustained but no time-avg | ✓ Same model, T_flame_thermite ~2300 K |
| Chunks (APCP basket sticks) | ✗ Not measured (APCP is propellant, not pyrogen) | ✓ Same model, T_flame_APCP ~3200 K |

**Change**: Single new field `Pyrogen.heat_delivery_mode: str =
'demar'` (default for pellet pyrogens) with values:

- **`'demar'`** — apply `heat_flux_cal_cm2_s` (already in YAML) as
  a constant time-averaged flux to ALL cartridge cells uniformly
  (analog of forward_plenum's cell-0 DeMar, distributed across
  the cartridge). Off when pyrogen is consumed.
  - Trusted for: BKNO3/MTV pellets where DeMar measured directly.
  - Defensibly pinned per-pyrogen YAML.
  - **Empirically grounded**, but doesn't capture distance falloff
    beyond the cartridge cells.

- **`'radiation'`** — Stefan-Boltzmann pellet radiation model:
  `q[j] = σ · ε · T_flame⁴ · F_ij · exp(-(x_j - x_i)/L_atten)`
  summed over emitting cartridge cells `i` and receiving cells
  `j`, where:
  - `F_ij = A_port[j] / (4·π·(x_j - x_i)² + A_port[j])` —
    geometric view factor (saturates to 1 for adjacent cells,
    falls as ~1/d² far field; user-corrected from my earlier
    exponential-falloff proposal).
  - `L_atten` is the radiation absorption length (captures
    aggregate gas + particle absorption without modeling each
    separately). New field `Pyrogen.radiation_absorption_length_m`
    defaults: clean pyrogen (BKNO3) ~1 m; aluminized exhaust
    ~0.1 m.
  - `ε` is `Pyrogen.pellet_emissivity` (new field; defaults 0.7,
    lit range 0.5-0.9).
  - Off when pellet `T_emitter` falls below ambient (pyrogen
    consumed).
  - Trusted for: powder (thermite) and chunks (APCP) where DeMar
    data is unavailable.
  - **Physically modeled** with one tunable per pyrogen
    (`pellet_emissivity`) and one optional attenuation length.
    Reuses already-pinned constants (Stefan-Boltzmann, T_flame).

- **`'none'`** — neither pathway applied. Recovers Phase A.2
  byte-for-byte. Backward-compat regression target.

**Default mapping by pyrogen form** (set in YAML, not
auto-derived):
| Form | Default mode | Rationale |
|---|---|---|
| Pellets (BPNV, MTV) | `'demar'` | DeMar measured; empirical wins |
| Powder (thermite) | `'radiation'` | No DeMar data; physical model only |
| Chunks (APCP) | `'radiation'` | No DeMar data; physical model only |

**Test plan (independent A/B per the user's framing)**:
- Verify `'none'` mode recovers Phase A.2 byte-for-byte (regression).
- Verify `'demar'` mode applied uniformly across cartridge cells
  delivers `heat_flux_cal_cm2_s` × A_burn_per_cell to local T_surf
  (per-cell DeMar distribution).
- Verify `'radiation'` mode: F_ij = 1 at adjacent cell, ~1/d² far,
  sum bounded by total emitted power. exp(-d/L_atten) reduces
  far-cell flux by expected factor.
- **Run ISP Super Loki with `'demar'` AND with `'radiation'`
  separately**; tabulate which one matches the 8.8 MPa
  experimental plateau better.
- **Run Hasegawa A aft_basket with both modes**; check if either
  reaches ignition cascade, and if so whether the simultaneous-
  ignition artifact persists under reversed topology.
- Cross-check forward_plenum motors are byte-for-byte regressed
  (this field doesn't apply to `forward_plenum` topology).

**Files**:
- `propellant.py` — `Pyrogen.heat_delivery_mode`,
  `pellet_emissivity`, `radiation_absorption_length_m`
- `motors/pyrogens/bpnv.yaml`, `mtv.yaml` — set `heat_delivery_mode: 'demar'` default
- `simulation.py` — new Numba kernel
  `_compute_pyrogen_radiation_flux(...)`; topology-branch
  enables the chosen mode for uncontained topologies
- `tests/test_pyrogen_heat_delivery.py` (new) — covers all three
  modes + mode-switching regression

**Scope**: ~200 LOC + ~10 new tests. Replaces the original 250
LOC of separate B.4 + B.5 with a unified mechanism that side-
steps the double-counting.

**Energy-balance double-counting safeguard**: When
`heat_delivery_mode == 'demar'`, the radiation model is OFF (and
vice-versa). Hard-coded mutual exclusivity at the topology-
branch level in `_run_time_loop`. If a future user wants to
combine them, they explicitly opt in via a separate (deferred,
not in scope) `'both'` enum value with a documented warning.

### B.5 — DROPPED (merged into B.4)

The original B.5 (view-factor radiative falloff) is now part of
B.4's `'radiation'` mode. Removing as a separate item.

### B.6 — re-run validation + recalibration

**Goal**: Verify the four fixes together let uncontained topologies
ignite and validate against Super Loki experimental.

**Steps**:
- Rerun `examples/ISP_Super_Loki.py` with head_basket + new
  defaults. Target: P_peak within 2× of 8.8 MPa experimental at
  default knobs (not a tight fit — that's a v0.7.4 LHS task).
- Rerun `examples/hasegawa_motor_a_aft_basket.py`. Target:
  reaches an ignition cascade (P > 1 MPa, T_surf > 850 K
  somewhere). The DIAGNOSTIC question — does the simultaneous-
  ignition artifact persist under reversed topology? — is finally
  answerable.
- Cross-check forward_plenum motors: `hasegawa_motor_a.py`,
  `Zerox_test.py`, `BALLSstick.py`. Confirm B.0/B.1/B.2 don't
  break the existing v0.7.2 calibrations or document the
  required recalibration.
- Update example artifacts + the v0.7.3 docs with the validated
  pressure traces.

### B.7 — docs + memory + tag

- Update `docs/v0_7_3/TASKS.md` Phase B section with the
  scoped/landed work.
- Append v0.7.3-phaseB DEVNOTES API-breaking-change log entry.
- Update CLAUDE.md banner + roadmap.
- Update `project_v0_7_3_phaseA_state` memory or create a new
  `project_v0_7_3_phaseB_state` memory.
- Tag `v0.7.3-phaseB`.

## Ordering / dependencies

```
B.0 IC fix (must be first; everything else depends on this)
   ↓
B.2 radiation_emitter gating fix (independent, single-line; can
                                    land alongside B.0)
   ↓
B.3 form archetypes (pyrogen YAML defaults; cheap)
   ↓
B.4 unified pyrogen-to-surface heat delivery
    ('demar' | 'radiation' | 'none' mode enum)
   ↓
B.6 validation: run BOTH 'demar' and 'radiation' independently;
                tabulate which fits Super Loki + answers the
                Hasegawa A aft_basket diagnostic question
   ↓
B.7 docs + tag
```

Total estimated scope: **~310 LOC** + ~17 new tests.
- B.1 dropped (already exists: `max(Nu, 3.66)` floor)
- B.2 reduced from 80 → 30 LOC (one-line gating fix)
- B.4 + B.5 merged into one mechanism (250 → 200 LOC, removes
  double-counting risk)

3-5 sessions including independent A/B validation of 'demar'
vs 'radiation' modes against ISP Super Loki experimental.

## Risk flags / open questions

1. **B.0 may regress forward_plenum calibration**. If the IC
   change shifts Hasegawa A's peak pressure or timing
   significantly, the v0.7.0 knobs (roughness=37.1µm, kappa=0.45)
   may need re-LHS. This is a calibration risk but not a code
   risk.
2. **B.2 gas emissivity default**. Enabling at 0.3 by default
   could over-predict ignition spike in motors that were
   calibrated without it. Safest: ship at 0.0 default initially
   and add the value to per-motor YAMLs as the calibration
   progresses.
3. **B.5 view-factor numerical stability**. The 1/distance² in
   `F_ij` could create stiff coupling if too many cells
   contribute. Numba kernel needs efficient bounded-radius cutoff
   (drop contributions where `F_ij < 1e-6`).
4. **`'pellets'` becoming the new default in B.3** could change
   existing motor traces if their pyrogen mass was tuned to the
   old single-sphere A_burn. This affects only the auto-sized
   defaults; any motor with explicit `pyrogen_burn_area=X` is
   unaffected.

## What to do at burst-diaphragm (deferred per user feedback)

Per `[[feedback_defaults_reflect_majority_use]]` — burst diaphragm
is industry-common but amateur-rare. Defer to v0.7.4 or v0.8 as an
opt-in `Nozzle.plug_burst_pressure_pa` / `plug_eject_time_s`
field. Not in Phase B scope.
