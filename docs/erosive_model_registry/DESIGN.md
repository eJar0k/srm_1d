# Erosive burn-rate model registry — DESIGN

**Goal:** make the erosive-burning closure a selectable, registered model —
mirroring the v0.8.0 solver-plugin pattern (`srm1d_plugin.py`) one level down.
Ships two models: `ma2020` (the current Ma heat-transfer closure — default,
byte-identical to today) and `mukunda_paul` (MP 1997 Eq. 12 universal law).
The registry is the extension point for future closures (Lenoir-Robillard,
Greatrix, …).

**Motivation:** the 2026-07/08 spike-reopen sessions
(`docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` §10–§11): the high-L/D startup
over-spike equals Ma's high-g excess over MP per-motor; Ma's own 2022 star-motor
paper (`docs/references/energies-15-01564.pdf`) runs BOTH models on one code as
a validation axis — precedent for exactly this feature. The A/B axis also lets
users trade Ma's no-tuning shipping default against MP's high-g behaviour
deliberately (never silently — default stays `ma2020`).

## Scope and non-goals

- Swappable piece = the **erosive closure only**. The Saint-Robert tab lookup
  (`saint_robert_from_tabs`) is propellant *data*, shared by all models — every
  model returns `r_total >= r0` built on the same `r0`.
- The **ignition-gate h_c path is untouched** in all models
  (`simulation.py` Goodman heating uses Gnielinski/Haaland independently).
  `roughness` / `kappa` keep their gate role under every model; under
  `mukunda_paul` they simply stop influencing the *burn rate*. Document this
  in the kwarg docstrings — it is what makes model A/B a clean isolation of
  the erosive closure (identical ignition sequence both arms).
- Model choice is **solver config, not propellant data** (matches the
  `srm1d_plugin.py` precedent note: "transport variant is a per-propellant
  property … not a run-control param" — the erosive model is the opposite).
  No propellant-editor / .ric-propellant change.
- Z-N (Phase Z) composes unchanged: it relaxes whatever quasi-steady
  `r_total` the active model produces.
- `erosive_mp_cap` remains a **`ma2020`-only option** (its whole point is to
  bound Ma by MP; under `mukunda_paul` it is meaningless).

## Architecture

### 1. Registry (Python level, `srm_1d/burn_rate.py`)

```python
@dataclass(frozen=True)
class ErosiveModel:
    name: str          # registry key, e.g. 'ma2020'
    code: int          # int crossing the @njit boundary
    label: str         # human string for GUI enum
    description: str   # one-liner + citation
    params: tuple      # run_simulation kwargs this model consumes

EROSIVE_MODELS = {
    'ma2020':       ErosiveModel('ma2020', 0, 'Ma 2020 (heat-transfer)',
                                 '...', ('roughness', 'kappa', 'erosive_mp_cap')),
    'mukunda_paul': ErosiveModel('mukunda_paul', 1, 'Mukunda-Paul 1997 (universal law)',
                                 '...', ()),
    'saint_robert': ErosiveModel('saint_robert', 2, 'Saint-Robert only (erosive off)',
                                 '...', ()),
}

def resolve_erosive_model(name) -> int   # ValueError lists registered names
def register_erosive_model(model)        # future extension hook
```

Numba cannot dispatch on objects — the **int `code` is the only thing that
crosses the JIT boundary** (precedent: `topology_code`,
`heat_delivery_mode_code` in `_run_time_loop`).

### 2. Kernel dispatch (`burn_rate_cell`)

Signature gains `erosive_model_code` (after `erosive_mp_cap`). Branch:

- `code == 0` — current Ma path verbatim (bisection + optional MP cap).
  MUST remain byte-identical at defaults (gate G1).
- `code == 1` — closed form, no iteration:
  `G = Re_local * mu_gas / D_hyd`;
  `r_total = r0 * mukunda_paul_eta(G, r0, D_hyd, mu_gas, rho_p)`;
  `r_erosive = r_total - r0`. (`mukunda_paul_eta` already exists and is
  PDF-verified.) The existing `Re_local < 100` early-return stays above the
  branch — for MP that regime is `g << 35 → eta = 1 → r0` anyway (identical
  result, cheaper).
- `code == 2` — **pure Saint-Robert**: `return r0, 0.0` immediately after
  `r0` is computed (before the Re guard — no flow dependence at all). This is
  the first-class "erosive off" mode. The existing
  `diagnostic_disable_erosive` flag is RETAINED unchanged (it post-subtracts
  in the time loop and reaches the same r0) — gate G2c pins the equivalence.
  Beyond diagnostics, `saint_robert` is the decomposition tool for the
  spike/plateau decoupling investigation (REOPENED §10.5/§11.4): plateau_SR
  vs plateau_Ma vs plateau_MP apportions how much plateau lift the erosive
  closure must supply per motor — i.e., whether Chunc's plateau preference
  for Ma is model fidelity or Saint-Robert a/n calibration debt.

`compute_burn_rates` passes the code through unchanged.

### 3. Threading (mechanical)

`_run_time_loop` gains `erosive_model_code` (int64) in the same argument
group as `erosive_mp_cap`; its single call site in `run_simulation` passes
`resolve_erosive_model(erosive_model)`.
`run_simulation(..., erosive_model='ma2020')` — string kwarg, resolved +
validated at the Python boundary:

- unknown name → `ValueError` listing registered names;
- `erosive_mp_cap > 0` with `erosive_model != 'ma2020'` → `ValueError`
  (the cap is a Ma option — fail loudly, never ignore silently).

`run_from_ric` needs nothing: `**sim_overrides` already forwards.

### 4. Provenance

Stamp `result['summary']['erosive_model'] = name` so traces, probes, and CSV
exports are self-describing (the cross-motor probes compare models — unlabeled
results are how mistakes happen).

### 5. openMotor GUI leg (separate phase)

`_transient_config_schema()` in `srm1d_plugin.py` gains
`EnumProperty('Erosive Model', ...)` (exists in `motorlib.properties`,
verified) with values = registry keys, default `ma2020`, plus
`FloatProperty('MP Erosive Ceiling', '', 0.0, 2.0)` seeded 0.0 for
`erosive_mp_cap`. Schema keys already pass straight through
`simulate() → simulate_motor(**overrides) → run_simulation`, and persist
per-motor via `Motor.solverConfigs['srm_1d-transient']`
(verified: `motorlib/motor.py` getDict/applyDict; `uilib/simulationManager.py`
per-motor override; `uilib/preferencesManager.py` global fallback;
`uilib/fileIO.py` setdefault on load). **Zero new plumbing** — but T3 must
verify the fork's config screen renders `EnumProperty` generically before
relying on it, and no unit conversion applies (string + dimensionless).

## Model contract (for future additions)

A new model must supply: (1) a registry descriptor (name/code/label/
description-with-citation/params); (2) an `@njit`-compatible kernel branch
returning `(r_total, r_erosive)` with `r_total >= r0 > 0` semantics;
(3) unit gate tests (closed-form values, off-state identity);
(4) a cross-motor validation entry via the
`docs/v0_7_4/probes/mp_cap_crossmotor.py` harness (model axis);
(5) a DEVNOTES API-break entry if signatures change.

## Compatibility

`burn_rate_cell` / `compute_burn_rates` / `_run_time_loop` signatures gain one
argument — log in DEVNOTES "API Breaking Changes Log" (hard breaks fine per
project convention). Default (`ma2020`, cap 0.0) must be result-identical to
pre-change (G1). Numba: clear `srm_1d/__pycache__` after kernel edits.

## Validation gates

| gate | check |
|---|---|
| G1 | default identity: Chunc canonical baseline P_peak 11.47 MPa @ 16.5 ms / plateau 8.22 unchanged |
| G2 | unit: MP mode `r_total == r0·eta_MP` on synthetic cells; `g<35 → r0`; `Re<100 → r0` |
| G2b | unit: `saint_robert` returns `(r0, 0.0)` at any Re/P combination |
| G2c | integration: `saint_robert` run == `ma2020 + diagnostic_disable_erosive=True` run (same P trace) on a short Chunc run |
| G3 | unit: on a cell with `eta_Ma > eta_MP`, MP mode == Ma mode with `erosive_mp_cap=1.0` |
| G4 | cross-motor probe, model axis {ma2020, ma2020+cap1.1, mukunda_paul, saint_robert}: Chunc ratio ≈ {1.395, 1.135, ~1.07, no-early-spike}; BALLSstick pk/tr ≈ {1.93, 1.70, ~1.52, —}; Hasegawa peak ≈ {6.14, 6.16, ~6.03, —} MPa. saint_robert rows are *recorded*, not gated (they are the decoupling-decomposition data) |
| G5 | full pytest green (408+ known-good; the 2 `BALLSstick VANILLA.ric` scan failures are pre-existing) |
| G6 | (GUI phase) EnumProperty renders; solverConfigs round-trips through save/load |

G4's `mukunda_paul` expectations derive from the measured cap=1.0 runs
(min(Ma,MP)=MP nearly everywhere active on these motors) — small deviations
where MP > Ma are expected and worth recording.

## Execution note

Tasks are chunked in `TASKS.md` for inexpensive code-writing agents: each task
names exact files, signatures, and acceptance gates, and none requires
architectural judgment. Review order: T1 → T2 → (T3, T4 parallel).
