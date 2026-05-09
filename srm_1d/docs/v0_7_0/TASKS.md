# srm_1d v0.7.0 — Task Breakdown

Concrete file-level tasks for implementing the [DESIGN.md](DESIGN.md)
hot-gas plenum igniter. Tasks are ordered for incremental commits — each
should leave 107/107 pytest passing.

## Phase 1 — Pyrogen plenum (standalone, no main-motor integration)

Goal: a `PyrogenChamber` class that integrates ODEs (1)-(4) from
DESIGN.md correctly against an analytical closed-bomb test case. No
main-motor coupling yet.

### Task 1.1: `srm_1d/propellant.py` — add `Pyrogen` dataclass

```python
@dataclass
class Pyrogen:
    """Pyrogen propellant properties. Mirrors Propellant for the
    igniter chamber. Single-tab (no multi-tab support; pyrogens burn
    over a narrow pressure range)."""
    name: str
    a: float           # Saint-Robert coefficient, m/s/Pa^n
    n: float           # exponent
    rho: float         # propellant density, kg/m³
    T_flame: float     # adiabatic flame temperature, K
    M: float           # molecular weight of products, kg/mol
    gamma: float       # specific heat ratio
    impetus_W: float = 0.0   # measured impetus, psi·in³/g (DeMar units, optional)
```

Acceptance: dataclass importable, no breakage of existing `Propellant`.

### Task 1.2: `srm_1d/motors/pyrogens/bpnv.yaml` + `mtv.yaml`

Sibling to `<motor>.transport.yaml`. Schema:

```yaml
name: BPNV
a: ...        # m/s/Pa^n
n: 0.5        # typical for BKNO3
rho: 1700     # kg/m³, dipped/molded BPNV
T_flame: 2800 # K, approximate
M: 0.030      # kg/mol, products avg
gamma: 1.25
impetus_W: 5000  # psi·in³/g, DeMar Table
```

Values from DeMar 2021 + literature. Document any guess explicitly in
header comments.

### Task 1.3: `srm_1d/igniter_plenum.py` — new module

```python
@dataclass
class PyrogenChamber:
    pyrogen: Pyrogen
    m_pyrogen_initial: float    # kg
    A_burn_initial: float       # m², initial burning surface
    A_throat: float             # m², plenum vent area
    V_plenum: float             # m³, free volume
    burn_law: str = '0d'        # '0d' (sphere-equivalent) | 'end_burning' | 'cylindrical'

@njit(cache=True)
def _step_plenum_ode(state, pyrogen_params, ..., dt, P_main):
    """RK4 step of (m_p, m_ig, T_ig). Returns new state +
    mdot_choke + mdot_into_main."""
    ...

@njit(cache=True)
def _choked_orifice_mdot(P_ig, T_ig, A_t, gamma, R, M, P_main):
    """Choked / subsonic-fallback. Returns mdot."""
    ...
```

Acceptance: in a unit test, integrate the plenum with `mdot_choke = 0`
(closed bomb) for 50 ms — total mass conserves to <0.1%, energy
conserves consistent with `(-dm_p/dt)·c_p·T_flame` integrated.

### Task 1.4: `srm_1d/tests/test_igniter_plenum.py`

Tests:
1. **Closed bomb** (no outflow): mass conservation, monotone P_ig rise.
2. **Steady choked outflow** (constant `A_burn`): asymptotic
   `P_ig_steady = (rho·a·A_burn·c*/A_t)^(1/(1-n))` matches analytic.
3. **Burnout**: `m_p → 0` triggers `mdot_choke → 0` smoothly.
4. **Subsonic fallback**: artificially lower P_ig + raise P_main below
   choke threshold, verify subsonic formula (6) activates.
5. **Sutton sizing default** (`m = 0.12·V_F^0.7`): given a representative
   motor, default mass is computed correctly.

Target: 5+ tests pass. No main-motor integration yet.

## Phase 2 — Goodman solid-conduction sub-solver

Goal: a per-cell `T_surf(t)` and `δ(t)` Goodman solver, validated
against analytical constant-flux solution.

### Task 2.1: `srm_1d/solid_thermal.py` — new module

```python
@njit(cache=True)
def _step_goodman_ode(delta, T_surf, h_c, T_gas, T_initial, alpha, k_solid, dt):
    """RK4 step of Eq. (7); recompute T_surf from Eq. (8)."""
    ...

@njit(cache=True)
def _compute_T_surf(delta, h_c, T_gas, T_initial, k_solid):
    """Algebraic relation (8)."""
    return (3*k_solid*T_initial + h_c*delta*T_gas) / (3*k_solid + h_c*delta)
```

Initial condition: `δ(0) = 1e-6` (1 µm) to avoid singularity. The
ODE has a `1/δ` term that's stiff initially; RK4 handles this if
`dt < ~1e-7 s` early on. May need adaptive substepping during
induction.

### Task 2.2: `srm_1d/tests/test_solid_thermal.py`

Test against analytical constant-q solution `T_s - T_i = q·2√(αt/π)/k`.
Goodman gives `T_s - T_i = q·2√(αt/3)/k`. Ratio is √(π/3) ≈ 1.023.
Verify Goodman matches the analytical to within ~3% over 0-100 ms for
typical APCP propellant (α ≈ 1e-7 m²/s, k ≈ 0.3 W/m·K, q = 1 MW/m²).

Test that `T_surf > T_ignition` triggers as expected.

Acceptance: 3+ tests pass.

## Phase 3 — Wire into the main solver

Goal: pyrogen + Goodman sub-solvers integrated into `_run_time_loop`,
replacing the placeholder. Pyrogen mass-flow injects at cell 0; per-cell
`T_surf` controls `is_burning[i]`.

### Task 3.1: `srm_1d/simulation.py` — replace igniter section

- Remove existing `igniter_mass`, `igniter_tau`, `ignition_ramp_tau`,
  `P_ignition` sim_kwargs (lines ~430-450)
- Add `pyrogen_chamber: PyrogenChamber`, `T_ignition: float` (default 850)
- Allocate per-cell `T_surf[N]` (init `T_initial`) and `delta[N]`
  (init `1e-6`) in setup
- In `_run_time_loop`:
  - Each step: integrate `_step_plenum_ode` to get `mdot_choke`, `T_ig`
  - Inject `mdot_choke` into cell 0 source (mass + enthalpy + momentum)
  - For each unignited cell: integrate `_step_goodman_ode` to update
    `delta[i]`, recompute `T_surf[i]`
  - Trigger ignition: `T_surf[i] > T_ignition` → `is_burning[i] = True`,
    `delta[i]` no longer integrated
- Update structured `summary` dict to include pyrogen stats (mass burned,
  duration, peak P_ig)

Acceptance: existing tests still pass when run with `pyrogen=None`
fallback or default Sutton-sized BPNV. New tests for pyrogen-driven runs.

### Task 3.2: `srm_1d/openmotor_adapter.py` — `pyrogen.yaml` parsing

Add `load_pyrogen(path)` mirroring `load_transport(path)`. Auto-discover
`<motor>.pyrogen.yaml` sibling next to `.ric` files. Override via
explicit `pyrogen=...` kwarg.

```python
result = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    pyrogen='bpnv',  # OR pyrogen=Pyrogen(...) explicit object
    pyrogen_mass=None,         # Sutton default if None
    T_ignition=850,
    ...
)
```

If `pyrogen=None` AND no sibling YAML found, raise informative error
suggesting either path.

### Task 3.3: Update `Zerox_test.py` and `hasegawa_motor_a.py` examples

Both currently use the v0.6.0 placeholder kwargs. Update to:

```python
result = run_from_ric(
    motor_path,
    pyrogen='bpnv',
    pyrogen_mass=None,    # use Sutton default
    T_ignition=850,
    # ... rest same as before, but igniter_* removed
)
```

The Zerox `igniter_tau = 127ms` FSI-proxy is GONE. Documentation in the
script's docstring should explicitly note v0.7.0 removed this knob.

## Phase 4 — Validation against Hasegawa A

Goal: spike overshoot drops to <10%, MSE to ≈ 0.10 MPa², no calibration
knobs left untreated.

### Task 4.1: Re-run `hasegawa_motor_a.py`

Compare new pressure trace against experimental. Expected behavior:
- Spike no longer overshoots by 25% — should be within 5-10%
- Plateau and tail unchanged from v0.6.0
- LHS rank-1 igniter knobs (mass, throat, T_ignition) settle on
  physical values

If spike is still over by >15%, re-examine the pyrogen mass/throat
defaults (likely needs Sutton's `0.12·V_F^0.7` to be slightly different
for amateur APCP — the original constant came from a 1971 industrial
data set).

### Task 4.2: Re-run Zerox LHS with v0.7.0 parameters

`zerox_lhs.py` bounds:
```python
BOUNDS = {
    'erosion_coeff_scale': (1.5, 3.5),    # locked in v0.6.0 calibration
    'a_scale':             (0.85, 1.10),
    'pyrogen_mass_scale':  (0.5, 2.0),    # multiplier on Sutton default
    'pyrogen_throat':      (1e-6, 5e-5),  # m², orifice/vent area
    'T_ignition':          (700, 950),    # K
    'kappa':               (0.30, 0.60),
}
```

`pyrogen_volume` and pyrogen choice remain fixed at defaults. Run the
same 280-main + 16×6-pinned LHS pattern as v0.6.0.

Verify P_ignition (a v0.6.0 knob) is gone. The new LHS should converge
to a tighter cluster since FSI-proxy degenerate combinations are no
longer in the search space.

### Task 4.3: DEVNOTES update

Replace the v0.6.0 "Igniter (exponential-decay model)" section with a
v0.7.0 section pointing to:
- `igniter_plenum.py` for the pyrogen chamber model
- `solid_thermal.py` for the Goodman conduction model
- `srm_1d/motors/pyrogens/` for pyrogen datasheets
- The Hasegawa A and Zerox calibration tables

Add the Sutton Eq. 15-4 default with units footnote.

## Phase 5 — Memory + git tag

### Task 5.1: Tag `v0.7.0`

After all 4 phases pass and pytest is green:

```
git tag -a v0.7.0 -m "Hot-gas plenum igniter model — replaces v0.6.0 exponential-decay placeholder"
git push origin v0.7.0
```

### Task 5.2: Memory updates

Update [project_hasegawa_calibration_state](../../) memory to reflect
v0.7.0 numbers. Update [project_zerox_calibration_state](../../) memory
with new pyrogen-based parameters.

Mark [project_v0_7_0_design](../../) memory as "implemented; pointer
preserved for v0.7.x extensions".

## What's deferred (NOT in v0.7.0)

- Squib stage (electric → BPNV ramp → pyrogen) — v0.7.2
- Lumped radiation `C_hc(x/L)` — v0.7.1
- Multi-species `Y_ig` passive scalar — v0.7.3
- Head-end primary motor architecture (the user's long-term goal) — v0.8.0
- Cavallini-style 6-species mixture + Godunov — v0.9.0+

## Estimated effort

- Phase 1 (plenum standalone): ~4 hours
- Phase 2 (Goodman): ~3 hours
- Phase 3 (integration): ~6 hours
- Phase 4 (validation): ~4 hours (Hasegawa) + ~2 hours wall (LHS)
- Phase 5 (docs): ~1 hour

**Total: ~20 hours** of focused dev time, plus ~2 hours LHS wall time.
Doable in a single multi-session dev push if Numba caches stay warm.
