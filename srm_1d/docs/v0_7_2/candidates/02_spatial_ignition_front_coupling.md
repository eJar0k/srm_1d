# Candidate 2: Spatial Ignition-Front Coupling

**One-line**: Replace independent per-cell Goodman ignition kernels
with a coupled formulation where each cell's convective heat-transfer
coefficient `h_c[i]` is augmented by the cumulative upstream mass
flux from already-burning cells, producing physically realistic
sequential bore ignition.

## Motivation

The structural artifact is that all bore cells currently cross
`T_surf > T_ignition` within a few ms of each other because each
cell's Goodman ODE sees only its own local h_c, which at ignition
onset is similar across cells (driven by the pyrogen plenum
near-equilibrium pressure rise). Real SRMs ignite via a propagating
front: cells closer to the igniter outlet light first, the upstream
combustion gases flow downstream adding mass flux, the local h_c at
the next unignited cell goes up, that cell lights, and so on.

Literature consensus (Salita 2001, Han 2017, Kashiwagi 1982,
SPINBALL) is that **the Peretz per-cell critical-T criterion is fine
— the load-bearing physics gap is the axial heat-flux distribution**.
This candidate fixes the gap by adding cell-to-cell coupling through
the convective h_c.

## Physical model

Replace the current per-cell local h_c:
```
h_c[i] = h_c_Bartz(rho[i], u[i], T[i], ...)   # local-state Dittus-Boelter
```
with mass-addition-augmented:
```
h_c[i] = h_c_local(state[i]) · f_blow(G_cum[i] / G_ref)

G_cum[i] = G_igniter + sum_{j < i, burning} (rho_p · r_b[j] · P_b[j] · dx[j] / A_p[i])

f_blow(G_ratio) = (1 + G_ratio)^0.8     # Dittus-Boelter Re^0.8 scaling
```

where:
- `G_igniter` is the upstream mass flux from the pyrogen plenum entering
  cell 0
- The sum is over upstream cells `j < i` that have already ignited
- `rho_p` is propellant density, `r_b[j]` is the burning cell's
  regression rate, `P_b[j]` is its perimeter, `A_p[i]` is the port
  area at cell i
- `G_ref` is a reference mass flux (use cell 0's mass flux from the
  igniter at t=0, ~1 kg/(m²·s) for typical pyrogens)

**At t=0** before any cell has ignited, only the pyrogen contributes
to G_cum and downstream cells see weaker h_c — they don't ignite
immediately.

**As cells light sequentially**, G_cum grows along the bore and
augments the next unignited cell's h_c, producing a propagating
ignition front.

**At fully developed flow**, the augmentation factor saturates and
the model returns to behavior similar to the current build for the
plateau phase.

No new fitted constants are introduced if `f_blow` is taken as
Dittus-Boelter (Re^0.8 scaling, established correlation). The
existing Bartz local-h_c kernel is preserved unmodified; the
augmentation is a multiplicative factor.

## Implementation interface

Affects:
- `srm_1d/simulation.py:_run_time_loop` — add cumulative-sum pass
  per timestep before the Goodman ODE step (~10 LOC, Numba-friendly)
- `srm_1d/solid_thermal.py` — `_step_goodman_ode` already takes h_c
  as an argument, no signature change needed
- `srm_1d/burn_rate.py` — careful: existing Ma 2020 erosive burn rate
  also uses local Re from G. **Need to NOT double-count** by using
  the local-Re-only path for the erosive enhancement and the
  G_cum-augmented path only for the Goodman pre-ignition heating
- `srm_1d/propellant.py:Propellant` — optionally add
  `flame_spread_enabled: bool = True` for diagnostics

Proposed kernel addition (Numba):

```python
@njit(cache=True)
def _compute_cumulative_mass_flux(G_igniter, rho_p, r_b, P_b, A_p,
                                  is_burning, dx, N):
    """Compute G_cum[i] = upstream mass-flux integral up to cell i.

    G_cum[i] = G_igniter + sum_{j < i, burning} rho_p * r_b[j] * P_b[j] * dx[j] / A_p[i]
    """
    G_cum = np.zeros(N)
    running_mdot = G_igniter  # mass flux entering cell 0 from pyrogen
    for i in range(N):
        G_cum[i] = running_mdot / max(A_p[i], 1e-12)
        if is_burning[i]:
            # This cell contributes to downstream mass flux
            running_mdot += rho_p * r_b[i] * P_b[i] * dx[i]
    return G_cum

@njit(cache=True)
def _blowing_augmentation(G_ratio):
    """Dittus-Boelter Re^0.8 scaling on cumulative mass flux."""
    return (1.0 + G_ratio) ** 0.8
```

Per-step call site (in `_run_time_loop` before Goodman):

```python
G_cum = _compute_cumulative_mass_flux(
    G_igniter_this_step, rho_propellant, r_b, P_burn, A_port,
    is_burning, dx, N
)
G_ratio = G_cum / G_ref  # G_ref is a constant from initial pyrogen state
for i in range(N):
    if not is_burning[i]:
        h_c_aug = h_c_local[i] * _blowing_augmentation(G_ratio[i])
        # Goodman step uses h_c_aug instead of h_c_local
        delta[i], T_surf[i] = _step_goodman_ode(
            delta[i], T_surf[i], h_c_aug, T_gas[i], T_initial,
            solid_alpha, k_solid, dt,
        )
```

`G_ref` is a one-time constant captured at simulation init from the
pyrogen plenum's expected mass flux at first orifice opening (a
known function of pyrogen mass + throat area + chamber pressure).

## Validation strategy

1. **Unit tests**: kernel-level checks that
   - `_compute_cumulative_mass_flux` returns G_igniter / A_port[0]
     when no cells are burning
   - G_cum[i] increases monotonically with i when upstream cells are
     burning
   - `_blowing_augmentation(0) == 1.0` and increases with G_ratio
2. **Hasegawa A**: re-run `hasegawa_motor_a.py`; target: ignition
   spike drops from 8.5 → ~6 MPa (matching experimental), ignition
   propagation visible in the snapshot history (`is_burning` array
   ignites cells sequentially over 50-200 ms rather than
   simultaneously).
3. **Cross-motor**: re-run `cross_motor_frozen_vs_effective.py`;
   target: all 4 motors show spike-to-plateau ratio < 1.5, and the
   over-prediction shrinks meaningfully for Zerox/BALLSstick/Chunc.
4. **No regression of plateau / erosive peak**: the v0.7.1
   effective-transport plateau match for Hasegawa A is good; this
   change should not regress it (the augmentation should saturate
   by the time the bore fully ignites).

## Risks / open questions

- **Double-counting risk with Ma 2020 erosive burn rate**: the
  erosive term ALSO uses Re_local from local G. The Kashiwagi
  augmentation applies to pre-ignition heating only — once a cell
  ignites, its burn rate goes through the normal Ma chain. Need to
  verify in code that augmentation is gated on `is_burning[i] == 0`.
- **Numerical sensitivity to G_ref**: literature says G_ref ~ Dittus-
  Boelter cancels but in practice the augmentation factor's
  magnitude depends on the choice. Need a single-cell smoke test
  with two G_ref values (e.g., 0.5x and 2x nominal) to confirm
  spike magnitude is robust within ±10%.
- **Burning-cell mass-flux double-count vs PISO**: PISO already
  advects gas mass through the bore. The cumulative-G calculation
  should NOT add to the PISO momentum balance — it's purely a
  heat-transfer-coefficient driver. Need explicit comment that
  G_cum is a diagnostic input to Goodman h_c only, not a
  conservation-law source.
- **Recirculation regions** (sudden expansions, fin slots): the
  cumulative-G argument breaks when the local flow recirculates
  (Mukunda-Paul JPP 2007). For srm_1d's BATES + FMM motor library
  this should be a minor concern since axial flow is well-defined,
  but BALLSstick or any motor with strong port-area changes may
  show non-physical behavior at the discontinuity.

## Literature

See [../references/02_spatial_ignition_front_coupling.md](../references/02_spatial_ignition_front_coupling.md)
(extended lit dive — populated by subagent 2026-05-23). Key
references:
- Salita 2001 (AIAA 2001-3443) — Modern SRM ignition transient
  modeling architecture (Peretz Q1D gas + Goodman + axial
  impingement-region).
- Han & Cai 2017 (JPP, DOI 10.2514/1.B36024) — cleanest 1D-code
  demonstration of cell-coupled flame spread.
- Kashiwagi 1982 (CST 28) — closed-form flame-spread velocity from
  cumulative-G coupling.
- Cavallini SPINBALL — Sapienza impingement-region implementation.

## Estimated implementation cost

Medium of the 4 candidates. ~100-200 LOC in `simulation.py` +
`solid_thermal.py` integration + ~5 new tests. 1-2 focused sessions
including validation runs. Smaller than the structural change to
`igniter_plenum.py` candidates 3/4 would require.

## Related candidates

- This candidate is the structural ignition-stage fix. **Z-N
  (candidate 1) is complementary** — Z-N lags the burn rate
  ramp-up after ignition; coupling lags the ignition timing
  itself. Stacking both is the published architecture
  (SPP/SPINBALL).
- **Candidate 3 (pyrogen spatial distribution)** addresses the same
  symptom from a different physical lever: distribute the pyrogen
  energy axially rather than concentrating in cell 0. If the
  pyrogen footprint is distributed over the first several cells,
  the cumulative-G coupling effect is partially redundant for the
  head-end region but still important for sustained flame spread
  along the bore. Stacking both is also published (SPINBALL
  impingement-region + downstream flame spread).
