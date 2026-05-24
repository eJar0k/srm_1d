# Candidate 1: Z-N Dynamic Burn Rate

**One-line**: Replace the instantaneous Saint-Robert burn rate
`r(P) = a·P^n` with a relaxation ODE that lags the steady-state
response, capturing the propellant's transient pyrolysis dynamics.

## Motivation

Solid propellant pyrolysis is not instantaneous — there's a thermal
lag in the burning surface boundary layer (~ms timescale) between a
pressure change and the burn rate reaching its new steady-state
value. During sharp pressure transients (ignition, throat erosion
recovery, end-of-burn knee), the lag matters: instantaneous-rate
models over-predict the magnitude of pressure spikes because they
let the burn rate jump up immediately with no resistance.

Zeldovich-Novozhilov (Z-N) theory gives the canonical form for this
lag as a relaxation ODE on the burn rate. Under the simplest
formulation the relaxation time is `τ_ZN = κ · α / r²` where α is the
propellant thermal diffusivity and κ is an order-1 constant — meaning
**no fitted parameters are needed** if α and r are computed from
existing motor data.

For the v0.7.2 cross-motor spike artifact, Z-N would smooth the
simultaneous-cell ignition transient by preventing all cells from
jumping to their pressure-dependent steady-state burn rate at once;
each cell's burn rate ramps over τ_ZN ≈ 1-10 ms, blunting the
pressure spike's leading edge.

## Physical model

For each bore cell:

```
dr/dt = (r_steady(P, ...) - r) / τ_ZN
τ_ZN = κ · α / r²
α = k_solid / (ρ_propellant · Cps)
```

where `r_steady` is the current Ma 2020 erosive-augmented burn rate
(everything `srm_1d/burn_rate.py` currently computes), `r` is the
new dynamic burn rate that actually drives mass production / regression,
and `κ` is the Z-N coupling constant.

**At steady state**, `dr/dt → 0` so `r → r_steady` and Z-N is
invisible — the model reduces to current behavior.

**During transients**, `r` lags `r_steady` by τ_ZN. For α ~ 0.15
mm²/s (AP/HTPB) and r ~ 5 mm/s, τ_ZN ~ 0.15 / 25 ≈ 6 ms with κ=1
— matches the experimental ignition-transient timescale.

The ODE is **stiff when r→0** (τ_ZN→∞ implies r locks at zero). Need
to clamp: `τ_ZN = min(τ_max, κ·α/max(r, r_floor)²)` with sensible
floors (e.g., r_floor = 0.1 mm/s, τ_max = 100 ms).

## Implementation interface

Affects:
- `srm_1d/burn_rate.py` — `r_steady` becomes the post-erosive output;
  new dynamic-rate state lives elsewhere
- `srm_1d/simulation.py:_run_time_loop` — per-cell `r_dynamic[N]`
  array advanced each step; replaces direct use of erosive burn rate
- `srm_1d/propellant.py:Propellant` — add `kappa_zn: float = 1.0`
  attribute (one knob, with literature prior)

Proposed API addition (Numba kernel):

```python
@njit(cache=True)
def _advance_zn_burn_rate(r_dyn, r_steady, alpha_solid, kappa_zn,
                          r_floor, tau_max, dt, N):
    """Advance per-cell Z-N dynamic burn rate by one timestep.

    Explicit Euler is fine because τ_ZN >> dt (gas-side dt is
    typically <10 μs; τ_ZN is 1-10 ms).
    """
    for i in range(N):
        r_eff = max(r_dyn[i], r_floor)
        tau = min(tau_max, kappa_zn * alpha_solid / (r_eff * r_eff))
        r_dyn[i] += dt * (r_steady[i] - r_dyn[i]) / tau
    return r_dyn
```

Mass conservation: regression rate driver becomes `r_dyn` instead of
`r_steady`. Mass-source-by-species writes use `r_dyn`.

## Validation strategy

1. **Unit tests**: kernel-level checks that steady-state limit
   (`r_steady` constant) recovers `r → r_steady`, that step changes
   in `r_steady` give first-order exponential response with measured
   time constant τ_ZN.
2. **Hasegawa A**: re-run `hasegawa_motor_a.py` after Z-N lands;
   target: ignition spike P_peak drops from 8.5 MPa toward
   experimental 6.5 MPa without regressing plateau / erosive peak /
   tail-off shape.
3. **Cross-motor**: re-run `cross_motor_frozen_vs_effective.py`;
   target: all 4 motors show spike-to-plateau ratio < 1.5 (vs
   current 1.3-10 across the matrix).
4. **No regression of Hasegawa A LHS path B rank-1**: re-score the
   v0.7.1 effective rank-1 calibration against the same MSE; should
   improve (lower MSE) or hold within 5%.

## Risks / open questions

- **κ value**: literature claims κ = 1 with no fit, but real
  propellants may have material-specific values. Need a literature
  range from the lit dive subagent.
- **r_floor / tau_max sensitivity**: numerical robustness floors may
  affect ignition timing. Need to verify floors are well below the
  range where the physics matters.
- **Coupling to erosive burn rate**: does Z-N lag both the base and
  erosive contributions, or only the base? Cavallini 2009 may have
  a specific recommendation.
- **Stiffness near cell ignition** (r=0 immediately before, r>>0
  immediately after): the explicit-Euler step may need a special
  case at the ignition transition.

## Literature

See [../references/01_z_n_dynamic_burn_rate.md](../references/01_z_n_dynamic_burn_rate.md)
(extended lit dive — populated by subagent 2026-05-23).

## Estimated implementation cost

Smallest of the 4 candidates. ~150-300 LOC in
`simulation.py` + `burn_rate.py` + `propellant.py` + ~5-8 new tests.
1-2 focused sessions including validation runs.

## Related candidates

- Z-N is **fundamentally a burn-rate physics fix**; candidates 2/3/4
  are ignition-stage fixes. Z-N can stack with any of them.
- Predicted synergy with (3) pyrogen spatial distribution: spreading
  the early pyrogen energy axially AND giving each cell's burn rate
  a lag would both attack the simultaneous-ignition artifact, and
  they don't conflict.
