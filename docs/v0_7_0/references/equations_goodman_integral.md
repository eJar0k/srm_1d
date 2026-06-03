# Goodman Cubic-Polynomial Integral Method — Derivation

The load-bearing kernel for v0.7.0's per-cell solid-phase conduction.
Per Peretz 1973 §III.G (Eqs. III-44 to III-50). One ODE per cell —
Numba-friendly, ~5% error vs. exact PDE.

This document is the full derivation for any implementing agent. The
end result (Eqs. 7 and 8 in [DESIGN.md](../DESIGN.md)) drops out
algebraically.

## Problem statement

For each unignited propellant cell, solve the 1D heat conduction
equation in the propellant slab normal to the burning surface:

```
∂T/∂t = α · ∂²T/∂y²        (PDE)
α = k / (ρ · c)
```

with:
- `y = 0` at the burning surface, `y > 0` into the propellant
- IC: `T(y, t=0) = T_initial` (propellant bulk temperature, e.g. 290 K)
- BC at `y = 0`: convective heat flux from gas
  `-k · ∂T/∂y|_{y=0} = h_c · (T_gas - T_surf)`
- BC at large `y`: `T → T_initial`

The naive approach is to discretize the PDE per cell with N_y nodes and
solve implicitly. Cost: per-cell N_y matrix solve every gas-side step.
For a typical 100-cell motor this is 100 × O(N_y) extra work per
timestep. Numba can do this but it's painful.

Goodman's heat-balance integral method (1958, ASME Trans. 80) reduces
the per-cell PDE to **one ODE for a "penetration depth" `δ(t)`**, with
the surface temperature `T_surf(t)` recoverable algebraically. No inner
spatial grid per cell — N_y = 0 effectively.

## Step 1 — Cubic polynomial profile

Assume the temperature profile inside the thermal boundary layer
`y ∈ [0, δ(t)]` is a cubic polynomial:

```
T(y, t) = T_initial + (T_surf(t) - T_initial) · (1 - y/δ(t))^3      (P)
```

Outside the boundary layer (`y > δ`), `T = T_initial` and `∂T/∂y = 0`.

Verify boundary conditions on (P):

| BC | Check |
|---|---|
| `T(0, t) = T_surf` | `T_initial + (T_surf - T_initial)·1^3 = T_surf` ✓ |
| `T(δ, t) = T_initial` | `T_initial + 0 = T_initial` ✓ |
| `∂T/∂y(δ, t) = 0` | `-3(T_surf - T_initial)/δ · (1 - δ/δ)^2 = 0` ✓ |
| `∂²T/∂y²(δ, t) = 0` | `6(T_surf - T_initial)/δ² · (1 - δ/δ) = 0` ✓ |

Four BCs satisfied by the four polynomial coefficients (cubic = 4 free
parameters); the cubic is uniquely determined.

## Step 2 — Convective surface BC closes T_surf(δ)

At `y = 0`, the convective BC is:

```
-k · ∂T/∂y|_{y=0} = h_c · (T_gas - T_surf)
```

From (P): `∂T/∂y(y) = -3 (T_surf - T_initial) (1 - y/δ)² / δ`, so at
`y = 0`:

```
∂T/∂y|_{y=0} = -3 (T_surf - T_initial) / δ
```

Substituting:

```
-k · [-3 (T_surf - T_initial) / δ] = h_c · (T_gas - T_surf)
3k (T_surf - T_initial) / δ = h_c (T_gas - T_surf)
```

Solving for `T_surf`:

```
T_surf = (3k · T_initial + h_c · δ · T_gas) / (3k + h_c · δ)        (8)
```

This is **DESIGN.md Eq. 8** — `T_surf` is an algebraic function of `δ`
given the gas-side state `(T_gas, h_c)` and the propellant properties
`(k, T_initial)`.

Sanity checks:
- `δ → 0` (initial condition): `T_surf → T_initial` ✓ (no thermal layer
  yet, surface is at bulk temp)
- `δ → ∞`: `T_surf → T_gas` ✓ (thermal layer fully developed, surface
  in equilibrium with gas)

## Step 3 — Heat-balance integral

Integrate the PDE from `y = 0` to `y = δ(t)`:

```
∫₀^δ ∂T/∂t dy = α · ∫₀^δ ∂²T/∂y² dy
              = α · [∂T/∂y]₀^δ
              = α · [0 - ∂T/∂y|_{y=0}]
              = α · 3(T_surf - T_initial) / δ
```

The LHS uses Leibniz's rule (since `δ` depends on `t`):

```
d/dt ∫₀^δ T dy = ∫₀^δ ∂T/∂t dy + T(δ) · dδ/dt
               = ∫₀^δ ∂T/∂t dy + T_initial · dδ/dt
```

So:

```
∫₀^δ ∂T/∂t dy = d/dt ∫₀^δ T dy - T_initial · dδ/dt
```

Compute `∫₀^δ T dy` using (P) — substitute `u = 1 - y/δ`,
`dy = -δ du`, `y = 0 → u = 1`, `y = δ → u = 0`:

```
∫₀^δ T dy = ∫₀^δ [T_initial + (T_surf - T_initial)(1 - y/δ)³] dy
          = T_initial · δ + (T_surf - T_initial) · ∫₀^δ (1 - y/δ)³ dy
          = T_initial · δ + (T_surf - T_initial) · δ · ∫₀¹ u³ du
          = T_initial · δ + (T_surf - T_initial) · δ / 4
          = T_initial · δ + δ(T_surf - T_initial)/4
```

Differentiate wrt time:

```
d/dt ∫₀^δ T dy = T_initial · dδ/dt
                + (1/4) · d/dt[δ(T_surf - T_initial)]
              = T_initial · dδ/dt
                + (1/4) · [dδ/dt · (T_surf - T_initial)
                          + δ · dT_surf/dt]
```

So:

```
∫₀^δ ∂T/∂t dy = (1/4) · [(T_surf - T_initial) · dδ/dt + δ · dT_surf/dt]
```

Setting LHS = RHS:

```
(1/4) · [(T_surf - T_initial) · dδ/dt + δ · dT_surf/dt]
    = α · 3(T_surf - T_initial) / δ
```

Multiply by 4:

```
(T_surf - T_initial) · dδ/dt + δ · dT_surf/dt = 12α(T_surf - T_initial)/δ
```

This is the **heat-balance integral equation** — one equation in two
unknowns `(δ, T_surf)`. Eq. (8) closes the system.

## Step 4 — Substitute (8) to get a single ODE in δ

From (8), differentiate wrt `t` (assuming `T_gas` and `h_c` are
quasi-static — i.e., they vary on the gas-side timescale, which is
slow compared to the per-cell thermal layer development):

```
dT_surf/dδ · dδ/dt
```

Using quotient rule on (8), let:
- Numerator `N = 3k·T_initial + h_c·δ·T_gas`, so `dN/dδ = h_c·T_gas`
- Denominator `D = 3k + h_c·δ`, so `dD/dδ = h_c`

```
dT_surf/dδ = (dN/dδ · D - N · dD/dδ) / D²
           = (h_c·T_gas · (3k + h_c·δ) - (3k·T_initial + h_c·δ·T_gas) · h_c) / (3k + h_c·δ)²
           = h_c · [T_gas·(3k + h_c·δ) - 3k·T_initial - h_c·δ·T_gas] / (3k + h_c·δ)²
           = h_c · [3k·T_gas - 3k·T_initial] / (3k + h_c·δ)²
           = 3k·h_c·(T_gas - T_initial) / (3k + h_c·δ)²
```

Also from (8):

```
T_surf - T_initial = [(3k·T_initial + h_c·δ·T_gas) - T_initial(3k + h_c·δ)] / (3k + h_c·δ)
                   = h_c·δ·(T_gas - T_initial) / (3k + h_c·δ)
```

Substitute into the heat-balance integral equation:

```
[h_c·δ·(T_gas - T_initial) / (3k + h_c·δ)] · dδ/dt
   + δ · [3k·h_c·(T_gas - T_initial) / (3k + h_c·δ)²] · dδ/dt
= 12α · [h_c·δ·(T_gas - T_initial) / (3k + h_c·δ)] / δ
```

Factor `h_c · (T_gas - T_initial)` out of every term and multiply both
sides by `(3k + h_c·δ)²`:

```
[h_c·(T_gas - T_initial)] · dδ/dt · [δ·(3k + h_c·δ) + 3k·δ]
   = 12α · h_c·(T_gas - T_initial) · (3k + h_c·δ)
```

Cancel `h_c·(T_gas - T_initial)`:

```
dδ/dt · [δ·(3k + h_c·δ) + 3k·δ] = 12α · (3k + h_c·δ)
dδ/dt · δ · [3k + h_c·δ + 3k] = 12α · (3k + h_c·δ)
dδ/dt · δ · (6k + h_c·δ) = 12α · (3k + h_c·δ)
```

Solve for `dδ/dt`:

```
dδ/dt = 12α · (3k + h_c·δ) / [δ · (6k + h_c·δ)]                     (7)
```

This is **DESIGN.md Eq. 7** — a single ODE for the penetration depth
`δ(t)`, given the gas-side state `(T_gas, h_c)` (which is treated as
quasi-static within each gas timestep).

## Step 5 — Verify against analytical limits

### Limit 1: `h_c · δ ≪ 3k` (early times, very thin thermal layer)

```
dδ/dt ≈ 12α · 3k / (δ · 6k) = 6α / δ
```

So `d(δ²)/dt = 12α`, giving `δ(t) = √(12αt)`. From Eq. (8) in this
limit:

```
T_surf ≈ T_initial + (h_c·δ/3k) · (T_gas - T_initial)
       = T_initial + (h_c · √(12αt) / 3k) · (T_gas - T_initial)
```

For comparison, the **exact** semi-infinite-slab solution with constant
convective BC `q = h_c(T_gas - T_surf)`:

```
T_surf(t) - T_initial = (T_gas - T_initial) · [1 - exp(β² t) · erfc(β·√t)]
β = h_c/(ρ·c·√α) = h_c/√(k·ρ·c) = h_c/(k/√α)
```

For small `β·√t` (early times): `T_surf - T_initial ≈ (T_gas - T_init) · 2β√(t/π)`

Goodman: `T_surf - T_initial ≈ (h_c · √(12αt)/(3k)) · (T_gas - T_init)`
       `= (h_c/√k · √(12αt/k))/3 · (T_gas - T_init)`
       `= 2 h_c √(αt/3) / k · (T_gas - T_init)`
       `= 2 β √(αt/3) · (T_gas - T_init) · √(α)/√(α)`
       `= 2β√(t/3)·√α · (T_gas - T_init)`

Wait — let me redo: `β = h_c/√(kρc)`. So `β √t = h_c √(t/(kρc))`. And
`T_surf-T_init` (exact, early) `= (T_gas-T_init) · 2β√(t/π)`.

Goodman gives `T_surf - T_init = 2 h_c √(αt/3) / k · (T_gas-T_init)`
where `α = k/(ρc)`. So `h_c √(αt/3)/k = h_c √(t/(3kρc)) = β √(t/3)`.

Ratio Goodman/exact = √(π/3) ≈ 1.023. **Goodman over-predicts the
surface temperature rise by ~2.3% in the early-time limit.**

### Limit 2: `h_c · δ ≫ 3k` (late times, thermal layer fully developed)

```
dδ/dt ≈ 12α / δ → δ² = 24αt → δ = √(24αt)
```

In this limit `T_surf → T_gas` per (8). The exact solution also
asymptotes to `T_gas` exponentially, so for ignition criterion purposes
both give the same answer.

### Worst-case error band

Per Peretz 1973, the cubic-polynomial Goodman method has been verified
against the exact PDE: **2-5% error on time-to-`T_ign`** at constant
heat flux. For the ignition transient where heat flux varies rapidly
(igniter shock arrival, then plateau, then erosive enhancement), the
error is similar, dominated by the assumption of self-similar cubic
profile at every instant.

## Step 6 — Numerical implementation notes

### Initial condition

`δ(0) = 1e-6 m` (1 µm). The ODE has a `1/δ` term that's stiff at
`δ = 0`. Starting at 1 µm: `dδ/dt ≈ 12α·3k/(1e-6 · 6k) = 6e6 · α`.
For typical α = 1e-7 m²/s: `dδ/dt ≈ 0.6 m/s` initially.

### Adaptive substepping during induction

For first ~10 µs after ignition, `δ` grows fast and `dδ/dt` is large.
The gas-side timestep (~1 µs) may be too coarse. Implementation
strategy: subdivide the gas-side timestep into N_sub Goodman substeps
when `dδ/dt · dt_gas / δ > 0.1` (CFL-like criterion on δ).

### When to stop integrating

Once `T_surf > T_ignition`, set `is_burning[i] = True` and freeze
`δ[i]`. The thermal layer is now combusting; the cell switches to the
existing Ma 2020 burn-rate kernel. Don't continue to integrate
`δ` after ignition — it's meaningless and will diverge.

### Numerical method

RK4 in the same time loop. Numba @njit-compatible. Per-cell scalar
state — no inner spatial loop. Roughly:

```python
@njit(cache=True)
def _step_goodman_rk4(delta, T_gas, h_c, alpha, k_solid, dt):
    """Integrate Eq. (7) with RK4. Returns new delta."""
    def rhs(d):
        return 12 * alpha * (3 * k_solid + h_c * d) / (d * (6 * k_solid + h_c * d))

    k1 = rhs(delta)
    k2 = rhs(delta + 0.5 * dt * k1)
    k3 = rhs(delta + 0.5 * dt * k2)
    k4 = rhs(delta + dt * k3)
    return delta + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
```

Compile-time properties:
- All scalar arithmetic — no allocations
- Branch-free
- Numba-compatible

### Stiffness watch

The `1/δ` term makes the ODE stiff at small δ. RK4 is conditionally
stable; the stability boundary scales like `dt < δ²/α`. For δ = 1 µm,
α = 1e-7: `dt_stable < 1e-12 / 1e-7 = 1e-5 s = 10 µs`. Gas timestep is
typically ~1 µs, so RK4 is fine. If a stiffer regime ever shows up,
swap to backward Euler.

## References

- Goodman, T. R. (1958). "The heat-balance integral and its application
  to problems involving a change of phase." Trans. ASME 80: 335-342.
  Cited as Ref. 49 in Peretz 1973.
- Peretz, A., Caveny, L. H., Kuo, K. K., Summerfield, M. (1973).
  "The Starting Transient of Solid-Propellant Rocket Motors with High
  Internal Gas Velocities." Princeton AMS Report 1100. §III.G,
  Eqs. III-44 to III-50, pp. 43-47.
  PDF: `19740005393.pdf` in repo root.
- Carslaw, H. S. & Jaeger, J. C. (1959). "Conduction of Heat in Solids,"
  2nd ed. Oxford. — for the exact analytical comparison solution.
