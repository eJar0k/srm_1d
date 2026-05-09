# Ma et al. 2019 — Igniter MFR History Model: Full Extraction

**Citation**: Ma, Y., Bao, F., Hui, W., Liu, Y., Gao, Y. "A Model for Igniter Mass Flow Rate History Evaluation for Solid Rocket Motors." *International Journal of Aerospace Engineering*, vol. 2019, Article ID 2593602, 12 pages. https://doi.org/10.1155/2019/2593602

**PDF in repo**: `International Journal of Aerospace Engineering - 2020 - Ma - A Model for Igniter Mass Flow Rate History Evaluation for.pdf`

**One-line summary**: This paper does **not** model the igniter chamber as a separate plenum producing mdot(t). It is an **inverse-problem** model: take an experimentally measured *main combustion chamber* pressure trace P(t) from an igniter-only firing (no main propellant — main grain replaced by inert aluminium dummy), and back-solve a 0D mass/energy/state system for the igniter mass-flow history dm_ig/dt(t). This distinction is central and changes how the model would map to srm_1d.

---

## 1. The 0D Combustion-Chamber Model — Equations

The control volume is the **main combustion chamber** (head end → throat). The igniter is just a black-box source delivering gas of known composition x_i and known flame temperature T_f at rate dm_ig/dt (the unknown).

### 1.1 Outflow boundary (Eq. 1, 2)

```
q_m = 0                    before slug rupture
q_m = P · A_t / c*         after slug rupture

c*  = √(R · T) / Γ
Γ   = √( k · (2/(k+1))^((k+1)/(k-1)) )
```

A_t is treated as **constant** (no throat erosion — short burn). Choking is **assumed** the moment the nozzle plug ruptures (assumption b); the paper acknowledges in §5.1 that this overestimates outflow late in tail-off when chamber P is too low to actually choke.

### 1.2 Mass conservation (Eq. 3, 4)

Total chamber mass:
```
dm/dt = dm_ig/dt - q_m
```

Per-species i (igniter gas has prescribed mass fractions x_i; outflow carries chamber-average composition m_i/m):
```
dm_i/dt = x_i · dm_ig/dt - q_m · m_i/m
```

### 1.3 Energy conservation (Eq. 5–7)

Starting from d(m C_v T)/dt:
```
d/dt(m·C_v·T) = C_p,f(T_f)·T_f · dm_ig/dt
              - q_m · C_p(T) · T
              - q_h
```

The three RHS terms are: enthalpy inflow at the **igniter flame temperature T_f**, enthalpy outflow at chamber T, and heat loss q_h to walls.

After expanding the LHS and substituting the species mass equation, the **energy ODE** becomes (Eq. 7):

```
dT/dt = [ (C_p,f(T_f)·T_f - C_v,f(T)·T) · dm_ig/dt
          - q_m · R · T
          - q_h ] / B

B ≡ T · Σ m_i · dC_v,i/dT + m · C_v(T)
```

Note: the printed Eq. 7 in the PDF has a slightly garbled term order from Wiley's typesetting; the form above is the consistent reconstruction from Eq. 5 + Eq. 6.

### 1.4 Equation of state (Eq. 20–22)

Ideal gas, differentiated in time:
```
V/R* · dP/dt = T · dn/dt + n · dT/dt
n = Σ m_i / M_i
```

Combined with Eqs. 4 and 7, the **closure for the unknown dm_ig/dt** (Eq. 22 — the workhorse of the paper):

```
dm_ig/dt = [ V·(dP/(R*·dt)) + q_m·(T·n·R + 1/M*) + (n/B)·q_h ]
           / [ T/M_f + (n/B)·(C_p,f(T_f)·T_f - C_v,f(T)·T) ]
```

(M* and M_f are mixture and igniter-gas molar masses respectively.) The numerator's dP/dt is taken by **numerical differentiation of the experimental pressure trace** — that is the input that makes the inverse problem well-posed.

### 1.5 Multi-species temperature-dependent C_p (Eq. 23)

Each pure species uses the 7-coefficient NASA-CEA polynomial:
```
C_p,i(T)/R* = a_1·T^(-2) + a_2·T^(-1) + a_3 + a_4·T + a_5·T² + a_6·T³ + a_7·T⁴
```

Coefficients pulled directly from the CEA database. Five species considered (Table 1): CO₂ (56.3%), CO (21.9%), N₂ (19.0%), SO₂ (2.8%), O₂ (0% — placeholder). Black powder igniter (75/15/10 KNO₃/S/C); CEA gives 43.7% condensed-phase mass which is dropped (assumption c).

### 1.6 Igniter burn-rate model

**There is none.** The paper deliberately avoids modelling the igniter chamber. Reason given in the introduction: small SRM igniters often use powders/granules/pellets/strips with no defined burning surface, and confined-igniter cases rupture during firing, so neither a Saint-Robert exponent nor a burning area is meaningful. dm_ig/dt is solved *for*, not solved *from* a burn law.

---

## 2. Coupling to the Main Motor

This is where the paper's framework diverges from what srm_1d's v0.7.0 needs. Two distinct couplings appear:

### 2.1 Inside the inverse-problem run

The "main motor" and the "igniter chamber" are **the same control volume**. The igniter's chemical output is delivered as a source term `x_i · (dm_ig/dt)` carrying enthalpy `C_p,f(T_f) · T_f`. Because the model is 0D, **spatial distribution does not exist** — head-end vs. anywhere-else is meaningless inside the model proper.

The chamber's rising P does *not* feed back to the igniter chamber's outflow, because the igniter chamber is not modelled. The plug-rupture transition (Eq. 1) is the only "back-pressure" effect: while P < P_rupture, q_m = 0; afterwards, choked flow at the main throat.

### 2.2 In the validation run (Fluent 2D axisymmetric, §5.3)

For external CFD validation only, the recovered dm_ig/dt(t) is fed into Fluent 6.3 as a **mass-flow-inlet boundary condition** at the **igniter outlet location** — the head end with 9 orifices (1 axial, 8 radial). Cubic-spline-interpolated in time. Single boundary patch, not multi-cell injection.

The Fluent run (k-ε, species transport with no reaction, T-dependent properties, 121,195 cells, dt = 1 µs, residuals < 10⁻³) recovered peak P = 3.06 MPa vs. experimental 3.16 MPa (3.13% error) — that is the validation.

---

## 3. Heat Transfer (q_h)

q_h is computed **only as a sink** in the chamber energy equation; it doesn't drive the igniter chamber because there isn't one.

### 3.1 Surface partitioning ("segment way")

The gas-solid interface is split into three segments (Eq. 8 applied to each, summed):
- gas–grain (aluminium-alloy dummy)
- gas–case (steel)
- gas–nozzle (steel + steel throat lining)

For each segment:
```
q_h = h · A · (T - T_w)
h   = h_c + h_r
```

### 3.2 Convective coefficient (Eq. 9)

Dittus-Boelter form:
```
h_c = 0.023 · Pr^(-2/3) · C_p · (ρu)^0.8 / D^0.2
```

with Pr = 4k/(9k − 5), and viscosity μ = 1.187×10⁻⁷ (1000 M)^0.5 T^0.6 (Huzel & Huang correlation). Note: NOT Gnielinski/Haaland — much simpler than srm_1d's Ma 2020 burn-rate stack uses for the *propellant* surface.

### 3.3 Velocity for Re (Eq. 10–13)

Although the chamber is 0D, h_c needs a velocity. **Two estimators** computed in parallel and the larger is used for grain and case:

1. **Area-ratio estimator**: assume mass flux at section A_p equals throat mass flux. Solve A_t/A_p = f(Ma, k) by Newton iteration; recover u from u = Ma·√(kRT_static).
2. **Port-throughput estimator**: assume all igniter gas passes through the port: dm_ig/dt = ρ_static · u · A_p.

Justification: estimator 1 overestimates u; estimator 2 underestimates u. Taking max(·) is an empirical compromise.

### 3.4 Radiative coefficient (Eq. 14)

```
h_r = C · σ · (T² + T_w²)·(T + T_w),   C = 0.25
```

C = 0.25 cited from prior literature; Ma's case (4) shows h_r is **small relative to h_c**. Lumps gas-emission and condensed-particle emission into one constant.

### 3.5 Wall temperature T_w — 1D heat conduction (§2.3, Eq. 15–19)

For each segment:
```
∂T_s/∂t = α · ∂²T_s/∂y²,   α = λ/(ρ_s · C_p,s)
```

Boundary conditions:
- Outer wall: Dirichlet T_s(t, L) = T_a
- Aluminium grain: half-thickness modelled with Neumann symmetry
- Inner gas-side wall: Robin, −λ_s · ∂T_s/∂y|_w = h(T − T_w)

Surface temperature recovered from a 2-cell ghost extrapolation:
```
T_w = [h·T + (λ_s/3Δy)·(9·T_0 - T_1)] / [h + 8·λ_s/(3Δy)]
```

40 cells per segment, central-difference space, first-order implicit time.

---

## 4. Numerical Method

- **Gas-phase ODEs (Eq. 7 + Eq. 22 + species, Eq. 4)**: 4th-order Runge-Kutta, **fixed Δt = 0.05 ms**.
- **Solid-phase 1D conduction**: FVM, central-difference space, first-order implicit time, 40 cells per segment.
- **Coupling**: segregated. Each Δt: solve gas RK4 → update h_c, h_r → solve solid conduction → recompute T_w → re-solve gas → iterate **until ΔT_chamber between successive sweeps < 10⁻⁹** (relative). Then advance time.
- No CFL constraint reported; the implicit conduction step is unconditionally stable, and 0.05 ms with peak |dP/dt| ~1.0 MPa/ms is comfortable.

---

## 5. Validation

### 5.1 Apparatus (§4)

Small SRM, real grain replaced by inert aluminium-alloy dummy of same geometry: L = 73.5 mm, r_in = 2.5 mm, r_out = 18 mm. Steel case (r_in = 19.25, r_out = 35 mm). Long-tail nozzle: 31 mm cylindrical inlet (r = 7 mm), throat r = 2.55 mm, 6.65 mm minimum thickness. Igniter: aluminium case, 9 orifices, 1.35 g black powder (75/15/10). Burn duration ~26.5 ms; peak P ~3.16 MPa at ~5 ms. Pressure transducer mid-case, sampled every 0.5 ms.

### 5.2 What was measured

**Only chamber P(t)**. The paper does NOT independently measure mdot_ig(t), T_0_ig(t), or igniter chamber pressure — those are recovered by the model. Validation is by:

1. **Total recovered mass**: 4 cases (different Cp/species assumptions) all within 2.6–4.0% of the experimental gas-phase mass 0.76 g (Table 3). Best is case 3 (species-resolved Cp(T), no radiation) at 2.57% high.
2. **2D Fluent re-run** with recovered dm_ig/dt(t) as inlet BC: peak P = 3.06 MPa vs. 3.16 MPa experimental (3.13% low). Pressure trace shape is in good visual agreement.

Quality summary: total-mass error 2–4%, peak-P error 3% on a single test article.

---

## 6. Required Inputs

For the **chamber** model:
- V (chamber volume, m³)
- A_t (throat area, m²) — assumed constant
- P_rupture (slug rupture pressure, Pa)
- A_p (port cross-section area, m²) for the velocity estimator
- L_grain, geometry of dummy grain
- Solid material properties (ρ_s, C_p,s, λ_s)
- T_a (ambient T)

For the **igniter gas chemistry**:
- Igniter mass fractions {x_i} per species (from CEA — Table 1)
- T_f (igniter flame temperature, K)
- Per-species NASA-CEA 7-coefficient C_p(T) polynomials (a_1…a_7)
- Per-species M_i

For the **measurement**:
- P(t) sampled finely enough that dP/dt can be differentiated cleanly.

---

## 7. Things Explicitly NOT Modelled

Stated by the paper itself:

- **Condensed-phase particles** in the igniter products (assumption c) — 43.7% of black-powder products by mass.
- **Secondary combustion** between igniter CO/H₂ and chamber O₂ (assumption d, "frozen flow"). Identified as the dominant remaining error source (~3% mass overprediction).
- **Subsonic outflow** in tail-off (assumption b) — choking is forced from rupture to end-of-burn.
- **Throat erosion** (A_t = const). Justified by short burn.
- **Igniter chamber internal dynamics** — no plenum, no igniter-charge burn-rate law, no choked-orifice from igniter to main chamber. The igniter is purely a prescribed-composition source term.
- **Igniter ignition delay** itself — t = 0 of the model is the start of the experimental P trace.
- **Plenum filling / ullage of the igniter case**.
- **Kinetic and potential energy** of the gas in the chamber energy balance.
- **Variable A_t** (no nozzle erosion).
- **Gradients in the chamber** (0D — explicitly L/D-limited).

---

## 8. Implementation Cost for srm_1d

### 8.1 Critical orientation

The model as published is **inverse**: P(t) experimental → mdot_ig(t). srm_1d wants the **forward** direction: an *a priori* mdot_ig(t) producing P(t). The paper cannot be used directly without either:

- **Path A — Pretabulated source**: run the Ma model offline against a measured P(t), bake mdot_ig(t) and (T_f, x_i) into a CSV/.npy, and have srm_1d read it as a table. Fast, JIT-compatible (just np.interp inside the @njit loop), but every new igniter needs a hot-fire P-trace. This sidesteps the paper's "we don't model the igniter chamber" stance entirely — you'd still have the same calibration tax srm_1d is trying to escape.

- **Path B — Forward 0D igniter chamber**: do what Ma explicitly chose *not* to do. Build a separate igniter-plenum control volume with its own mass/energy/state, a Saint-Robert burn law and known burning area for the igniter charge, and a choked-orifice coupling to the main chamber. **This is what srm_1d v0.7.0 does** (see DESIGN.md). Ma 2019 supplies useful sub-pieces (the multi-species temperature-dependent C_p with NASA polynomial, ideal-gas closure, choked outflow form) but the igniter plenum itself must come from elsewhere.

### 8.2 Numba-JIT viability

The forward 0D plenum (Path B) is JIT-friendly: it is a small ODE system (1 + N_species states + a choked-orifice algebraic relation). RK4 at Δt = 0.05 ms — same order as Ma's gas-side step — could run as a sub-cycler inside the main @njit time loop. The NASA polynomial C_p(T) is a 7-term series, also @njit-friendly. The 1D wall-conduction sub-solver is what would force a separate, slower path — but heat loss into the plenum walls is a second-order correction that v0.7.0 can defer (Ma's case 3-vs-case 4 comparison shows radiation only changes the answer by 0.2% on total mass).

### 8.3 Pieces of Ma 2019 that **do** transplant cleanly

Independent of which path is chosen, these are publication-ready building blocks:
- NASA-CEA 7-coefficient C_p(T) polynomial (Eq. 23) and per-species C_v = C_p − R/M.
- Mixture rules: C_p,mix = Σ x_i C_p,i, M_mix = (Σ x_i/M_i)⁻¹.
- B-function (the heat-capacity-derivative term in Eq. 6) — only matters if T-dependent C_p is wanted.
- Choked-throat mass-flow form q_m = P A_t / c* (Eq. 1, 2).
- Robin BC + 2-cell extrapolation for T_w (Eq. 18, 19).

---

## 9. Files referenced

- Paper: `Rocketry/Code Stuff/Erosive Burning Solver/srm_1d/International Journal of Aerospace Engineering - 2020 - Ma - A Model for Igniter Mass Flow Rate History Evaluation for.pdf`
- Target replacement site: `srm_1d/srm_1d/simulation.py` lines ~430-450 (igniter sim_kwargs in `run_simulation`)
