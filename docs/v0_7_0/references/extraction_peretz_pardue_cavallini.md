# Foundational References for srm_1d v0.7.0 Igniter Model

Three foundational SRM ignition transient documents extracted: Peretz et al. 1973 (the foundational 1D HVT model), Pardue & Han 1992 (two-phase Shuttle SRM extension), and Cavallini 2009 (SPINBALL/SPIT — Vega-program internal ballistics).

PDFs in repo:
- `19740005393.pdf` — Peretz 1973 (Princeton AMS Report 1100)
- `pardue1992.pdf` — Pardue & Han 1992 (AIAA 92-3277)
- `74323997.pdf` — Cavallini 2009 (Sapienza PhD thesis)

---

## 1. Peretz, Caveny, Kuo, Summerfield (1973) — Princeton AMS-1100

**Full title:** "Starting Transient of Solid-Propellant Rocket Motors with High Internal Gas Velocities" (the AIAA Journal Vol. 11, 1719, 1973 paper is a condensed version; this is the full 250-page report).

**Scope:** "HVT" = high-velocity-transient = low port-to-throat area ratio (Ap/At = 1.06–2.0), where erosive burning during ignition matters. They built a laboratory motor with a head-end pyrogen igniter discharging into a rectangular slab grain (composite PBAA-AP, 80% AP).

### 1.1 Mathematical formulation

Conservation equations are 1D, **single-phase gas only** (no two-phase). After order-of-magnitude analysis, the surviving forms are (their Eqs. III-7 through III-9):

```
∂u/∂t + u·∂u/∂x + (1/ρ)·∂p/∂x = -(f·u²)/(2·d_h) · (P_w/P_b)     (momentum)
∂p/∂t + u·∂p/∂x + γ·p·∂u/∂x = (γ-1)·[q_g·P_w/A_p + ρ_p·r_b·P_b·h_f/A_p]
                              − γ·p·r_b·P_b/A_p                  (energy)
∂T/∂t + u·∂T/∂x = ...                                            (coupled)
```

The **rate of mass accumulation in the free volume created by surface regression** is dropped after order-of-magnitude analysis (gas density is 0.002–0.02 of propellant density). Viscous dissipation, axial heat conduction, and axial viscous stress are dropped. Friction is retained. Thermophysical props (cp, W, γ) **assumed constant** during the entire transient — propellant combustion gas and igniter gas are forced to share γ, cp, W.

The system is **fully hyperbolic**: three real distinct eigenvalues (u, u+c, u-c) — the canonical 1D Euler structure.

### 1.2 Igniter coupling — head-end "entrance section" boundary condition

This is **the single most directly relevant section for v0.7.0**. They do NOT inject mass into a discrete cell — instead they treat the head-end volume upstream of the propellant as a 0D plenum, governed by ODEs (Section III.C.4, Eqs. III-12 through III-15):

**Mass conservation in entrance section:**
```
d/dt(ρ_e·V_e) = mdot_ig(t) − mdot_e→port
```

**Energy conservation in entrance section** (kinetic energy term dropped after OOM analysis):
```
d/dt(ρ_e·V_e·c_p·T_e) = mdot_ig(t)·c_p·T_ig − mdot_e→port·c_p·T_e + heat_loss
```

After substituting the equation of state and rearranging, this becomes a **pair of coupled ODEs** for (p_e, T_e) at the head end (their III-14, III-15). These two ODEs supply two of the three head-end boundary conditions. The third comes from the left-running characteristic compatibility relation Eq. IV-13.

**Critical detail:** mdot_ig(t) is **prescribed from experiment** (a measured trace, not modeled). T_ig is taken **constant during the transient** — and importantly, it's NOT the igniter adiabatic flame temperature; it's an **effective mean temperature** that is ~50–65% of the adiabatic value, accounting for "incomplete combustion (~5%)" + "shock-pattern losses + entrance heat transfer (25–35%)".

The three head-end BCs are stitched to the interior gasdynamics through a single point of coupling: the left-running characteristic at x=0 provides the third compatibility relation. The pair of plenum ODEs + the characteristic form a closed 3×3 algebraic system solved by Cramer's rule each timestep.

### 1.3 Ignition criterion — critical surface temperature

From Section III.G:

> "The ignition criterion adopted in this study is that a point on the propellant surface ignites when it attains some critical ignition (or autoignition) temperature, denoted here by T_ps,ig, at which rapid runaway reactions leading to ignition take place."

Note: this is **not Frazer-Hicks (heat-flux integral) and not Arrhenius rate** — it is literally just T_s(x,t) ≥ T_ps,ig, threshold-of-temperature. They cite that this "leads to results often indistinguishable from those obtained using other criteria."

**Solid-phase heat conduction** is 1D normal to the surface (their Eq. III-36):
```
ρ_p·c_p·∂T_p/∂t = k_p·∂²T_p/∂y²
```

with surface BC III-39:
```
-k_p·∂T_p/∂y|_{y=0} = h_c(t)·[T_gas(t) - T_s(t)]
```

**Numerical method for the solid:** Integral method with **cubic temperature polynomial (Goodman 1958)**, reducing the PDE to an ODE for T_s(t) and the penetration distance δ(t) at each cell. **This is the load-bearing kernel for v0.7.0** — see [equations_goodman_integral.md](equations_goodman_integral.md) for the full derivation.

Solved by 4th-order Runge-Kutta. They verified vs. exact 1D conduction PDE: 5% error on time-to-700K at constant heat flux, 2% for cubic polynomial under constant q. **This is computationally much cheaper than per-cell PDE solves** — a key portability point for srm_1d.

After local ignition: heat-transfer coefficient over propellant surface set to zero (h_cp=0, only inert wall h_cw remains), friction over the burning perimeter set to zero (large blowing).

### 1.4 Multi-species treatment

**Absent.** Igniter gas and propellant gas share γ, cp, W (Assumption 4). They argue this matters most during induction (when igniter gas dominates) and the difference is negligible there.

### 1.5 Numerical method

**Six-point implicit centered finite difference** (their Eq. IV-1) with weighting parameter θ (typically 0.6, in [0.5, 1.0] for unconditional stability). All three governing equations cast as:
```
∂φ/∂t + F(φ)·∂φ/∂x = I(φ)
```
with quasi-linearization on the inhomogeneous I term and **predictor-corrector** on the F coefficients. The resulting linear tridiagonal system is solved by matrix methods.

**Boundary conditions** use compatibility relations along characteristics (right-running at right boundary, left-running at left, plus particle-path at right) — Eqs. IV-9 through IV-14.

**Aft-end nozzle BC:** Quasi-steady isentropic flow between port end and throat. Choking check each step.

### 1.6 Comparison to experiment

Test motor: rectangular slab grain, head-end pyrogen igniter via central sonic nozzle, replaceable rectangular nozzle at aft end. Five pressure stations along port. Ap/At ∈ {1.06, 1.2, 1.5, 2.0}.

Agreement:
- **<10% error** on time-to-first-ignition and time-to-peak head-end pressure for Ap/At ∈ {1.2, 1.5, 2.0}
- **<10% (mostly <5%) error** on peak pressures
- Only one parameter adjusted: erosive burning exponent β (Lenoir-Robillard form).
- For Ap/At = 1.06 (lowest, most erosive): time-axis error up to 15%

### 1.7 Things Peretz is silent on

- **No two-phase / particle treatment.** Aluminum is excluded explicitly.
- **No igniter chamber thermodynamics.** mdot_ig(t) is just a measured curve, T_ig is a fitted constant.
- **No dynamic burning** (Zeldovich-Novozhilov). Discussed and dismissed in Appendix C.
- **No flame-spreading sub-model** beyond per-cell ignition criterion. Flame spread emerges from the heat-conduction + critical-T criterion.

---

## 2. Pardue & Han (1992) — AIAA 92-3277, Two-Fluid Shuttle SRM

**Full title:** "Ignition Transient Analysis of a Solid Rocket Motor Using a One-Dimensional Two-Fluid Model" (Tennessee Tech, NASA/MSFC sponsorship).

**Scope:** The Shuttle SRM (large L/D, low Ap/At, ~17% Al loading) — Peretz's framework cannot represent the Al2O3 particulate phase, so this is a direct extension.

### 2.1 Two-fluid governing equations

Subscript i ∈ {1=gas, 2=Al2O3 condensed phase}. With r_i the volume fraction:

**Mass (Eq. 2):**
```
∂/∂t(r_i·ρ_i·A) + ∂/∂x(r_i·ρ_i·u_i·A) = (w_i·ρ_PR − r_i·ρ_i)·r_b·b + mdot_ig,i
```

**Momentum (Eq. 5):** Each phase has its own momentum equation, with inter-phase drag A·F_ij as coupling.

**Energy (Eq. 11):** Each phase has its own internal-energy equation with inter-phase heat transfer A·Q_ij.

Constraints: Σr_i = 1, Σw_i = 1. Gas alone provides the thermodynamic pressure: P = ρ_1·R_1·T_1·f(ρ_1).

### 2.2 Inter-phase coupling

**Drag (Eq. 9, particle Reynolds form):**
```
F_21 = (18·r_2/d_2²)·μ_1·(u_2 − u_1)·(1 + Re^(2/3)/6)
Re   = ρ_1·|u_2 − u_1|·d_2/μ_1
```

**Heat transfer (Eq. 12):**
```
Q_21 = (6·k_1·r_2/d_2²)·(T_2 − T_1)·(0.58·Re^0.7·Pr^0.3)
```

Particle diameter d_2 = 10 µm assumed throughout.

### 2.3 Burning rate

Identical structure to Peretz/Lenoir-Robillard:
```
r = r_ref·(P/P_ref)^n + α_e·G^0.8·D_h^(-0.2)·exp(−β_e·r·ρ_PR/G)
```
with G = |ρ_1·u_1|.

**Star perimeter ratio α** (Eq. 1, defined 0 ≤ α ≤ 6.7): a fudge factor multiplying actual burning perimeter / cylinder circumference, to model the head-end star segment in 1D. Their Figure 11 shows this α(t) was *empirically tuned* to match experimental data — explicitly admitted as arbitrary.

### 2.4 Igniter and ignition criterion

**Igniter:** Same approach as Peretz — prescribed mdot_ig(t) trace plus prescribed igniter gas adiabatic flame temperature **T_ig = 2450 K** and propellant **T_flame = 3361 K**. The mdot_ig source is split between phases by mass fraction.

**Ignition criterion:** T_surf ≥ 850 K — autoignition temperature. Same form as Peretz (critical T_s).

### 2.5 Numerical method

**SIMPLE-derived IPSA** (Inter-Phase-Slip-Algorithm). Sequence:
1. Solve mass eq. for r_2 (volume fraction of particles); r_1 = 1 - r_2.
2. With assumed pressure field, solve momentum equations for u_i each phase.
3. Pressure-correction equation, update ρ, u, p.
4. Solve energy equations to update T_i.
5. Iterate to convergence within timestep.

Grids tested: 61, 76, 94 control volumes; CPU time on a VAX 8800 was 21–40 minutes for a 1.0-second simulation.

### 2.6 Comparison to experiment

Shuttle SRM data. Particle density swept (1500, 3000 kg/m³), loading swept (10%, 17%). Two-fluid model is **markedly better than one-fluid** on thrust (Figure 9 shows error reduced by ~50%). Head-end pressure gradient is mostly insensitive to particle density.

The takeaway: **two-fluid matters for thrust and chamber filling, but the head-end pressure gradient (the critical IT trace) is dominated by the gas-phase model.** For srm_1d's IT problem this implies **single-phase is acceptable for first-pass** unless thrust prediction matters.

---

## 3. Cavallini (2009) — Sapienza PhD, SPINBALL/SPIT

**Full title:** "Modeling and Numerical Simulation of Solid Rocket Motors Internal Ballistics" — Tutor: Maurizio Di Giacinto. ~330 pp.

### 3.1 Framework: SPINBALL = SPIT-extended

**SPIT** (1990s heritage, by Di Giacinto / Favini / Serraglia at Sapienza) was for ignition transient only, with three-gas formulation (separate species: igniter, pressurizing gas, propellant combustion products). **SPINBALL** is the doctoral-research extension covering the **entire combustion time** (IT + QSS + tail-off) using an "infinite-gases" mixture formulation.

Architecture (Section 1.4.2, Eq. 1.3):
```
∂(ρ_i·A_p)/∂t + ∂(ρ_i·u·A_p)/∂x =
    r_b·P_b·ρ_p·δ_{i,prop} + (mdot_s·A_p/V)·δ_{i,cav} + (mdot_ig·A_p/V)·δ_{i,ig}
                                                                  i=1,...,6

∂(ρ·u·A_p)/∂t + ∂[(ρ·u² + p)·A_p]/∂x − p·∂A_p/∂x =
    (mdot_ig·A_p·v̄_inj/V) + ½·ρ·u²·c_F·P_w

∂(E·A_p)/∂t + ∂[(E+p)·u·A_p]/∂x =
    r_b·P_b·ρ_p·H_f + (mdot_ig·A_p·H_ig/V) + (mdot_s·A_p·H_s/V)
```

Six mass-conservation equations track six gas species independently.

**Sub-models:**
- igniter sub-model
- ignition criterion
- heat-transfer model (convection + radiation)
- conduction in solid propellant
- cavity model (slots, submergence) — 0D ODE bolt-on per cell
- burning rate (APN + erosive)
- grain burnback (GREG = Level Set 3D)
- nozzle seal rupture, nozzle throat ablation

### 3.2 Igniter sub-model

**Cavallini does not give the SPIT igniter equations in detail in this thesis** — he refers to Refs from the Sapienza group. What he describes:

- mdot_ig(t) is **prescribed** from igniter design / experimental measurement
- Distributed over a "control volume defined by the igniter and impingement region sub-models" — i.e., **not just cell zero**: it's a multi-cell injection region driven by the geometry of the igniter jets and where they impinge on the grain
- The **impingement region** is a separate sub-model treating the radial igniter nozzle jets, their impingement angle, and the local enhanced heat transfer + erosive burning where they hit
- Igniter gas has its own enthalpy H_ig which enters the energy source term
- Igniter jets have an "average axial velocity" v̄_inj which enters the momentum source — this is **the modern improvement over Peretz**: igniter jets have momentum, not just mass+enthalpy
- The "head-end and nozzle throat wall BCs" hold initially (igniter gas pressurizes a sealed chamber). Nozzle seal rupture at a prescribed differential pressure releases the gas

### 3.3 Ignition criterion

Cavallini explicitly defers the IT criterion to SPIT references:

> "An Ignition Criterion for the solid propellant based on a temperature pressure dependent of combustion phenomena activation"

**No equations given in this thesis.** This is "temperature-pressure dependent" — implying T_ign(p) rather than constant T_ign as in Peretz/Pardue. Likely Frazer-Hicks-style or Vilyunov-Zarko (the standard European IT criteria).

### 3.4 Burning-rate model — relevant to srm_1d's Ma 2020

**APN (Saint-Robert/De Vieille, Eq. 2.28):**
```
r_b = a·(p/p_ref)^n
a   = a_ref · exp(σ_p·(T_i − T_ref))
```

**Lenoir-Robillard erosive (Eq. 2.31):**
```
r_be = α·(ρu)^0.8 / L^0.2 · exp(−β·r_b/u)
```

**Lawrence modification (Eq. 2.33):** Replace L with D_h (hydraulic diameter).

This is **fundamentally different from Ma 2020** which uses Haaland → Gnielinski friction-factor closure with a roughness Reynolds number. Lenoir-Robillard is a 1962-vintage semi-empirical correlation; Ma 2020 has more physics. **Not directly portable.** Both, however, share the additive split: r = r_APN + r_erosive.

**Dynamic burning (Zeldovich-Novozhilov, Eq. 2.42):** Discussed in Section 2.2.3 but **NOT implemented** in SPINBALL.

### 3.5 Numerical method — Godunov

**Quasi-1D Euler in conservative form, finite-volume Godunov scheme** with **exact Riemann solver modified for gas-mixture interfaces**. First or second order in space (piecewise constant or piecewise linear with MinMod limiter). Second order in time via Heun's method. CFL condition: Δt ≤ Δx/λ_max.

**Riemann problem for mixture-of-gases:** Modified iterative algorithm. Iterates until p_2 = p_3 across the interface. **Each cell carries its own thermophysical state** — γ, Cp, R, MW vary in space and time, evaluated from local mixture concentrations.

**Critical implementation notes:**
- Grid resolution: 100–1000 cells for motor lengths 5–40 m
- ~25,000 timesteps per simulated second of motor operation (at CFL ~0.5)
- Strong unsteady discontinuity capture is mandatory for IT — ENO/Godunov chosen explicitly because PISO-class methods can't handle strong nozzle-seal rupture shocks

### 3.6 Validation against Vega — Zefiro 23

Z23 is the second solid stage of Vega, ~6 m long, HTPB1912, 23 t propellant, finocyl + aft-star, submerged nozzle, 77 s burn, MEOP 106 bar.

**IT comparison (Figure 7.5a):** SPINBALL Q1D matches the experimental head-end pressure trace during IT very well "as strong heritage of the SPIT model" (Cavallini's words).

**Quote:** "the flame spreading and the igniter jets impinging directly on the motor burning surface have some effects in terms of non-uniform regression of the motor grain geometry during the motor start-up" — Figure 7.9b shows the **erosive contribution** to burning rate is comparable to or larger than the APN term in the impingement region during IT, decaying rapidly outside it.

**Quote:** "Considering only the simulation of the IT, a blocked geometry approach for the grain geometry during the motor start-up is, as expected, an acceptable assumption, as the burning surface evolution has some effects only in the final part of the IT" — port area changes <few% during typical IT.

### 3.7 Things Cavallini is silent on

- **Detailed igniter equations** — referred to SPIT papers
- **Detailed ignition criterion equation** — referred to SPIT papers
- **Pre-ignition transient (igniter chamber thermodynamics)** — out of scope, considers mdot_ig(t) prescribed
- **Two-phase / particles** — explicitly assumes single-phase non-reacting mixture of perfect gases.

---

## 4. Cross-document synthesis

### 4.1 The 1973 → 1992 → 2009 timeline

**Common architecture (preserved across all three):**
1. **Q1D Euler in the bore** with mass/momentum/energy conservation, area-varying duct
2. **Source terms for grain combustion** (ρ_p × P_b × r_b for mass; corresponding enthalpy and zero-axial-momentum injection)
3. **Source term for igniter** (prescribed mdot_ig(t) + prescribed T_ig or H_ig)
4. **Per-cell ignition criterion** based on local solid-phase surface temperature
5. **1D heat conduction in the propellant slab** normal to the surface
6. **Choked-nozzle aft-end BC** with quasi-steady isentropic relations
7. **Saint-Robert (APN) + Lenoir-Robillard erosive burning** as the burning-rate decomposition

**Flame spread is emergent in all three.**

**What evolved:**

| Aspect | Peretz 1973 | Pardue 1992 | Cavallini 2009 |
|---|---|---|---|
| Phases | Gas only | Gas + Al2O3 | Gas mixture (6 species) |
| Numerical scheme | Implicit centered FD + predictor-corrector | SIMPLE/IPSA finite volume | Godunov + exact Riemann + ENO |
| BC at head-end | 0D plenum ODEs + characteristic | Solid wall + cell-zero source | Wall BC until seal rupture; impingement region for mdot_ig |
| Igniter model | Prescribed mdot(t), constant T_ig | Prescribed mdot(t), constant T_ig | Prescribed mdot(t), H_ig(t), **+ jet axial momentum v̄_inj** |
| Ignition criterion | Constant T_ign | Constant T_ign = 850 K | T_ign(p) pressure-dependent |
| Solid phase | Cubic-polynomial integral method (ODE) | 1D radial PDE | 1D Fourier PDE normal to surface |
| Erosive burn | Lenoir-Robillard (Eq. III-35) | Lenoir-Robillard | LR + Lawrence + Beddini corrections |
| γ, Cp, R | Constant | Constant | Variable in space/time per cell |

### 4.2 Most directly portable formulation for srm_1d

**Recommended architecture for v0.7.0** — synthesized from the three sources:

**(a) Igniter as a 0D plenum upstream of cell 0** — use the **Peretz formulation** (his Eqs. III-12 to III-15), modernized:

```
V_e · dρ_e/dt = mdot_ig(t) - mdot_e→0
V_e · c_v · d(ρ_e·T_e)/dt = mdot_ig(t)·c_p·T_ig(t) - mdot_e→0·c_p·T_e - q_loss
```

Cell-0 boundary flux is computed from (ρ_e, T_e, p_e). The two ODEs are **trivially Numba-friendly** — two scalar state variables, RK4 within the time loop.

**(b) Igniter inputs:** mdot_ig(t) as a tabulated profile, T_ig as a constant (with optional 50–65% efficiency factor à la Peretz). For phase 1, model T_ig constant.

**(c) Ignition criterion:** Start with **constant T_ign** (Peretz/Pardue form). It's robust, Numba-friendly, and matches the simplicity philosophy. 850 K as default per Pardue.

**(d) Solid-phase conduction:** **Use Peretz's cubic-polynomial integral method** (Eqs. III-44 to III-50). One ODE for δ(t) per cell. T_s(t) algebraic. Solved with RK4 in the same time loop. **Dramatically more Numba-friendly than a per-cell PDE.**

**(e) Burning rate during IT:** Keep the existing Ma 2020 erosive model.

**(f) Friction & heat transfer over un-ignited propellant:** h_c from a Bartz/Dittus-Boelter correlation. Friction over un-ignited surface from Colebrook + entrance correction. After ignition: blowing kills both h_c → 0 and f → 0 over propellant perimeter.

**(g) Multi-species:** **Skip for v0.7.0.** Treat igniter gas as same γ, Cp, MW as propellant gas (Peretz/Pardue assumption).

### 4.3 Equations to use as-is

**From Peretz 1973:**
- **Eqs. III-13, III-14, III-15** — head-end plenum ODEs (igniter coupling). Direct port.
- **Eqs. III-44, III-50** — cubic-polynomial integral method for surface temperature. Direct port. See [equations_goodman_integral.md](equations_goodman_integral.md).
- **Eq. III-31** — convective heat transfer coefficient.
- **Eq. III-33** — friction coefficient with entrance correction.

**From Pardue 1992:**
- **Eq. 15** — Bartz form for h_c. Direct port if Peretz-style isn't already implemented.
- **Eq. 14** — convective heat flux to unignited propellant: q_1 = h_c·(T_1 − T_surf).

### 4.4 Numerical-stability notes

- All three references warn about stiffness during induction → ignition transition. Peretz uses implicit centered FD (θ=0.6). Pardue uses SIMPLE (which is implicit). Cavallini uses CFL-limited Godunov.
- srm_1d's adaptive CFL should already handle the stiff cell-by-cell ignition. Watch for **non-monotone pressure gradients near the ignition front** — this is where Godunov shines and PISO/SIMPLE can ring.
- **Igniter-mass-flow ramp rate**: a step function in mdot_ig(t) is the worst case. For v0.7.0 use a finite (not zero) rise time — 5–20 ms is physically realistic and numerically benign.
- **Solid-phase ODE timestep**: RK4 on Peretz's Eq. III-50 is stable for Δt < O(δ²/α_thermal) during induction. For δ ~ 100 µm and propellant α ~ 10⁻⁷ m²/s, this gives Δt ~ 100 µs — much longer than the gas-phase CFL timestep, so the solid ODE is essentially free.
