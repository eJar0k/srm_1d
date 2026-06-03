# SRM Ignition Transient Modeling: Three-Paper Extraction (2001)

Three 2001 conference papers extracted in detail. Salita 2001 establishes the modern framework; Wang 2001 is a CFD outlook; d'Agostino 2001 is an independent, validated quasi-1D solver.

PDFs in repo:
- `salita2001.pdf` — AIAA 2001-3443
- `wang2001.pdf` — AIAA 2001-3447
- `10.2514@6.2001-3449.pdf` — AIAA 2001-3449

## 1. Salita 2001 (AIAA 2001-3443) — "Modern SRM Ignition Transient Modeling (Part 1): Introduction and Physical Models"

**Author**: Mark Salita (TRW Space and Missile Systems, San Bernardino, CA).

### 1.1 The 5-part series structure

The series, all presented at the 37th AIAA/ASME/SAE/ASEE Joint Propulsion Conference (Salt Lake City, July 2001):

- **Part 1 (Salita, AIAA 2001-3443)** — History, physical sub-models, volume-filling and 1-D equations, ignition criteria, grain preheating, water films, viscoelastic propellant response, grain overhang.
- **Part 2 (Lundgreen, AIAA 2001-3444)** — "Analyses Routinely Conducted at Alliant Techsystems." NOT IN REPO.
- **Part 3 (Luke & Meyer, AIAA 2001-3445)** — "3D Flame Spreading in Aft-end Fins" (Atlas 5). NOT IN REPO.
- **Part 4 (Graham, AIAA 2001-3446)** — "SBAS" (Structural Ballistic Analysis System). NOT IN REPO.
- **Part 5 (Wang, AIAA 2001-3447)** — "Prospective Developments in CFD Simulation." Available.

### 1.2 What is "SHARPTT"?

The title-page graph caption reads "SRM IGNITION TRANSIENT PREDICTED BY 'SHARPTT' RSRM FOR TCRIT=685K, FSR=12 m/s WITH RADIATION." This is **SHARP1DIT** (or SHARP1D-IT) — a 1-D ignition-transient solver developed at Thiokol (Rozanski, ref. 13: TWR-40265, 1990), specialised from Thiokol's 2-D CFD code SHARP. It uses Salita's ignition criterion and wet-grain modelling, allows igniter injection at any axial station, and is the workhorse used in Salita's RSRM sample case (Figs. 13–15).

### 1.3 Physical-models framework — subsystems treated

Salita's framework breaks the ignition transient into modules:

1. **Flow modelling** — choice between zero-D ("volume-filling") and 1-D conservation form.
2. **Igniter mass-flow boundary condition** — either prescribed `mdot_ig(t)` or derived from a coupled igniter-chamber pressurisation equation (Eq. 5).
3. **Grain preheating** — convective + radiative heat flux feeding a 1-D heat-conduction equation in the propellant solid (Eq. 7).
4. **Ignition criterion** — critical surface temperature, surface heat flux, in-depth temperature, or Baer/Ryan time-to-ignition correlation.
5. **Flame spreading** — sequential ignition of N axial sub-segments, each tracked independently.
6. **Wet-grain (water film) handling** — three-phase model (heat-up, vaporisation, then propellant preheat).
7. **Grain-overhang flow** — quasi-steady ΔP around forward-facing steps.
8. **Case expansion / joint rotation** — `β = dV_void/dp` feedback into volume-filling equation.
9. **Viscoelastic propellant response** — spring-dashpot whose damper coefficient scales with `dp/dt`.

### 1.4 Equations and submodel choices

**Volume-filling (VOLFIL4)** — coupled mass + energy:
```
dp_c/dt = (RT_c/V_c) (mdot_in − mdot_out)
        + (p_c/T_c) dT_c/dt
        − (p_c/V_c) dV_c/dt
        all divided by (1 + β·p_c/V_c)

dT_c/dt = (RT_c)/(p_c·c_v·V_c) · [
            mdot_in·c_p·T_in
          − mdot_out·c_p·T_c
          − (mdot_in − mdot_out)·(c_p − R)·T_c
          ]
        − Q_w / c_v
```
with `mdot_in = mdot_ig + ρ_p·S·r` and `mdot_out = α·p_c·A*/C*`.

**Igniter chamber** (Eq. 5), if not prescribed, follows the same lumped equation with V_ig and T_ig assumed constant.

**1-D conservation form (SHARP1DIT)** — Eq. 6:
```
∂E/∂t + ∂F/∂x = W + ∂R/∂x
```
with `E = [ρ, ρU, e]^T · A`, `F = [ρU, ρU² + p, (e+p)U]^T · A`. The source `W` contains:
- mass: `ρ_p·S·r + mdot_ig`
- momentum: `p · dA − ½ρU² · f · A_w + mdot_ig²/(ρ·A)`
- energy: `ρ_p·S·r·c_p·T_f − Q_w·A_w + mdot_ig·c_p^ig·T_ig`

Solved with a flux-vector-splitting TVD scheme borrowed from multi-D CFD.

**Grain-preheat equation** (Eq. 7) is the 1-D heat equation in the solid with an in-depth radiative source and Beer's-law absorption (Eqs. 8–10). For opaque propellants Salita gives the closed-form:
```
T_s(t) = T_p0 + (Q_0/k) · √(4αt/π)        (Eq. 16)
```

**Ignition criterion** — Baer/Ryan (1965, 1966) constant-flux fit:
```
t_ig^(1/2) = B / Q_0^C                                   (15)
T_ig = T_0 + √(4α/π) · (B/k) · Q_0^(1−C)                (17)
```
For AP-composite propellants: B ≈ 5.5, C ≈ 0.92 (Q_0 in cal cm⁻² s⁻¹). Because C ≈ 1, T_ig is nearly flux-independent. Predicted T_ig ≈ 646–773 K depending on flux — much lower than the historical 850 K.

**Convective heat flux** — modified Dittus–Boelter with jet-impingement correction (Eq. 14):
```
h_conv = [C_1 + C_2·F(X)] · Re^0.8 · (A_p/A*)^0.4 · (k/D_p) · G
```

**Radiative heat flux** — Salita's reduced two-parameter Mie correlation (Eq. 13):
```
ε = 0.79 · [1 - 0.215·(4.05·T_D - T_D²)] · κ / (1 + 0.84·κ)
κ = 1.5·(R_w/D_m)·(ρ_mix/ρ_d)·Φ·f_s
T_D = T_d / 1000
```
with ρ_d(Al₂O₃) = 5.632 - 0.001127·T_d g/cc. For RSRM, ε ≈ 0.4–0.5, giving Q_r ≈ 100 cal cm⁻² s⁻¹ at 3500 K — comparable to or larger than convective flux early in the transient.

### 1.5 Pyrogen vs pyrotechnic igniters

Salita addresses igniter types in one paragraph (p. 3): boosters use **head-end pyrogen** (axial or canted), space motors use **aft-end**, **head-end tangential**, or **consumable wafer** igniters; tactical motors use **bag igniters** that spew hot particles. Pyrotechnic-only igniters are not modelled distinctly — they fold into the prescribed `mdot_ig(t)` or the lumped V_ig, T_ig pressurisation. The paper does **not** give separate equations for pyrotechnic chemistry.

### 1.6 Recommendations for 1-D codes

The paper closes with four explicit guidelines (p. 10):

1. If the igniter shock does not appear in the numerical solution, something is wrong.
2. Radiative grain preheating is important; do not absorb it into a fudged convective coefficient — its geometric and time dependencies matter.
3. Use Eq. 17 with Baer/Ryan B ≈ 5.5, C ≈ 0.92 for T_ig in AP-composites; for low-pressure conditions Eqs. 15 and 17 must be modified.
4. Structural/ballistic interaction must use a viscoelastic (time-dependent) propellant model.

The RSRM sample case (Fig. 13) achieves excellent match (including the two "knees" from igniter-shock reflections) using T_crit = 685 K and a 12 m/s flame-spread rate.

---

## 2. Wang 2001 (AIAA 2001-3447) — "Prospective Developments in CFD Simulation"

**Author**: J.C.T. Wang (Aerospace Corporation).

### 2.1 What CFD developments does Wang anticipate as essential?

Wang's argument is that 3-D CFD is now tractable but expensive, and that several enabling pieces are needed:

- **TVD shock-capturing schemes** — without them the igniter shock and subsequent compression-wave reflections cannot be resolved.
- **Parallel computing with domain decomposition** — Wang cites a Delta II 7925 plume-on simulation (Intel Paragon, 506 nodes, 48 hours) vs. ~3 months on a vector Cray-YMP.
- **Quasi-3-D as a stop-gap** — axisymmetric in the bore, 2-D Cartesian in the fins, joined at block interfaces. He claims this exposes 3-D phenomena (reverse flow, fin-induced pressure oscillations) that 1-D and 2-D miss, at 36 CPU-hours on Cray-SV1 vs. far more for full 3-D.

### 2.2 Multi-species, particles, turbulence

Wang is **silent** on multi-species transport in his own simulations — Eqs. (5)–(9) are single-species ideal-gas Euler. Particle phases / two-phase flow are not addressed.

Turbulence: he reviews two prior efforts but Wang's own runs are **inviscid Euler with empirical heat-transfer correlations bolted on** — he explicitly states "the flow will be assumed inviscid and the simplest heat transfer and ignition models will be applied" (p. 3).

### 2.3 Boundary condition for the igniter

Wang is sparse here. His simulations are "driven by a specified input for igniter mass flow rate and enthalpy" (p. 7). No igniter chamber is modelled — `mdot_ig(t)` and `h_ig(t)` are prescribed at the head-end inflow patch.

### 2.4 Why are simpler models insufficient?

Wang is moderate, not dismissive. He explicitly notes (p. 1) that the volume-filling method, applied as **two interconnected volumes**, predicted the RSRM head-end pressure rise rate "in good agreement" with 1-D CFD. His case for 3-D rests on:

- Star/fin propellant grains create flow features (large reverse-flow zones, fin-tip pressure oscillations) that 1-D and 2-D cannot capture.
- Configuration changes can materially alter the transient pressure history at fixed L/D and equivalent fin area.
- For propellant-crack effects on performance, full 3-D is required.

### 2.5 Heat-transfer and ignition models reviewed

Wang reproduces the Churchill / Duhamel integral form of the surface-temperature equation (Eq. 1):
```
T_s(t,x) = T_p + √(α/π) · ∫₀^t (q(τ,x)/k) · dτ/√(t-τ)
```
He notes this is memory-expensive in multi-D.

Convective flux (Johnston's form, Eq. 3):
```
q_c(t,x) = C_c · Pr^(-2/3) · c_p · (μ/L)^0.2 · (ρ|v|)^0.8 · (T_ad − T_s)
```

Radiative flux (Eq. 4):
```
q_r = (1−η) · C_r1 · σ · (T⁴ − T_s⁴)
    + η     · C_r2 · σ · (T_flame⁴ − T_s⁴)
```
with `η = ΣA_lm/(ΣA_lm + A_ij)` a weighting that ramps from "view-factor to general gas" to "view-factor to adjacent ignited cells" as flame spread proceeds.

Burning rate is plain Saint-Robert `r = a·p^κ` (Eq. 10), no erosive term.

### 2.6 Useful for srm_1d

Wang's paper is most useful as **rationale for what srm_1d intentionally does not do**: no turbulence model, no 3-D, no fin-resolved geometry. The radiative-flux blending function `η` is a clean idea that could inform how srm_1d's ignition front transitions a cell from "igniter heating" to "neighbour-flame heating." Otherwise this paper is largely orthogonal to srm_1d's amateur 1-D scope.

---

## 3. d'Agostino, Biagioni & Lamberti 2001 (AIAA 2001-3449)

**Authors**: L. d'Agostino (University of Pisa), L. Biagioni and G. Lamberti (Centrospazio).

### 3.1 What's different from Salita

Where Salita's Part 1 is a survey + checklist of physical models, d'Agostino et al. present a **complete, validated, single-coordinate solver** for the entire ignition transient — and they validate against full Ariane 4 and Ariane 5 static-firing pressure traces with sub-1% error. Major differences:

- **Erosive burning is built in** as a Lenoir–Robillard term (Salita does not include erosive burning in his ignition equations).
- **Dimensionality**: pure quasi-1-D throughout; no volume-filling fallback.
- **Igniter chamber**: not modelled — `mdot_ig(t)` is prescribed (Fig. 6 shows the experimentally measured Ariane 4 igniter mass-flow input).
- **Numerics**: modified Lax–Friedrichs on a finite-volume variable-step grid, with method-of-characteristics compatibility relations at boundaries.
- **Radiative preheating**: lumped into a single empirical correction factor `C_hc` on the convective coefficient (their Fig. 4). Salita explicitly warns against this approach.

### 3.2 The math — equations

**Gas-phase quasi-1-D conservation laws** (with A_p port area, P_b burning perimeter, r_b burning rate, ρ_pr grain density, h_f gas enthalpy):

Continuity:
```
∂ρ/∂t + ∂(ρu)/∂x + (ρu/A_p)·∂A_p/∂x = (r_b·P_b/A_p)·ρ_pr
```

Momentum:
```
∂(ρu)/∂t + ∂(ρu²)/∂x + (ρu²/A_p)·∂A_p/∂x = -∂p/∂x − (r_b·P_b/A_p)·ρu
```

Energy:
```
∂(ρe_T)/∂t + ∂(ρu·e_T)/∂x = -∂(p·u)/∂x + (r_b·P_b/A_p)·ρ_pr·h_f
```

Terms explicitly dropped: convective heat loss to walls, viscous stresses, axial heat conduction, mass accumulation in the volume opened by surface regression.

**Solid-phase surface temperature** (their integrated-balance "BI" method, ODE per axial node):
```
dT_ps/dt = [4·α_pr·h_c²·(T − T_ps)³] / [3·k_pr²·(T_ps − T_pi)·(2T − T_ps − T_pi)]
```
The t = 0 singularity is sidestepped by initialising T_ps(0) = T_pi + ε_T. **Ignition criterion** is critical surface temperature — same as Salita.

**Heat-transfer correlation** — Dittus–Boelter with one empirical lumped factor:
```
h_c = [0.026·(ρu)^0.8·(μ/D_H)^0.2·Pr^0.6·c_p] · C_hc
```
C_hc absorbs radiation, conduction, geometry, and jet impingement. For Ariane 4 it varies from ~2.5 at the head end to ~0.5 at the aft end (Fig. 4).

**Burning rate** — modified Lenoir–Robillard:
```
r_b = a·p^n + k_eb·h_c·exp(-β·r_b·ρ_pr/(ρ·u))
```
β (erosive exponent) is calibrated against water-quench tests + experimental p-t curves.

**System form**:
```
∂U/∂t + ∂F/∂x = S
U = [ρ, ρu, ρe_T]^T
F = [ρu, ρu² + p, u·(ρe_T + p)]^T
```

### 3.3 Validation cases

- **Ariane 4 SRM** (segmented; Fig. 5 shows a single-port grain) — used to **calibrate** C_hc and β. Figs. 13 vs. 14: without erosive burning the predicted peak underestimates and lags; with erosive burning the match is essentially exact.
- **Ariane 5 SRM** (segmented, larger scale) — used as a **blind validation**. Pressure trace shows clear discontinuities at the inter-segment slots. Match with experimental data is also very good (Fig. 16). **Stated overall error: <1%** for both motors.
- The same code is then applied to the entire burn life with an empirical "BRAF" (Burning Rate Anomaly Factor; Fig. 17) correction over web fraction.

### 3.4 Single-species or multi-species

**Single-species perfect gas** throughout. Stated explicitly in assumption 4: "the chamber gas obeys the perfect gas law." The igniter and main-propellant gases are distinct only in their prescribed enthalpies h_f; no species transport equation is solved. **This matches srm_1d's current frozen / effective single-gas approach.**

### 3.5 Numba-JIT tractability

This paper's formulation is **the most directly translatable to a Numba-JIT 1-D code** of the three:

- Identical conservation form to srm_1d's Euler-type 1-D solver, with three state variables [ρ, ρu, ρe_T]^T and a source vector that already isolates the burning-surface injection terms.
- The solid-phase surface-temperature ODE is a per-cell scalar update — trivially `@njit` with one array per axial node.
- The Lenoir–Robillard erosive term shares structure with srm_1d's existing Ma-2020 implementation (both are surface-velocity-dependent corrections); only the functional form differs.
- The Lax–Friedrichs scheme is simpler than the PISO+TDMA already in srm_1d. (srm_1d's current scheme is more accurate; Lax–Friedrichs is mentioned only because it's what they used.)
- The igniter is **prescribed** mdot_ig(t), h_f^ig — matching exactly the v0.6.0 placeholder structure and the v0.7.0 hot-gas-plenum target.

The single tractability concern: the C_hc(x/L) profile is not generic — it's calibrated per motor from experimental data. For srm_1d's amateur use case this becomes a tuning knob, not a predictive input. Salita's separately-modelled radiation (his ε from Eq. 13) is more physically grounded but requires Al₂O₃ size distribution data not available for amateur propellants.

---

## 4. Cross-paper synthesis

### 4.1 Compatibility / contradictions

The three papers are largely **compatible at the conservation-law level** and **diverge at the heat-flux and igniter levels**:

| Aspect | Salita (2001-3443) | Wang (2001-3447) | d'Agostino (2001-3449) |
|---|---|---|---|
| Dimensionality | 0-D and 1-D | quasi-3-D, axisymmetric | quasi-1-D |
| Gas model | single calorically-perfect | single calorically-perfect | single calorically-perfect |
| Gas equations | conservation form, TVD flux-split | Euler conservation, Wang–Widhopf TVD | conservation form, Lax–Friedrichs |
| Igniter | prescribed OR coupled V_ig, T_ig pressurisation | prescribed mdot, h | prescribed mdot, prescribed h_f^ig |
| Radiative preheat | explicit ε·σ·(T_d⁴ − T_s⁴) via Mie/Salita correlation | explicit, with cell-view-factor blending η | folded into empirical C_hc |
| Convective preheat | Dittus–Boelter + jet impingement | Dittus–Boelter form (Johnston) | Dittus–Boelter × C_hc |
| Ignition criterion | Baer/Ryan T_ig ≈ 646–773 K | Churchill integral; uncoupled | critical T_ps via BI ODE |
| Erosive burning | not in main equations | not modelled (r = ap^κ) | Lenoir–Robillard built in |
| Flame spreading | sequential ignition of N axial sub-segments | per-cell ignition state | sequential node ignition |

**Direct contradiction**: Salita's Recommendation 2 ("don't lump radiation into convection") is explicitly violated by d'Agostino's C_hc approach. d'Agostino's <1% match on Ariane 4 and 5 is the empirical counter-argument.

**Direct contradiction**: Salita argues the historical 850 K ignition temperature is too high; his RSRM match works at 685 K.

### 4.2 What all three agree on (the "must-have" list)

1. **1-D unsteady Euler conservation form** for the gas, with explicit injection source terms tied to burning rate.
2. **Critical surface temperature ignition criterion**, derived from a 1-D heat-conduction equation in the solid.
3. **Quasi-steady burning rate** at any instant (no dynamic-burning ODE).
4. **Sequential ignition of axial nodes** — flame spreading is each cell's surface temperature crossing T_ig.
5. **Prescribed (not modelled) igniter mass-flow profile** is the production approach; coupled igniter chambers are optional.
6. **Single calorically-perfect gas** is sufficient at this fidelity.
7. **Igniter shock propagation must appear naturally** in the solution.

### 4.3 Most defensible subset for srm_1d (amateur 1-D)

- **Adopt d'Agostino's quasi-1-D conservation form wholesale** — structurally identical to what srm_1d already has and the validation is the strongest of the three (sub-1% on two real motors).
- **Adopt Salita's recommendation on T_ig**: use Baer/Ryan Eq. 17 with B ≈ 5.5, C ≈ 0.92 rather than a hard-coded 850 K.
- **For radiation, the honest answer for amateurs is d'Agostino's C_hc lumping** — Salita's ε formula needs Al₂O₃ size distributions not available for hobby propellants. Document as known limitation.
- **For the igniter (the v0.7.0 target)**: adopt the prescribed-mdot_ig(t), h_f^ig approach used by both d'Agostino and Wang. Salita's optional coupled-igniter formulation (Eq. 5) is the v0.7.0 target since it provides pressure feedback.
- **Erosive burning** — srm_1d already has Ma 2020, which is more sophisticated than Lenoir–Robillard. Keep it.

### 4.4 Relation to NASA SP-8051 (m = 0.12·V_F^0.7)

None of the three papers cite SP-8051 (1971) by name. SP-8051's empirical sizing rule is a **boundary-condition input**: it tells the designer "for chamber free volume V_F, deliver an igniter mass flow of approximately 0.12·V_F^0.7." It says nothing about the time-history shape or the gas temperature.

**All three papers consume that boundary condition** — Salita prescribes mdot_ig(t) or pressurises a lumped igniter volume; Wang prescribes mdot_ig(t) and h; d'Agostino prescribes mdot_ig(t) measured experimentally.

These models therefore **make SP-8051 deeper, not obsolete**. The two are complementary.
