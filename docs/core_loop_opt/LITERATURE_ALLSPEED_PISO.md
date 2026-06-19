# Literature review: numerical methods for 1D transient SRM internal ballistics

**Prepared 2026-06-18** by a background research agent, for the
core-loop-optimization thread (see [`README.md`](README.md)).
**Scope:** Performance (low-Mach CFL escape), Correctness (ignition
fill), Domain practice (established SRM codes).

> Solver context given to the agent: pressure-based PISO, non-iterated,
> 2 pressure correctors, staggered grid, TDMA/Thomas pressure solve,
> adaptive CFL; compressible ideal gas, quasi-1D port with distributed
> mass + heat addition, choked-nozzle BC; operating point M ≈ 0.3–0.4.

---

## GOAL 1 — PERFORMANCE: escaping the acoustic CFL at low Mach

An explicit-acoustic PISO has `dt = CFL·dx/(u+a)`, with `a ≈ 1000 m/s`
dominating at M ≈ 0.3–0.4 (`u ≈ 300–400 m/s`). The acoustic CFL is
`(1 + 1/M) ≈ 3.5×` more restrictive than the convective CFL
`dt = CFL·dx/u`. Recovering that factor is the whole perf question.

### 1A. Semi-implicit / IMEX acoustic–convective splitting  ★ recommended
Split into slow (convective/material) + fast (acoustic) subsystems;
integrate material waves explicitly (upwind/WENO), acoustic pressure-wave
subsystem implicitly via a linear elliptic solve. Klein (1995, §4) gives
the exact 1D staggered form: acoustic flux (`dp/dx` in momentum, `u·a²`
in density) → implicit side; all convective terms (`u·du/dx`, mass
injection, enthalpy advection) → explicit side. Degond–Tang (2011) /
Haack–Jin–Liu (2012) give compact asymptotic-preserving (AP) versions:
one scalar elliptic (density/pressure) solve per step. Boscheri–Pareschi
(2021) extend to 3rd order, no iterative nonlinear solve.
- **CFL limit:** convective only (`dt ≤ CFL·dx/u`) → ~3× larger steps at M=0.35.
- **Cost/step:** one explicit hyperbolic sweep + one tridiagonal solve
  (O(N), same as the existing TDMA corrector).
- **Net speedup vs our PISO:** ~3× fewer steps, ~1.3–1.5× cost/step →
  **~2× wall-clock.** Parameter-free.
- **Implementation cost (ours):** MODERATE. Staggered grid + TDMA already
  present. Reclassify the acoustic flux as implicit (linearized about
  current state), solve the p–ρ elliptic eq (one TDMA pass, existing
  structure), update velocity explicitly. ~50–150 LOC in `solver.py` +
  careful source-term classification (combustion mass addition explicit;
  the pressure rise it drives is implicit). Refactors the JIT
  `_run_time_loop` into split explicit/implicit stages.

### 1B. Low-Mach preconditioning (density-based route) — NOT recommended
Turkel-style preconditioner rescales the time-derivative so eigenvalues
are O(u). Great for **steady-state** density-based convergence, but
corrupts **unsteady** acoustic timing — wrong for our transient/ignition
physics. Also requires converting to a density-based solver. Skip.

### 1C. All-speed pressure-based (modified Rhie–Chow / AUSM⁺-up)  ★ alt
Bauer–Zeidan (2013) modify Rhie–Chow to retain a dt-independent acoustic
correction (standard RC coupling vanishes as dt→0, killing it at large
dt/low M). Liou (2006) AUSM⁺-up splits convective+pressure fluxes with
Mach-weighted functions that become one-sided for M>1 — fixes low-M
decoupling AND high-M choking at once. Miettinen (2015): AUSM⁺-up is
time-step-independent, robust at large dt.
- **CFL limit:** convective. **Cost/step:** same as current PISO (only
  the face-flux formula changes). **Speedup:** ~3× at M=0.35.
- **Implementation cost (ours):** LOW–MODERATE — modify the momentum
  face-flux assembly to AUSM⁺-up pressure splitting / acoustic-correction
  term. Localized to face-flux assembly in `solver.py`. **Serves Goal 2 too.**

### 1D. Multiple Pressure Variables (MPV) — NOT recommended
Klein/Munz scale separation for wide atmospheric/combustion domains.
Poor fit: our closed domain has physically present acoustics coupled by
combustion feedback + a choked nozzle reflector. HIGH cost.

**GOAL 1 bottom line.** Implement **IMEX acoustic–convective splitting
(1A, Klein/Degond–Tang)** on the existing staggered+TDMA structure:
convective-CFL stepping (~3× fewer steps), ~1.3–1.5× cost/step → **~2×
net**, zero tuning. Lower-effort alternative with the same payoff that
*also* fixes Goal 2: **AUSM⁺-up face flux (1C, Liou 2006).**

---

## GOAL 2 — CORRECTNESS/STABILITY: the high-Mach ignition fill

Igniter forces hot gas into a near-empty, cold, not-yet-choked port;
pressure-based PISO has no interior choking limit and is not
shock-capturing → contact velocity reaches grid-divergent M≈3–12.
Physically a constant/expanding-area duct with distributed mass addition
cannot exceed M=1 without a C-D passage (Fanno–Rayleigh–Shapiro choking).

- **Physics (2A):** distributed mass addition drives both sub- and
  supersonic flow toward M=1 (Shapiro 1953; Anderson 1982). Real
  phenomenon in high-L/D motors (Peretz 1973; Laubacher 2000) — but our
  artifact is **grid-divergent** (→ primarily numerical), while the
  pressure spike is **grid-convergent** (the real erosive response) and
  **decoupled** from the velocity artifact (per
  `IGNITION_SPIKE_REOPENED.md`).
- **Pressure-based behavior at M~1–12 (2B):** SIMPLE/PISO lose fidelity
  above M≈0.5 (linearized continuity assumes low-M acoustics; RC becomes
  very diffusive). AUSM⁺-up extends faithfully through/above M=1. SPINBALL
  uses density-based Godunov + Riemann + minmod → naturally shock-capturing.
- **Parameter-free closures (2C):**
  1. **`port_mach_cap` velocity clip (already shipped, commit `16bc527`)**
     — clips face velocity to `≤ a_local`, parameter-free, grid-converges
     the Mach field, does NOT change P_peak. Sufficient as solver hygiene.
  2. **Sonic-cell upwinding** — when a face detects M→1, switch its
     pressure flux to one-sided/upwinded (local Riemann BC). ~20 LOC,
     zero parameters.
  3. **AUSM⁺-up face flux** — shock-consistent through M=1, parameter-free,
     also serves Goal 1.
  4. Density-based Godunov rewrite — most robust, NOT warranted.

**GOAL 2 bottom line.** Existing **`port_mach_cap` is sufficient** for the
decoupled artifact. Principled long-term fix = **AUSM⁺-up (also serves
Goal 1)**; Goal-2-only = **sonic-cell upwinding (~20 LOC)**. No
preconditioning/relaxation (those need tuning).

---

## GOAL 3 — DOMAIN PRACTICE: what established SRM codes do

- **SPINBALL / Q1D (Sapienza — Di Giacinto, Favini, Cavallini,
  Serraglia):** Godunov-type FV, exact Riemann solver + minmod, 1st/2nd
  order, **density-based, conservative, shock-capturing**. Acoustic CFL
  fully explicit (not flagged as a bottleneck). Riemann solver handles
  M→1 and supersonic fill with no special closure; one model ignition→burnout.
- **SPP (Solid Performance Program, JANNAF standard):** IGT + TAB 1D
  modules; heritage = density-based quasi-1D explicit FD (MoC or MacCormack
  predictor-corrector), per Peretz 1973. Not pressure-based. Details
  distribution-controlled.
- **NAWC / heritage (Peretz–Kuo–Caveny–Summerfield 1973):** foundational
  Q1D ignition-transient paper; explicit density-based Lax–Wendroff-type
  predictor-corrector; explicitly treats high internal gas velocity
  (M~0.5–1, aft-end) and pressure overpeak at high L/D and loading density.
- **Modern SRM CFD:** fully implicit compressible **density-based**
  (Roe/HLLC, RANS k-ω/SST).

**Key finding:** *Every* established SRM Q1D code is **density-based** and
handles the full Mach range natively. srm_1d's **pressure-based PISO is
the domain outlier** — well-suited to M≈0.3 quasi-steady, structurally
ill-matched to the fill. `port_mach_cap` + the Goal-1 IMEX upgrade bring
srm_1d into conformance **without a density-based rewrite.**

---

## CONSOLIDATED REFERENCES (DOIs recorded; access flagged)

| # | Citation | DOI / ID | Access | Priority |
|---|----------|----------|--------|----------|
| 1 | Klein, R. "Semi-implicit extension of a Godunov-type scheme based on low Mach number asymptotics I." *J. Comput. Phys.* 121(2), 213–237, 1995. | 10.1016/0021-9991(95)90034-9 | PAYWALLED | **HIGH (Goal 1 core)** |
| 2 | Degond, P. & Tang, M. "All speed scheme for the low Mach number limit of the isentropic Euler equation." *Commun. Comput. Phys.* 10(1), 1–31, 2011. | arXiv:0908.1929; CiCP DOI unverified | OPEN (arXiv) | **HIGH (Goal 1 scheme)** |
| 3 | Haack, Jin, Liu. "An All-Speed Asymptotic-Preserving Method for the Isentropic Euler and Navier–Stokes Equations." *Commun. Comput. Phys.* 12(4), 955–980, 2012. | CiCP v12 DOI unverified | OPEN (author PDF: sites.math.duke.edu/~jliu) | **HIGH (Goal 1 AP)** |
| 4 | Boscheri, W. & Pareschi, L. "High order pressure-based semi-implicit IMEX schemes for the 3D NS equations at all Mach numbers." *J. Comput. Phys.* 434, 110206, 2021. | 10.1016/j.jcp.2021.110206; arXiv:2008.01789 | OPEN (arXiv) | Medium |
| 5 | Modesti, D. & Pirozzoli, S. "An efficient semi-implicit solver for DNS of compressible flows at all speeds." 2016. | arXiv:1608.08513 | OPEN | Medium |
| 6 | Liou, M.-S. "A sequel to AUSM, Part II: AUSM⁺-up for all speeds." *J. Comput. Phys.* 214(1), 137–170, 2006. | 10.1016/j.jcp.2005.09.020 | PAYWALLED | **HIGH (Goals 1 & 2)** |
| 7 | Turkel, E. "Convergence acceleration for computing steady-state compressible flow at low Mach numbers." *Comput. Fluids* 27(5–8), 385–404, 1998. | 10.1016/S0045-7930(98)00058-9 | PAYWALLED | Low (not our path) |
| 8 | Guillard, H. & Viozat, C. "On the behavior of upwind schemes in the low Mach number limit." *Comput. Fluids* 28(1), 63–86, 1999. | DOI unverified | PAYWALLED | Medium (background) |
| 9 | Miettinen, A. "Application of pressure- and density-based methods for different flow speeds." *Int. J. Numer. Methods Fluids* 79(5), 243–267, 2015. | 10.1002/fld.4051 | OPEN (Wiley) | Medium |
| 10 | Cavallini, Favini, Di Giacinto, Serraglia. "SRM Internal Ballistic Numerical Simulation by SPINBALL Model." AIAA 2009-5512. | 10.2514/6.2009-5512 | PAYWALLED | **HIGH (SPINBALL)** |
| 11 | Cavallini, E. "Modeling and Numerical Simulation of Solid Rocket Motors Internal Ballistics." PhD thesis, Sapienza, 2010. | handle 11573/918117 (iris.uniroma1.it) | OPEN (institutional) | **HIGH (SPINBALL detail)** |
| 12 | Di Giacinto, Favini, Cavallini, Serraglia. "An Ignition-to-Burn Out Analysis of SRM Internal Ballistic and Performances." AIAA 2008-4749. | 10.2514/6.2008-4749 | PAYWALLED | High |
| 13 | Peretz, Kuo, Caveny, Summerfield. "Starting Transient of Solid Propellant Rocket Motors with High Internal Gas Velocities." *AIAA J.* 11(12), 1719–1727, 1973. | 10.2514/3.50676; NTRS 19740005393 | OPEN (NTRS) | **HIGH (foundational Q1D)** |
| 14 | Laubacher, B.A. "Internal Flow Analysis of Large L/D Solid Rocket Motors." AIAA-2000-3803. | NTRS 20000064699 | OPEN (NTRS) | Medium |
| 15 | Issa, R.I. "Solution of the implicitly discretised fluid flow equations by operator-splitting." *J. Comput. Phys.* 62(1), 40–65, 1986. | 10.1016/0021-9991(86)90099-9 | PAYWALLED | Medium (PISO origin) |
| 16 | Moukalled, Mangani, Darwish. *The Finite Volume Method in CFD (OpenFOAM/Matlab).* Springer, 2015. | ISBN 978-3-319-16873-9 | Textbook | Low |
| 17 | Chalons, Girardin, Kokh. "Large Time Step and AP Numerical Schemes for the Gas Dynamics Equations with Source Terms." *SIAM J. Sci. Comput.* 35(6), A2874–A2902, 2013. | 10.1137/130908671 | PAYWALLED | Medium |
| 18 | Thornber, Mosedale, Drikakis. "An improved reconstruction method for compressible flows with low Mach number features." *J. Comput. Phys.* 227(10), 4873–4894, 2008. | 10.1016/j.jcp.2008.01.036 | PAYWALLED | Low |
| 19 | Cavallini et al. "Internal Ballistics Simulation of NAWC Tactical Motors with SPINBALL Model." AIAA 2010-7163. | 10.2514/6.2010-7163 | PAYWALLED | Medium |
| 20 | "A new simulation strategy for solid rocket motor ignition (CFD + 1D boundary flame)." AIAA 2021-3695. | 10.2514/6.2021-3695 | PAYWALLED | Low |

**Open-access (direct download):** #2, 3, 4, 5, 11, 13, 14, 9.
**Top paywalled to obtain:** #1 (Klein), #6 (Liou AUSM⁺-up), #10 (SPINBALL),
#13 (Peretz — also on NTRS, so effectively open).
**Fabrication guard:** entries marked "DOI unverified" were not confirmed
by the agent — verify before citing.
