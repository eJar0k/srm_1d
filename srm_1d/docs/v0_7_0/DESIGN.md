# srm_1d v0.7.0 — Hot-gas Plenum Igniter Model

**Status**: Phase 4 complete on branch `v0.7.0-phase4`; ready to tag
v0.7.0. The Hasegawa A v0.7.0 calibration LHS rank-1 produces
`mse_all = 0.0968 MPa²` (better than v0.6.0's 0.24 MPa² with the
now-removed `igniter_tau = 127 ms` FSI proxy). See the per-step
journey in `audits/2026-05-20_radiation_collapse_localT.md` and
`audits/2026-05-21_hasegawa_a_lhs_v0_7_0.md`.

**Target**: replace the v0.6.0 single-knob exponential-decay igniter
with a physically grounded forward model that removes `igniter_tau` as
an FSI-cushioning calibration knob. **Done.**

This document is self-contained. A fresh coding agent should be able to
check out the repo, read this doc + [TASKS.md](TASKS.md) + the
[references/](references/) directory, understand the implemented
v0.7.0 design, and continue Phase 4 validation without re-deriving the
literature.

## Context — why this is needed

Both Hasegawa A and Zerox calibrations (committed in v0.6.0) hit the
same residual: ~25% pressure-spike overshoot at ignition. The LHS in
[zerox_lhs.py](../../examples/zerox_lhs.py) showed this is a **structural
artifact**, not a parametric one — no combination of `igniter_mass`,
`igniter_tau`, `ignition_ramp_tau`, `P_ignition`, `kappa` brings the
spike to within experimental. The current placeholder

```
mdot_igniter(t) = (m_init/τ) · exp(-t/τ)
```

distributed full-grain has neither pressure feedback nor temperature, so
LHS-tuned `igniter_tau ≈ 127 ms` for Hasegawa A is a numerical
FSI-cushioning proxy, not a physical timescale (per
[project_hasegawa_calibration_state](https://github.com/eJar0k/srm_1d) memory).
v0.7.0 replaces this with the canonical pyrogen architecture.

## Architectural decisions

After reading 7 papers (Ma 2019, Salita 2001, Wang 2001, d'Agostino 2001,
Peretz 1973, Pardue & Han 1992, Cavallini 2009) plus Sutton 9e and the
DeMar 2021 amateur deck, the literature converges on the following
"must-haves" for a 1D ignition transient model. v0.7.0 adopts all of
them:

1. **Forward 0D pyrogen plenum** with its own state `(P_ig, T_ig, m_ig)`,
   integrated alongside the main loop.
2. **Choked sonic orifice** boundary condition coupling plenum to cell 0
   of the main motor. Subsonic fallback (rare; pyrogens designed to stay
   choked).
3. **Per-cell critical surface temperature ignition criterion**.
4. **Goodman cubic-polynomial integral method** for the solid-phase
   conduction (one ODE per cell — Numba-friendly, ~5% error vs. exact
   PDE; Peretz 1973 Eqs. III-44 to III-50). See
   [equations_goodman_integral.md](references/equations_goodman_integral.md).
5. **Single calorically-perfect gas** — igniter and main propellant share
   γ, MW, T_flame at the bulk-flow level. Multi-species transport
   deferred to v0.7.x.
6. **Sequential per-cell ignition**: flame spread is *emergent* from the
   per-cell criterion, not separately modeled.
7. **No squib stage in v0.7.0**: pyrogen ignites instantly at t=0
   (instantaneous full burning surface, propellant at T_flame). v0.7.x
   adds a squib (electric → squib → pyrogen → main).

8. **Robust ignition-transient nozzle boundary**: the aft throat is a
   signed open boundary, not a permanent choked sink. It supports
   subsonic outflow, choked outflow, reverse ambient inflow, and
   un-choking during failed ignition or taildown.
9. **Energy-audited surface heating and spread**: DeMar pyrogen
   surface heating and adjacent-burning-cell radiation feed the Goodman
   surface update and are recorded in energy/momentum diagnostics.

## Pyrogen chamber model

State: `m_p` (pyrogen propellant remaining mass, kg), `m_ig` (gas mass
in plenum, kg), `T_ig` (plenum temperature, K). Plenum volume `V_plenum`
is fixed — pyrogen burns inside it; no cavity expansion modeled.

```
dm_p/dt    = -ρ_pyro · A_burn(m_p) · r_pyro                  (1)
dm_ig/dt   = (-dm_p/dt) - mdot_choke                          (2)
d(m_ig·c_v·T_ig)/dt = (-dm_p/dt) · c_p · T_flame
                     - mdot_choke · c_p · T_ig                (3)
P_ig       = m_ig · R / (M_ig · V_plenum)                     (4)
```

where:
- `r_pyro = a_pyro · P_ig^n_pyro` — Saint-Robert burn rate
- `A_burn(m_p)` — burning surface as a function of remaining mass.
  For amateur **0D pyrogen geometries** (loose flakes, BKNO3/MTV pellets,
  dipped squibs), treat `A_burn = (m_p/m_p,0)^(2/3) · A_burn,0` (sphere-equivalent
  surface area regression). For COTS hobby igniter shapes (long
  cylindrical 10:1 dipped squibs), use `A_burn = constant ≈ A_burn,0`
  until burnout (end-burning approximation).

**Choked outflow** (when `P_main / P_ig < (2/(γ+1))^(γ/(γ-1))`,
i.e., critical pressure ratio):

```
mdot_choke = P_ig · A_t,pyro · Γ / √(R · T_ig / M_ig)        (5)
Γ = √(γ_pyro · (2/(γ_pyro+1))^((γ_pyro+1)/(γ_pyro-1)))
```

For amateur hobby pyrogens which often don't have a discrete throat
(e.g., loose pyrogen in a paper tube), `A_t,pyro` is a conceptual
quantity — the user specifies it as the effective vent area through the
container wall + igniter leads + any deliberate orifice.

**Subsonic fallback** (`P_main/P_ig ≥` critical ratio):

```
mdot_subsonic = P_ig · A_t,pyro · √(2γ/(γ-1) · M_ig/(R·T_ig)) ·
                √[(P_main/P_ig)^(2/γ) - (P_main/P_ig)^((γ+1)/γ)]   (6)
```

This case is rare in practice; pyrogens are designed to stay choked
throughout. Implement as a fallback with a warning.

## Coupling to cell 0

Inject `mdot_choke` and enthalpy `h_inject = c_p · T_ig` into cell 0 as
source terms. The implementation uses separate per-cell `mass_source`
and `thermal_source` arrays: propellant/end-face sources contribute at
`T_flame`, and pyrogen contributes at `T_ig`.

Igniter momentum (`mdot_choke · v_inject`) is now an explicit
face-centered momentum source. The default projects the pyrogen orifice
momentum into face 1 of the bore as
`mdot_ig*v_exit/(A_face*dx)`. Diagnostics record expected force,
deposited force, and residual so the source can be verified without a
tuning multiplier.

Direct pyrogen-to-propellant surface heating uses DeMar measured heat
flux for the pyrogen composition. It applies only to the first unignited
grain cell, feeds Goodman through an equivalent heat-transfer
coefficient/driver temperature, and subtracts the same delivered power
from the gas temperature-source ledger.

Adjacent-burning-cell radiation is implemented as a geometry-local
opt-in: only unignited grain cells adjacent to burning cells receive
the radiation term, and the emitter temperature is the local gas
`T[neighbor]` (not the constant adiabatic `T_flame`, which overstated
the cold-start transient). The model uses Stefan-Boltzmann transfer
with `Propellant.radiation_emissivity` as a material property and
debits the same emitted power from the burning-cell gas energy ledger.

The default `radiation_emissivity` is `0.0` for all propellants
(adapter no longer auto-defaults aluminized .ric files to 0.45). A
2026-05-14 sweep on Hasegawa A showed that turning radiation on under
ambient initial gas drives an unphysical ignition chain (~1 ms/cell)
that pushes interior flow supersonic before the signed-throat PISO
boundary can vent, producing a non-monotonic discrete resonance over
both `radiation_emissivity` and any `tau_establishment` ramp meant to
slow it. Sutton 9e Section 15.3 also documents pyrogen ignition as
primarily convective. Set `radiation_emissivity` explicitly in the
.ric (or on `Propellant`) to opt back in once the spread-rate /
PISO-stability interaction is understood.

Multi-cell impingement-region distribution (Cavallini SPINBALL) remains
deferred unless a motor provides physical igniter basket/jet geometry.

## Ignition criterion

Replace the current global `is_burning` boolean (which fires when
`P_head > P_ignition`) with a **per-cell `T_surf[i]` field**. Cell `i`
ignites when `T_surf[i] > T_ignition`. Default `T_ignition = 850 K`
per Pardue & Han 1992; for AP-composite propellants, Salita 2001's
Baer/Ryan correlation gives a more nuanced 646-773 K depending on local
heat flux. Treat as a single user input for v0.7.0, default 850 K, with
documentation that this is conservative.

Once `is_burning[i] = True`:
- The existing Ma 2020 erosive burn-rate kernel takes over for cell `i`
- Convective coefficient `h_c[i] → 0` over the burning surface (large
  blowing kills convection — Peretz, Pardue, d'Agostino all agree)
- `T_surf[i]` is no longer integrated (the thermal layer is now
  combusting; integration is irrelevant)

## Solid-phase conduction (Goodman cubic-polynomial integral)

Per cell, one ODE:

```
dδ/dt = 12α · (3k + h_c·δ) / [δ · (6k + h_c·δ)]              (7)
T_surf = (3k·T_initial + h_c·δ·T_gas) / (3k + h_c·δ)         (8)
```

where:
- `α = k/(ρ·c)` — propellant thermal diffusivity (m²/s)
- `k` — propellant thermal conductivity (W/m·K)
- `h_c[i]` — convective coefficient at cell `i` (already in srm_1d
  via Gnielinski/Haaland — see [burn_rate.py](../../burn_rate.py))
- `T_gas[i]` — local gas temperature
- `T_initial` — propellant bulk temperature (room temp ~290 K)
- `δ[i]` — Goodman penetration depth (m), per-cell state variable

Initial condition: `δ[i](t=0) = ε` (small positive, ~1 µm) to avoid the
ODE singularity at δ=0. RK4 integration in the same time loop. See
[equations_goodman_integral.md](references/equations_goodman_integral.md)
for the full derivation including the cubic-polynomial profile, BCs,
heat-balance integral, and asymptotic behaviors.

## Pyrogen sizing (default behavior)

When the user doesn't specify `pyrogen_mass`, default to **Sutton Eq.
15-4**:

```
m_pyrogen [grams] = 0.12 · V_F^0.7  [V_F in cubic inches]    (9)
```

where `V_F` is the motor free volume — the void in the case not
occupied by propellant. For amateur motors with negligible head-end
ullage, `V_F ≈ V_bore + V_inter_segment_gaps`.

Rationale: this is industry standard since 1971 (NASA SP-8051 lineage).
DeMar's amateur formula `m = V·P/W` is a useful sanity check — if the
user supplies a specific pyrogen with measured impetus W, the model
should produce a consistent ignition pressure.

User overrides:
```python
result = run_from_ric(
    motor_path,
    pyrogen='bpnv',                # selects srm_1d/motors/pyrogens/bpnv.yaml
    pyrogen_mass=None,              # None → use Sutton Eq. 15-4 default
    pyrogen_throat_area=None,       # None → designed for choked flow at peak Pc
    pyrogen_volume=None,            # None → 1.5× pyrogen volume rule of thumb
    T_ignition=850,                 # K, propellant ignition criterion
    ...
)
```

## Pyrogen datasheet conventions

Pyrogen properties (a, n, ρ, T_flame, MW, γ, default geometry) live in
sibling YAML files under `srm_1d/motors/pyrogens/`. Two reference
pyrogens to ship with v0.7.0 (per DeMar 2021 measurements):

- **`bpnv.yaml`** — Boron + KNO3 + Viton (25:60:15 by mass). Industry
  benchmark. W = 5000 psi·in³/g, q ≈ 69 cal/cm²/s, T_flame ≈ 2800 K.
- **`mtv.yaml`** — Mg + Teflon + Viton. Higher heat flux (q ≈ 110
  cal/cm²/s), faster action, electrostatic-sensitivity caution.

Schema mirrors `<motor>.transport.yaml` for consistency. See
[primary_sources_summary.md](references/primary_sources_summary.md) for
the full DeMar measurements table.

## What v0.7.0 deliberately does NOT model

Per the **no-unfounded-smoothing principle**
(see [feedback_no_unfounded_smoothing](https://github.com/eJar0k/srm_1d) memory),
we exclude:

- **Multi-species gas transport**. Igniter and main propellant share γ,
  MW, T_flame in the bulk flow. Defer to v0.7.x.
- **Full particle-radiation physics**. The Phase 4 adjacent-cell model
  uses a material emissivity and local burning neighbors. It does not
  implement Salita's full Mie-correlation or Al₂O₃ size-distribution
  physics, and it does not use d'Agostino's lumped `C_hc(x/L)` tuning
  multiplier.
- **Squib stage** (electric match → pyrogen). Pyrogen ignites instantly
  at t=0 with full burning surface. Squib added in v0.7.x.
- **Impingement-region multi-cell injection**. Deferred unless tied to
  physical igniter basket/jet geometry. Do not use cell count as a
  numerical flow-field fix.
- **Two-phase flow** (Pardue & Han's Al₂O₃ condensed phase). Only
  matters for highly metallized propellants (>15% Al). Most amateur
  APCP is non-metallized or low-Al.
- **Dynamic burning rate** (Zeldovich-Novozhilov). Cavallini explicitly
  didn't implement and notes minor effect at IT and tail-off.

## Validation strategy

Primary target: **Hasegawa Motor A** (single-segment BATES, no erosive
nozzle). The motor's calibration was locked in v0.6.0
(`igniter_tau = 127 ms` FSI-proxy, MSE = 0.24 MPa²). The original
v0.7.0 success criterion was:

- **Spike overshoot drops from ~25% to <10%** — gross structural fix
- **MSE drops from 0.24 → ~0.10 MPa²** — comparable improvement
- **`igniter_tau` is removed from the parameter set** — replaced by
  pyrogen mass, throat, volume, and T_ignition (all physically grounded)
- **All tests still pass**

Current Phase 4 work added the diagnostic startup workflow, signed
open-throat boundary, pyrogen surface heating, momentum/energy ledgers,
and adjacent-cell radiation. Latest Hasegawa A smoke results show:

- `ambient_initial_gas` spreads over a finite window
  (`t10-t90 = 0.760146 s`) with startup-window peak
  `0.511 MPa @ 0.348159 s`.
- Removing either direct pyrogen surface heating or adjacent radiation
  restores the degenerate/no-spread ambient behavior.
- Baseline and `no_momentum` traces are nearly identical, so momentum is
  implemented and auditable but not the Hasegawa driver.

The historical hot-fill baseline still ignites essentially instantly,
so post-ignition burn establishment may still be needed for calibration.
Inspect the new energy/momentum audit outputs before adding that model.

Secondary target: Zerox (forward-Finocyl + aft-BATES, *with* erosive
nozzle — calibrated in v0.6.0 with `erosionCoeff` 2.34× openMotor
default). The v0.7.0 LHS will rerun on Zerox using the new pyrogen
parameters.

The user is searching for additional **high-quality static-fire data
without erosive nozzles** for further validation. Add as available.

## Roadmap (post v0.7.0)

- **v0.7.1**: post-ignition burn establishment / participation ramp
  after Goodman surface ignition, if energy-audit review still shows the
  hot-fill baseline activates too abruptly
- **v0.7.2**: radiation refinement if better material data become
  available; avoid lumped `C_hc(x/L)` as a default calibration knob
- **v0.7.3**: squib stage (electric → BPNV ramp → pyrogen → main)
- **v0.7.4**: optional multi-species via passive scalar `Y_ig[i]`
  (Cavallini-style mass-fraction-weighted thermo)
- **v0.8.0**: head-end primary motor with own nozzle (the user's
  long-term goal — Shuttle-SRB-style architecture). Pyrogen chamber
  becomes a `PrimaryMotor` with its own grain, multi-tab propellant,
  nozzle erosion. Multi-species needed at this point.
- **v0.9.0+**: full Cavallini-style mixture transport (6+ species,
  Godunov + exact Riemann mixture solver) — required for dual-pulse
  SRMs.

## File map (where the implementation lives)

| File | v0.7.0 change |
|---|---|
| New: `srm_1d/igniter_plenum.py` | `Pyrogen` dataclass + `PyrogenChamber` state + `_step_plenum_ode` Numba kernel (RK4) + `_choked_orifice_mdot` |
| New: `srm_1d/solid_thermal.py` | Goodman cubic-polynomial integral solver: `_step_goodman_ode` per-cell + `_compute_T_surf` |
| New: `srm_1d/motors/pyrogens/bpnv.yaml` | Reference Boron-KNO3-Viton pyrogen |
| New: `srm_1d/motors/pyrogens/mtv.yaml` | Reference Mg-Teflon-Viton pyrogen |
| [propellant.py](../../propellant.py) | `Pyrogen` dataclass plus `Propellant.k_solid` for Goodman ignition and `Propellant.radiation_emissivity` for adjacent-cell radiation |
| [simulation.py](../../simulation.py) | Replaced legacy igniter sim kwargs; integrates plenum step into `_run_time_loop`; adds per-cell `T_surf` and `δ`; switches ignition criterion to per-cell `T_surf > T_ignition`; records pyrogen state plus energy/momentum audit histories |
| [solver.py](../../solver.py) | `piso_step` consumes separate `mass_source`, `thermal_source`, and momentum-source arrays; nozzle boundary uses signed isentropic throat helper |
| [openmotor_adapter.py](../../openmotor_adapter.py) | Adds pyrogen loading, sibling `<motor>.pyrogen.yaml` discovery, and default `PyrogenChamber` builder |
| [tools/sensitivity.py](../../tools/sensitivity.py) | Adds segmented pressure metrics and quiet LHS progress controls for Phase 4 analysis |
| [tools/ignition_diagnostics.py](../../tools/ignition_diagnostics.py) | Startup reducer/classifier plus source, energy, momentum, ignition-time, overview, and x-t diagnostic artifacts |
| [DEVNOTES.md](../../DEVNOTES.md) | Current gotchas, validation state, and v0.7.0 breaking-change notes |

See [TASKS.md](TASKS.md) for the concrete file-level task breakdown.

## References

All cited papers are in the repo root as PDFs. Detailed extractions in
[references/](references/):

- **Ma 2019** ([extraction](references/extraction_ma2019.md)) — INVERSE
  problem (back-solves mdot_ig from measured P). NOT directly usable for
  forward prediction but provides multi-species Cp(T) machinery.
- **Salita 2001** ([extraction](references/deep_extractions.md#salita-2001))
  — physical-models framework, Baer/Ryan T_ign correlation,
  recommendation against radiation lumping.
- **Wang 2001** — CFD outlook; rationale for what we're not doing.
- **d'Agostino 2001** — quasi-1D Euler with Lenoir-Robillard erosive,
  prescribed mdot_ig, sub-1% match Ariane 4&5. Directly portable
  formulation.
- **Peretz 1973** ([extraction](references/deep_extractions.md#peretz-1973))
  — foundational 1D paper. Goodman cubic-polynomial integral method
  (Eqs. III-44 to III-50). Head-end plenum ODEs (Eqs. III-12 to III-15).
- **Pardue & Han 1992** — two-fluid extension for Shuttle SRM (Al₂O₃
  particles). Not used (single-phase sufficient for amateur APCP).
- **Cavallini 2009** — SPINBALL/SPIT Sapienza thesis. 6-species
  mixture, Godunov, Vega validation. Reference for v0.9.0+ direction.
- **Sutton 9e** — Eq. 15-4 default sizing `m = 0.12·V_F^0.7`.
- **DeMar 2021** — amateur empirical impetus measurements for 8
  pyrogen compositions; sizing formula `M = V·P/W`.

External references:
- NASA SP-8051 (Barrett 1971) — design-criteria document, available at
  [Internet Archive](https://archive.org/stream/nasa_techdoc_19710020870/19710020870_djvu.txt)
- Goodman 1958 — heat-balance integral method, ASME Trans. 80 cited in
  Peretz Ref. 49.
