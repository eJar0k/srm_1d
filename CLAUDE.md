# srm_1d — Claude Code orientation

A 1D transient finite-volume solid rocket motor internal ballistics
simulator with the Ma et al. (2020) erosive burning model. Numba-JIT
compiled time loop hits ~45-90k steps/s.

**v0.8.1 SHIPS (branch `main`, tag `v0.8.1`, public)**: the post-v0.8.0
`openmotor-frontend` work merged to `main` for public use. New since v0.8.0:
**parametric grain taper** (bore + OD/end, srm_1d transient + openMotor QS +
GUI) + longitudinal slice viewer / per-station axial viz; **core-loop
performance** — `fastmath` on the hot `@njit` kernels (~+30%, result-identical)
plus the acoustic-CFL "Lever B" design package (`docs/core_loop_opt/`,
implementation deferred); ignition-transient fixes (ProPep-anchored igniter gas
via `Pyrogen.gas_mass_fraction`; `forward_plenum` DeMar flux gated to the flame
front) + opt-in `port_mach_cap` (default off); and a CFD-light **contributor
guide** (`docs/contributor_guide/`, 6 chunks). Calibration defaults unchanged
from the v0.8.0 tag gate (roughness 32µm / kappa 0.44 / T_ignition 756K /
k_solid 0.271). See `docs/core_loop_opt/` + `docs/contributor_guide/` +
DEVNOTES "API Breaking Changes Log" (v0.8.1: OD-taper `cell_D_outer` @njit
signature change). The v0.8.0 narrative below is retained as historical.

**v0.8.0 (tag `v0.8.0`)**: openMotor frontend
integration (data-model + channels backbone + `motorlib` solver-plugin
contract; the 1-D PISO transient solver runs in openMotor's GUI alongside
its quasi-steady solver) + per-station axial viz + igniter-as-data. **Tag
gate satisfied**: v0.7.5 cross-motor re-LHS folded in and its rank-1
shared physical optimum applied to the canonical defaults — **roughness
32µm, kappa 0.44, T_ignition 756K, k_solid 0.271** (was 50µm/0.45/850/0.30;
`docs/v0_7_5/RESULT.md`). Re-validated: Hasegawa A 6.14 MPa (0.95× exp);
Zerox/Chunc residual over-prediction remains the documented high-L/D
ignition-transient QS-erosive limitation. The old nested-layout v0.7.x
line (`v0.7.0-phase4` + phase1/2/3) is retired. Full v0.8.0 narrative in
`docs/v0_8_0/` (DESIGN/TASKS/CLOSEOUT + STATION_VIZ_DESIGN). Older
per-version notes below are historical.

**v0.7.4 (tag `v0.7.4`, historical)**:
ignition-transient spike work — two opt-in (default-OFF) features plus
one always-on energy-balance bug fix. **Phase F** flame-spread front
(`Propellant.flame_front_enabled`, `flame_front_velocity≈3 m/s`):
ignition propagates as a front from the igniter; cells ahead are not
bulk-heated. **Phase Z** Zeldovich-Novozhilov dynamic burn-rate
relaxation (`zn_enabled`, `kappa_zn≈1`; `τ=κ·α_s/r²`). Each cut the
Chunc ignition spike ~35% but they DON'T stack (F+Z≈F); residual ~1.25×
overshoot remains. **Energy-balance bug fix (always on)**: the bore gas
now loses heat to unignited walls — a convective wall-loss sink AND the
previously DOUBLE-COUNTED pyrogen radiation (full enthalpy injected as
gas *and* radiated to walls) is now debited. Effect negligible (wall
heat-loss power is small during the low-Re fill) but correct.
**Audit conclusion**: the Goodman kernel, the `T_ceiling` clip, and the
wall sinks are NOT the spike cause; the gas reaching flame temperature
is physically correct; **the spike is the EROSIVE burn-rate
over-response (Root B)** — Ma's quasi-steady erosive firing instantly
off the genuine peak-G at the smallest-bore condition. Next lever:
transient/unsteady erosive closure. 291/291 pytest green. Full narrative
in `docs/v0_7_4/` (README = research synthesis, TASKS = outcome).
**Post-audit work (2026-06-01, uncommitted→committing)**: (1) **MTV
burn-rate recal** — `mtv.yaml` `a=3e-5→4.4e-5, n=0.5→0.35` (Kubota 1987,
the old seed was 6-8× too fast); old seed kept as `mtv_fast.yaml` for A/B.
Cut Chunc head_basket spike 2.02×→1.55× at Sutton 0.9 g. (2) **Realistic
basket geometry** — `PyrogenChamber.basket_fill_fraction=0.5` +
`pellet_packing_fraction=0.60` replace the solid-puck `L_cart`; distribution
proven NEUTRAL on spike (correctness only); particle diameter is the unifying
specific-surface knob. (3) **Finding**: at the as-fired 6 g charge (Sutton
gives 0.9 g — Eq.15-4 is a central AP/Al ±2× fit) the spike re-inflates to
3.19× via ignition-simultaneity; bed flame-spread is fast (τ_bed≈5-10 ms) so
an igniter ramp is 2nd-order. **INVESTIGATION CLOSED (2026-06-01, user decision
= document as known limitation).** Exhaustive differential diagnosis proved the
spike is NOT: igniter mass/topology/IC, our ignition kernel (Hasegawa A/B/C —
Ma's exact validation motors — are FAITHFUL: 1.15×/0.68×/no-spike), burn-rate
magnitude (fixed), ignition sequencing (flame-front can't suppress at physical
speeds), a derivable erosive lag (Beddini gives NO parameter-free τ → would need
tuning), or numerical resolution (every-step burn+geom + cfl 0.15 identical). It
IS the genuine, faithful **Ma quasi-steady erosive response to transient mass-
flux during fast ignition of a HIGH-L/D motor** — a regime Ma EXCLUDED from every
error figure and never benchmarked. Under the no-tuning dogma there is NO
literature closure; the sim is faithful to Ma, not broken. Full record +
elimination table: `docs/v0_7_4/IGNITION_SPIKE_CLOSEOUT.md`. (`bpnv` same-check
moot — Hasegawa is fine.)

**v0.7.3-phaseB ships (branch `v0.7.0-phase4`, tag `v0.7.3-phaseB`)**:
heat-flux completeness for uncontained ignition. Four fixes close
the Phase A.3 gap (uncontained topologies stalling at atmospheric P):
**B.0** IC fix (T_initial_gas = T_ambient instead of T_flame —
realistic physics, side effect of larger ignition spikes on
calibrated motors); **B.2** radiation_emitter gating extension
(pyrogen-hot cells now emit, no-op when emissivity=0);
**B.3** pyrogen form archetypes (powder/pellets/chunks with ×20/×5/×1
A_burn multipliers; pellets is the new default); **B.4** unified
pyrogen-to-surface heat delivery enum (`demar`/`radiation`/`none`).
**Empirical finding**: all three B.4 modes give identical P_peak on
Super Loki head_basket — the load-bearing fixes are B.0 + B.3, B.4
is diagnostic refinement. Hasegawa A aft_basket stalls under all
modes because the cartridge is too close to the nozzle (deferred
`aft_fore_firing` topology needed). **Provenance correction**:
the Super Loki "experimental" overlay was actually mis-labeled
Chunc data; removed. 272/272 pytest green (test windows widened
to ±150% on Hasegawa A baseline gates pending v0.7.4 Phase C
re-LHS). See `docs/v0_7_3/` for the full narrative.

**v0.7.3-phaseA baseline** (carried forward): uncontained-pyrogen
topology architecture (`head_basket` + `aft_basket`) wired into
the time loop via a shared `PyrogenChamber.injection_topology`
field. Each pyrogen pellet burns at its host cell's LOCAL bore
pressure (no plenum, no orifice). Diagnostic visualization helpers
(`plot_flow_snapshots`, `plot_field_heatmap`, sign-banded
`u_cell` panel in `plot_flow_snapshot`).

**v0.7.2-phaseA baseline** (carried forward): pyrogen axial
distribution (Phase A) shipped as a real Zerox win (P_peak
10.20→9.69 MPa, t_peak 0.035s→0.27s); Phase B (h_c augmentation)
shipped but DISABLED by default after both v1 (cumulative-G) and
v2 (flame-front) formulations amplified rather than smoothed the
spike. Structural ignition-kernel artifact persists for Hasegawa
A / BALLSstick / Chunc — PISO's local-Re tracking already
captures upstream-mass-flux contributions, so the Kashiwagi/Han
augmentation double-counts. See `docs/v0_7_2/`.

**v0.7.1.1 (historical)**: the N-species bore-gas refactor (SPINBALL-style
"infinite-gases mixture", per-cell γ/R/Cp) — still the current architecture.
NOTE its two now-SUPERSEDED claims: the effective-RPA-default and the old
`hasegawa_motor_a.py` knobs (37.1µm/850) are gone — transport is now embedded
per-propellant-tab (frozen won the v0.7.3.2 re-LHS) and the canonical defaults
are the v0.7.5 optimum (roughness 32µm / kappa 0.44 / T_ignition 756K /
k_solid 0.271).

This file is loaded on every session — keep it tight. Pointers to
deeper docs at the bottom.

## Quick start

**v0.8.0 flat layout (run everything from the repo root).** The package
is `srm_1d/` at the repo root; `tests/`, `examples/`, `tools` (now
`srm_1d/tools/`), `docs/`, `motors/`, `static_fire_data/` are siblings.
`pip install -e .` is set up in the pyenv 3.10.5 env (a root `conftest.py`
also puts the repo root on `sys.path`, so pytest works without the install).

```bash
# Tests (pyenv 3.10.5 -- has numba, pytest, scikit-fmm installed)
"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" -m pytest tests/

# Hasegawa A example (loads motors/hasegawa_a.ric) -- run as a module from repo root
"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" -m examples.hasegawa_motor_a
```

System Python (3.14 on PATH) does NOT have these deps; always use the
pyenv 3.10.5 path explicitly.

## Module map (one-liners)

```
srm_1d/                  <- importable package (core, ships); everything below
│                            tools/ is repo-root sibling, dev-only (not shipped)
├── solver.py            PISO + TDMA + adaptive CFL (pure numerics, no project deps)
├── burn_rate.py         Ma 2020: Haaland → Gnielinski → bisection. Multi-tab Saint-Robert lookup.
├── propellant.py        Propellant + tabs (PropellantTab) + GasProperties + GasSpecies + Pyrogen + thermo utilities
├── grain_geometry.py    GrainSegment / MotorGeometry / build_snapped_geometry; per-cell regress[i]
├── nozzle.py            openMotor-aligned Nozzle: thrust, Isp, CF, throat erosion. Adjusted-CF formula.
├── fmm_grain.py         Bridge to local openMotor checkout; FmmTable extraction + Numba lookup
├── igniter_plenum.py    Pyrogen chamber, choked/subsonic venting, Sutton sizing defaults
├── solid_thermal.py     Goodman integral solid-heating ignition subsolver
├── simulation.py        run_simulation wrapper + @njit _run_time_loop (pyrogen + Goodman); v0.7.1: _advect_species, _compute_mixture_cell, _refresh_mixture_arrays, _compute_T_ceiling_arr
├── plotting.py          matplotlib plots (pressure, thrust, flow snapshots, summary)
├── openmotor_adapter.py .ric reader, transport YAML loader, convert_propellant/_geometry/_nozzle, CSV export
├── pyrogens/            Runtime pyrogen material library (bpnv/mtv/...; package-data, ships)
├── srm1d_plugin.py      openMotor SolverPlugin adapter (registry boundary)
└── tools/sensitivity.py Latin Hypercube sweeps + tools/ignition_diagnostics.py (ship as srm_1d.tools)

repo root (siblings, dev-only — NOT in the package/wheel):
├── motors/              Motor data: <motor>.ric (transport embedded per-tab; sidecars retired)
├── examples/            hasegawa_motor_a, bates_4seg, hasegawa_a_lhs, machbusterNew, ... (run: python -m examples.X)
├── tests/               pytest suite (run: python -m pytest tests/)
├── docs/, static_fire_data/, pyproject.toml, conftest.py
```

## Dev workflow

- **Versioning is git tags**, not folder names. Latest tag: `v0.8.1`
  (2026-06-19, on `main`, public). Active dev on branch `openmotor-frontend`
  (== main at release). Bump on hard API breaks; document each break in
  DEVNOTES "API Breaking Changes Log."
- **Hard API breaks are fine** — refactor cleanly, no backward-compat
  shims. (See `feedback_api_breaks` memory.)
- **Defer to openMotor's architecture** when adding data structures
  (field names, semantics). UNITS are the documented exception:
  srm_1d keeps human-readable engineering units (μm/(s·MPa) for
  erosion_coeff, etc.); adapter converts at the boundary. (See
  `feedback_openmotor_alignment` memory.)
- **Named motors live as data**, not Python factories: add a
  `<motor>.ric` to `motors/` (gas transport is embedded per-propellant-tab
  in the .ric — the `.transport.yaml` sidecars are retired), load via
  `run_from_ric`. Parametric geometry uses `build_snapped_geometry` directly.
- **Repo: github.com/eJar0k/srm_1d** (private).
- **Pytest before commit**; clear `__pycache__/` and `.nbi`/`.nbc`
  after edits to @njit functions (Numba cache persistence is the #1
  source of "the fix didn't work" bugs — see DEVNOTES gotchas).

## Critical gotchas (concentrated, full list in DEVNOTES)

1. **Numba cache** — delete `srm_1d/__pycache__/` after any @njit
   edit, or you'll run the old compiled code.
2. **Mass conservation** — the burnout ramp `f_active` MUST multiply
   both `C_burn` AND the regression rate. Either alone causes 3–40%
   mass error.
3. **End-face injection (v0.6.0+)** uses a partition-of-unity hat
   function: each face's mass splits over 2 adjacent cells with weights
   summing to 1.0. Coupled to snapping (which puts faces on cell edges).
   Gated by `tests/test_endface_conservation.py`.
4. **Igniter API hard-broke in v0.7.0** -- the exponential knobs
   (`igniter_mass`, `igniter_tau`, `ignition_ramp_tau`, `P_ignition`)
   are gone. Use `pyrogen_chamber` directly or `run_from_ric(...,
   pyrogen='bpnv')`. The plenum injects mass + enthalpy + axial
   momentum (via `_orifice_exit_velocity`) into cell 0 / face 1.
   The `igniter_axial_momentum_fraction` kwarg scales the momentum
   contribution (default 1.0 = pure axial head-end jet).
   **v0.7.3 Phase A** adds two new optional fields:
   `PyrogenChamber.injection_topology` (default `'forward_plenum'`,
   alternatives `'head_basket'` | `'aft_basket'` for uncontained
   pellets that burn at LOCAL bore P with no plenum / no orifice /
   no DeMar surface heat flux), and `cartridge_length_m` (default
   `-1.0` sentinel = derive from pyrogen mass via
   `L_cart = m_pyrogen / (rho_p * A_port_avg)`). Both `run_from_ric`
   and `build_pyrogen_chamber` accept these as kwargs.
5. **Frozen vs effective gas transport — v0.7.1.1 ships EFFECTIVE as
   default for ALL fired motors**. Hasegawa A (v0.7.1), then Zerox /
   Chunc (machbusterNew) / BALLSstick (v0.7.1.1 patch) all flipped to
   their RPA-effective pair; frozen siblings preserved at
   `<motor>.frozen.transport.yaml` for diagnostic reference (load via
   `transport_path=...` explicitly). The v0.7.1 Phase 5 effective LHS
   for Hasegawa A landed k_solid at the literature center 0.331 W/(m·K)
   (vs frozen pegging lower bound 0.20 — free-parameter compensation
   for under-heat-transfer through the gas film). `k_thermal` and
   `mu_gas` remain scalar (per-cell array deferred to v0.7.2 if benefit
   shown). **Unfired motors (ChaseRed BATES, L3035, ivanO25k) still
   ship frozen** — no validation signal to motivate the switch.

   **IMPORTANT post-tag finding (v0.7.1.1 cross-motor cleanup
   2026-05-23)**: at default knobs (k_solid=0.4, roughness=35µm,
   kappa=0.45, T_ign=900K, Sutton pyrogen), effective transport
   AMPLIFIES the ignition spike for every fired motor by +30-55% vs
   frozen at the same knobs (Hasegawa A 5.84→8.27, Zerox 7.85→10.20,
   BALLSstick 9.33→14.48, Chunc 13.14→20.27 MPa). The pre-existing
   Zerox YAML comment ("effective amplifies the spike") is confirmed
   universally. **`hasegawa_motor_a.py` with the v0.7.1 effective
   default + v0.7.0 knobs now over-predicts the ignition spike by
   ~31%** (P_peak 8.5 MPa @ t≈0.05 s vs experimental 6.5 MPa @ t=1.1 s),
   while matching the plateau + erosive peak shape better than the
   v0.7.0 frozen baseline. This is the structural ignition-kernel
   artifact manifesting in the canonized example — it is NOT a
   regression, but anyone running `hasegawa_motor_a.py` post-v0.7.1
   will see this spike overshoot. v0.7.2 structural work targets it.
6. **v0.7.1 thermal_source units changed to W/m (Phase 3 step 1)**.
   Previously kg·K/(s·m) — multiply by Cp_gas to convert legacy
   external builds. Each source site multiplies its mdot·T contribution
   by Cp_gas; the PISO energy equation now treats the input as direct
   enthalpy injection per unit length. `_pyrogen_surface_thermal_sink`
   and `_thermal_source_power` signatures changed accordingly.
7. **v0.7.1 PISO + post-PISO take per-cell arrays**.
   `_piso_step_with_energy_diagnostics` and `piso_step` take
   `gamma_arr / R_arr / Cp_arr / T_ceiling_arr` instead of scalar gas
   thermo. Energy advection is sensible-enthalpy (Cp·T) — face fluxes
   carry upwind Cp·T to conserve energy across cells with different Cp.
   Nozzle BC uses cell-N-1 mixture. T_ceiling is strict DESIGN §5
   (per-cell Y > 0.05 filter) with an IC guard at T_initial_gas · 1.01
   to preserve the v0.7.0 IC (T = T_flame_prop while Y = 100% ambient).
   See `_compute_T_ceiling_arr` docstring.
8. **v0.7.1 per-species Cp at source sites (Phase 3.5)**. Each
   combustion source multiplies its `mdot · T_source` by its OWN
   species's Cp: propellant grain → `Cp_propellant`; pyrogen plenum →
   `Cp_pyrogen`. The `_pyrogen_surface_heat_power` sensible-power cap
   uses Cp_pyrogen. This is ~33% lower than the prior scalar Cp_gas
   for BPNV-class pyrogens, which suppresses the v0.7.0 Hasegawa A
   ignition spike (P_peak shifts from t=0.041s to t=3.36s). Expected
   regime change; Phase 5 LHS will recalibrate.

## External deps

- **openMotor checkout** sits as sibling: `Erosive Burning Solver/openMotor/openMotor/`.
  `srm_1d.fmm_grain` walks upward to find it; override with
  `SRM1D_OPENMOTOR_PATH` env var if needed.
- pip: `numba`, `numpy`, `scipy`, `pyyaml`, `matplotlib`, `pytest`,
  `scikit-fmm` (for FMM grains), `scikit-image` (for openMotor's
  Custom grain). The Cython `mathlib._find_perimeter_cy` is replaced
  at import-time with a Numba marching-squares shim, so MSVC build
  tools aren't needed.

## Where to look for more

**Start here (current):**
- `docs/contributor_guide/` -- CFD-light onboarding walkthrough of the sim
  core + newcomer code map + opt-in add-ins catalog (README + 01_SIM_CORE
  through 05_IO_AND_OPENMOTOR). Best first read for a new contributor.
- `srm_1d/README.md` -- public API, motor designation, validated parameters.
- `srm_1d/ARCHITECTURE.md` -- function-level map of every module.
- `srm_1d/DEVNOTES.md` -- gotchas, calibration state, API breaking-change
  log per minor version, performance profile.
- `srm_1d/tools/README.md` -- the LHS-sweep + ignition-diagnostics tooling.
- `docs/core_loop_opt/` -- core-loop profiling + the deferred acoustic-CFL
  "Lever B" design (IMEX) + the all-speed PISO literature review (with DOIs).

**Design packages by version (historical — read the one matching the area
you're changing):**
- `docs/v0_8_0/` -- openMotor frontend integration + per-station axial-viz
  design (DESIGN/TASKS/CLOSEOUT + STATION_VIZ_DESIGN + UPSTREAM_TAPER_PR_SCOPE).
- `docs/v0_7_4/` -- ignition-spike investigation + **CLOSE-OUT**: the high-L/D
  ignition over-spike is a documented, accepted Ma quasi-steady-erosive
  limitation — **do not re-open** (see IGNITION_SPIKE_CLOSEOUT / _REOPENED and
  the `project_ignition_model_audit` memory).
- `docs/v0_7_0..v0_7_3/` + `docs/post_v0_7_0/` -- pyrogen plenum + Goodman
  ignition (v0.7.0), N-species mixture (v0.7.1), ignition-model candidates +
  Phase A/B findings (v0.7.2), uncontained topologies (v0.7.3), SPINBALL
  research. Each has DESIGN/TASKS/references.

## Open roadmap (priority order)

1. **Code review → proper upstream PRs (near-term north star).** Chunk up and
   review ALL the code over several sessions so areilley (upstream openMotor
   maintainer) + the user can build proper upstream PRs. The CFD-light
   contributor guide (`docs/contributor_guide/`) is the first deliverable
   (DONE). This **pauses solver-core churn** (the core-loop items below) until
   the review cadence is set. See `project_code_review_for_prs` memory.
2. **Upstream openMotor taper PRs.** The generic tapering is prepped as three
   fork branches (`taper/core`, `taper/gui`, `taper/mainwindow-layout`) off
   `upstream/staging`, gated (oM unit + offscreen GUI smoke) — **local-only,
   HELD** pending the areilley discussion + a `/code-review ultra` pass on each
   (user-triggered; the agent can't launch it). Scope:
   `docs/v0_8_0/UPSTREAM_TAPER_PR_SCOPE.md`.
3. **Core-loop performance (deferred).** `fastmath` shipped (+~30%,
   result-identical). The big structural win — **Lever B: escape the acoustic
   CFL via IMEX** (~2×, parameter-free) — is DEFERRED: research-grade surgery on
   the co-reviewed core; the plan is the user reads the literature
   (Klein 1995 / Degond–Tang 2011) then builds an isolated git-worktree
   prototype BEFORE touching the real solver. `A2` burn-rate bisection
   warm-start also deferred. Full record: `docs/core_loop_opt/`.
4. **Ignition over-spike: CLOSED — do not re-open.** The high-L/D ignition
   pressure over-prediction is the faithful Ma quasi-steady-erosive response in
   a regime Ma never benchmarked — a documented, accepted limitation with no
   non-tuned closure (`docs/v0_7_4/`, `project_ignition_model_audit` memory).
   **Marked for later:** extend the ProPep igniter-gas fix (`bfc2f3f`) to
   MTV/thermite via CEA.
5. **Smaller open items.** RodTube grain support (small `PerforatedGrain`
   extension alongside `FmmGrain` in `from_openmotor`); the OD-tapered
   end-cell slice-render sliver at burnout; the un-instrumented per-station
   FMM-resolve setup phase (no progress bar); BALLSStick CAD QS-vs-transient
   validation.
