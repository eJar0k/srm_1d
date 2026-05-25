# srm_1d — Claude Code orientation

A 1D transient finite-volume solid rocket motor internal ballistics
simulator with the Ma et al. (2020) erosive burning model. Numba-JIT
compiled time loop hits ~45-90k steps/s.

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
re-LHS). See `srm_1d/docs/v0_7_3/` for the full narrative.

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
augmentation double-counts. See `srm_1d/docs/v0_7_2/`.

**v0.7.1.1 baseline** (carried forward): N-species bore-gas refactor
(SPINBALL-style "infinite-gases mixture") + EFFECTIVE RPA transport
default for ALL 4 fired motors (Hasegawa A, Zerox, Chunc, BALLSstick).
`hasegawa_motor_a.py` retains v0.7.0 knobs (roughness=37.1µm,
kappa=0.45, T_ign=850, k_solid=0.3 default).

This file is loaded on every session — keep it tight. Pointers to
deeper docs at the bottom.

## Quick start

```bash
# Tests (pyenv 3.10.5 -- has numba, pytest, scikit-fmm installed; 206 tests)
"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" -m pytest srm_1d/tests/

# Hasegawa A example (loads srm_1d/motors/hasegawa_a.ric)
"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" -m srm_1d.examples.hasegawa_motor_a
```

System Python (3.14 on PATH) does NOT have these deps; always use the
pyenv 3.10.5 path explicitly.

## Module map (one-liners)

```
srm_1d/
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
├── motors/              Canonical motor data: <motor>.ric + <motor>.transport.yaml pairs
├── tools/sensitivity.py Latin Hypercube parameter sweeps with parallel execution
├── examples/            hasegawa_motor_a, bates_4seg, hasegawa_a_lhs, Zerox_test, ZeroxOptimizer
└── tests/               15 files, 206 tests
```

## Dev workflow

- **Versioning is git tags**, not folder names. Latest tag: `v0.7.0`
  (2026-05-21), still on branch `v0.7.0-phase4` (no merge to main yet).
  Bump on hard API breaks; document each break in DEVNOTES "API
  Breaking Changes Log."
- **Hard API breaks are fine** — refactor cleanly, no backward-compat
  shims. (See `feedback_api_breaks` memory.)
- **Defer to openMotor's architecture** when adding data structures
  (field names, semantics). UNITS are the documented exception:
  srm_1d keeps human-readable engineering units (μm/(s·MPa) for
  erosion_coeff, etc.); adapter converts at the boundary. (See
  `feedback_openmotor_alignment` memory.)
- **Named motors live as data**, not Python factories: add a
  `<motor>.ric` + sibling `<motor>.transport.yaml` to `srm_1d/motors/`,
  load via `run_from_ric`. Parametric geometry uses
  `build_snapped_geometry` directly.
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

- `srm_1d/README.md` -- public API, motor designation, validated parameters
- `srm_1d/ARCHITECTURE.md` -- function-level map of every module
- `srm_1d/DEVNOTES.md` -- full gotchas, calibration state, complete
  API breaking-change log per minor version, performance profile
- `srm_1d/docs/v0_7_0/` -- v0.7.0 hot-gas plenum design package:
  `DESIGN.md` (implemented architecture), `TASKS.md` (phase status),
  `references/` (extracted papers + Goodman integral derivation +
  Sutton/DeMar summaries).
- `srm_1d/docs/v0_7_1/` -- v0.7.1 N-species design package:
  `DESIGN.md` (mixture architecture, Phase plan, ambient species
  decision), `TASKS.md` (Phases 1+2+3+3.5+4+5 complete, tagged).
- `srm_1d/docs/v0_7_2/` -- v0.7.2 ignition-model rework design
  package: `README.md` (problem statement + decision criteria),
  `candidates/01..04_*.md` (4 original candidate design docs: Z-N
  dynamic burn rate, spatial ignition-front coupling, pyrogen axial
  distribution, submerged pyrogen modes), `references/01..04_*.md`
  (extended literature digests), `TASKS.md` (Phase A complete +
  Phase B negative findings + Phase C close-out — tagged
  `v0.7.2-phaseA`). Candidate 3 (pyrogen distribution) shipped as
  Phase A. Candidate 2 (spatial coupling via h_c augmentation)
  attempted twice (v1 cumulative-G, v2 flame-front gating) and
  shipped DISABLED by default after both amplified rather than
  smoothed the spike.
- `srm_1d/docs/v0_7_3/` -- v0.7.3 uncontained-pyrogen topology
  package: `README.md` (problem statement + Phase A architecture +
  Phase B candidate space), `TASKS.md` (Phase A.1/A.1.1/A.2/A.3
  complete — tagged `v0.7.3-phaseA`). Validation findings: ISP
  Super Loki head_basket (P_peak 0.12 MPa vs ~8.8 MPa experimental)
  and Hasegawa A aft_basket (P_peak 0.10 MPa) both stall at
  atmospheric P, exposing the ignition-initiation gap. Cross-
  version architectural anchor is `docs/v0_7_2/candidates_post_phaseA.md`.
- `srm_1d/docs/post_v0_7_0/references/` -- SPINBALL research that
  motivated v0.7.1: Cavallini 2009 + DiGiacinto 2008 extractions plus
  the `spinball_walkthrough.md` decision document (recommends Z-N
  dynamic burn rate as the spike-taildown candidate; N-species is the
  prerequisite infrastructure being delivered first).
- `gemini summary.md` (repo root) -- historical record of the v0.6.0
  development cycle that originated build_snapped_geometry, the new
  end-face kernel, and the exponential-decay igniter
- `generic agent instructions.md` -- short current handoff for future
  coding agents. Older external agent-memory references are historical;
  this repo's committed Markdown is the source of truth.

## Open roadmap (priority order)

1. **v0.7.3 Phase B — uncontained-pyrogen ignition initiation** --
   v0.7.3-phaseA shipped the uncontained-pyrogen topology
   architecture (`head_basket` + `aft_basket`) but Phase A.3
   validation surfaced a structural gap: pellets at atmospheric
   bore P burn slowly enough that the bore never pressurizes
   without an external thermal kick. forward_plenum hides this
   via the choked-orifice startup transient; uncontained has no
   equivalent. The natural Phase B is to add an ignition-
   initiation pathway. Recommended options pending user
   decision (smallest → largest scope):
   - **Initial thermal pulse** (~150 LOC) — kick `T_surf` and/or
     local bore T at t=0 in cartridge cells. Models the
     e-match heat dump.
   - **Per-pellet surface heat flux** (~200 LOC) — re-enable
     DeMar-style surface flux PER CARTRIDGE CELL (not cell 0
     only). Most direct mapping from forward_plenum.
   - **Coupled e-match dataclass** (~400 LOC) — `Igniter` with
     explicit `t_kick`, `Q_kick`, `tau_kick`. Opens door to
     candidate 6 (plenum-as-option refactor).

   Once Phase B lands, re-run **both** validation examples to test
   (a) ISP Super Loki head_basket fit vs ~8.8 MPa experimental,
   and (b) Hasegawa A aft_basket diagnostic — does the
   simultaneous-ignition artifact persist under reversed mass-
   injection topology?

   Broader v0.7.3+ candidate space still on the table (independent
   of Phase B):
   - **Z-N dynamic burn rate** — burn-rate-ramp lag, stacks
     cleanly with any topology.
   - **Per-cell coupling alternatives** — reverse-polarity Phase B
     damping, solid-phase axial conduction.
   - **Pardue 1992 Al2O3 condensation** — secondary spike-taildown
     candidate.
   - **Plenum-as-option refactor** (candidate 6) — best done AFTER
     ignition-initiation work lands so the abstraction is informed
     by real requirements.

   Cross-version anchor: `srm_1d/docs/v0_7_2/candidates_post_phaseA.md`.
   v0.7.3 Phase B narrative: `srm_1d/docs/v0_7_3/TASKS.md`.
2. **Cross-motor effective-transport recalibration** -- v0.7.1 only
   flipped Hasegawa A to effective; Zerox/BALLSstick/machbusterNew/
   ChaseRed/L3035/ivanO25k YAMLs are still frozen. After v0.7.2's
   structural ignition fix lands, re-run cross-motor LHS using
   effective + the new kernel.
3. **ε = 0.05 single-cell ignition spike** (radiation-collapse residual)
   -- the only remaining outlier in the 27-variant radiation matrix
   trips on the last-grain-cell ignition transition. Would require
   source sub-stepping (split-operator within one PISO step).
4. **RodTube grain support** -- small extension (PerforatedGrain in
   addition to FmmGrain in `from_openmotor`).
5. **Al2O3 two-phase thermal lag** -- Pardue 1992 form, the secondary
   spike-taildown candidate from the SPINBALL walkthrough. Higher
   implementation cost than Z-N (full re-cal needed). v0.7.3+ if Z-N
   alone doesn't close the residual.
