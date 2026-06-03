# srm_1d v0.7.3 — Uncontained Pyrogen Topologies

**Status (2026-05-24)**: v0.7.3 Phase A landed as the
**uncontained-pyrogen architecture** (`head_basket` + `aft_basket`
topologies on a shared kernel). Tag `v0.7.3-phaseA` ships as an
intermediate milestone — the topology infrastructure validates
in tests but the diagnostic-against-fired-motors validation
surfaced a **structural gap**: uncontained models lack an
ignition-initiation pathway and stall at atmospheric pressure
without a thermal kick. See [TASKS.md](TASKS.md) for the full
narrative.

## Context

v0.7.2-phaseA shipped Phase A (pyrogen axial distribution) as a
real Zerox win but Phase B (h_c augmentation via cumulative-G then
flame-front gating) both empirically amplified the simultaneous-
ignition spike. The post-Phase-B candidate breakdown
([../v0_7_2/candidates_post_phaseA.md](../v0_7_2/candidates_post_phaseA.md))
identified the **submerged pyrogen 4b — aft-inserted impinging
cartridge** as the user-flagged cleanest diagnostic test: reverse
the mass-injection topology and see whether the simultaneous-
ignition artifact persists.

v0.7.3 unifies candidates 4a (head-end basket) and 4b
(aft-inserted) under a single **uncontained-pyrogen** model
(distinct from the v0.7.0+ **forward_plenum** with its choked
orifice). The naming pivoted during A.1.1 from "submerged" to
"uncontained" to clarify the physics: pellets sit physically
inside the bore with no plenum wall, orifice, or pressure
separation — each pellet burns at its host cell's local bore
pressure.

## What shipped (Phase A.1 + A.1.1 + A.2)

- **`_compute_uniform_band_weights`** Numba kernel — mass-
  conservative top-hat axial weights over arbitrary
  `[i_start, i_end]` cell range.
- **`_compute_uncontained_pyrogen_mdot`** Numba kernel — per-cell
  pyrogen mdot from local bore P, with mass-conservation cap.
- **`PyrogenChamber.injection_topology`** field with three
  values: `'forward_plenum' | 'head_basket' | 'aft_basket'`.
  Default preserves v0.7.0+ behavior byte-for-byte.
- **`PyrogenChamber.cartridge_length_m`** field with `-1.0`
  sentinel that derives cartridge length from pyrogen mass.
- **`_run_time_loop` topology branch** — uncontained topologies
  use per-cell mass / species / enthalpy delivery, skip momentum
  injection, skip DeMar surface heat flux, and use volume-avg
  bore pressure over cartridge cells as the P_ig diagnostic.
- **574 lines of new tests** covering the kernels and integration
  paths (`test_uniform_band_weights.py`,
  `test_uncontained_pyrogen.py`, `test_submerged_topology.py`).

## What shipped (Phase A.3 — this session)

- **Diagnostic-visualization helpers** in `plotting.py`:
  - `plot_flow_snapshot` now renders a 3x2 grid that adds
    velocity (with sign-banding for at-a-glance reverse-flow
    diagnosis) and gas temperature.
  - `plot_flow_snapshots` — multi-time subplot grid (rows = time,
    cols = field).
  - `plot_field_heatmap` — `pcolormesh(x, t, field)` for
    visualizing back→front cascades and reverse-flow bands.
- **`ISP_SUPER_LOKI_EXPERIMENTAL`** dataset added to `plotting.py`
  with proper labeling (was previously mis-labeled commented-out
  block in the example script).
- **Two validation example scripts** wired to the new topologies:
  - `examples/ISP_Super_Loki.py` — head_basket fit
  - `examples/hasegawa_motor_a_aft_basket.py` — aft_basket
    diagnostic
- **API extension**: `run_from_ric` and `build_pyrogen_chamber`
  accept `injection_topology=` and `cartridge_length_m=` kwargs.

## What didn't validate

- **ISP Super Loki head_basket**: P_peak=0.12 MPa vs experimental
  ~8.8 MPa. The MTV pyrogen burns to completion in 298 ms at
  atmospheric pressure (Saint-Robert at P_atm gives tiny r_b)
  without ever pressurizing the bore enough to ignite the main
  grain. Mass-flux pathway alone is too slow.
- **Hasegawa A aft_basket diagnostic**: P_peak=0.10 MPa, never
  reaches ignition cascade. Same failure mode — pyrogen burns
  out at atmospheric P without lighting the main grain.

**The diagnostic question** — "does the simultaneous-ignition
artifact persist under reversed topology?" — is therefore
**inconclusive** because the aft_basket run never reaches an
ignition state to compare against forward_plenum. The structural
finding is itself informative: the uncontained model correctly
captures "pellets at P_atm burn slowly" but exposes a gap that
forward_plenum hides via its choked-orifice startup transient.

## v0.7.3 Phase B candidate space

The natural Phase B is to **add an ignition-initiation pathway**
for uncontained topologies. Options:

1. **Initial thermal pulse** — kick the cartridge cells'
   `T_surf` and/or local bore T at t=0 (small Gaussian or
   step pulse). Models the e-match / squib heat dump that
   physically initiates pyrogen pellets.
2. **Initial pressure pulse** — small initial mass injection
   (separately from pyrogen pellet burn) at t=0 that primes the
   bore to a low ignition pressure. Models the same physical
   mechanism via a different mathematical lever.
3. **Per-pellet surface heat flux** — a `head_basket_heat_flux`
   parameter that re-enables the DeMar-style surface heat flux
   path PER CARTRIDGE CELL (not at cell 0 only as in
   forward_plenum). Most direct mapping from the forward_plenum
   ignition behavior to the uncontained world.
4. **Coupled e-match model** — add a small `Igniter` dataclass
   for the electrical ignition transient with a user-specified
   `t_kick`, `Q_kick`, `tau_kick`. Most physical but largest
   scope.

Once Phase B lands, the v0.7.3 aft_basket diagnostic can finally
be answered.

## Cross-version pointers

- v0.7.2-phaseA close-out: [../v0_7_2/TASKS.md](../v0_7_2/TASKS.md)
- v0.7.3+ candidate analysis (cross-version architectural anchor):
  [../v0_7_2/candidates_post_phaseA.md](../v0_7_2/candidates_post_phaseA.md)
- Architectural source-of-truth for uncontained vs plenum split:
  [PyrogenChamber docstring](../../igniter_plenum.py) L52-L120
