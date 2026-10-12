# Pyrogen realism — IG1 diagnosis (agent-executable)

**User priority (2026-10-11): this comes BEFORE the LHS re-sweep and the
Mukunda-Paul experiment.** Two symptoms, reported from openMotor GUI runs of
Chunc (`motors/machbusterNew.ric`, BKNO3 = library `bpnv`, head basket):

1. The sim ignites motors with ~10× less pyrogen than reality. The real Chunc
   was lit with ~0.01 lb (≈4.5 g) BKNO3 in a head basket, and its trace has
   NO startup spike (flat ~8.8 MPa; `static_fire_data/thomas_chunc_firing.csv`).
2. Realistic charges spike the sim: 0.001 lb → ~1440 psi, 0.01 lb → >3000 psi.

**IG1 is DIAGNOSIS ONLY** — find and rank the causes with evidence; propose
physically grounded fix options; change no defaults. The user decides fixes.

## Reviewer's reproduction (main `daf367c`, current defaults, committed Chunc
`.ric`, roughness 32 µm, kappa 0.44, T_ign 756 K, library BPNV, t_max 0.1 s)

| case | peak | first ignition |
|---|---|---|
| head_basket 0.0001 lb (0.045 g) | does not ignite (0.11 MPa) | — |
| head_basket 0.001 lb (0.45 g) | 9.91 MPa (1437 psi) @ 63.7 ms | 28.5 ms |
| head_basket 0.01 lb (4.5 g) | 21.66 MPa (3141 psi) @ 9.6 ms | 0.68 ms |
| same, propellant cannot ignite (T_ign 1e5) | 2.03 MPa (igniter gas alone) | — |
| same, erosion off (`saint_robert`) | 7.27 MPa (no spike) | 0.68 ms |
| forward_plenum 0.01 lb | 15.97 MPa @ 13.7 ms | 0.64 ms |

Reading: the 4.5 g charge lights the grain within ~1 ms; the pyrogen gas
itself is only ~2 MPa; the spike is the erosive response to that
near-simultaneous early ignition. Two families of cause are live:
**propellant ignition too easy** and/or **pyrogen output too fast**.

## Prior evidence to read first (do not re-derive)

- `docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` §9.2 (Keller-Baer-Ryan 1966: the
  Gnielinski h_c used for ignition is ~7× Keller's measured h at matched G;
  q″ ~31 MW/m² vs Keller ≤ 6.7; a transient h_c knock-down only POSTPONED the
  spike in the old tree), §8 thread A (ms snap-on; endothermic / integrated-
  energy ignition criterion = "candidate 2"), §12–13 (T8 / P5 / T11).
- Memory-level facts recorded in the repo: DeMar 69.4 cal/cm²/s is a
  TIME-AVERAGED igniter flux; igniter pellet-bed flame spread τ_bed ≈ 5–10 ms
  (canonical igniter models assume simultaneous pellet ignition); BPNV gas
  properties are ProPep-anchored (`bfc2f3f`).
- `docs/v0_7_3/` (head_basket / uncontained topology design),
  `srm_1d/igniter_plenum.py`, `srm_1d/solid_thermal.py` (Goodman ignition),
  `_goodman_ignition_sources_and_mass` in `srm_1d/simulation.py`.
- `docs/references/keller1966.pdf`, `cain2006.pdf`, and the Ma 2020 igniter
  paper (`docs/references/International Journal of Aerospace Engineering -
  2020 - Ma - A Model for Igniter Mass Flow Rate History Evaluation for.pdf`).

## Hypotheses to discriminate

- **H1 heat transfer to the propellant is too high** during the igniter phase
  (convective h_c from the developed-pipe Gnielinski correlation applied to a
  transient fill; plus DeMar surface flux / pyrogen radiation, for the
  topologies that use them).
- **H2 the ignition criterion is too permissive** (surface T ≥ 756 K with no
  finite ignition energy / no endothermic delay).
- **H3 pyrogen gas production is too fast** (all pellet surface burning from
  t = 0, i.e. no bed flame spread; r = a·Pⁿ at LOCAL bore pressure in the
  head basket → P↑ → r↑ feedback; BPNV a/n and pellet size vs BKNO3
  literature).
- **H4 the cold-bore fill delivers the heat too fast** (hot-gas front arrival
  → snap-on), now with the T11 solver.
- (Not a cause by itself but the amplifier: the erosive response, H5.)

## Measurements (no default changes; probe-only instrumentation allowed on
your branch)

1. **Ignition-threshold sweep:** head_basket and forward_plenum, pyrogen
   mass log-spaced 0.00003 → 0.03 lb (~10 points), t_max 0.2 s: ignites?,
   first-ignition time, time to 90% of grain cells ignited
   (`ignition_time_by_cell`), peak, t_peak, pyrogen burn duration. Report the
   sim's minimum ignition mass per topology.
2. **Time histories at 0.001 and 0.01 lb (head_basket):** pyrogen mass rate
   and cumulative gas; bore P at the basket; per-cell surface heat flux by
   channel (convective h_c·ΔT, DeMar/surface, radiation) and T_surf for the
   first cells to ignite; the hot-gas front position vs time; aft-cell G and
   erosive fraction near the peak. If the convective flux is not exposed,
   add probe-only instrumentation on your branch.
3. **H1 check vs Keller:** the sim's convective h (and q″) at the cells that
   ignite first vs Keller 1966's measured h at matched G / conditions —
   update §9.2's 7× figure for the current tree.
4. **H3 checks:** pyrogen depletion time and peak gas rate vs (a) a
   bed-flame-spread-limited estimate (τ_bed 5–10 ms) and (b) BKNO3 burn-rate
   literature (you may do a targeted literature lookup; cite sources).
5. **Isolating sensitivities (diagnostic only, NOT fixes):** at 0.01 lb
   head_basket, one factor at a time — (i) convective ignition heat flux ×0.5
   / ×0.25 (probe flag), (ii) pyrogen output stretched (e.g. A_burn ramped
   over 5 / 10 ms), (iii) erosive model `mukunda_paul`. For each: first
   ignition, t90, peak. These sizes are probes, not proposals.

## Deliverables

`docs/pyrogen_realism/IG1_REPORT.md`: the tables; a ranked list of causes
with the evidence for each; and 2–4 fix OPTIONS, each physically grounded
(no tuned fudge factors — the project rule is "no unfounded smoothing";
a parameter must come from literature or first principles), with what each
would change and how to validate it against the real Chunc trace (and
BALLSstick / Hasegawa A). Probes + raw results committed on your branch.