# Ignition-spike investigation probes

Empirical measurement scripts from the ignition-spike investigation (v0.7.4 +
the 2026-06-16..18 re-open). **These are investigation artifacts, not shipped
tooling** — scratch quality, most carry **machine-absolute paths** (this
developer's box) and some need an external experimental trace (see caveats).
They are committed so a fresh session can *reproduce the exact measurements*
that produced the recorded findings instead of re-deriving them — but they are
**not gospel**: re-run, question, and modify them freely (see
`../../v0_7_4/` and the reopen brief).

Run with the pyenv 3.10.5 interpreter from the repo root, e.g.
`"C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" docs/v0_7_4/probes/<name>.py`.
Chunc = `motors/machbusterNew.ric` (the cleanest simultaneous-ignition
diagnostic — real motor has ZERO spike → flat 8.8 MPa; sim over-spikes ~2×).

## Mass-flux G — the crux measurement
The erosion-driving aft-port mass flux `G = ρ|u|`. Key finding: **sim G is only
~1.29× experimental, with NO transient gas-dynamic inflation** — the over-spike
is (mostly) a faithful erosive *response*, not an inflated *driver*.
- `chunc_G_reconstruction.py` — runs Chunc, dumps fine snapshots + meta; computes
  G three ways (exp-QS / sim-QS / sim-actual `ρ|u|`). **Needs the high-res
  experimental xlsx** (`ThomasMach5_edited.xlsx`, this box's Downloads).
- `chunc_G_analyze.py` — reconstructs exp-vs-sim aft-port G from the saved
  `_chunc_snaps.npy`. **Run `chunc_G_reconstruction.py` first.**
- `chunc_G_plot.py` — overlay plot of exp vs sim G + pressure.

## Supersonic fill — artifact vs physical
Finding: the Mach 1.3–12 fill is a **grid-divergent PISO artifact** but
**DECOUPLED from P_peak** (`port_mach_cap` bounds the Mach without changing the
spike). Relevant to the user's "sonic propagation / shock modeling" suspect.
- `chunc_mach_audit.py` — is the fill Mach>1 continuity-consistent, or a PISO
  blowup? Flags the structural misuse of feeding Mach>1 blowdown velocity into
  the steady-pipe Gnielinski correlation.
- `chunc_mach_convergence.py` — refine cells + tighten CFL; a physical blowdown
  converges, a numerical over-acceleration drifts (it drifts → artifact).

## Erosive-closure levers (what's been tried)
- `chunc_K_relaminarization.py` — favorable-pressure-gradient relaminarization
  gate `K=(ν/u²)(du/dx)` (crit ~3e-6). **TESTED DEAD**: K≈1e-7 at erosive cells,
  so it can't gate the turbulent Ma term without wrongly killing steady erosion.
- `chunc_mukunda_compare.py` — sim Ma-2020 `η=r/r₀` vs the Mukunda-Paul (1997)
  universal law per cell. Findings: MP threshold g=35 is NOT tripped-suppressed
  (g≫35 at the spike); **Ma over-predicts ~33% vs MP at high g**. The main
  physics frontier (but MP-vs-Ma is a fidelity call, not a free lunch).
- `chunc_flux_vs_keller.py` — sim convective ignition flux q″ vs the
  Keller-Baer-Ryan 1966 convective-AP ignition regime (validates the surface-T
  criterion over Fs=20–160 cal/cm²·s).
- `chunc_hc_decomp.py` — decompose the ignition flux into `h_c` and ΔT at each
  cell's ignition; compare `h_c` to Keller's measured convective range.
- `chunc_hc_probe_sweep.py` — knock down the shared convective `h_c` (ignition
  gate + Ma erosive) for `t<window`. **Probe-FALSIFIED as a fix**: it only
  POSTPONES the spike (t_peak tracks the window), doesn't remove it.

## Ignition fill / gate physics
- `chunc_fill_physics.py` — fill-transient audit: port/throat area ratios (is
  high bore Mach even expected?), fill mass balance, nozzle choke state.
- `chunc_fill_probe.py` — how does the surface reach T_ign across the whole ~1 m
  grain in <1.5 ms? Dumps T_gas/T_surf/P/is_burning at 0.25 ms cadence.
- `chunc_spike_probe.py` — general structural spike probe. **Its own docstring
  captures the reopen ethos: "treat the spike as a candidate STRUCTURAL
  artifact, not assumed-faithful Ma physics."**

## Igniter gas generation (ProPep fix)
Finding: correcting bpnv gas generation (ProPep impetus + condensed-phase gas
fraction) cut the Chunc spike 1.50→1.36× (Hasegawa unchanged); shipped as
`bfc2f3f`. **Marked for later:** the same fix for MTV/thermite via CEA.
- `igniter_impetus_test.py` — impetus match (M 30→41 g/mol) + condensed-phase
  gas-fraction A/B.
- `igniter_survey.py` — Chunc igniter sizing / impetus / gas-generation vs lit.

## Sonic / acoustic-CFL behavior (core-loop-opt, sonic-propagation relevant)
- `cfl_probe.py` — Courant-headroom sweep on Hasegawa A (the current scheme's
  acoustic-Courant ceiling ~0.5–0.6; divergence at 0.7 is in the *fill*).
- `cfl_crossmotor.py` — cfl 0.3 vs 0.5 across chunc/zerox/BALLSstick (BALLSstick
  diverges at 0.5 in the fill → no safe blanket bump). See `../../core_loop_opt/`.

## Numerical-collapse forensics
- `bisect_lhs_collapse.py` — bisects why some LHS calls collapse where the
  canonical example runs clean (v0.7.3.2 Kn-throat / cfl / source_cfl fix era).

---

**Caveats for a fresh session:** (1) paths are absolute for this box — adjust if
moved. (2) The G-reconstruction chain needs `ThomasMach5_edited.xlsx` (the
high-res Chunc static-fire trace) — the single most valuable validation asset;
consider committing it under `static_fire_data/`. (3) Constants like
`MU=9.125e-5`, `RHO_P=1700` are per-probe stand-ins — verify against the actual
propellant/gas before trusting absolute numbers. (4) These encode *what was
measured and found*; the **methods** are the durable value, not the exact code.
