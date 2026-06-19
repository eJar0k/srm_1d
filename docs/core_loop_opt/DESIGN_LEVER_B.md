# Lever B — acoustic-CFL escape (design + mechanism + revert reference)

**Status:** DESIGNED / **DEFERRED — not started on the real solver**
(2026-06-19). The probes (§2) revealed the full fix is research-grade
surgery on the co-reviewed solver core — bigger than originally greenlit.
**Decision: do NOT touch the real solver yet.** Next steps before any
implementation: (1) user reads the primary literature (Klein 1995,
Degond–Tang 2011 — possibly with research help), (2) build an isolated
**git-worktree prototype** to measure real tractability + speedup, (3)
only then decide on a real implementation, coordinated with the areilley
code review. This doc is the spec for that work. Companion to
[`README.md`](README.md) (profiling) and
[`LITERATURE_ALLSPEED_PISO.md`](LITERATURE_ALLSPEED_PISO.md) (lit + DOIs).

This is the **review + revert reference doc**: it states the exact
mechanism, the diagnosis it rests on, the planned code changes, the
validation gates, and how to back the change out. The "AS-BUILT" section
at the end is filled in as code lands.

---

## 1. Problem (from profiling)

The core loop spends **99.9 % of steps in the steady plateau**, and the
plateau time step is **acoustic-CFL-limited**:

```
dt = cfl_target · dx / (|u|_max + a_max)        # simulation.py:1617, compute_dt_cfl
```

`a_max ≈ 1000 m/s` (hot-gas sound speed) is added to every cell's wave
speed, while the flow is only M ≈ 0.36 (u ≈ 360 m/s). We therefore take
**~(1 + 1/M) ≈ 3.8× more steps than the convective/flow timescale needs.**
Removing the acoustic term from the step limit is the structural win
(~2–3× net after the extra per-step cost), and is parameter-free —
aligned with the project's no-tuning dogma.

---

## 2. Diagnosis — two probes (2026-06-19)

Probes are non-invasive sweeps of `cfl_target` via `run_from_ric`
(`c:/tmp/cfl_probe.py`, `c:/tmp/cfl_crossmotor.py`).

### 2.1 Hasegawa A Courant headroom (`cfl_probe.py`)

| cfl | P_peak (MPa) | ΔP | n_steps | status |
|---|---|---|---|---|
| 0.30 (current default) | 6.1439 | — | 2,335,770 | ok |
| 0.50 | 6.1439 | +0.00 % | 1,401,549 | ok (1.67× fewer steps) |
| 0.70 | 0.64 | −90 % | 991 | **diverges in the fill** |
| ≥0.70 | — | — | — | blows up |

- The current scheme's **acoustic-Courant ceiling is ~0.5–0.6.**
- cfl 0.3→0.5 is bit-identical P_peak (0.3 is over-conservative
  post-v0.7.3.2 throat fix — the historical collapse was the
  cold-IC/throat bug, since fixed, not the CFL number itself).
- The hard divergence at 0.7 occurs in the **first ~1000 steps =
  ignition fill**, i.e. the ceiling is set by the high-Mach fill
  transient, not the plateau.

### 2.2 Cross-motor at cfl 0.5 (`cfl_crossmotor.py`)

| motor | P_peak(0.3) | P_peak(0.5) | ratio | verdict |
|---|---|---|---|---|
| chunc (machbusterNew) | 11.4733 | 11.4527 | 0.998 | ✓ stable |
| zerox | 7.7836 | 7.7785 | 0.999 | ✓ stable |
| BALLSstick | 8.3864 | 45.16 / 20 steps | — | **✗ diverges in the fill** |

**BALLSstick blows up at cfl 0.5 in the ignition fill.** (BALLSstick is a
known high-spike outlier and currently has uncommitted CAD-based edits in
the working tree; regardless, it is a real motor and must not break.)

### 2.3 Conclusions

1. **No safe blanket CFL bump.** 0.5 is fine for 3/4 fired motors but
   breaks BALLSstick in the fill.
2. **A phase-aware CFL (conservative fill, aggressive plateau) is
   REJECTED** as the primary mechanism: the phase trigger requires
   hand-picked Mach/pressure-rate thresholds (tuning), and it only
   reaches the same ~0.5 acoustic-Courant ceiling (~1.67×), not the full
   convective-CFL win. Fragile + tuned + small → fails the project's
   standards.
3. **The acoustic-Courant ceiling is real**, set by the fill transient,
   and the only parameter-free way past it (and to the full ~3.8×, for
   ALL motors incl. BALLSstick) is to make the **acoustic coupling
   implicit** so `dt` is limited by the convective CFL only.

---

## 3. Why the scheme is acoustic-CFL-limited (root cause)

The PISO step (`solver.py:_piso_step_with_energy_diagnostics`,
276–707) is a **non-iterated, 2-corrector segregated** scheme:

1. **Momentum predictor** (380–428): `u* = u + dt/ρ_f·(−∂P/∂x|ⁿ + conv +
   friction + S_mom)`. The pressure gradient uses the **old** pressure.
2. **Pressure correction ×2** (430–563): TDMA solve for `P'`, correct
   `u` and `P`. The transient density coefficient
   `a_t = A·dx/(R·T·dt)` (454, 521) uses the **old** temperature and the
   **isothermal** compressibility `∂ρ/∂P = 1/(RT)`.
3. **Energy** (602–682): advect sensible enthalpy `Cp·T` with the **new**
   velocities but **old** face temperatures; `T_new` is then derived.
4. **Density** (684–692): `ρ = P_new/(R·T_new)`.

The acoustic speed is `a = √(γRT)`. Two explicit/lagged couplings cap the
acoustic Courant number at ~O(1):
- the **energy/temperature lag** (T is explicit; `a` depends on T), and
- the pressure equation uses the **isothermal** compressibility
  `1/(RT)`, not the **isentropic/acoustic** one `1/a² = 1/(γRT)`.

A correct acoustic-implicit (IMEX) scheme couples `p`–`ρ`–`(T or ρe)`
implicitly so the acoustic eigenmode is unconditionally stable, leaving
only the convective CFL.

---

## 4. Mechanism — IMEX semi-implicit acoustic coupling

Following Klein (1995, §4; DOI 10.1016/0021-9991(95)90034-9) and the
Degond–Tang (2011, arXiv:0908.1929) all-speed pattern, adapted to our
staggered + TDMA structure. The **convective terms stay explicit** (so
the existing upwind momentum/enthalpy advection and the combustion
mass/heat sources are unchanged); the **acoustic sub-operator becomes
implicit**:

- **Momentum:** pressure gradient on the **new** pressure
  (`−∂Pⁿ⁺¹/∂x`) instead of `Pⁿ`. In the staggered PISO this is already
  what the pressure-*correction* delivers; the change is to form the
  predictor/corrector so the implicit acoustic operator is consistent.
- **Pressure/continuity:** use the **acoustic (isentropic)
  compressibility** `∂ρ/∂P|ₛ = 1/(γRT)` in the transient diagonal so the
  implicit operator is the true acoustic Helmholtz operator, and couple
  the temperature/energy response to the pressure change implicitly
  (the crux — the energy lag must be removed or linearized into the
  pressure equation).
- **Time step:** `compute_dt_cfl` switches to the **convective** wave
  speed `|u|` (drop `a`) when the implicit-acoustic mode is active,
  capped by a convective Courant ≈ the scheme's convective stability
  (to be measured; expected ~0.5–1.0). Net `dt` gain ≈ `(|u|+a)/|u| ≈
  3.8×` at the plateau.

**Crux / risk:** the energy–pressure implicit coupling (removing the T
lag) is the hard part and the main correctness risk. Increasing PISO
correctors alone does NOT reach the convective CFL (the probe shows a
hard ceiling, not a soft one), so the energy coupling must be addressed.

---

## 5. Implementation strategy — OPT-IN, default OFF (revert-safe)

Mirror the project's established pattern for new solver mechanisms
(`port_mach_cap`, `flame_front_enabled`, `zn_enabled` all shipped
default-off, validated, then canonized):

- Add a boolean flag, e.g. **`implicit_acoustic=False`**, plumbed
  `run_simulation` → `_run_time_loop` → `piso_step` /
  `_piso_step_with_energy_diagnostics`, and the convective-`dt` branch in
  `compute_dt_cfl`. Default OFF ⇒ **byte-identical current behavior**;
  the whole feature is a no-op until explicitly enabled.
- A/B harness (`c:/tmp/`) compares OFF vs ON: P_peak, full trace, steps/s,
  cross-motor stability.
- **Canonize (flip default) ONLY after** the validation gates (§6) pass
  on all four fired motors. Until then it ships off — zero risk to the
  validated v0.8.0 physics.

This is the "in case of revert" guarantee: with the flag off the change
is inert; reverting is either `implicit_acoustic=False` (runtime) or
`git revert` of the feature commit (source). No entanglement with the
validated default path.

---

## 6. Validation gates (must all pass before any default flip)

1. **Result fidelity:** P_peak within ±1 % and full head-pressure trace
   visually matching the OFF baseline on **Hasegawa A, Zerox, Chunc** (the
   experimentally-validated set), and a stable, non-diverging BALLSstick.
2. **Stability:** no `verify_run_health` failure; no fill divergence on
   any of the four fired motors at the new convective `dt`.
3. **Speedup:** measurable net wall-clock reduction (target ≥1.5×, ceiling
   ~3.8× minus extra per-step cost) at N=100, confirmed by benchmark.
4. **Regression:** full pytest green (currently 406 pass; the 2
   motor-file-conformance failures are pre-existing/unrelated).
5. **Grid:** P_peak grid-convergence at N=50/100/200 unchanged in
   character vs OFF.

---

## 7. Honest scope / risk assessment

- This is a **research-grade change to a validated solver core**, not the
  "50–150 LOC" the lit survey optimistically suggested. The energy–pressure
  implicit coupling is genuinely intricate and is where correctness risk
  concentrates.
- Mitigations: opt-in/default-off (above); incremental landing
  (plumbing → convective-dt branch → implicit pressure operator →
  energy coupling), validating at each step; A/B against the frozen
  baseline; no default flip without §6.
- Coordinated with the ongoing areilley codebase review — this doc is the
  review artifact.

---

## 8. Rejected alternatives (recorded so they aren't re-litigated)

- **Blanket cfl 0.3→0.5** — breaks BALLSstick in the fill (§2.2).
- **Phase-aware CFL** — needs tuned phase thresholds; only ~1.67×; fragile.
- **AUSM⁺-up face flux** — paradigm mismatch: it is a density-based /
  collocated FV interface-flux scheme; grafting it onto a staggered
  pressure-correction PISO is more invasive than the lit survey implied.
  Keep as a fallback if the IMEX energy coupling proves intractable.
- **Density-based Godunov rewrite** — the SPINBALL/SPP approach; correct
  but a full solver replacement, out of scope.
- **More PISO correctors** — does not reach convective CFL (hard ceiling
  at ~0.6, per probe); would only nudge the acoustic-Courant limit.

---

## 9. AS-BUILT (filled in as code lands)

> _Pending implementation. Will record: exact files/functions/lines
> changed, the flag name + plumbing path, the convective-dt branch, the
> implicit operator form actually used, measured speedup + P_peak on each
> fired motor, and the exact revert command._

- Flag: _TBD_
- Files touched: _TBD_
- Measured speedup: _TBD_
- Cross-motor P_peak (OFF vs ON): _TBD_
- Revert: `implicit_acoustic=False` (runtime no-op) / `git revert <sha>`.
