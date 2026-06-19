# Core-loop optimization + PISO grid analysis

**Thread opened 2026-06-18** (branch `openmotor-frontend`). Goal:
speed up the core PISO time loop, analyze the PISO grid bugs
(especially the core-Mach flow divergence), and enable denser grids.

This README holds the **profiling findings + decisions**. The companion
[`LITERATURE_ALLSPEED_PISO.md`](LITERATURE_ALLSPEED_PISO.md) holds the
full numerical-methods literature review (with DOIs) that informs the
structural performance work.

---

## 1. Profiling (macro, non-invasive)

Harness: `c:/tmp/bench_srm1d.py` (not committed) — runs canonical
Hasegawa A (`motors/hasegawa_a.ric`, v0.7.5 knobs) to `P_cutoff` at
`target_propellant_cells` = 50 / 100 / 200, JIT warmed.

| metric | N=56 | N=106 | N=206 | scaling |
|---|---|---|---|---|
| wall (10 s sim) | 15.4 s | 44.0 s | 146.8 s | **~N¹·⁷³** |
| n_steps | 1.17 M | 2.34 M | 4.67 M | ~N¹·⁰⁶ |
| steps/s | 75.7 k | 53.1 k | 31.8 k | ~N⁻⁰·⁶⁷ |
| fill steps (% of total) | 0.1 % | 0.1 % | 0.1 % | — |

*(steps/s here is from the first warmed bench; a later quiet-machine
A/B measured baseline ~35–39 k at N=100, i.e. ±10 % machine variance.
The scaling exponents are the durable result, not the absolute steps/s.)*

### Key findings

1. **The ignition fill — where the Mach divergence lives — is 0.1 % of
   the step budget.** The smallest `dt` is in the *plateau*, not the
   fill (`dt_min` global < `dt_min`-in-fill). **The high-Mach fill
   artifact costs essentially nothing.** This OVERTURNS the prior
   working assumption (shared in the design discussion and implied by
   `IGNITION_SPIKE_REOPENED.md` §6) that the Mach divergence is the
   performance tax / densification-cost gate. It is not.

2. **Cost is 99.9 % the steady plateau, and the plateau is
   acoustic-CFL-limited.** `dt` halves *exactly* as `dx` halves
   (8415 → 4207 → 2103 ns), and `CFL·dx/dt ≈ 1200 m/s ≈ sound speed`.
   The flow itself is only M ≈ 0.36 (u ≈ 360 m/s), so we resolve
   acoustic waves we do not care about — **~3.3× more steps than the
   flow/convective timescale requires.** This is the real structural
   lever (see Lever B).

3. **Densification cost exponent ≈ N¹·⁷.** Doubling resolution → ~3.3×
   wall. `n_steps ~ N` (acoustic CFL → linear step growth);
   `steps/s ~ N⁻⁰·⁷` (per-step cost is sub-linear because of a
   ~6 µs/step N-independent component).

### Allocation microbench (A1 de-risk)

Harness: `c:/tmp/alloc_microbench.py`. The PISO hot path does **17 heap
allocations/step** (11 in `_piso_step_with_energy_diagnostics` + 3×2 in
the two `thomas_solve` calls).

```
17 np.zeros @ N=100:  ~950 ns/step   (~56 ns/array alloc+free)
real per-step @ N=100: ~18.9 µs
A1 preallocation ceiling: ~0.95 µs  →  ~5 % of per-step
```

**Allocation is only ~1 µs of the ~6 µs fixed cost** — NOT the
dominant fixed cost as first hypothesized. **A1 (preallocate PISO +
Thomas scratch buffers) is DEPRIORITIZED**: ~5 % for a ~17-buffer
signature break across multiple @njit functions + cache management +
test churn is a poor trade.

### fastmath A/B (the cheap keeper)

Flipping `@njit(cache=True)` → `@njit(cache=True, fastmath=True)` on the
hot kernels in `solver.py` / `burn_rate.py` / `simulation.py`:

| build | steps/s @ N=100 | n_steps | P_peak |
|---|---|---|---|
| baseline | 35.4 k / 39.0 k | 2,335,770 | 6.1439 MPa |
| **fastmath** | **48.5 k / 51.2 k** | **2,335,770** | **6.1439 MPa** |

**~+30 % steps/s, with `n_steps` and `P_peak` bit-identical** — the
trajectory is unperturbed, so fastmath is physically safe here.
**VALIDATED + COMMITTED 2026-06-19:** full suite 406 passed; the 2
failures (`test_all_repo_motors_have_parseable_igniter`,
`test_repo_motors_all_migrated_and_loadable`) are pre-existing
working-tree motor-file conformance issues (`BALLSstick VANILLA.ric`
has no igniter block; a `'Chase - Energy'` motor has unmigrated frozen
transport) — unrelated to fastmath (which touches only @njit numeric
kernels, not YAML parsing).

---

## 2. Lever map (revised after profiling)

**Lever A — per-step cost (steps/s); helps every resolution:**
- **A1 preallocate PISO/Thomas scratch** — DEPRIORITIZED (~5 %, see above).
- **A3 `fastmath=True`** — ✓ ~+30 %, result-identical. The cheap keeper.
- **A2 burn-rate bisection warm-start** — `compute_burn_rates` does a
  per-cell Haaland→Gnielinski→bisection root-find each step; seed from
  the previous step's root. Likely the dominant O(N) per-cell cost;
  needs a per-function attribution to size. OPEN.

**Lever B — step budget (n_steps); the structural densification win:**
- **Escape the acoustic CFL** so `dt` is limited by the flow CFL, not the
  acoustic CFL. Theoretical ceiling ≈ 3.3× fewer steps across the whole
  plateau. **Lit-backed (see companion doc): IMEX acoustic–convective
  splitting (Klein 1995 / Degond–Tang 2011) OR AUSM⁺-up face flux
  (Liou 2006).** ~2× net wall-clock after the extra per-step cost; both
  parameter-free. This is the real prize. NOT started — needs explicit
  go-ahead (moderate but well-defined refactor of the JIT loop into
  split explicit/implicit stages).

**Mach divergence (correctness / densification stability, NOT perf):**
the existing default-off `port_mach_cap` limiter (commit `16bc527`) is
confirmed by the lit review as the correct, sufficient, parameter-free
fix. AUSM⁺-up would resolve it more principally *and* serve Lever B.

---

## 3. Domain-practice context (for the upstream / areilley discussion)

Every established SRM Q1D internal-ballistics code (SPINBALL, SPP, NAWC
heritage) is **density-based** (Godunov / Roe / MacCormack), which
handles the full Mach range — including the fill — natively.
srm_1d's **pressure-based PISO is the architectural outlier**: ideal for
the M≈0.3 plateau, structurally mismatched to the fill transient. Lever B
(IMEX) + `port_mach_cap` bring our behavior into conformance without a
density-based rewrite.

---

## 4. Status / next steps

- [x] Macro profiling (overturned the Mach-as-perf-tax assumption).
- [x] A1 de-risk → deprioritized (~5 %).
- [x] fastmath → ~+30 %, result-identical; **validated (406 passed) +
      committed** 2026-06-19.
- [ ] **A2** burn-rate bisection attribution + warm-start (cheap O(N) win).
- [~] **Lever B** acoustic-CFL escape (IMEX or AUSM⁺-up) — the 2× prize;
      **IN PROGRESS** (greenlit 2026-06-19). Design + mechanism in
      [`DESIGN_LEVER_B.md`](DESIGN_LEVER_B.md); companion lit doc has the
      recommended path and the paywalled papers to obtain.

Benchmark harnesses live in `c:/tmp/` (`bench_srm1d.py`,
`alloc_microbench.py`, `quick_bench.py`) — promote to `srm_1d/tools/`
if/when this work is committed.
