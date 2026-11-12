# NB1: the nozzle boundary uses the true stagnation state

**User decision (2026-11-10):** "Prefer a physical fix, so taking the isentropic
relations to find actual stagnation properties at that choked boundary."

## The defect

`_nozzle_boundary_flow` (`srm_1d/solver.py`, about line 246) feeds the last
cell's STATIC P and T into the throat flow function as if they were the
stagnation (reservoir) values. The gas leaving the port is moving, so the
true stagnation pressure is higher:

    P0 = P (1 + (γ-1)/2 M²)^(γ/(γ-1))

At Chunc's grain-exit Mach (about 0.35) the boundary understates P0 by about
7 %. For a given gas generation, the throat then passes too little mass and
the chamber's static pressure settles too high. The error grows with port
Mach, so it is largest during the startup transient. That makes it relevant
to the ignition over-spike, the project's foremost open problem
(CLAUDE.md, roadmap item 1).

G1 currently works around this assumption. `grain_geometry.bc_reservoir_area_ratio`
and `BC_RESERVOIR_P0_TOL` (1 %) size the optional `aft_end='cone'` cut so the
boundary sits where static ≈ stagnation. Review them after the fix (step 5).

## Physics to implement

The port-exit face N has flow area `A_face[N] = A_port[N-1]`, called A_e
here. `u[N]` is the velocity at that area, not at the throat. Quasi-steady
isentropic flow runs from the port exit (static P, T of cell N-1, Mach M_e)
through the throat A_t.

**Solve the boundary state self-consistently. Do not lag a velocity.**

- **Choked:** M_e is the SUBSONIC root of the area-Mach relation for
  A_e/A_t. Then:
  - P0 = P · g(M_e), with g as in the formula above;
  - T0 = T (1 + (γ-1)/2 M_e²);
  - ṁ = Γ(γ) A_t P0 / sqrt(R T0).

  For fixed M_e, P0 is linear in P, so `dmdp` follows directly. Derive it,
  and check it against a finite difference.
- **Subsonic outflow:** M_e and the throat Mach come from continuity plus
  isentropic relations with the throat static pressure equal to P_ambient.
  This is a 1-D root solve: a few Newton steps, safeguarded (bisection
  fallback, bounded iterations), @njit-safe. Derive its `dmdp` too.
  - The choked/unchoked switch must be continuous: compare P_ambient with
    the critical ratio of P0, not of the static P.
- **A_e ≤ A_t** (port exit no larger than the throat): there is no subsonic
  area-Mach root. Handle it physically: the port exit itself is the sonic
  section. Document and test it; it must not crash or NaN.
- **Reverse inflow** from the ambient (reservoir at rest): unchanged.

**Energy flux, a trap to avoid.** The callers build the outflow energy flux
as `mdot * (Cp_i * T_upstream + ke_old_i)`, i.e. static T plus kinetic
energy, which already equals h0. If the function starts returning T0 as
`T_upstream`, kinetic energy is double-counted. Keep the enthalpy at static
T + KE, or change the return signature so the mass-flow stagnation state and
the advected enthalpy are separate. Audit every call site in `solver.py`
(about 10) and `simulation.py:~2764`.

Use one consistent boundary state in all of them: the predictor, correctors
1 and 2, the energy update, the total-energy path, the accounting call
(~1129) and the ledger.

`grep -n _nozzle_boundary_flow srm_1d/*.py` lists the call sites.

## Toggle and default

- Add `nozzle_bc_stagnation: str = 'isentropic'` as the DEFAULT (user
  decision). `'static'` keeps the old behaviour for A/B comparisons.
- Thread it through `run_simulation`, the plugin's solver-config defaults
  (generated from the signature) and `run_from_ric` via `sim_overrides`.
- Log the @njit signature change in the DEVNOTES "API Breaking Changes Log".
  Clear `srm_1d/__pycache__/` after every @njit edit.

## Verification (exact answers first)

1. **Unit tests of the boundary function:**
   - choked P0/P against g(M_sub(A_e/A_t)) to 1e-10;
   - ṁ against the textbook choked flow;
   - continuity of ṁ and of dmdp across the choked/unchoked switch;
   - dmdp against a central finite difference on both branches;
   - the A_e ≤ A_t case;
   - reverse flow unchanged.
2. **Steady duct test.** A constant-area duct with a steady head-end mass
   injection (or a uniform wall source) feeding a choked throat, run to steady
   state. Compare the converged head stagnation pressure with the analytic
   value (continuity + isentropic throat + the 1-D flow with mass addition, or
   the simplest case you can solve exactly). The new BC should match to
   discretization error; the old one should miss by the predicted g(M_e)
   factor. Grid N = 25/50/100.
3. **Conservation.** The mass and energy ledgers (the existing EN1 harness)
   still close to their current tolerance.
4. **Canonical motors under the guard**, at N = 100, `'static'` against
   `'isentropic'`:
   - Chunc (`machbusterNew.ric`): peak, plateau, peak/plateau;
   - BALLSstick: peak/trough, peak/mid;
   - Hasegawa A: early and late peak against the 6.44 MPa experiment;
   - for each: the grain-exit Mach and P0/P at the plateau.
   - Use the committed `.ric` files (`git show HEAD:motors/<f>.ric`); the main
     checkout has user edits.
5. **G1's reservoir assumption.** With the boundary using true stagnation
   values, does `bc_reservoir_area_ratio` / `BC_RESERVOIR_P0_TOL` still have a
   purpose? Recommend keep, retire or repurpose, with reasons. Do not delete it
   without the user's decision; it is only used by the opt-in `aft_end='cone'`.
6. **Report only, do not implement: the aft dump loss.**
   - In most motors the port exit dumps into a larger aft plenum or convergent
     entrance (about the case ID) before the nozzle converges. A sudden
     expansion loses stagnation pressure (Borda-Carnot, from the 1-D momentum
     balance, no free constants).
   - For each canonical motor's aft geometry, estimate the plateau loss as a
     fraction of P, next to the isentropic correction. The user wants to know
     whether it matters.
7. **Thrust and performance.** Check what pressure `nozzle.compute_motor_performance`
   and the thrust path use as chamber stagnation pressure. If they reuse the
   static-as-stagnation assumption, report it and propose the consistent fix.

## Gates

- The full pytest suite, under the guard with a 60-min timeout (it takes
  about 45 min). Set
  `SRM1D_OPENMOTOR_PATH="C:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/openMotor/openMotor"`;
  the worktree is outside the repo tree and the FMM tests need it.
- **B0** (`docs/solver_verification/probes/baseline_gate.py`): values WILL move
  with the new default. Re-record, keep the old values in comments, and give
  the before/after table with an explanation of each move.
- Another agent's branch, G2 (`g2/snap-lump`), will probably merge before you
  finish. When it lands, merge `openmotor-frontend` into your branch (merge,
  not rebase) and re-run B0 on top.

## Rules

- **Compute:** every job through
  `python -m srm_1d.tools.guarded_run --mem-gb 8 --timeout-min 60 -- <cmd>`.
  At most 4 simulation processes. Warm the Numba cache with one small run
  before parallel jobs. Probes are memory-bounded and have max-step stops.
- No tuned constants. Do not hand-edit `motors/*.ric`. Do not touch the main
  checkout.
- Commit on `nb1/stagnation-bc` in `C:/tmp/sv/nb1`. Do not merge into
  `openmotor-frontend` and do not push.
- Write `docs/solver_verification/worklog/NB1.md`: what changed, the
  verification results and the tables above.
- Report:
  - the before/after table;
  - the step 5 and step 6 findings;
  - the step 7 finding;
  - anything surprising.
