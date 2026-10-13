# EN1 — conservative total-energy formulation (opt-in)

**User decision 2026-10-13.** This runs before the head-end model (HE1), because
the gas-side fill flow feeds ignition heat transfer. Thread F
(`docs/v0_7_4/IGNITION_SPIKE_REOPENED.md` §8, §13) and the T6 energy audit
(`worklog/T6.md`) are the background.

## The problem

`solver.py` STEP 3b transports **sensible enthalpy only**:

```
m_new·Cp·T_new = m_old·Cp·T_old + dt·(Σ face mdot·Cp·T_upwind + thermal_source·dx)
```

For a fixed control volume the stored quantity should be the internal energy
m·Cv·T plus the kinetic energy m·u²/2, and the face fluxes should carry the
total enthalpy mdot·(Cp·T + u²/2). Storing m·Cp·T leaves out the pressure work
V·dp/dt and the kinetic energy.

Consequences:
- Gas that accelerates does not cool. At Mach 1.4 with γ ≈ 1.2, static T should be
  ~16 % below stagnation.
- Shocks cannot satisfy the Rankine-Hugoniot jump conditions, so their position,
  strength and post-shock state are wrong.
- Compression heating during the fill is missing.
- The between-step EOS reconciliation jump (T6: ~30 % of the fill's source
  energy) is plausibly the same root cause.

The user observed Chunc fill at ~Mach 1.4 at 5 ms (T6 measured Mach ~2.9 in
earlier fill snapshots, with the KE term up to 80 % of local enthalpy).

## Part 1 — formulation (opt-in kwarg, e.g. `energy_formulation='sensible'` default | `'total'`)

Conservative total energy per cell:
- storage: m·(e + u²/2) with e = Cv·T, using the per-cell mixture Cv = Cp − R;
- face fluxes: mdot·(Cp·T + u²/2), upwind, consistent with the existing species
  and Cp treatment;
- sources: injected gas brings its own enthalpy, as `thermal_source` does now,
  with zero axial kinetic energy; wall heat-loss sinks unchanged.

Design choices to make and justify:
- how to build a cell kinetic energy from the staggered face velocities;
- how T is recovered;
- how the energy update couples to the PISO pressure correction and to the P5
  energy-coupled iteration (`piso_energy_coupling_iters=1`).

Keep the T ceiling and floor clips and their ledgers. Read the T6 and P5
worklogs first. Also read the core-loop notes on the compressibility term in the
pressure operator (`docs/core_loop_opt/`).

## Part 2 — verification against exact solutions (both formulations)

Use `solver.py` directly. Its module is pure numerics: build minimal test
harnesses for:

1. **Sod shock tube** (ideal gas, γ = 1.4): exact Riemann solution for density,
   velocity, pressure and T. Report shock position and speed, plateau values, and
   L1 errors vs grid. The sensible form should show wrong jump conditions; the
   total form should converge to the exact solution, apart from first-order
   smearing.
2. **Steady isentropic quasi-1-D nozzle flow**: T/T0 = 1/(1 + (γ−1)/2·M²) along a
   converging-diverging duct.
3. **Normal shock in a converging-diverging duct** at a back pressure that puts it
   in the divergent section: exact quasi-1-D shock location and the post-shock
   state.
4. **Closed-vessel compression**: gas injected into a closed volume. The energy
   balance gives the exact pressure and temperature rise, including the
   pressure-work heating the sensible form leaves out.

State the boundary conditions you use and how they map onto the motor's
head-closure and nozzle BCs.

## Part 3 — effect on the motors

Run both formulations, with the current defaults otherwise:
- Chunc fill: Mach map vs time (snapshots every 1 ms), supersonic extent, shock
  position, bore-gas T at the first-lit cells;
- ignition timing: first light and t90;
- startup peak, plateau and ratios for Chunc, BALLSstick and Hasegawa A
  (canonical B0 rows);
- the T6 energy ledger: does the EOS jump shrink?
- the head-end dead zones at fixed charges (BALLSstick 4.83 g, Hasegawa A
  14.98 g, forward plenum, DeMar off): cells never lit;
- wall-time cost.

## Gates and deliverables

- B0 PASS with defaults (`'sensible'`, bit-identical). Full suite green.
- `docs/solver_verification/worklog/EN1.md` with the verification tables, the
  motor tables, the cost, and a recommendation. The user decides the default; a
  default flip re-records B0.
- Opt-in implementation, tests (exact-solution checks; default identity), probes.
