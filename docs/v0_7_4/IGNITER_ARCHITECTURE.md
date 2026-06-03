# Igniter / Pyrogen Code — Architectural Map (v0.7.4)

Map of the igniter code paths + the Chunc head_basket audit that
verified the igniter is NOT over-injecting and NOT vented through a
choked throat (user questions, 2026-05-31).

## Two topology paths

```
run_from_ric (openmotor_adapter.py)
  └─ build_pyrogen_chamber(pyrogen, geo, nozzle, injection_topology, ...)
        sizing (if not user-given):
          m_pyrogen      = sutton_pyrogen_mass(free_volume_in3)   # 0.12·V_F^0.7 g
          A_burn_initial = particle geometry:                     # v0.7.4 Phase C.1
                             sphere   6·m/(ρ·d)
                             cylinder m·(4λ+2)/(ρ·λ·d)
          A_throat       = A_burn/Kn(=100), clamped [1,100] mm²   # vestigial for baskets
          V_plenum       = 1.5·m/ρ                                # vestigial for baskets
        → PyrogenChamber(injection_topology, cartridge_length_m, ...)

run_simulation → _run_time_loop  (simulation.py), per step, STEP 3:
  if topology_code == 0  (forward_plenum):
      _step_plenum_ode(...)              # igniter_plenum.py
        ├─ _plenum_pressure  P_ig = m_gas·R·T/(M·V_plenum)   (ENCLOSED volume)
        ├─ _burn_area·a·P_ig^n           burn at PLENUM P_ig
        └─ _choked_orifice_mdot(P_ig,T,A_throat,...)  ← CHOKED THROAT vent into cell 0
      + axial momentum injected at face 1 (igniter_axial_momentum_fraction)
      + mass/enthalpy distributed by Phase A exponential axial weights

  else  (head_basket / aft_basket — UNCONTAINED):
      _compute_uncontained_pyrogen_mdot(...)   # simulation.py
        r_b[i]  = a · max(P_bore[i], 0)^n           ← burns at LOCAL BORE P
        mdot[i] = ρ_p · r_b[i] · A_burn_per_cell     A_burn_per_cell = A_burn_initial/n_cart
        (cartridge cells [cart_i_start, cart_i_end] from resolve_injection_cells)
      → mass + enthalpy injected at T_flame into the cartridge cells
      → NO choked orifice, NO plenum pressure, NO momentum injection
        (PISO handles axial flow via the pressure gradient)
```

## Chunc head_basket audit (mtv, machbusterNew.ric)

| quantity | value |
|---|---|
| mtv a, n | 3e-5, 0.5 (fast pyrogen) |
| ρ, T_flame, γ, M | 1800 kg/m³, 3000 K, 1.22, 0.032 kg/mol → Cp≈1441 J/kgK |
| particle | 3.2 mm sphere |
| m_pyrogen (Sutton) | **0.90 g** |
| A_burn_initial | 9.4 cm² (A/m = 1.0 m²/kg) |
| cartridge cells | **[0,0] = 1 cell** (L_cart 3.94 mm < dx 8.13 mm → snaps to 1) |
| mdot vs local bore P | 16 g/s @0.1 MPa → 51 @1 → **203 @16 MPa** |
| enthalpy power | 0.07 → 0.22 → 0.88 MW |
| depletion time | 56 → 18 → 4.4 ms |

**Verified answers to the user's questions:**
1. **Not choked-throat vented.** head_basket uses `_compute_uncontained_pyrogen_mdot` (local-bore-P burn); the choked-orifice path (`_choked_orifice_mdot`) and `A_throat`/`V_plenum` belong to `forward_plenum` and are vestigial for baskets.
2. **Not over-injecting.** The igniter is modest — 0.9 g (Sutton), concentrated in cell 0, mdot 16–203 g/s, ~0.07–0.88 MW. It is NOT the driver of the fast fill; it ignites the head region and the propellant cascade does the rest.
3. **Burns at LOCAL bore P** → mdot ∝ P^0.5 surges with the chamber spike (positive feedback P↑→mdot↑), but the 0.9 g mass makes it a brief kick (depletes in 4–18 ms), so it does not dominate the spike.

## How the fast supersonic fill actually arises (not the igniter)
1. Pyrogen in **cell 0** ignites; its radiation lights the near cells (cell 3 ≈ 51,000 W/m², falling as 1/d² to 31 W/m² at the aft — so radiation lights only the head region).
2. Those cells produce propellant gas → the hot products expand into the **cold-ambient bore** (Phase B.0 IC: `T_initial_gas = T_ambient`). This is a shock-tube driver/driven setup.
3. The hot/cold contact accelerates to ~sonic **at the interface** (M≈1–1.5) and sweeps the bore at ~800 m/s, igniting cells by advection as it passes; the surface ignites ≈ as fast as the gas fills.
4. A brief M≈3.7 PISO velocity spike occurs at the last bore cell when the sharp contact vents into the nozzle (t≈2.76 ms) — a transient at the discontinuity, *before* the 8 ms pressure spike; the nozzle BC itself (`_nozzle_boundary_flow`) is a correct quasi-steady choked/subsonic mass boundary.

## Conclusion (cross-checked from 5 angles)
Goodman kernel, T_ceiling clip, wall heat-loss sinks, radiation reach,
advection over-speed, AND the igniter are all either physically correct
or quantitatively negligible. The gas genuinely fills fast (acoustic /
shock-tube into the cold-ambient bore), matching the real motors'
"snap-on." So the **simultaneous ignition is largely physical, and the
2× spike is the erosive quasi-steady over-response (Root B)** firing off
the genuine peak-G at the smallest-bore condition — not an
ignition/igniter/advection artifact. The literature offers no
off-the-shelf transient erosive closure; the one physically-grounded
(non-smoothing) lever is the Beddini turbulent-BL-development factor on
the erosive Nu (erosion requires developed core turbulence; during the
un-developed ignition fill it is below the quasi-steady Ma value).
