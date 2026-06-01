# v0.7.4 — Transient/Unsteady Erosive-Closure Lit Dive

**Question** (2026-05-31): the ignition-model audit concluded the Chunc
spike is the *erosive* burn-rate over-response (Root B) — Ma's quasi-steady
erosive firing instantly off the genuine peak-G at the smallest-bore
condition (erosive fraction 0.13 → 0.61 → 0.13). Is there a physically-
grounded **transient / development-lag erosive closure** that ramps the
erosive augmentation over a real timescale (so it doesn't surge at
ignition) while keeping Ma 2020 as the steady limit?

**Bottom line: NO off-the-shelf transient erosive closure exists.** Every
erosive model in the literature is quasi-steady; unsteadiness is handled
by a Z-N dynamic burn rate on the *total* rate (which we already shipped
as Phase Z and found self-attenuating at the spike). The one defensible
new lever is a turbulent-boundary-layer *development* factor on the
erosive enhancement (Beddini framework), but its ~ms timescale likely
gives only modest additional suppression.

## Findings

1. **All erosive closures are quasi-steady** — instantaneous functions of
   local mass flux G / Reynolds number, no time lag:
   - Lenoir-Robillard (1957); Ma 2020 (our model); Mukunda-Paul "universal"
     (steady BL); King 1993 review (JPP 9(6)); Beddini turbulent-BL
     (AIAA J, DOI 10.2514/3.61303). WebSearch + Hasegawa/Ma refs.
2. **Unsteadiness → Z-N dynamic burn rate on the TOTAL rate**, not an
   erosive-specific lag. Confirmed in the Cavallini SPINBALL thesis
   (`docs/references/74323997.pdf`): burn-rate model = §2.2.1 quasi-steady
   APN + §2.2.2 **erosive (Lenoir-Robillard + Lawrence + Beddini),
   quasi-steady** + §2.2.3 **dynamic (Z-N), a separate term**. This is
   exactly the architecture we already have (Phase Z = Z-N on the total).
   Phase Z's τ = κ·α_s/r² **self-attenuates** (high r at the spike → tiny
   τ → negligible lag) — measured: 16.93 → 11.50 MPa, insufficient.
3. **Velocity-coupled response (King/Lengellé/Beddini)** — the only body
   of "dynamic erosive" work — is built for **combustion instability**
   (response to *oscillatory* cross-flow), on the same condensed-phase +
   gas-phase response functions. It is not an ignition-development lag.
   (Lengellé 1975, "Erosive Combustion and Velocity Response," AIAA J
   13(3), DOI 10.2514/3.49697 = `lengelle1975.pdf`.)
4. **The literature treats ignition-transient erosion as physically REAL**
   — "severe during the ignition transient, due to the small port area"
   for high-L:D motors. So the field would call our spike *expected*, not
   a bug. The real Chunc's *lack* of a spike is therefore geometry-specific
   (fore-tapering conical port → G is LOWER toward the aft, not
   concentrated) plus its erosive lift being a slow, seconds-scale rise
   (bore opening), per `srm-1d-chunc-erosive-spike-diagnostic`.
5. **Energy-side corroboration:** SPINBALL's energy equation (thesis Eq.
   1.3) sources combustion enthalpy `r_b·P_b·ρ_p·h_f` with **no convective
   wall-loss sink** — the wall heating for the ignition criterion is a
   separate conduction sub-model (their analog of our Goodman). So our
   v0.7.4 energy books — now *with* the convective + pyrogen-radiation
   sinks — are *more* complete than SPINBALL's standard practice, and the
   audit's conclusion (energy side is not the spike driver) is consistent
   with SPINBALL omitting the sink yet not spiking on VEGA (different
   regime: low bore mass flux, erosion "totally negligible").

## The one defensible new lever
A **turbulent-BL-development factor on the erosive Nu/h** (Beddini): the
erosive enhancement comes from core-flow turbulence penetrating the
combustion zone, which requires a *developed* turbulent boundary layer.
During the ignition transient the flow is still establishing, so the
actual erosive enhancement is **below** the quasi-steady (developed-BL)
Ma value. Ramping the erosive component by a development factor
`f_dev(t_since_flow_arrival)` rising over the BL-development time
(~L_e/u ~ ms) would suppress the ignition surge while leaving the steady
lift intact. Caveat: that timescale is ~ms (the same short scale that made
the convective wall sink and Z-N only partially effective), so expect
*modest* additional suppression — not a full fix to 8.5 MPa.

## Options for the user
- **(A) Implement the BL-development erosive factor** (Beddini-grounded;
  the shared h_c/Nu already exists). Test on Chunc; expect partial.
- **(B) Accept the F+Z residual** (2× → ~1.25×). The lit dive shows there
  is no literature-backed transient erosive closure to fully close it
  without ad-hoc tuning, and the literature considers ignition-transient
  erosion physically real; document Chunc's residual as a
  quasi-steady-erosive + geometry/simultaneity limitation.
- **(C) Revisit simultaneity** — but the ignition audit showed the gas
  fill is physically fast and slowing it (flame-front) is unphysical.

**Recommendation:** try (A) as the last physically-grounded lever, with
clear expectations that it's likely partial; if it underwhelms, take (B)
and ship F+Z (opt-in) + the energy-balance fix as the v0.7.4 result,
documenting the residual honestly.
