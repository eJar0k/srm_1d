# Head-end heat transfer: Peretz Eqs. III-28 .. III-30 (note for a follow-up task; no code)

Origin: SM1 Part 4. The user's point: a 1-D model cannot represent the turbulent mixing and recirculation at a real
head end, which likely deliver hot gas to the surface much faster than the 1-D bore flow does. SM1 shows (SM1_REPORT.md,
Part 3.1) that the head-end dead zones are not a solid-model artifact but a gas-side result: the local h collapses there
because the 1-D head end stagnates. This note records what the only head-end heat-transfer data behind Peretz's closure
say and what a follow-up would need.

## 1. The equations, read from the page image (thesis `peretz1974.pdf`, p. 40, Eqs. III-28 .. III-30)

Cohen, Derr and Price (Peretz ref. 158), "examining carefully" the Rocketdyne igniter-duct data of Carlson, Seader and
Wrubel (Peretz refs. 50-53), "proposed the following correlations for the convective heat transfer from an axial igniter
with a sonic nozzle":

```
At the region of igniter-jet impingement (maximum heat transfer):
  Nu_d = 0.0615 Re_d^0.8 (Ap/At)^0.4                                              (III-28)

Upstream of the impingement point:
  Nu_d = 0.7 Re_d^0.8 (Ap/At)^0.4 (0.022 + 0.066 x/x_im)                          (III-29)

Downstream of the impingement point:
  Nu_d = 0.7 Re_d^0.8 (Ap/At)^0.4 [0.031 + 0.057 exp(-0.5 x_r/d)]                 (III-30)
```

Consistency checks (done here, arithmetic only): at x = x_im, III-29 gives 0.7 x 0.088 = 0.0616 = III-28; at x_r = 0,
III-30 gives 0.7 x 0.088 = 0.0616, so the three pieces join continuously. Far downstream III-30 tends to 0.7 x 0.031 =
0.0217, close to the Dittus-Boelter constant 0.023. Relative to Dittus-Boelter the multiplier on Nu is therefore about 0.67
at x = 0, rising linearly to 2.7 at the impingement point and decaying back to ~0.94 downstream (Peretz: the maximum is
"up to 6 times" the Dittus-Boelter value of III-21; the factor 6 is not reproduced by III-28 and this is not resolved here).

Peretz's text after the equations (p. 40): "the data in Refs. 50 through 53 extend only to x/d = 10 approximately. Eq. (III-30)
gives almost a constant value for Nu_d after that point." He then fixed a6 = 2.56e-2 in his own bore correlation III-22/III-31 so
that III-31 equals III-30 "over the propellant leading edge and up to x/d_h = 10". The impingement position x_im is
"found by the jet-boundary calculations presented by Love, et al. (ref. 147); for most of the experiments ... approximately
at one port diameter downstream of the igniter nozzle exit, i.e., x_im ~ d_h".

## 2. How Peretz used them

1. Only indirectly for the bore heat transfer: III-30 calibrated the constant a6 and the (Ap/At)^0.4 factor of III-31, the
   correlation actually used in the surface-heating ODE (and in this repo, `_peretz_heat_transfer_coeff`, IG2/IG3).
2. Directly for the **effective mean igniter-gas temperature T_ig** of his entrance (head-end) section: "Equations (III-28 to
   III-30) are used to calculate the effective mean temperature of the igniter gas, T_ig." The igniter gas loses 5 % to
   incomplete combustion and heat transfer inside the igniter, then a further **25-35 %** (by Ap/At) in the high-heat-transfer
   entrance section, so T_ig is the adiabatic flame temperature minus those losses (his run D-9, Ap/At = 1.5: T_ig = 2340 K).
   T_ig is a prescribed constant during the transient and drives the entrance-section pressure/temperature ODEs (III-14,
   III-15); sensitivity: +10 % T_ig shortens the ignition time 24 % and the flame-spreading time 14 %.
3. Inputs the correlations need: the igniter-jet quantities, not only bore flow -- Re_d (see section 4), Ap/At, the position
   x/x_im (x_im ~ one port diameter, a jet-boundary calculation), and x_r/d (distance from the point of reattachment of the
   igniter jet, in port diameters). An igniter nozzle diameter does not appear in III-28..30 as quoted (the Carlson & Seader
   figures carry it as the parameter Dp/Di, see below).

## 3. Set against Carlson & Seader 1967 (summary in `R9_SOURCES.md`, commit 358d83a; PDF not re-read for this note)

Carlson & Seader: hot N2 at 570-620 psig, 561-783 K, uniform copper ducts Dp = 1, 2, 4 in, Re_Dp 5e4-6e5, Ap/At 1.1, 2, 3
(Ap = "port or test-section flow area", Ap/At a "secondary variable" through the jet structure). Their head-end correlation is
Dittus-Boelter, Nu = 0.023 Re^0.8 Pr^0.4 with Re built from the **igniter mass flow** and Dp, multiplied by a **Nu ratio read
from graphs** (Figs. 8-10, vs X/Dp, parameter Dp/Di, one figure per Ap/At, two times), valid with or without impingement; the head
end is either jet impingement (maximum h at impingement) or a disintegrating free jet with periodic recirculation zones.
Nothing in the paper covers the head end after the igniter burns out.

Are III-28..30 "the analytical fit of these figures"? Plausibly their reduction, but it is not proven here:
- Same structure: Dittus-Boelter-like Re^0.8 scaling, a multiplier that rises to a peak at impingement and relaxes downstream
  (III-29 / III-30 have the shape of a Nu-ratio-vs-X/Dp curve), an Ap/At dependence (Peretz's 0.4 power; Carlson & Seader call Ap/At secondary).
- Differences: III-28..30 contain no Dp/Di parameter and a single Ap/At power law, whereas the figures are plotted per Ap/At and
  per Dp/Di; the figures are for two times, III-28..30 for none. So III-28..30 look like Cohen/Derr/Price's collapse of the
  figure family onto one law, not a transcription. Peretz himself calls them proposals made after "examining carefully these
  data" and says the data extend to x/d ~ 10 only.
- To settle it, the figures (Carlson & Seader Figs. 8-10) must be read against III-29/III-30 at the same Ap/At, Dp/Di and X/Dp, and
  Cohen/Derr/Price (CPIA 220 vol. II, 1971; not found, R9) is needed for the Re_d and x_r definitions. Peretz's own text does not define
  Re_d beyond "the dimensionless parameters ... calculated using the local flow properties" for his bore law III-22; whether Re_d in
  III-28..30 uses igniter or local mass flow is NOT established from Peretz. By Carlson & Seader's description (igniter mass flow) it
  would be the igniter's.

## 4. Why this matters for the dead zones (observation, not a result)

- BALLSstick cells 3-11 (x = 36-117 mm, D_port ~ 13 mm) are at x/d ~ 3-9: inside the x/d <= 10 range of the duct data and
  downstream of the impingement point (x_im ~ d, i.e. cell 1), so III-30 is the regime of those cells. Hasegawa A cells 3-10 sit
  at x = 59-176 mm. Peretz's III-31 (what the sim uses) was fixed to equal III-30 over that range **but with Re built from the LOCAL bore
  flow**. In the 1-D sim the local flow at those cells collapses (u ~ 0.3-6 m/s after the downstream grain lights; h falls 5-15x), so
  III-31 gives a small h. If III-30's Re_d is built from the **igniter** flow (the Carlson & Seader statement), the same cells would see
  a head-end h set by the igniter jet and the recirculation, independent of the 1-D bore velocity -- which is what a real head end does and
  what the 1-D flow cannot represent. That is the gas-side gap the user pointed at, and the only data-based route to a physically
  sourced fix.
- Limits: the data are during igniter operation only (nothing after burnout), Ap/At 1.1-3 (Chunc 0.56, BALLSstick 0.27-0.58 are below
  it), N2 / uniform ducts, Dp 1-4 in.

## 5. Suggested follow-up task

1. Obtain Cohen/Derr/Price (CPIA 220 II) or AFRPL-TR-67-267 for the Re_d / x_r / x_im definitions; compare III-28..30 with Carlson & Seader
   Figs. 8-10 digitised at Ap/At 1.1, 2, 3 (decides the fit question above).
2. Add an opt-in head-end h term for cells x/d <= ~10 (and x/x_im), with Re_d from the igniter mass flow while the igniter burns, switched
   off (or blended to the local-flow law) after burnout; no new tuned constants -- only III-28..30 as published. Test: does it light the
   BALLSstick / Hasegawa A dead zones, and does it change the Chunc ignition front (it will raise h near the head by up to ~2.7x DB)?
3. Decide T_ig: Peretz's effective T_ig (adiabatic minus 5 % and 25-35 %) is a driving-temperature input the sim does not have (it uses local gas T).
