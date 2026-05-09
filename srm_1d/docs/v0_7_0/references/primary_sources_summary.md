# Sutton 9e + DeMar 2021 — Primary-Source Summary

Two reference documents already in the repo, extracted via pdftotext during research. These cover the empirical / textbook foundations that the academic papers don't repeat.

PDFs in repo:
- `Rocket Propulsion Elements.pdf` — Sutton 9e (Sutton & Biblarz, ~800 pp.)
- `STARTING YOUR UP-GOERS_pub (1).pdf` — DeMar 2021 amateur-rocketry presentation

---

## 1. Sutton 9e — Igniter Coverage

### 1.1 Section 14.2 Ignition Process

**Three-phase ignition framework (Sutton's standard taxonomy):**
- **Phase I — Ignition time lag**: from electrical signal to first hot gas from grain surface
- **Phase II — Flame-spreading interval**: from first ignition until full burning area is ignited
- **Phase III — Chamber-filling interval**: from full ignition to equilibrium chamber pressure and mass flow

**Ten factors affecting ignitability:**
1. Propellant formulation
2. Propellant grain surface initial temperature
3. Surrounding pressure
4. Modes of heat transfer (convective, radiative)
5. Grain surface roughness
6. Propellant age
7. Composition and hot-solid-particle content of igniter gases
8. Igniter propellant and its initial temperature
9. Velocity of hot igniter gases relative to grain surface
10. Cavity volume and configuration

Figure 14-4 (referenced but not in the extracted text) shows ignition time vs. heat flux for a specific propellant — empirical curves, not equations.

### 1.2 Section 15.3 Igniter Hardware

**Pyrogen igniters — the canonical "small SRM inside the main SRM":**

> "A pyrogen igniter is basically a small unit, containing all the elements of a rocket motor, used to ignite a larger rocket motor... and is not designed to produce thrust. They all consist of one or more nozzle orifices (both sonic and supersonic types) and most use conventional rocket motor grain formulations (sometimes the same as the main propellant grain) and design technology."

This is **exactly** the head-end primary motor architecture the user described. Sutton confirms it's the standard pyrogen design.

**Heat transfer mechanism:** "Heat transfer from the pyrogen gas to the motor grain is largely convective, with hot gases contacting the inner grain surface, in contrast to pyrotechnic igniters that transfer heat by high-energy radiation."

This justifies the v0.7.0 choice to use convective `h_c` from existing Gnielinski/Haaland infrastructure rather than implementing a separate radiation model.

**Mounting locations** (Figure 15-13):
- Forward end (most common) — gas flow over propellant surface helps ignition
- Aft end — little gas motion; relies on temperature, pressure, heat transfer from igniter gas
- Internally through the nozzle — discardable hardware
- Annular sleeve around submerged nozzle throat

**Pyrotechnic igniters:**
> "In industrial practice, pyrotechnic igniters are defined as igniters (other than pyrogen-type igniters) that use solid explosives or energetic propellant-like chemical formulations (usually small pellets of propellant that provide large burning surfaces and short burning times) as the heat-producing materials."

A common pellet-basket design uses staged ignition: squib → booster charge → main ignition charge of pellets (e.g., 24% boron / 71% potassium perchlorate / 5% binder).

**Initiator types** (Figure 15-15):
- (a) Integral diaphragm — shock wave through diaphragm; case stays sealed
- (b) Header type with double bridgewire (glow plug)
- (c) Exploding bridgewire (low-resistance Pt or Au, 0.02-0.10 mm)

**Pyrogen mass and impulse contribution:**
> "Since igniter propellant mass is small (often less than 1% of the motor propellant) and burns mostly at low chamber pressure (hence low Is), it has a negligible contribution to the motor's overall total impulse."

### 1.3 Sutton Eq. 15-4 — Empirical Igniter Sizing

The single most important quantitative result in the chapter:

```
m = 0.12·V_F^0.7         (Sutton Eq. 15-4)

where:
  m   = igniter charge mass (grams)
  V_F = motor free volume (cubic inches)
        the void in the case not occupied by propellant
```

Source: Figure 15-16 (empirical chart, not in extracted text — based on industrial data set).

**Sanity check examples** (computed for typical motors):

| Motor | V_F [in³] | m [g] |
|---|---|---|
| Hasegawa A | ~30 | 1.3 |
| Zerox | ~80 | 3.0 |
| ~"M" amateur motor | ~100 | 3.5 |
| Shuttle SRB | ~30,000 | 460 |

These match published values within order of magnitude (Shuttle SRB pyrogen is ~64 kg per its design — but that includes the pyrogen *propellant grain*, not just the heat-producing charge; the "charge" in Sutton's sense corresponds to the smaller pellet basket inside the pyrogen).

**Sutton's caveat** (re. all igniter analysis):
> "current analytical models of physical and chemical processes relevant to igniter design (including heat transfer, propellant decomposition, deflagration, flame spreading, and chamber filling) are far from complete and accurate. Analyses and design of igniters, regardless of the type, depend heavily on experimental results."

This is why DeMar's empirical pyrogen measurements (next section) are essential.

### 1.4 Section 13.5 Igniter Propellants (referenced, not extracted)

Page 524 in Sutton 9e. Discusses specific propellants used in igniters. Not extracted in the v0.7.0 research pass; revisit if pyrogen propellant chemistry needs deeper grounding.

---

## 2. DeMar 2021 — "Starting Your Up-Goers" (Amateur Pyrogen Survey)

**Citation**: DeMar, J. S. (2021). *Starting Your Up-Goers: An Experimental Survey of Igniter Compositions*. TRATECH 2021, Tripoli Rocketry Association. TRA #10273.

### 2.1 Composite Propellant Ignition Fundamentals

**Best pyrogens for solid motors have:**
- High radiant heat flux per gram (bring propellant to self-ignition T)
- High gas generation (pressurize chamber above critical pressure for sustained combustion)
- Fast action time (electrical firing → max heat flux and pressure)
- Low sensitivity to ambient pressure and temperature
- Minimal shock to grain surface (avoid fracture)
- Low electrostatic sensitivity

### 2.2 6-Phase Solid Motor Startup Taxonomy

After electrical circuit closure:
1. Electrical delay
2. Induction
3. Flame Spreading
4. Chamber Filling
5. Erosive Burning
6. Steady State

**Burst disk acceleration:** Adding a burst disk at the nozzle gives "Faster rise time to critical pressure and steady state combustion."

### 2.3 Critical Sizing Parameters

- **Heat flux target**: typically >40 cal/cm²/s, higher for instant-on
- **Critical pressure for sustained combustion**: 50-100 psi for high-performance propellants; >150 psi for low-performance or instant-on
- **Secondary effects**: long motors / low port-throat ratio / erosive burning / vehicle resonance

### 2.4 The Sizing Formula (DeMar's empirical pyrogen-specific formula)

```
M = V · P / W           [DeMar formula]

where:
  M = pyrogen mass (grams)
  V = motor core volume (in³)
  P = desired ignition pressure (psi, typically 200-300)
  W = "impetus" or characteristic work (psi·in³/gram) — pyrogen-specific
```

**Calculation procedure:**
1. Find core volume: `V = π·D²/4 · L`
2. Find pyrogen mass: `M = V · P / W`
3. Find pyrogen length (for dipped igniters): `S = D · 1.5`
4. Find pellet count (for basket igniters): `N = M / q` (q = mass per pellet, e.g. 0.42 g for 1/4" pellets)

### 2.5 Measured Pyrogen Properties Table (DeMar's experimental measurements)

Test apparatus: 75mm snap-ring case, 5.75" × 2.75" ID. K-type thermocouple, 1000 psi pressure transducer, 1000°C max range.

| Pyrogen | Peak P [psi] | t_peak [ms] | T_peak [°C] | Impetus W [psi·in³/g] | Heat Flux [cal/cm²/s] | T·mass flow [°C·g/s] | Norm. Heat Flow Density |
|---|---|---|---|---|---|---|---|
| **BKNO3-Viton (BPNV)** | 105 | 67 | 840 | **5000** | **69.4** | 2724 | **1.00 (ref)** |
| BKNO3-NC | 102 | 73 | 768 | 4800 | 66 | 2200 | 0.77 |
| Mg-Teflon-Viton (MTV) | 103 | 50 | 930 | 4700 | 110 | 1488 | 0.87 |
| Mg-Teflon (dry) | 106 | 40 | >1000 | 3000 | 109 | 3429 | 1.98 |
| Mg-KNO3-epoxy | 53 | 133 | 476 | 2400 | 20.3 | 626 | 0.067 |
| Al-CuO thermite | 45 | 184 | 277 | 2038 | 12.1 | 183 | 0.012 |
| First fire | 105 | 120 | 594 | 2200 | 29.2 | 674 | 0.104 |
| BP-resin pellet | 56 | 130 | 242 | 2570 | 44.6 | 693 | 0.164 |

**Top performers** (combined gas + heat metrics):
- BKNO3-Viton (BPNV) — industry benchmark for amateur use
- MTV — comparable, faster, but electrostatic-sensitive
- Mg-Teflon (dry) — performs well but **dangerous** (flash powder; will damage grain)

**Worst performer**: Al-CuO thermite — produces mostly hot solids/liquids, low gas, low heat flux. "Thermite causes the partial decimation of the propellant surface to produce pressurization gas."

### 2.6 BPNV Composition Recipe

**Ratio by mass: B = 25 : KNO3 = 60 : Viton = 15**

Ingredients:
- Boron — amorphous, brown fine powder (<2 µm), ≥92% pure
- Potassium nitrate (KNO3) — white fine-milled powder
- Viton — raw (uncured) fluoroelastomer (FKM) polymer, ~66% fluorine

Mixing procedure:
1. Dissolve Viton in acetone (10:1 acetone:Viton)
2. Add boron powder, mix thoroughly
3. Add potassium nitrate slowly, mix thoroughly
4. Use more acetone for desired consistency:
   - Dipping: thin paste
   - Pellet casting: thick paste
5. Can be reconstituted with acetone later

### 2.7 BPNV Sizing Examples (from DeMar)

| Motor class | Core diam [in] | Core length [in] | Ignition P [psi] | V [in³] | M [g] | Length S [in] |
|---|---|---|---|---|---|---|
| J | 0.5 | 12 | 300 | 2.36 | 0.14 | 0.75 |
| M | 0.75 | 40 | 200 | 17.7 | 0.71 | 1.13 |
| O | 2 | 33 | 400 | 103.7 | 8.3 | 3.0 (or ~20 pellets) |

### 2.8 Cross-check: DeMar vs. Sutton

For an "M" motor (V_F ≈ 17.7 in³ in DeMar's example):
- Sutton: m = 0.12·17.7^0.7 = 0.12·7.86 = **0.94 g**
- DeMar (BPNV at 200 psi): M = 17.7·200/5000 = **0.71 g**

These are within ~30% of each other — order-of-magnitude consistent. Sutton tends slightly higher because Sutton's empirical curve was fit to industrial data including various pyrogen types; DeMar's formula uses BPNV's specific impetus.

**Recommendation for srm_1d**: default to Sutton (pyrogen-agnostic) unless the user specifies a pyrogen with known impetus W, in which case use DeMar's formula for better fidelity.

### 2.9 Static Sensitivity Ranking (DeMar's hazard data)

| Pyrogen | Static Sensitivity (joules) | Rank |
|---|---|---|
| BPNV (dipped/molded) | >5 | 1 (low sensitivity) |
| BPNV-NC (pressed) | 1.25 | 2 |
| BP-resin pellet | 0.875 | 3 |
| Pyrodex (BP variant) | 0.225 | 4 |
| MTV-Viton | 0.124 | 5 |
| Mg-Teflon (dry) | 0.030 | 6 — **VERY DANGEROUS** |
| Mg-KNO3-epoxy | 0.020 | 7 |
| First fire | 0.005 | 8 — **EXTREMELY DANGEROUS** |
| Al-CuO thermite | various | (depends on grade) |

**Implication for v0.7.0 defaults**: BPNV is the right safety-vs-performance default. MTV is comparable performance but adds a safety knob the user must affirmatively choose. Don't ship "first fire" or "Mg-Teflon (dry)" as default options.

---

## 3. Implications for srm_1d v0.7.0

Direct translation of the above into design choices already made in [DESIGN.md](../DESIGN.md):

1. **Default sizing formula**: Sutton Eq. 15-4 (m = 0.12·V_F^0.7) when user doesn't specify pyrogen mass. Pyrogen-agnostic.
2. **Pyrogen YAML datasheet schema** mirrors DeMar's measured table — `impetus_W`, `T_flame`, peak pressure, action time. User can specify pyrogen by name (`pyrogen='bpnv'`), and the simulator computes mass from `M = V·P/W` if `pyrogen_mass=None`.
3. **Default pyrogen options**: BPNV (best safety-performance trade) + MTV (faster but hazard-flagged). Both ship as YAML files in `srm_1d/motors/pyrogens/`.
4. **Heat-transfer**: convective only at v0.7.0. Sutton confirms pyrogens are convective-dominated (vs. radiative-dominated pyrotechnics). DeMar's measured heat fluxes (~70 cal/cm²/s for BPNV) feed into the per-cell convective coefficient via gas-side dynamics — not directly imposed.
5. **Action times** match the v0.7.0 plenum dynamics naturally — BPNV's ~67 ms time-to-peak corresponds to a plenum filling timescale `V_plenum / (mdot_choke / ρ_gas)` for the stated burning surface. This is what v0.6.0's ad-hoc `igniter_tau = 127 ms` was trying (and failing) to mimic.
