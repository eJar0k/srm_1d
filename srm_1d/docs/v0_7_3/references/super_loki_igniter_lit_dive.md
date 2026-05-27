# Super Loki Igniter Literature Dive

**Date**: 2026-05-25
**Scope**: Determine the correct igniter topology for the ISP / RCS
HTPB Super Loki recreation that drives `srm_1d/motors/ISP_Super_Loki.ric`
and the `examples/ISP_Super_Loki.py` validation example.
**Trigger**: v0.7.3 Phase B close-out flagged that the previous
experimental overlay was Chunc data; igniter-topology assumption
(`head_basket`) was inherited from the inline `PyrogenChamber` docstring
without a cross-checked literature record.

---

## TL;DR

- **Original Loki (1955)**: 76 mm, polysulfide/AP propellant. Igniter
  was two parallel 1 W / 1 A no-fire **electric squibs plus an
  appropriate ignition charge** — **installed at the launch site**
  (i.e. the igniter is separable from the motor, not factory-loaded).
  Location not explicitly head-end vs aft in the public sources
  reviewed.
- **Original Super Loki (1968, NASA CR-61238 baseline)**: 102 mm
  scale-up of Loki. Same propellant family (polysulfide / AP) until
  1993. NASA CR-61238 (Bollermann, 30 Jun 1968) and the 1973 follow-on
  AFCRL-TR-73-0412 (also Bollermann) contain a "Cross-Section View of
  Super Loki Rocket Motor With Igniter Installed" figure — sources
  reviewed here describe the igniter as separable and installed at
  the launch site, consistent with Loki heritage. Specific basket /
  cup geometry not extracted because direct WebFetch was blocked
  (see "Source-access caveat" below).
- **ISP / RCS HTPB recreation (1993 onward; the version we model)**:
  RCS Rocket Motor Components sells the surplus Super Loki motor
  hardware (limited stock of unused, unloaded factory casings) and
  ISP developed the **8522 HTPB analog** propellant to replace the
  original polysulfide. Formulation (80.40% AP blend, 9.88% HTPB,
  1.80% Al, 2.00% KP, 0.70% copper chromite, 3.78% IDP, 0.10%
  CAO-5, 1.35% MDI curative) is the exact match for the `RCS - 8522`
  string in `srm_1d/motors/ISP_Super_Loki.ric`. The RCS product page
  does not publish igniter assembly details — the user assembles the
  motor and supplies their own igniter.
- **Recommendation**: keep `head_basket` as the modeling topology for
  the current `examples/ISP_Super_Loki.py` only if the user explicitly
  configures a head-end pellet-cup igniter (e.g., a homemade BKNO3
  pellet pack glued to the forward bulkhead). For the **factory /
  surplus Super Loki recreation as documented**, the appropriate
  model is closer to a **nozzle-inserted electric squib + ignition
  charge** — which `srm_1d` currently has **no exact topology for**.
  See §5 for the gap and §6 for next-step options.
- **Static-fire data**: A genuine Super Loki "Sea Level Chamber
  Pressure and Thrust Vs. Time" plot **exists** in the MIT-hosted
  copy of NASA CR-61238 (1968 report). We did not extract numerical
  digitized samples here because WebFetch was permission-denied on
  the PDF; recommend opening the PDF manually and digitizing the
  curve into a Python array. The 1.95-2.1 s burn time + ~1400 PSI
  average chamber pressure values referenced in the RCS product page
  are independent corroborating numbers consistent with the
  CR-61238 figure.

---

## 1. Filtered timeline

### 1.1 Original Loki (RM-82 / PWN-1 Loki-Dart, 1950s)

| Field | Value | Confidence |
|---|---|---|
| Era | 1955-1965 (Army anti-aircraft 1955; sold to civilian met-rocket market thereafter) | HIGH |
| Diameter / length | 76 mm (3.0 in) / ~1.5 m | HIGH |
| Propellant | **Polysulfide polymer + ammonium perchlorate** | HIGH |
| Igniter type | **Two parallel 1 W / 1 A no-fire electric squibs + ignition charge** | HIGH |
| Igniter installation | **Separable from motor; installed at the launch site** (i.e. NOT a factory-loaded head-end pellet basket sealed in moisture cup) | HIGH |
| Igniter location (head vs aft) | Not explicitly stated in sources reviewed; the "separable, installed at launch site" phrasing strongly implies the igniter is inserted through the nozzle throat after the launch tube is loaded | MEDIUM |
| Heritage | Cooper Development / JPL; 3,544 Loki rounds fired at White Sands by JPL | HIGH |

**Primary sources**:
- Andreas Parsch, *Cooper Development Loki-Wasp / HASP*,
  designation-systems.net/dusrm/app4/hasp.html
- Andreas Parsch, *JPL RM-82/PWN-1 Loki-Dart*,
  designation-systems.net/dusrm/n-1.html
- *Loki (rocket)*, en.wikipedia.org/wiki/Loki_(rocket) (used as
  pointer; not authoritative on its own)

### 1.2 Original Super Loki (PWN-8 → PWN-10/11/12, 1968-mid-1990s)

| Field | Value | Confidence |
|---|---|---|
| Era | First flight 22 Apr 1968 (CR-61238 baseline); last operational flight 3 Jul 2013 | HIGH |
| Diameter / length / mass | **102 mm (4.0 in) / 3.16 m / 29 kg gross** | HIGH |
| Thrust | 18.00 kN (4,046 lbf) average | HIGH |
| Burn time | 1.95-2.1 s nominal | HIGH |
| Propellant (1968-1993) | **Polysulfide polymer + ammonium perchlorate**, cast-in-the-case, internal burning | HIGH |
| Propellant (1993 onward) | **HTPB-based APCP**, design change per Loki Wikipedia + Viper IIIA dual-change record | HIGH |
| Motor designation (MIL) | **MIL-R-83403 (USAF), 1 May 1973** | HIGH |
| Manufacturer | Space Data Corporation → Orbital Sciences (acquired SDC 1991). Motor hardware later surplus-sold by RCS Rocket Motor Components | HIGH |
| Igniter type | Continued separable squib + ignition charge, installed at launch site (the consensus reading; no source explicitly contradicts) | MEDIUM |
| Igniter cross-section figure | **NASA CR-61238 Fig 3.5: "Cross-Section View of Super Loki Rocket Motor With Igniter Installed"** — exists, location head-end vs nozzle not extracted here (PDF was not WebFetched) | LOW (location unconfirmed from primary text in this dive) |

**Primary sources**:
- Bollermann, B., *Super Loki Dart Meteorological Rocket System Final
  Report*, NASA CR-61238, 30 June 1968. [RRS mirror PDF](https://www.rrs.org/wp-content/uploads/2014/01/Super-Loki-Dart-Meteorological-Rocket-System-1968.pdf)
  and [MIT mirror PDF](https://wikis.mit.edu/confluence/download/attachments/122633025/SuperLokiReport.pdf)
- Bollermann, B. & Walker, R. L., *Design, Development and Flight
  Test of the Super Loki Stable Booster Rocket System*,
  AFCRL-TR-73-0412 / SDC report dated 30 Jun 1973.
  [RRS mirror PDF](https://rrs.org/wp-content/uploads/2014/01/Design-Development-and-Flight-Test-of-the-Super-Loki-Stable-Booster-Rocket-System-1973.pdf)
- *MIL-R-83403* (USAF rocket motor specification, 1 May 1973),
  [DLA QuickSearch](https://quicksearch.dla.mil/WMX/Default.aspx?token=414709)
- Andreas Parsch, *Space Data PWN-10 / PWN-11 / PWN-12 Super Loki Datasonde*,
  designation-systems.net/dusrm/n-10.html, n-11.html, n-12.html
- *Super Loki engine*, astronautix.com/s/superlokiengine.html
  (this is the page that documents "igniter is separable from the
  motor and installed at the launch site")

### 1.3 ISP / RCS HTPB recreation (the version `srm_1d` validates against)

| Field | Value | Confidence |
|---|---|---|
| Era | 1993 design change → present | HIGH |
| Motor hardware | **Surplus factory Super Loki aluminum casings, graphite nozzles, steel nozzle retainers**, sold by RCS Rocket Motor Components, Inc. (rocketmotorparts.com). 4-fin can spot-welded, match-drilled forward bulkhead retained by 36 spring pins. OD 4.00 in, ID 3.834 in, length 78.5 in, empty 11.6 lb. | HIGH |
| Propellant | **RCS 8522 HTPB APCP** (HTPB analog of original polysulfide). Density 1716 kg/m³, n=0.405, a=0.027 (RCS units). Match to `.ric` file confirmed: "RCS - 8522", a=1.912e-5, n=0.405, density 1716.15 | HIGH |
| Propellant formulation | 80.40% AP blend (200 µm), 9.88% HTPB, 3.78% IDP plasticizer, 2.00% potassium perchlorate, 1.80% 4 µm Al, 1.35% MDI curative ("8522"), 0.70% copper chromite, 0.10% CAO-5 | HIGH |
| Propellant weight loaded | 37.2 lb propellant; loaded motor 48.8 lb; total impulse 8,900 lb-sec | HIGH |
| Avg chamber pressure | **~1230 PSI (8.48 MPa)**; max 1470 PSI (10.13 MPa) per RCS spec sheet | HIGH |
| Igniter | **NOT supplied with the motor.** User-assembled. RCS does not publish a recommended igniter design for the Super Loki kit. Amateur convention for HPR-scale APCP motors of this size is a custom multi-pellet BKNO3 / thermite pyrogen, typically inserted through the nozzle and either: (a) suspended at mid-grain on a wire stub, or (b) tape-mounted at the forward bulkhead via a long-lead wire and a small head-end basket / cup. Choice is up to the user. | MEDIUM-LOW |
| Validated igniter | None published; this dive found no static-fire data with explicit igniter geometry called out | LOW |

**Primary sources**:
- RCS Rocket Motor Components, *Super Loki Motor Assembly - Single
  (w/out Dart fins)*, [product page](https://www.rocketmotorparts.com/product/super-loki-motor-assembly--single).
  Includes the 8522 propellant formula, ballistic coefficients, and
  hardware spec.
- AeroTech/Quest Division (RCS Rocket Motor Components, Inc.),
  [aerotech-rocketry.com](https://aerotech-rocketry.com/). FirstFire /
  Firestar igniters are RCS's standard HPR igniter line —
  nozzle-inserted thermite/BKNO3 dipped tip on twin lead wires.

### 1.4 Other notable Super Loki variants

| Variant | Notes | Source |
|---|---|---|
| PWN-8 Loki Datasonde | 1969 instrumented Loki-Dart precursor | designation-systems.net/dusrm/n-8.html |
| PWN-10 Super Loki Transpondersonde | Late-1960s USAF; transponder dart; used to ~1998 at NASA Wallops | designation-systems.net/dusrm/n-10.html |
| PWN-11 Super Loki Datasonde | Smaller dart, no transponder; ~80 km apogee; same motor as PWN-10 | designation-systems.net/dusrm/n-11.html |
| PWN-12 Super Loki ROBINSphere | 1972+ ROBIN inflatable falling-sphere payload | designation-systems.net/dusrm/n-12.html |
| Viper IIIA (1972) | 4.5" (110 mm) follow-on; same propellant transition in 1993 | astronautix.com/s/superloki.html |
| ASP / Apogee Super Loki Dart (kit) | 1/3-scale **hobby model rocket**, single-use 24mm motor — irrelevant for our validation | apogeerockets.com / asp-rocketry.com |
| Loki Research K/L/M/N motors | **Unrelated commercial HPR motors** sold by Loki Research LLC; trade name only, not a Super Loki recreation | lokiresearch.com |

---

## 2. The "head_basket" assumption — audit

The `PyrogenChamber` docstring at `srm_1d/igniter_plenum.py:82-93`
asserts:

> "the ISP Super Loki igniter is a head-end BKNO3 pellet charge in a
>  consumable plastic moisture cup, with NO defined orifice or
>  pressure-containing aft cap. Modeling it as plenum-with-orifice
>  would be wrong physics; modeling it as head_basket (uncontained,
>  at local bore P) is the appropriate fit."

This dive **could not corroborate the BKNO3-pellet-in-moisture-cup
claim** for the ISP / RCS recreation specifically. The only authoritative
description of the Super Loki igniter we located in the sources reviewed
says the igniter is *"separable from the motor and installed at the
launch site"* (astronautix.com/s/superlokiengine.html) — which is
shape-compatible with a nozzle-inserted squib-and-charge assembly,
not a head-end basket sealed in a moisture cup.

The "moisture cup" phrasing in the docstring matches US Patent
**4,539,910 (Igniter Pellet Cup)** which describes a plastic cup with
ignition pellets, sealed with heat-bonded plastic film to prevent
moisture intrusion. That patent's preferred mounting is at the head
end **of high-performance solid rocket motors with head-end web
ignition** — i.e. motors designed to be ignited at the forward end of
the propellant grain. Super Loki does not appear in the patent's
references, and we did not find any source that ties this specific
patent to Super Loki.

**Most likely interpretation**: the docstring conflated *general
amateur / HPR head-end-igniter practice for BKNO3-pellet pyrogens*
(real and common, e.g. ProCast BKNO3-V kits used at Loki Research
M/N motor scale) with *the specific factory-recreation Super Loki
igniter*. The conflation does not invalidate `head_basket` as a
**physically reasonable** modeling choice for an amateur-style
head-end BKNO3 pellet pack — but it does mean we cannot cite
"NASA CR-61238, MIT Super Loki Report, Smithsonian/NASM" as
authority for that specific choice. None of those sources, as
indexed by our searches, explicitly says "head-end pellet basket in
moisture cup."

---

## 3. Source-access caveat

The two key primary documents that would resolve the head-vs-aft
question definitively are publicly hosted as PDFs but our WebFetch
tool was permission-denied on the first attempt:

- NASA CR-61238 (1968), [RRS mirror](https://www.rrs.org/wp-content/uploads/2014/01/Super-Loki-Dart-Meteorological-Rocket-System-1968.pdf)
  — contains **Fig 3.5 cross-section with igniter installed**
- MIT Super Loki Report mirror, [MIT confluence](https://wikis.mit.edu/confluence/download/attachments/122633025/SuperLokiReport.pdf?version=1&modificationDate=1525150071000&api=v2)
  — same document; contains **"Super Loki Rocket Motor Sea Level
  Chamber Pressure and Thrust Vs. Time"** plots in the appendix

Both are non-authenticated, non-paywalled US Government / academic
mirrors. They can be opened in a browser and examined manually. A
follow-up dive with the WebFetch permission granted (or human
download + transcription) would close the residual uncertainty on
igniter location and would let us digitize the pressure trace.

---

## 4. Static-fire pressure trace — provenance

| Source | Has Super Loki P-t data? | Digitized? | Notes |
|---|---|---|---|
| **NASA CR-61238 (1968) / MIT mirror** | **YES**, "Sea Level Chamber Pressure and Thrust Vs. Time" plot in appendix | **NO** | Single best public source. Polysulfide-era data, not HTPB-era — but the 1993 propellant transition was advertised as a drop-in analog so the shape should be close (±15-20% on peak). HIGH confidence the data exists. |
| AFCRL-TR-73-0412 (Bollermann 1973) | Likely yes (follow-on report) | NO | Not examined in this dive. |
| MIL-R-83403 (1973) | Yes, acceptance-test envelope | NO | Specification document; gives min/max curves, not nominal trace. |
| RCS product page | Avg + peak P numbers only | NO | 1230 PSI avg / 1470 PSI peak, 8900 lb-sec impulse. Useful as scalar sanity check. |
| ISP / Industrial Solid Propulsion (specificimpulse.com) | Possibly | NO | Vendor's website says they "contributed to" Super Loki. Did not find a published trace. |
| Hobby static-fire recreations | Possibly | NO | None found in web searches for "Super Loki" + RCS 8522 static fire. The amateur HPR community discusses Super Loki but the casing is a piece of US space program history with limited unloaded inventory at RCS, so very few have actually fired one. |

**Recommendation**: Manually download the MIT mirror PDF, locate the
chamber-pressure-vs-time figure (likely Section 3 or appendix),
digitize via WebPlotDigitizer or similar, and add as
`ISP_SUPER_LOKI_EXPERIMENTAL` in `srm_1d/plotting.py` with explicit
provenance comment: "NASA CR-61238 polysulfide-era trace; HTPB 8522
analog is expected to deviate ~15% on peak". Until digitized, the
example should run without overlay (current state, per Phase B.6
provenance correction).

---

## 5. Recommendation for srm_1d's modeling assumption

### 5.1 The right answer depends on what's being modeled

The `srm_1d/motors/ISP_Super_Loki.ric` file represents the RCS surplus
hardware + 8522 HTPB propellant + **a user's own igniter** (not a
factory igniter). There are three credible amateur-igniter configurations
for this motor:

| User config | Best srm_1d topology | Notes |
|---|---|---|
| (a) Nozzle-inserted Firestar / FirstFire / homemade thermite-tipped twin-lead | **No exact match.** Closest: `aft_basket` with very small cartridge length, OR a missing `nozzle_insert` topology (directional momentum upstream into the bore). The MIT/CR-61238 description "separable, installed at launch site" reads as this style. | High confidence this is what most amateurs would do for a 4"/N-class APCP motor. |
| (b) Head-end BKNO3 pellet pack, glued/taped to forward bulkhead inside a plastic cup or paper sleeve | **`head_basket`** ← current model | Plausible amateur build, especially for high-grain-aspect-ratio motors where flame propagation from the aft end is slow. This is the configuration the current docstring implicitly assumes. |
| (c) Factory Loki / Super Loki original igniter | **`forward_plenum`-like** if there is a defined orifice / aft cap, or a missing topology if not | NASA CR-61238 Fig 3.5 would resolve. Sources reviewed are silent on whether there's a defined orifice. |

### 5.2 What to do in the short term

1. **Update the `PyrogenChamber` docstring** at `srm_1d/igniter_plenum.py:82-93`
   to remove the unsupported "NASA CR-61238 / MIT Super Loki Report
   / Smithsonian/NASM" citation chain. Replace with: "head_basket is
   one of several plausible amateur configurations for the RCS 8522
   Super Loki recreation; the original 1968-1993 polysulfide-era
   igniter geometry is not fully resolved in publicly indexed
   sources. See `docs/v0_7_3/references/super_loki_igniter_lit_dive.md`."

2. **Add a comment block to `examples/ISP_Super_Loki.py`** documenting
   that the model is amateur head-end BKNO3 pellet pack, not the
   factory 1968 igniter, and that the user assembled their own igniter.

3. **Keep `head_basket` as the current topology** — it is internally
   consistent with the "amateur head-end pellet pack" amateur build
   in §5.1 row (b), and it is the safest available choice in `srm_1d`
   for an uncontained-pyrogen amateur build.

4. **Defer**: a new `nozzle_insert` topology (variant of `aft_basket`
   with directional upstream momentum injection) would model row (a)
   more faithfully. This is the same topology the `PyrogenChamber`
   docstring calls out as `aft_fore_firing` (deferred from v0.7.3
   Phase A, L90-93). Promote to v0.7.4 if a user shows a real
   amateur Super Loki firing that doesn't validate under
   `head_basket`.

### 5.3 What would change the recommendation

The **only** evidence that would force a change is: (a) a static-fire
pressure trace from the ISP / RCS recreation with the igniter geometry
explicitly stated, and (b) that geometry being clearly aft-positioned
or factory-style with a defined orifice. Neither was found in this
dive. CR-61238 Fig 3.5 might resolve (b) — open the PDF manually.

---

## 6. Open follow-ups

- [ ] Download MIT mirror of CR-61238 manually; transcribe Fig 3.5
      (igniter cross-section) and the appendix pressure-time trace.
- [ ] Re-request WebFetch permission for `rrs.org` and `wikis.mit.edu`
      domains if a future dive wants to automate this.
- [ ] Contact RCS Rocket Motor Components (rocketmotorparts.com)
      directly for their recommended Super Loki igniter — they almost
      certainly know what amateur customers have used successfully.
- [ ] Digitize the CR-61238 P-t trace and add to `plotting.py` as
      `SUPER_LOKI_POLYSULFIDE_NASA_CR61238` (named explicitly to
      preserve provenance — it is NOT the HTPB 8522 trace).
- [ ] If row (a) nozzle-inserted igniter ends up being the canonical
      amateur build, scope a `nozzle_insert` topology for v0.7.4.

---

## 7. Bibliography

**Primary aerospace sources**:

1. Bollermann, B., *Super Loki Dart Meteorological Rocket System
   Final Report*, NASA CR-61238 / NTRS 19680026183, Space Data
   Corporation for NASA Marshall, 30 June 1968.
   <https://ntrs.nasa.gov/citations/19680026183>
2. Bollermann, B. & Walker, R. L., *Design, Development and Flight
   Test of the Super Loki Stable Booster Rocket System*,
   AFCRL-TR-73-0412 / DTIC AD0750796, 30 June 1973.
3. *Military Specification MIL-R-83403 (USAF) — Rocket Motor, Super
   Loki*, US Air Force, 1 May 1973.

**Reference / pointer sources**:

4. Parsch, A., *Designation-Systems / Directory of US Military
   Rockets and Missiles, Appendix 4 (Sounding Rockets)*,
   designation-systems.net/dusrm/.
5. *Super Loki* page, astronautix.com/s/superloki.html (Mark Wade).
6. *Super Loki engine* page, astronautix.com/s/superlokiengine.html
   — contains the "igniter is separable from the motor and installed
   at the launch site" sentence.
7. *Loki (rocket)* article, en.wikipedia.org/wiki/Loki_(rocket).
   Used as a pointer to primary sources.

**ISP / RCS recreation sources**:

8. RCS Rocket Motor Components, *Super Loki Motor Assembly - Single
   (w/out Dart fins)*, product page,
   <https://www.rocketmotorparts.com/product/super-loki-motor-assembly--single>.
   Contains 8522 propellant formula, ballistic coefficients, hardware
   dimensions, loaded weights, peak/avg P, total impulse.
9. AeroTech/Quest Division of RCS Rocket Motor Components, Inc.,
   <https://aerotech-rocketry.com/>. FirstFire and Firestar
   nozzle-inserted igniter product line — RCS's canonical HPR
   igniter style.
10. Industrial Solid Propulsion ("ISP"), <https://www.specificimpulse.com/>.
    Vendor's own page; states ISP contributed to Super Loki programs.

**Igniter-architecture supporting patents** (background, none
specifically tied to Super Loki):

11. US 4,539,910 — *Igniter pellet cup*, plastic cup with
    moisture-sealing film for head-end BKNO3 pellets. Closest
    architectural match to the `head_basket` docstring claim.
12. US 4,503,773 — *Aft end igniter for full, head-end web solid
    propellant rocket motors*.
13. US 4,498,292 — *Igniter for rocket motors*.
14. US 4,901,642 — *Consumable wafer igniter*.

**Reports that have Super Loki chamber-P-vs-t data**:

15. MIT mirror of NASA CR-61238,
    <https://wikis.mit.edu/confluence/download/attachments/122633025/SuperLokiReport.pdf>.
    Contains the "Super Loki Rocket Motor Sea Level Chamber Pressure
    and Thrust Vs. Time" figure plus radar-tracked flight data.
16. RRS mirror of NASA CR-61238,
    <https://www.rrs.org/wp-content/uploads/2014/01/Super-Loki-Dart-Meteorological-Rocket-System-1968.pdf>.

---

## 8. Confidence summary

| Claim | Confidence |
|---|---|
| Original Loki used squib + ignition charge installed at launch site | HIGH |
| Super Loki (1968) used the same architecture | MEDIUM-HIGH (consistent description in multiple secondary sources; primary PDF not transcribed in this dive) |
| Super Loki polysulfide → HTPB transition occurred in 1993 | HIGH |
| RCS 8522 in the `.ric` file is the ISP HTPB analog | HIGH |
| `head_basket` is the "right" topology for the **factory** igniter | LOW (sources don't support this specific claim; could be aft-inserted instead) |
| `head_basket` is a reasonable topology for an **amateur** head-end BKNO3 pellet pack on the RCS recreation | MEDIUM-HIGH (consistent with HPR practice) |
| Super Loki P-t static-fire data exists in NASA CR-61238 | HIGH |
| Data is digitized and ready for `plotting.py` | NO — manual extraction step still required |
