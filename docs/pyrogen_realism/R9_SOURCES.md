# R9 — where the Rocketdyne igniter-duct data can be found (search 2026-10-11)

Peretz's ignition h constant (a6 = 2.56e-2) and its (Ap/At)^0.4 factor come from the
Carlson / Seader / Wrubel (Rocketdyne) hot-gas igniter-duct data, through the Cohen, Derr &
Price reduction (Peretz refs 50-53, 158). Delegated web search; DTIC Public Search was
offline ("offline while we enhance the capability"), arc.aiaa.org and apps.dtic.mil PDFs
returned 403 to the tool.

| Peretz ref | Document | Identifiers | Status |
|---|---|---|---|
| 52 | Carlson & Seader, "Heat-Transfer Characteristics of Hot-Gas Ignition", AIAA J. 5(7):1272-1279, July 1967 | **DOI 10.2514/3.4183 — verified** (Crossref record matches; doi.org redirects to arc.aiaa.org) | Paywalled (AIAA ARC). **Fetch this first.** |
| 51-ish | Carlson & Seader, "A study of heat transfer characteristics of hot-gas ignition", AIAA Paper 66-67, Jan 1966 | DOI 10.2514/6.1966-67 (Crossref only; unverified on the AIAA page) | Paywalled. Possibly the same content as CPIA 105 (unverified) |
| — | Jungowski W. M., Comment on "Heat-Transfer Characteristics of Hot-Gas Ignition", AIAA J. 1969, and the authors' reply | DOI 10.2514/3.55496 and 10.2514/3.55497 (Crossref/Semantic Scholar; unverified) | Paywalled; may discuss the correlation |
| 50 | Carlson & Seader, AFRPL-TR-65-158, June 1965 (final report, Contract AF04(611)-9884; the circular-duct "Ref. 1" of the follow-on) | No AD number found | Not online. When DTIC is back, search by contract number AF04(611)-9884 |
| 53 | Wrubel & Carlson, AFRPL-TR-67-267, July 1967 (final report, Contract AF04(611)-11613) | No AD number found | Not online. Search DTIC by contract AF04(611)-11613 |
| 51 | Carlson & Seader, ICRPG 2nd Combustion Conf., CPIA Publ. 105 Vol. 1, May 1966, pp. 563-597 | none | Not found (CPIA holdings) |
| 158 | Cohen, Derr & Price, "Study of Ignition and Extinguishment in Multiple Stop-Restart Duty Cycles", 8th JANNAF Combustion Mtg., CPIA Publ. 220 Vol. II, Nov 1971 | none | Not found. This is where the correlation Peretz used was derived |

Open but not useful for the correlation (same contract AF04(611)-11613, progress reports;
no h values, Nu fits or Ap/At exponents, checked by grep of the OCR text):
- Wrubel, 2nd Quarterly, Dec 1966, AFRPL-TR-66-358, AD0804893 — https://archive.org/details/DTIC_AD0804893
- Wrubel, 3rd Quarterly, Apr 1967, AFRPL-TR-67-83 (Rocketdyne R-6856-3), AD0811154, DOI 10.21236/ad0811154 — https://archive.org/details/DTIC_AD0811154

They state the correlations were "developed for cylindrical geometry in Ref. 1"
(AFRPL-TR-65-158) and that the follow-on extends them to **star and conocyl ports**, so
AFRPL-TR-67-267 is the most likely place where Ap is defined for a non-uniform port.

Checked locally: `Fundamentals_of_Solid-Propellant_Combustion.pdf` (Kuo & Summerfield 1984,
Ch. 6 Kumar & Kuo) only summarizes Peretz ("a modified Dittus-Boelter equation"); it
reproduces no duct data.

**Ap definition, as far as the sources go:** Ap = "cross-sectional area of the port"
(Peretz nomenclature). In every data set the port was uniform (Rocketdyne: a constant-area
copper-lined duct; Peretz: 2.54 x 0.635 cm along the whole motor, Ap/At set by four
throats), so local, head-end and aft-end Ap are the same there. Peretz attributes the
effect to igniter-jet expansion, a head-end mechanism; that leans toward the head-end
(jet-region) port area for a tapered bore, but no source states it.
