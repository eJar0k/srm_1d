# static_fire_data/ — experimental validation traces

Measured static-fire traces used to validate the simulator's pressure/thrust
history. **Primary consumer right now: the ignition-spike investigation**
(`docs/v0_7_4/SPIKE_REOPEN_BRIEF.md`), whose goal is to match empirical data
**across the full short→long L/D range** — so the value of this folder is
*coverage across motors*, not just one trace. Also feeds the example overlays.

---

## Formats

Three forms coexist. Prefer **(1) CSV** for new drops.

**(1) Cleaned CSV — the drop-in format.** Loaded by
`srm_1d.plotting.load_experimental_csv(path, ...)` (does unit conversion).
Header + columns as in `pathfin54.csv`:
```
time(s),force(N),pressure(Pa)
0.0000,80.1472,90148.5969
...
```
Note: **pressure in Pa, force in N, time in s** here (the loader converts).
`force` is optional if you only have pressure — keep the column, leave blank or
0 if absent.

**(2) In-code digitized dict** — `plotting.py` `<MOTOR>_EXPERIMENTAL`, used
directly as a plot overlay:
```python
{'label': 'Experimental (…)',
 'time': np.array([...]),      # seconds
 'pressure': np.array([...]),  # MPa  (note: MPa here, not Pa)
 'time_offset': -0.3}          # optional s — shifts exp ignition to sim t=0
```
These are hand-digitized (coarse, ~30–60 pts). Good enough for eyeballing a
trace overlay; **not** ideal for the quantitative L/D-matching work — prefer a
high-resolution CSV where one exists.

**(3) Raw source (Excel / transducer dumps)** — multi-column instrument data
(e.g. `Zerox Data.xlsx`). Must be cleaned/extracted to a CSV (form 1) before
use. Put raw, unprocessed dumps in [`raw/`](raw/); commit the cleaned CSV at the
top level.

---

## Inventory (what exists + state)

| Motor | L/D | Data here | Form | State | Notes |
|---|---|---|---|---|---|
| **Hasegawa A** | ~42 (high) | `HASEGAWA_MOTOR_A_EXPERIMENTAL` (plotting.py) | dict, 36 pt, MPa | clean (digitized) | Ma 2020 Fig. 10; the canonical validation. Its peak is the *late progressive* peak, not an ignition spike. |
| **Zerox** | mid | `Zerox Data.xlsx` + `ZEROX_EXPERIMENTAL` (plotting.py, 58 pt, `time_offset=-0.3`) | Excel (raw) + dict | raw + digitized | Risky Batman V3 (fwd-Finocyl + aft-BATES). Excel is the full trace; extract a hi-res CSV from it. |
| **Chunc / machbusterNew** | high | `CHUNC_EXPERIMENTAL` (plotting.py, 59 pt, ~8.8 MPa plateau) | dict | clean (digitized) | The clean spike diagnostic (real = **zero** spike). **Hi-res 346-pt trace `ThomasMach5_edited.xlsx` is NOT in the repo yet (user Downloads)** — it made the G-reconstruction possible; **committing it here (as CSV) is the single highest-value add** for the spike work. |
| **Pathfinder 54** | — | `pathfin54.csv` (409 pt, force+pressure) | CSV | clean | The CSV format exemplar. Run via `examples/run_pathfin54.py`. |
| **54-2800 ST2.0** | — | *(expected `…54-2800 ST2.0.csv`, MISSING)* | CSV | absent | `examples/run_pathfin54st2_0.py` + `motors/2025.12.25 54-2800 ST2.0.ric` exist but the referenced CSV isn't here. |
| **BALLSStick** | high | *(none)* | — | **not in repo** | **No 3″-scale static fire.** A **2″ subscale** firing exists as **raw pressure-transducer data** (user-held) that needs (a) cleaning and (b) **subscale→3″ scaling** before it can validate the 3″ `BALLSstick.ric`. Drop the raw dump in `raw/`, document the scaling. |

**L/D coverage note:** the spike is ~L/D-proportional, so the set needs *both*
ends. High-L/D is well-covered (Hasegawa A, Chunc, BALLSStick-when-added);
**short/low-L/D coverage is thin** — worth adding a short motor's trace so a
"fix" can be shown not to *break* the low-L/D (near-zero-spike) case.

---

## Adding a trace

1. **Cleaned data →** drop a `<motor>.csv` (form 1) at the top level.
2. **Raw dump →** put it in [`raw/`](raw/), then clean → CSV. Record any
   processing (filtering, scaling, transducer calibration, subscale factors) in
   a sibling `<motor>.notes.md` so the provenance is auditable — especially the
   **BALLSStick 2″→3″ scaling**, which must be justified, not eyeballed
   (see `feedback_no_unfounded_smoothing`).
3. **Quick overlay only →** a digitized `<MOTOR>_EXPERIMENTAL` dict in
   `plotting.py` is fine, but note it in the table above.

## Known issues (pre-existing, flagged not fixed)
- `examples/run_pathfin54.py` points at `…/pathfin54.csv/` (**trailing slash** —
  won't load as a file).
- `examples/run_pathfin54st2_0.py` points at a nested `…/srm_1d/srm_1d/
  static_fire_data/2026.01.03 54-2800 ST2.0.csv` (wrong depth + missing file).

*(Do not hand-edit `.ric` motor files — those are openMotor's format;
this README and the data files are fair game.)*
