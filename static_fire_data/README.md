# static_fire_data/ — experimental validation traces

Measured static-fire traces used to validate the simulator's pressure/thrust
history. **Primary consumer: the ignition-spike investigation**
(`docs/v0_7_4/SPIKE_REOPEN_BRIEF.md`), whose goal is to match empirical data
**across the full short→long L/D range** — so the value here is *coverage across
motors*, not just one trace. Also feeds the example plot overlays.

---

## How to load these (there is no single schema — read this)

`srm_1d.plotting.load_experimental_csv` is **generic**: you tell it which
columns hold time and pressure and what unit the pressure is in, and it returns
`{'time' (s), 'pressure' (MPa), 'label'}`.

```python
load_experimental_csv(path, time_col=0, pressure_col=1,
                      pressure_unit='MPa',   # 'Pa' | 'kPa' | 'MPa' | 'psi'
                      skip_header=1, label=None)
```

**The files here do NOT share a column layout.** Use these exact recipes:

| File | Recipe |
|---|---|
| `thomas_chunc_firing.csv` | `time_col=0, pressure_col=1, pressure_unit='MPa'` |
| `zerox_data.csv` | **`time_col=1, pressure_col=0`**, `pressure_unit='MPa'` — ⚠️ **columns are REVERSED** (pressure first). Loading with the defaults silently swaps time and pressure. |
| `pathfin54.csv` | `time_col=0, pressure_col=2, pressure_unit='Pa'` (col 1 is force in N) |
| `raw/ballsstick_subscale_raw.csv` | **do not load directly** — raw DAQ, needs cleaning (see below) |

A third form exists in code: hand-digitized dicts in `plotting.py`
(`HASEGAWA_MOTOR_A_EXPERIMENTAL`, `ZEROX_EXPERIMENTAL`, `CHUNC_EXPERIMENTAL`) —
`{'label', 'time' (s), 'pressure' (MPa), 'time_offset'?}`, where `time_offset`
shifts the experimental ignition to sim `t=0`. These are coarse (30–60 pt);
**prefer the CSVs above where one exists.**

---

## Inventory (verified 2026-06-19)

| Motor | L/D | File / dict | Points | State | Notes |
|---|---|---|---|---|---|
| **Chunc / machbusterNew** | high | **`thomas_chunc_firing.csv`** (+ `.xlsx` source) | **346** | ✅ clean | t=[0.010, 2.353] s; plateau ≈ **8.77 MPa**; P_peak 8.894 @ 0.397 s → **spike ratio 1.015 = essentially NO ignition spike.** The high-res trace behind the G-reconstruction (REOPENED §9.1). **The key spike-diagnostic asset.** Supersedes the 59-pt `CHUNC_EXPERIMENTAL` dict. |
| **Zerox** | mid | `zerox_data.csv` + `Zerox Data.xlsx` (raw) + `ZEROX_EXPERIMENTAL` dict | 58 | ✅ clean (coarse) | t=[0, 8.113] s, P_peak 3.994 MPa @ 0.43 s. Same 58 pts as the dict (no extra resolution). The `.xlsx` is the full DAQ source — a higher-res trace could still be extracted from it. |
| **Hasegawa A** | ~42 (high) | `HASEGAWA_MOTOR_A_EXPERIMENTAL` (plotting.py) | 36 | ✅ clean (digitized) | Ma 2020 Fig. 10. Its peak is the *late progressive* peak, not an ignition spike. |
| **Pathfinder 54** | — | `pathfin54.csv` | 409 | ✅ clean | Has force (N) + pressure (Pa). |
| **BALLSStick (2″ subscale)** | high | **`raw/ballsstick_subscale_raw.csv`** | 417 | ⚠️ RAW: absolute P unusable — **but the SHAPE/ratio IS usable** | **Shows a REAL startup over-pressure of ~1.35–1.50×** (peak 1330 psi @ +0.070 s → trough 885 @ +0.300 s → mid-burn 985 @ +1.03 s) — vs Chunc's 1.015×. Ratio is offset/gain-robust, so it stands despite the calibration problem. **Absolute P not trustworthy:** no quiescent baseline (PSI ramps from sample 1 while volts stay flat → PSI isn't a linear map of the volts); **do not apply a constant offset**; load-cell impulse 0.85× expected; time on the DAQ power-on clock. **No 3″ firing exists**; needs justified 2″→3″ scaling. See `raw/ballsstick_subscale_raw.notes.md`. |
| **54-2800 ST2.0** | — | *(CSV missing)* | — | ❌ absent | `examples/run_pathfin54st2_0.py` + `motors/2025.12.25 54-2800 ST2.0.ric` exist; the CSV does not. |
| **short / low-L/D** | low | *(none)* | — | ❌ absent | Held for now (no clean data). **Still the key coverage gap** — a fix must be shown not to break the near-zero-spike low-L/D case. |

**L/D coverage:** high-L/D is now well covered (Chunc hi-res, Hasegawa A,
BALLSStick-once-cleaned). **Low-L/D remains uncovered.**

---

## Adding a trace

1. **Cleaned →** drop `<motor>.csv` at the top level; add its load recipe to the
   table above (columns/units are per-file, so this is mandatory).
2. **Raw →** put it in [`raw/`](raw/) and write a sibling
   `<name>.notes.md` recording provenance + every processing step (calibration,
   offsets, filtering, scaling). Justify scalings physically — do not eyeball
   (`feedback_no_unfounded_smoothing`).

## Known issues (pre-existing, flagged not fixed)
- `examples/run_pathfin54.py` points at `…/pathfin54.csv/` (**trailing slash** —
  won't load).
- `examples/run_pathfin54st2_0.py` points at a nested `…/srm_1d/srm_1d/
  static_fire_data/…` path (wrong depth + missing file).

*(`.ric` motor files are openMotor's format — never hand-edit. Data files and
this README are fair game.)*
