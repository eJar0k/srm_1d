# static_fire_data/raw/ — unprocessed dumps

Landing spot for **raw, unprocessed** static-fire data: pressure-transducer
CSVs/dumps, instrument exports, oscilloscope traces, etc. — anything that still
needs cleaning (filtering, calibration, zeroing, resampling) or **scaling**
before it can validate a motor.

Workflow: raw here → clean/scale → commit a `<motor>.csv` at the parent level
(`time(s),force(N),pressure(Pa)`) → record the processing steps + provenance in
`../<motor>.notes.md`.

**Expected first tenant:** the **BALLSStick 2″ subscale** pressure-transducer
data — needs cleaning *and* a justified subscale→3″ scaling (Kn / c\* / L\*
reasoning, not eyeballed) before it can validate `motors/BALLSstick.ric`.
Document the scaling in the notes file.
