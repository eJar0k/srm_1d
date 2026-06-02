# v0.8.0 — openMotor Frontend Integration: Phase Plan (SKELETON)

> **Status: skeleton — phases are placeholders to be turned into real specs
> once [`DESIGN.md`](DESIGN.md) open questions are resolved.** Companion to
> DESIGN.md (scoping). Branch: `openmotor-frontend`, parallel to the v0.7.x
> re-LHS line on `v0.7.0-phase4`.

## Sequencing intent

Land the **data-model refactor (channels + units)** before the **GUI hook**,
because the GUI consumes the channel objects. Keep each phase independently
testable (the project's per-phase gate discipline). Default-OFF / additive
where possible until the GUI hook makes the break unavoidable.

## Phases (TODO — scope each before coding)

- **Phase 0 — Scoping/design.** Resolve DESIGN.md open questions: invocation
  model, interchange schema, channel-object shape, units boundary, what lives
  where. Output: a filled DESIGN.md + an interface sketch. *(This package is
  the start of Phase 0.)*
- **Phase 1 — Channel-object refactor.** Replace the bare results dict with
  named, unit-tagged channel objects (aligned to openMotor's results model).
  Map the `*_hist` diagnostics to channels (curated subset vs all — decide in
  Phase 0). Gate: existing examples/tests read results via channels;
  byte-for-byte numerics unchanged.
- **Phase 2 — Unit-aware + generic plotting.** `plot_channels(channels, ...)`
  replaces the bespoke per-figure functions; units carried on channels,
  converted at display. Gate: current pressure/thrust/summary figures
  reproduced via the generic path.
- **Phase 3 — openMotor GUI hook.** Invoke the srm_1d transient backend from
  openMotor (model per Phase 0) and render the result channels in openMotor's
  plot widgets alongside its quasi-steady output. Gate: a motor designed in
  openMotor produces an srm_1d transient trace in-GUI.
- **Phase 4 — Validation + docs.** End-to-end check on a canonical motor;
  DEVNOTES API-break log entries; close-out doc.

## Cross-line sync (important)

When **v0.7.5 (cross-motor re-LHS)** tags on `v0.7.0-phase4`/main, **merge or
rebase it into `openmotor-frontend`** so the frontend inherits the
recalibrated model. Cut **v0.8.0** only from a base that already contains
v0.7.5. Do not let this branch's API churn block the calibration line, and do
not let it drift from the validated physics.

## Gates / discipline (carried from prior packages)

- Delete `srm_1d/__pycache__/` (+ `.nbi/.nbc`) after any `@njit` edit (gotcha #1).
- Hard API breaks OK (v0.8.0) — log each in DEVNOTES "API Breaking Changes Log."
- Defer to openMotor data-structure conventions; units convert at the boundary.
- Never hand-edit `.ric`.
- `pytest srm_1d/tests/` green before each phase close.
