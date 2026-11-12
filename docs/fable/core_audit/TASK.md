# Task: blind audit of a solid-rocket-motor simulator's numerical core

You are auditing the numerical core of a 1-D transient simulator of
solid-rocket-motor internal ballistics. You get only three things:
- the code (`srm_1d/`);
- a statement of the intended physics (`SPEC.md`);
- the primary references for the closures it implements (`references/`).

That is deliberate. The developers have a long history with this code. You
should find what is wrong without inheriting their assumptions, their
explanations, or their sense of what has already been checked.

**Every comment, docstring, option name and changelog note in the code is an
unverified claim.** Many describe past bugs as fixed, or say a term is
"exact" or "conservative". Verify each claim from the code itself, the
references and first principles. A claim of correctness is not evidence.

## What the developers want from you

Defects that change computed pressure histories, ranked by how much. The
startup pressure transient matters most. After that, anything that makes
results depend on the grid or the time step, or that breaks conservation.
Correct-but-ugly code is low priority.

Equally valuable: **a list of what you checked and found correct, with the
evidence.** That lets the developers stop doubting it.

## Environment

- Windows. Python with numba, numpy, scipy, pyyaml and matplotlib:
  `C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe`.
  The system `python` on PATH lacks these packages; do not use it.
- Run everything from this folder (the bundle root) as modules, so that
  `import srm_1d` resolves to the copy here:
  ```
  <python> -m probes.my_probe
  ```
- **Numba cache:** the core is Numba-compiled and caches to `__pycache__`.
  - If you edit any `@njit` function, delete `srm_1d/__pycache__/` first, or
    you will run the old compiled code.
  - The first run of a session compiles for 1 to 3 minutes. Later runs are
    fast.
- **Example motors in `motors/`:**
  - `hasegawa_a.ric`: one BATES grain;
  - `machbusterNew.ric`: BATES plus a Conical grain;
  - `example_bates.ric`: four BATES segments.
- To run a motor:
  ```python
  from srm_1d.openmotor_adapter import run_from_ric
  result, perf, nozzle, geo, prop = run_from_ric(
      'motors/hasegawa_a.ric', target_propellant_cells=50, t_max=0.3,
      verbose=False)
  ```
  - `result` is a dict of time histories plus a `summary` dict.
  - Lower-level entry points are in `srm_1d/simulation.py`
    (`run_simulation`) and `srm_1d/grain_geometry.py`
    (`build_snapped_geometry`).
  - Build minimal cases yourself: a sealed vessel, a constant-area duct, a
    nozzle blowdown. Analytic answers beat motor runs for verification.

## Compute budget (the machine's owner is working on it at the same time)

- **Run every simulation through the guard:**
  ```
  <python> -m srm_1d.tools.guarded_run --mem-gb 4 --timeout-min 30 -- <python> -m probes.my_probe
  ```
  It caps memory for the whole process tree and enforces the timeout.
- Limits: at most 2 simulation processes at once, 4 GB each, 30 min per job.
- If you need more, **stop and ask**. Measure one small instance first; the
  guard prints `peak_job_gb` and `wall_s`.
- **Probes must be memory-bounded:**
  - no per-step `list.append` in long loops; preallocate or decimate;
  - every time-marching probe needs a step-count or dt-floor stop.
- Keep grids small (20 to 100 cells) and `t_max` short unless the question
  needs more.
- Before you finish, check that none of your processes are still running.

## Rules

- **Do not change `srm_1d/` on the `bundle` branch.** Write probes in
  `probes/`.
  - To demonstrate a fix, make it on a branch `fix/<slug>`, show the before
    and after numbers, and cite the commit in the report.
- Do not tune numbers to match anything. This is verification against physics
  and references, not calibration.
- **References:** use the PDFs in `references/`.
  - Extract text with `pdftotext file.pdf -` and grep it; view single pages
    as images only when an equation's typesetting matters.
  - If you need a paper that is not here, list it (citation, DOI marked
    verified or unverified, what you need from it). Do not reconstruct a
    paywalled paper's content from abstracts or secondary sources.
- Commit your work in this folder's git repo as you go.

## Usage-limit discipline

This session may hit usage limits and be resumed later, possibly by a
different session.
- Create `AUDIT_REPORT.md` within your first few steps.
- Keep a status checklist at its top; update it and commit after each area.
- A resumed session should be able to pick up from the checklist alone.
- Read code by targeted sections: `grep -n "^def "`, then specific line
  ranges. `srm_1d/simulation.py` is about 4800 lines; do not read it whole.

## Suggested order (highest leverage first)

1. **Conservation.**
   - Build ledgers: mass, species and total energy, integrated over the
     domain plus the integrated sources minus the nozzle outflow. Check them
     on:
     - a sealed vessel with sources;
     - a regressing port;
     - the lumped head and aft volumes;
     - the end-face injection.
   - Derive the moving-boundary source terms (SPEC section 2) and compare
     them with the code.
2. **The PISO step.**
   - Derive the compressible pressure-correction equation for this staggered
     grid and compare every coefficient, including the storage terms, both
     correctors and the energy-pressure coupling.
   - Check against exact solutions: a Sod shock tube in a constant-area duct,
     and a choked-vessel blowdown.
   - Check temporal and spatial convergence orders.
3. **The nozzle boundary.**
   - The stagnation state used, the choked/unchoked switch and reverse flow.
   - Whether the dṁ/dP in the pressure equation matches the flow function.
4. **Burn rate against Ma et al. 2020.**
   - Equation by equation: definitions, units and property evaluation.
   - The iterative solve: its bracket, convergence and failure modes.
   - Continuity of the tab lookup.
5. **Ignition and solid heating.**
   - Transcription of the Peretz correlations and how they are evaluated.
   - The recovery temperature.
   - Verify the conduction solver against analytic solutions: constant h
     (erfc / exp-erfc), and a step heat flux.
   - The ignition switch.
6. **Geometry.** Analytic BATES and Conical perimeter, area and web against
   closed forms; regression bookkeeping; burnout; snapping of lengths to the
   grid.
7. **Species and mixture.** Boundedness of the advection, the mixture
   properties, and any temperature or pressure clipping.
8. **Every clamp, floor, ceiling, cap and epsilon in the core.**
   - List each one.
   - Judge whether it can bind in a normal run, and what it violates when it
     does: conservation, physics or consistency.

## Deliverable: `AUDIT_REPORT.md`

1. **Status checklist:** each area marked done, partial or not started.
2. **Findings**, ranked by severity:
   - **S1:** wrong physics or a conservation violation that materially
     changes results on the default path;
   - **S2:** wrong only in some regimes (startup, transonic, burnout, reverse
     flow) or on a non-default option;
   - **S3:** latent or robustness issues;
   - **S4:** clarity.

   For each finding give:
   - `file:line`;
   - what the code does;
   - what it should do, with a derivation or the reference equation and page;
   - the demonstration: probe path and numbers;
   - a suggested fix;
   - the estimated effect on the head-end pressure history.
3. **Verified correct:** what you checked, how, and the evidence.
4. **Open questions:** places where the spec, the references or physics leave
   the right answer unclear.
5. **References you still need** (if any).

End your session with a short summary: the S1 findings, then the S2 findings,
then what remains unchecked.
