# Fable task 1 — independent audit of the ignition over-spike

Paste this file's path into a NEW Claude Code session, started in the srm_1d repo
root, with the instruction: "Read docs/fable/SPIKE_AUDIT_PROMPT.md and do it."
The repo's CLAUDE.md and the project memory load automatically. This task is
meant to have that context.

---

You are auditing the foremost open problem in srm_1d, a 1-D transient
solid-rocket-motor simulator: **the startup (ignition) pressure spike is
over-predicted against static-fire data.** Months of work have eliminated many
causes and changed many defaults. Your job is to step back, judge the whole
record with fresh eyes, and produce a **ranked list of testable hypotheses,
each with the cheapest experiment that would separate it from the others.**

This is an audit, not an implementation task:
- Do not change solver code.
- Small probes are allowed, inside the compute budget below.
- Cheaper agents will run the experiments you design.

## Ground rules (the user's, non-negotiable)

- **No tuned constants and no smoothing added to make a residual vanish.**
  Every parameter must come from a cited source or from first principles.
- **Ma et al. 2020 stays the erosive-burning model of record.** Alternatives
  such as Mukunda-Paul are opt-in only, and the user decides any flip.
- **Paywalled literature:** never infer a paywalled paper's numbers from
  abstracts or secondary mentions. List what you need (citation, DOI marked
  verified or unverified, what to extract); the user fetches the PDFs into
  `docs/references/`. Check that folder first; it holds many papers.
- **Compute budget:** follow the CLAUDE.md "Compute budget" section. Run every
  simulation through `python -m srm_1d.tools.guarded_run --mem-gb 4
  --timeout-min 30 -- <cmd>`. Use at most 2 processes. Ask before anything
  heavier.
- Never hand-edit `motors/*.ric`. Use the pyenv 3.10.5 Python (see CLAUDE.md).

## Usage-limit discipline

Your session may hit usage limits, so work in phases.

- Create `docs/fable/SPIKE_AUDIT_REPORT.md` early, and save it after every phase.
- Read selectively: skim first, then go deep only where it matters.
- Prefer `grep` and targeted section reads over whole files.
- PDFs: extract text with `pdftotext` and grep it. Never read a whole PDF as
  images.

## Phase 1 — establish the current state (cheap)

1. Read, in this order (skim, then target):
   - `CLAUDE.md`, roadmap item 1 (current defaults and findings);
   - `docs/v0_7_4/SPIKE_REOPEN_BRIEF.md`;
   - `docs/v0_7_4/IGNITION_SPIKE_REOPENED.md`, sections 8 and 10-13.
     Section 10 holds the Ma vs Mukunda-Paul cross-motor result, measured
     before the recent default changes.
   - the reports, for their conclusions and tables (skim):
     - `docs/pyrogen_realism/`: `IG1_REPORT`, `IG2_REPORT`, `OB1_REPORT`,
       `SM1_REPORT`, `SM2_REPORT`, `HE1_REPORT`, `GS1_REPORT`;
     - `docs/solver_verification/worklog/`: `T11.md`, `EN1.md`, `RT1.md`,
       `MU1.md` (on branch `mu1/muscl`; `git show
       mu1/muscl:docs/solver_verification/worklog/MU1.md`), `D1.md`, `G1.md`,
       and `AS1.md` / `G2.md` if present.
2. Re-measure the current defaults on the canonical motors yourself, with short
   guarded runs. The defaults moved recently and some report numbers are stale.
   - Use `docs/solver_verification/probes/baseline_gate.py`, or
     `run_from_ric(...)` on `motors/machbusterNew.ric` (Chunc),
     `motors/BALLSstick.ric` and `motors/hasegawa_a.ric`.
   - Use the COMMITTED `.ric` versions (`git show HEAD:motors/<f>.ric`): the
     main checkout has the user's uncommitted edits to two of them.
   - Record startup peak, plateau, peak/plateau (Chunc), peak/trough and
     peak/mid (BALLSstick), Hasegawa A's early and late peaks, first light, and
     the auto charges.

## Phase 2 — scrutinize the evidence, including the data

Known facts to weigh:
- **Chunc static fire (`static_fire_data/thomas_chunc_firing.csv`):**
  - sampled every 7 ms (about 143 Hz);
  - the first sample is at 10 ms, already 8.47 MPa (97 % of plateau), so the
    first 10 ms is unrecorded;
  - experiment peak/plateau is 1.027;
  - the real igniter was 4.5-6 g of about 3 mm cubes of BKNO3 or MTV in a head
    basket, at a zero head gap.
- **Hasegawa 2006 (`docs/references/hasegawa2006.pdf`):** the igniter was
  "minimized"; the main motor lights after igniter burnout. The measured
  full-scale traces DID show ignition peaks, which the authors attributed to
  unsteady burning plus erosive burning and excluded from their erosive
  evaluation. Ma 2020 also excluded the startup from its error figures.
- **BALLSstick:** the experimental ratios (~1.50 / ~1.35) come from a 2-inch
  subscale motor. The full-scale motor is not fired.
- **The sim at the current defaults:**
  - the startup is grid-independent (G1);
  - with the physical zero head gap it is violent: Chunc 2.62 g peak/plateau
    about 2.0;
  - the over-prediction tracks each motor's Ma-minus-Mukunda-Paul excess
    (REOPENED section 10, pre-change).

Ask whether the comparison itself is sound: sampling, the window, the igniter,
transducer location, and what "the spike" even means in each data set.

Recent changes (2026-11-11), now the defaults. Account for them, but don't
spend the session rediscovering them:
- **NB1, nozzle boundary** (`worklog/NB1.md`). The outlet takes the true
  isentropic stagnation state; it used to take the last cell's static P and T.
  - The plateau grain-exit Mach is 0.11-0.13, so the plateau moves only
    0.7-1 %.
  - The startup moves 4-10 %: Chunc peak -5.6 %, BALLSstick -4 %, Hasegawa A
    early peak -9.5 %.
  - TP1 (`worklog/TP1.md`) moved thrust onto the nozzle-entrance P0.
    `P_head` stays the pressure trace. During the startup `P_head`/P0
    reaches 4.4 (BALLSstick), a sign of how steep the startup port-pressure
    gradient is.
- **AS1 + AW1, auto charge** (`worklog/AS1.md`, `results/aw1/`). The charge
  must give a prompt start: at least 95 % of grain cells lit within 0.1 s.
  - It steps up ×1.25 from a DeMar 300 psi guess, which the user treats as
    the amateur floor.
  - Ignition comes out all-or-nothing: 99-100 % lit by 9-17 ms, or 0-7 % lit.
  - Auto charges: Chunc 1.27 g, BALLSstick 2.03 g, Hasegawa A 9.83 g.
- **G2, grid** (`worklog/G2.md`). Grain snapping preserves mass, gap volumes
  are lumped, and lumps are spread across cells. Open: Aeroblast is 2 % off,
  unexplained.
- **B0 at the current defaults** (N = 100):
  - Chunc peak 11.09 MPa at 23 ms, plateau 8.16, pk/plateau **1.36**
    (exp 1.027);
  - BALLSstick peak 8.41, trough 4.33, mid 4.85: pk/trough 1.94,
    pk/mid 1.74 (exp ~1.50 / ~1.35);
  - Hasegawa A late peak 5.957 MPa (0.925× exp), early peak 5.33.
- **MU1, tabled by the user.** A MUSCL reconstruction branch exists. The user
  judged it a symptomatic fix, so it is not on the default path.

## Phase 3 — hypotheses and experiments

Cover at least these, and add your own:
- **Unsteady / transient erosive burning.** Boundary-layer development time
  versus the ignition timescale; Ma's quasi-steady closure applied to a
  transient mass flux. What closures exist without tuning? Give DOIs.
- **Ma's model outside its validated envelope.** Chunc's head port is smaller
  than its nozzle throat (Ap/At about 0.56 at the head), so the port chokes
  during startup and is near Mach 1. What does the literature say about
  erosive burning and heat transfer at near-sonic port flow?
- **Pressure-coupled transient burn rate** (Zeldovich-Novozhilov). It was
  tested opt-in in v0.7.4 under the OLD defaults; does that conclusion survive
  the new ones?
- **Igniter realism:** the real charges against the auto-sized ones; the basket
  model; pellet ejection.
- **Ignition model calibration:** T_ignition 756 K and k_solid 0.271 were
  fitted under old defaults; the LHS re-fit is pending.
- **1-D limits:** the head-end jet and recirculation; 3-D flow at the fore end.
- **Anything in the record that looks inconsistent or wrong.**

## Deliverable: `docs/fable/SPIKE_AUDIT_REPORT.md`

1. **Current state:** one page, with your re-measured numbers and the data
   caveats.
2. **Ranked hypotheses table**, one row each:
   - mechanism;
   - evidence for and against, with file and section references;
   - what would falsify it;
   - the cheapest discriminating experiment (configuration, motors, metrics,
     estimated compute);
   - the expected outcome under each competing hypothesis.
3. **Recommended next three experiments,** in order.
4. **Literature needed:** DOIs, verified or unverified, and what to extract.
5. **Inconsistencies or errors** you found in the record or the code.

Commit the report on a new branch, `fable/spike-audit`. Do not push. Finish
with a short summary of the top three hypotheses and the first experiment.
