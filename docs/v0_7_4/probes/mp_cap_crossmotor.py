"""Cross-motor differential test + Mukunda-Paul erosive-ceiling validation
(2026-07-16 session; REOPENED §10) + erosive-model-registry validation
(2026-08-08 session; docs/erosive_model_registry/, T2).

THE measurement behind the erosive_mp_cap probe flag AND the erosive-model
registry's cross-motor validation gate (DESIGN.md G4). Three parts:

(1) DIFFERENTIAL TEST (no code change needed): the per-motor startup
    over-pressure ratio, sim vs real, across the two high-L/D motors with
    static-fire data. Finding: the sim gets the motor-to-motor ORDERING
    right (fore-tapered Chunc lowest; BALLSstick genuinely spikes) and
    over-predicts BOTH by a near-uniform factor (1.27-1.36x) that matches
    each motor's own Ma-vs-Mukunda-Paul excess at its spike g (1.26-1.33),
    while Ma/MP is only 1.03-1.05 at trough/mid-burn.

(2) CAP VALIDATION: run all three motors with erosive_mp_cap in
    {0.0, 1.0, 1.1} (erosive_model='ma2020' throughout) and compare startup
    ratios + plateau absolutes vs experiment. See REOPENED §10 for the
    results table.

(3) MODEL-AXIS VALIDATION (erosive_model_registry T2, DESIGN.md G4): the
    same three motors run across the full erosive-model registry —
    (erosive_model, erosive_mp_cap) in {('ma2020', 0.0), ('ma2020', 1.0),
    ('ma2020', 1.1), ('mukunda_paul', 0.0), ('saint_robert', 0.0)} — via
    MODEL_AXIS below. `saint_robert` rows are the erosive-off decomposition
    tool for the spike/plateau decoupling investigation (REOPENED §10.5 /
    §11.4): they are RECORDED, not gated. `erosive_mp_cap > 0` is a
    ma2020-only option (run_simulation raises otherwise) — MODEL_AXIS never
    pairs a nonzero cap with a non-ma2020 model.

CAVEATS: (a) run against the COMMITTED motors/*.ric. As of the 2026-08-08
T2 session the working tree still has `machbusterNew.ric` re-saved with
transportVariant: effective (spike-inflating; likely accidental) and a
re-geometried `BALLSstick.ric` (both `M` in `git status`, same as the
2026-07-16 session that first flagged this). `_committed_ric_path()` below
detects this automatically via `git status --porcelain` and extracts the
committed HEAD copy to `C:/tmp/<stem>_committed.ric` when needed, printing
which source was used — no manual `git show` step required. `hasegawa_a.ric`
is unmodified and is routed through the same helper for a uniform header
note (it resolves to the plain `motors/` path).
(b) BALLSstick's real ratios are from the 2" SUBSCALE raw DAQ (ratio-only
usable, absolute P not trustworthy); sim is the 3" motor. (c) real ratios:
Chunc 1.027 (peak/plateau0.4-0.9 of thomas_chunc_firing.csv), BALLSstick
~1.50 peak/trough, ~1.35 peak/mid-burn(1.03s).

Run from repo root with the pyenv 3.10.5 interpreter. Each row is one full
transient sim (Chunc ~2-4 min, BALLSstick ~3-6 min, Hasegawa A ~4-8 min);
the full 5-row model-axis sweep across all three motors is 15 runs,
~45-75 min wall clock — see docs/erosive_model_registry/WORKLOG.md T2 for
how that was chunked into separate interpreter invocations with a JSON
accumulator (not part of this script; ad hoc driver, not committed).
"""
import subprocess
import numpy as np
from pathlib import Path

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import load_experimental_csv

ROOT = Path(__file__).resolve().parents[3]
MOTORS = ROOT / "motors"
TMP_DIR = Path("C:/tmp")

CANON = dict(roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None,
             T_ignition=756.0, print_interval=0.0, verbose=False)

# Model axis for the erosive-model registry cross-motor gate (DESIGN.md G4).
# (erosive_model, erosive_mp_cap). erosive_mp_cap is ma2020-only —
# run_simulation raises ValueError if paired with a non-ma2020 model, so no
# row here does that. The old (cap-only) sweep is the ma2020 sub-axis of
# this list; kept as separate constant for callers that only want the A/B
# cap comparison without the model axis.
CAP_AXIS = (0.0, 1.0, 1.1)
MODEL_AXIS = (
    ('ma2020', 0.0),
    ('ma2020', 1.0),
    ('ma2020', 1.1),
    ('mukunda_paul', 0.0),
    ('saint_robert', 0.0),
)


def _committed_ric_path(name):
    """Resolve `motors/<name>` to a path safe for validation runs.

    If `git status --porcelain` shows the file modified in the working
    tree, extracts the committed `HEAD` copy to
    `C:/tmp/<stem>_committed.ric` and returns that path (re-extracting
    every call is cheap and keeps this session-independent). Otherwise
    returns the plain `motors/<name>` path. Either way, prints a one-line
    note identifying which source was used, so probe output is
    self-documenting about provenance.
    """
    rel = f"motors/{name}"
    status = subprocess.run(
        ['git', 'status', '--porcelain', '--', rel],
        cwd=str(ROOT), capture_output=True, text=True, check=True,
    )
    if status.stdout.strip():
        TMP_DIR.mkdir(parents=True, exist_ok=True)
        dest = TMP_DIR / f"{Path(name).stem}_committed.ric"
        with open(dest, 'w', encoding='utf-8', newline='') as f:
            subprocess.run(['git', 'show', f'HEAD:{rel}'], cwd=str(ROOT),
                            check=True, stdout=f)
        print(f"[ric-source] {rel}: MODIFIED in working tree -> using "
              f"committed HEAD extract at {dest}")
        return dest
    print(f"[ric-source] {rel}: unmodified -> using working-tree copy")
    return MOTORS / name


def run(ric, erosive_model, cap, t_max, p_cutoff=0.01e6):
    res = run_from_ric(str(ric), erosive_model=erosive_model,
                        erosive_mp_cap=cap, t_max=t_max,
                        P_cutoff=p_cutoff, **CANON)[0]
    summary = res['summary']
    return (np.asarray(res['time']), np.asarray(res['P_head']) / 1e6, summary)


def chunc_row(erosive_model='ma2020', cap=0.0):
    ric = _committed_ric_path("machbusterNew.ric")
    t, P, summary = run(ric, erosive_model, cap, 1.0)
    k = int(P.argmax())
    pl = float(np.median(P[(t >= 0.4) & (t <= 0.9)]))
    return dict(peak=float(P[k]), t_peak=float(t[k]), plateau=pl,
                ratio=float(P[k]) / pl,
                erosive_model_stamp=summary['erosive_model'],
                termination=summary.get('termination'),
                n_steps=summary.get('steps'))


def balls_row(erosive_model='ma2020', cap=0.0):
    ric = _committed_ric_path("BALLSstick.ric")
    t, P, summary = run(ric, erosive_model, cap, 2.0)
    ms = t < 0.3
    kp = int(np.argmax(P[ms])); Ppk, tpk = float(P[ms][kp]), float(t[ms][kp])
    after = (t > tpk) & (t < 1.2)
    Ptr = float(np.min(P[after]))
    Pmid = float(P[int(np.argmin(np.abs(t - 1.03)))])
    return dict(peak=Ppk, trough=Ptr, mid=Pmid,
                pk_tr=Ppk / Ptr, pk_mid=Ppk / Pmid,
                erosive_model_stamp=summary['erosive_model'],
                termination=summary.get('termination'),
                n_steps=summary.get('steps'))


def haseg_row(erosive_model='ma2020', cap=0.0):
    ric = _committed_ric_path("hasegawa_a.ric")
    t, P, summary = run(ric, erosive_model, cap, 10.0, p_cutoff=0.05e6)
    k = int(P.argmax())
    return dict(peak=float(P[k]), t_peak=float(t[k]),
                erosive_model_stamp=summary['erosive_model'],
                termination=summary.get('termination'),
                n_steps=summary.get('steps'))


if __name__ == '__main__':
    exp = load_experimental_csv(
        str(ROOT / 'static_fire_data' / 'thomas_chunc_firing.csv'),
        time_col=0, pressure_col=1, pressure_unit='MPa', skip_header=1)
    te, Pe = np.asarray(exp['time']), np.asarray(exp['pressure'])
    pl_exp = float(np.median(Pe[(te >= 0.4) & (te <= 0.9)]))
    print(f"Chunc experimental: peak {Pe.max():.2f}, plateau {pl_exp:.2f}, "
          f"ratio {Pe.max()/pl_exp:.3f}")

    for erosive_model, cap in MODEL_AXIS:
        c = chunc_row(erosive_model, cap)
        print(f"\n--- erosive_model = {erosive_model!r}, erosive_mp_cap = {cap} "
              f"(stamp: {c['erosive_model_stamp']!r}) ---")
        print(f"  Chunc: peak {c['peak']:.2f} @ {c['t_peak']*1e3:.0f} ms, "
              f"plateau {c['plateau']:.2f} ({c['plateau']/pl_exp:.2f}x exp), "
              f"ratio {c['ratio']:.3f} (real 1.027) "
              f"[termination={c['termination']}, n_steps={c['n_steps']}]")
        b = balls_row(erosive_model, cap)
        print(f"  BALLSstick: peak {b['peak']:.2f}, pk/trough {b['pk_tr']:.3f} "
              f"(real ~1.50), pk/mid {b['pk_mid']:.3f} (real ~1.35) "
              f"[termination={b['termination']}, n_steps={b['n_steps']}]")
        h = haseg_row(erosive_model, cap)
        print(f"  Hasegawa A: global peak {h['peak']:.2f} @ {h['t_peak']:.2f} s "
              f"(exp ~6.5 late progressive) "
              f"[termination={h['termination']}, n_steps={h['n_steps']}]")
