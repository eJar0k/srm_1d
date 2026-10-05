"""Zerox discriminator: cross-propellant test of the T2 plateau-decomposition
degeneracy (docs/erosive_model_registry/WORKLOG.md, T2 + "T2 review note").

T2 found that Chunc's saint_robert / mukunda_paul / ma2020 plateau
decomposition cannot distinguish "Ma is right at the plateau" from "Chunc's
Saint-Robert a/n calibration is ~10% low and Mukunda-Paul is right" -- the two
worlds are observationally IDENTICAL on Chunc alone (a ~10% a/n error and a
band-edge erosive lift are degenerate on one motor/one propellant).

Zerox is a DIFFERENT propellant (n=0.456 per its .ric propellant tab, vs
Chunc's; mid L/D; has its own static-fire trace,
`static_fire_data/zerox_data.csv`) run through the same erosive-model
registry axis:
  - if mukunda_paul's plateau deficit RECURS on Zerox too -> the deficit is
    CLOSURE-SYSTEMATIC (favors "Ma is right at plateau" / MP under-predicts
    universally, not just for Chunc's propellant);
  - if Zerox's plateau is fine (or over-predicted) under mukunda_paul while
    Chunc's plateau is short -> the Chunc deficit is CHUNC-SPECIFIC a/n debt
    (favors "Chunc's a/n needs a ~10% correction, then MP is the right
    closure").

Zerox's own calibration has known-stale igniter-era history
(`project_zerox_calibration_state` memory) -- part of this measurement is
establishing whether Zerox is even usable as a validation motor at all (see
`compute_verdict` / part C of the task brief: unfit if the ma2020 row's
mid-burn reference pressure or burn-duration proxy miss experiment by more
than +/-25%).

Four-row model axis (a trimmed sub-axis of `mp_cap_crossmotor.py`'s 5-row
MODEL_AXIS -- `(ma2020, cap=1.0)` is dropped here since T2 already
established it is bit-identical to `mukunda_paul` on Chunc/BALLSstick and
near-identical on Hasegawa A, so it carries no discriminating information
this sweep doesn't already get from the `mukunda_paul` row):
  (erosive_model, erosive_mp_cap) in
  {('ma2020', 0.0), ('ma2020', 1.1), ('mukunda_paul', 0.0),
   ('saint_robert', 0.0)}

CAVEATS: (a) run against the COMMITTED motors/*.ric. `motors/zerox.ric` is
unmodified in the working tree (confirmed via `git status --porcelain`) but
is still routed through `_committed_ric_path()` below for a uniform
provenance note across every run printed to stdout -- self-documenting,
mirrors `mp_cap_crossmotor.py`'s convention exactly (it resolves to the
plain `motors/` path, same as that script's `hasegawa_a.ric` row).
`motors/machbusterNew.ric` (Chunc, used for the part-D cross-check plots) IS
modified in the working tree (re-saved with `transportVariant: effective` --
the same anomaly WORKLOG T2 flagged) so its row is always served from the
`git show HEAD:...` extract.
(b) `zerox_data.csv` has REVERSED columns -- pressure is column 0, time is
column 1 (`static_fire_data/README.md`). Loading with `load_experimental_csv`
defaults silently swaps time and pressure; every load site below passes
`time_col=1, pressure_col=0` explicitly.

Run from repo root with the pyenv 3.10.5 interpreter:
  "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" \\
      docs/v0_7_4/probes/zerox_model_axis.py [--row MOTOR MODEL CAP] [--plots-only]

Zerox rows are ~1-3 min each (t_max=8.2s against an ~8.1s real burn); Chunc
rows (t_max=1.0s) are ~10-15s each warm. The real 8-row sweep (4 Zerox + 4
Chunc) was executed row-by-row via a throwaway driver in C:/tmp (not
committed, same pattern as WORKLOG T2's execution harness) using this
script's `--row` flag, to keep each interpreter invocation well inside the
tool timeout. See docs/erosive_model_registry/WORKLOG.md for the actual
runtimes and full results table.
"""
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import load_experimental_csv

ROOT = Path(__file__).resolve().parents[3]
MOTORS = ROOT / "motors"
TMP_DIR = Path("C:/tmp/zerox_axis")
ART_DIR = ROOT / "artifacts" / "erosive_model_axis"
RESULTS_JSONL = TMP_DIR / "results.jsonl"

CANON = dict(roughness=32e-6, kappa=0.44, pyrogen='bpnv', pyrogen_mass=None,
             T_ignition=756.0, print_interval=0.0, verbose=False)

# The 4-row model axis this measurement needs -- see module docstring for why
# (ma2020, cap=1.0) is dropped relative to mp_cap_crossmotor.py's 5-row
# MODEL_AXIS. erosive_mp_cap is a ma2020-only option (run_simulation raises
# otherwise), so it is only ever nonzero paired with 'ma2020' here.
ROWS = (
    ('ma2020', 0.0),
    ('ma2020', 1.1),
    ('mukunda_paul', 0.0),
    ('saint_robert', 0.0),
)

# Experimental comparators (Zerox) -- given in the task brief; cross-checked
# in-script against static_fire_data/zerox_data.csv by load_zerox_experimental
# + zerox_metrics (see __main__ / print_verdict).
ZEROX_EXP_PEAK = 3.994        # MPa @ 0.43 s
ZEROX_EXP_REF_1_2 = 3.355     # MPa, median t in [1.0, 2.0] s
ZEROX_EXP_RATIO = 1.19        # peak / ref_1_2
ZEROX_EXP_BURN_DURATION_BAND = (6.4, 6.6)   # s
ZEROX_EXP_BURN_DURATION_CENTER = 6.5        # s -- the part-C verdict target

# Chunc experimental comparator (plateau, median P over t in [0.4, 0.9] s) --
# from WORKLOG T2 / thomas_chunc_firing.csv.
CHUNC_EXP_PLATEAU = 8.6611  # MPa


def _committed_ric_path(name):
    """Resolve `motors/<name>` to a path safe for validation runs.

    Identical helper to `mp_cap_crossmotor.py`'s `_committed_ric_path` --
    duplicated here (not imported) so this probe is self-contained per the
    existing probes' convention. If `git status --porcelain` shows the file
    modified/untracked in the working tree, extracts the committed `HEAD`
    copy to `C:/tmp/zerox_axis/<stem>_committed.ric` and returns that path;
    otherwise returns the plain `motors/<name>` path. Either way prints a
    one-line provenance note.
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


def _run(ric, erosive_model, cap, t_max, p_cutoff):
    res = run_from_ric(str(ric), erosive_model=erosive_model,
                        erosive_mp_cap=cap, t_max=t_max,
                        P_cutoff=p_cutoff, **CANON)[0]
    summary = res['summary']
    return (np.asarray(res['time']), np.asarray(res['P_head']) / 1e6, summary)


def _tag(erosive_model, cap):
    return f"{erosive_model}_cap{cap:g}"


def zerox_metrics(t, P):
    """Startup peak (t<0.8s) + its time; ref_1_2 = median P, t in [1,2]s;
    startup ratio = peak/ref_1_2; burn-duration proxy = last t where
    P > 0.5 MPa."""
    startup = t < 0.8
    k = int(np.argmax(P[startup]))
    peak, t_peak = float(P[startup][k]), float(t[startup][k])
    ref_mask = (t >= 1.0) & (t <= 2.0)
    ref_1_2 = float(np.median(P[ref_mask]))
    burn = t[P > 0.5]
    burn_duration = float(burn.max()) if burn.size else float('nan')
    return dict(peak=peak, t_peak=t_peak, ref_1_2=ref_1_2,
                ratio=peak / ref_1_2, burn_duration=burn_duration)


def chunc_metrics(t, P):
    """Global peak + its time; plateau = median P, t in [0.4,0.9]s; ratio =
    peak/plateau. Mirrors mp_cap_crossmotor.py's chunc_row metric shape."""
    k = int(np.argmax(P))
    peak, t_peak = float(P[k]), float(t[k])
    plateau_mask = (t >= 0.4) & (t <= 0.9)
    plateau = float(np.median(P[plateau_mask]))
    return dict(peak=peak, t_peak=t_peak, plateau=plateau, ratio=peak / plateau)


def zerox_row(erosive_model='ma2020', cap=0.0):
    ric = _committed_ric_path("zerox.ric")
    t, P, summary = _run(ric, erosive_model, cap, t_max=8.2, p_cutoff=0.05e6)
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    np.save(TMP_DIR / f"zerox_{_tag(erosive_model, cap)}.npy",
            np.stack([t, P], axis=1))
    m = zerox_metrics(t, P)
    m.update(erosive_model_stamp=summary['erosive_model'],
              termination=summary.get('termination'),
              n_steps=summary.get('steps'))
    return m


def chunc_row(erosive_model='ma2020', cap=0.0):
    ric = _committed_ric_path("machbusterNew.ric")
    t, P, summary = _run(ric, erosive_model, cap, t_max=1.0, p_cutoff=0.01e6)
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    np.save(TMP_DIR / f"chunc_{_tag(erosive_model, cap)}.npy",
            np.stack([t, P], axis=1))
    m = chunc_metrics(t, P)
    m.update(erosive_model_stamp=summary['erosive_model'],
              termination=summary.get('termination'),
              n_steps=summary.get('steps'))
    return m


def load_zerox_experimental():
    exp = load_experimental_csv(
        str(ROOT / 'static_fire_data' / 'zerox_data.csv'),
        time_col=1, pressure_col=0, pressure_unit='MPa', skip_header=1)
    return np.asarray(exp['time']), np.asarray(exp['pressure'])


def load_chunc_experimental():
    exp = load_experimental_csv(
        str(ROOT / 'static_fire_data' / 'thomas_chunc_firing.csv'),
        time_col=0, pressure_col=1, pressure_unit='MPa', skip_header=1)
    return np.asarray(exp['time']), np.asarray(exp['pressure'])


def _append_jsonl(record):
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_JSONL, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record) + '\n')


def run_all_rows():
    """Run all 4 Zerox + 4 Chunc rows sequentially, appending each result to
    RESULTS_JSONL as it completes. Reference sequential path; the actual
    session execution instead drove this row-by-row via `--row` from a
    throwaway C:/tmp chunking script (see module docstring)."""
    for erosive_model, cap in ROWS:
        r = zerox_row(erosive_model, cap)
        r.update(motor='zerox', erosive_model=erosive_model, cap=cap)
        _append_jsonl(r)
        print('zerox', erosive_model, cap, r)
    for erosive_model, cap in ROWS:
        r = chunc_row(erosive_model, cap)
        r.update(motor='chunc', erosive_model=erosive_model, cap=cap)
        _append_jsonl(r)
        print('chunc', erosive_model, cap, r)


def load_all_zerox_metrics():
    """Recompute metrics fresh from the cached per-row .npy traces (single
    source of truth -- avoids staleness against a possibly-re-appended
    RESULTS_JSONL)."""
    return {key: zerox_metrics(*_load_npy('zerox', *key)) for key in ROWS}


def load_all_chunc_metrics():
    return {key: chunc_metrics(*_load_npy('chunc', *key)) for key in ROWS}


def compute_verdict(zerox_rows):
    """zerox_rows: dict {(erosive_model, cap): metrics}. Verdict is computed
    off the ma2020(cap=0.0) row only, per task brief part C: Zerox counts as
    'radically divergent / unfit' if that row's ref_1_2 misses the
    experimental 3.355 MPa by more than +/-25%, OR its burn-duration proxy
    misses ~6.5 s by more than +/-25%."""
    ma = zerox_rows[('ma2020', 0.0)]
    ref_err = abs(ma['ref_1_2'] - ZEROX_EXP_REF_1_2) / ZEROX_EXP_REF_1_2
    dur_err = (abs(ma['burn_duration'] - ZEROX_EXP_BURN_DURATION_CENTER)
               / ZEROX_EXP_BURN_DURATION_CENTER)
    unfit = (ref_err > 0.25) or (dur_err > 0.25)
    return dict(ref_1_2=ma['ref_1_2'], ref_err=ref_err,
                burn_duration=ma['burn_duration'], dur_err=dur_err,
                unfit=unfit)


def print_verdict():
    v = compute_verdict(load_all_zerox_metrics())
    print("\n=== FITNESS VERDICT (task brief part C) ===")
    print(f"ma2020 zerox ref_1_2       = {v['ref_1_2']:.4f} MPa vs exp "
          f"{ZEROX_EXP_REF_1_2} MPa -> {v['ref_err']*100:.1f}% error "
          f"(threshold +/-25%)")
    print(f"ma2020 zerox burn_duration = {v['burn_duration']:.4f} s vs exp "
          f"~{ZEROX_EXP_BURN_DURATION_CENTER} s -> {v['dur_err']*100:.1f}% "
          f"error (threshold +/-25%)")
    print(f"VERDICT: {'RADICALLY DIVERGENT / UNFIT' if v['unfit'] else 'USABLE'}")
    return v


# ---------------------------------------------------------------- plotting -

_COLORS = {
    ('ma2020', 0.0): '#2a78d6',        # categorical slot 1 (blue)
    ('ma2020', 1.1): '#eb6834',        # categorical slot 2 (orange)
    ('mukunda_paul', 0.0): '#1baf7a',  # categorical slot 3 (aqua)
    ('saint_robert', 0.0): '#eda100',  # categorical slot 4 (yellow)
}
_LABELS = {
    ('ma2020', 0.0): 'ma2020',
    ('ma2020', 1.1): 'ma2020 (cap=1.1)',
    ('mukunda_paul', 0.0): 'mukunda_paul',
    ('saint_robert', 0.0): 'saint_robert',
}


def _load_npy(motor, erosive_model, cap):
    arr = np.load(TMP_DIR / f"{motor}_{_tag(erosive_model, cap)}.npy")
    return arr[:, 0], arr[:, 1]


def plot_zerox_full():
    te, Pe = load_zerox_experimental()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(te, Pe, 'k-', lw=1.8,
            label=f'experimental (peak {Pe.max():.2f} MPa)')
    for key in ROWS:
        t, P = _load_npy('zerox', *key)
        m = zerox_metrics(t, P)
        ax.plot(t, P, '-', lw=1.3, color=_COLORS[key],
                label=f"{_LABELS[key]} (ref[1,2] {m['ref_1_2']:.2f} MPa)")
    ax.set_xlim(0, 8.2)
    ax.set_xlabel('time (s)'); ax.set_ylabel('P_head (MPa)')
    ax.set_title('Zerox -- erosive-model axis vs experimental (full trace)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    ART_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(ART_DIR / 'zerox_models_full.png', dpi=130)
    plt.close(fig)


def plot_zerox_startup():
    te, Pe = load_zerox_experimental()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(te, Pe, 'k-', lw=1.8,
            label=f'experimental (peak {Pe.max():.2f} MPa @ '
                  f'{te[np.argmax(Pe)]:.2f} s)')
    for key in ROWS:
        t, P = _load_npy('zerox', *key)
        m = zerox_metrics(t, P)
        ax.plot(t, P, '-', lw=1.3, color=_COLORS[key],
                label=f"{_LABELS[key]} (peak {m['peak']:.2f} @ "
                      f"{m['t_peak']:.2f} s)")
    ax.set_xlim(0, 1.5)
    ax.set_xlabel('time (s)'); ax.set_ylabel('P_head (MPa)')
    ax.set_title('Zerox -- erosive-model axis vs experimental (startup, t<1.5s)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    ART_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(ART_DIR / 'zerox_models_startup.png', dpi=130)
    plt.close(fig)


def plot_chunc():
    te, Pe = load_chunc_experimental()
    pl_exp = float(np.median(Pe[(te >= 0.4) & (te <= 0.9)]))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(te, Pe, 'k-', lw=1.8,
            label=f'experimental (346 pt, plateau {pl_exp:.2f} MPa)')
    for key in ROWS:
        t, P = _load_npy('chunc', *key)
        m = chunc_metrics(t, P)
        ax.plot(t, P, '-', lw=1.3, color=_COLORS[key],
                label=f"{_LABELS[key]} (plateau {m['plateau']:.2f} MPa)")
    ax.set_xlim(0, 1.0)
    ax.set_xlabel('time (s)'); ax.set_ylabel('P_head (MPa)')
    ax.set_title('Chunc -- erosive-model axis vs experimental (cross-check)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    ART_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(ART_DIR / 'chunc_models.png', dpi=130)
    plt.close(fig)


def plot_plateau_decomposition():
    """Grouped bar: two groups (Chunc, Zerox) x three bars (saint_robert,
    mukunda_paul, ma2020) -- the a/n-vs-closure discriminator view (WORKLOG
    T2 review note)."""
    order = [('saint_robert', 0.0), ('mukunda_paul', 0.0), ('ma2020', 0.0)]
    chunc_vals = [chunc_metrics(*_load_npy('chunc', *k))['plateau'] for k in order]
    zerox_vals = [zerox_metrics(*_load_npy('zerox', *k))['ref_1_2'] for k in order]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    group_x = np.array([0.0, 1.0])
    width = 0.22
    offsets = (-width, 0.0, width)

    for i, key in enumerate(order):
        heights = [chunc_vals[i], zerox_vals[i]]
        xs = group_x + offsets[i]
        bars = ax.bar(xs, heights, width, color=_COLORS[key],
                       edgecolor='black', linewidth=0.6, label=_LABELS[key])
        for rect, val in zip(bars, heights):
            ax.annotate(f'{val:.3f}',
                        (rect.get_x() + rect.get_width() / 2, rect.get_height()),
                        ha='center', va='bottom', fontsize=8)

    span = 1.5 * width
    ax.hlines(CHUNC_EXP_PLATEAU, group_x[0] - span, group_x[0] + span,
               color='black', lw=1.5, ls='--')
    ax.text(group_x[0], CHUNC_EXP_PLATEAU, f' exp {CHUNC_EXP_PLATEAU:.3f}',
            va='bottom', ha='left', fontsize=8)
    ax.hlines(ZEROX_EXP_REF_1_2, group_x[1] - span, group_x[1] + span,
               color='black', lw=1.5, ls='--')
    ax.text(group_x[1], ZEROX_EXP_REF_1_2, f' exp {ZEROX_EXP_REF_1_2:.3f}',
            va='bottom', ha='left', fontsize=8)

    ax.set_xticks(group_x)
    ax.set_xticklabels(['Chunc\n(plateau, t in [0.4,0.9]s)',
                         'Zerox\n(ref, t in [1,2]s)'])
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel('P_head (MPa)')
    ax.set_title('a/n-vs-closure discriminator: plateau decomposition, Chunc vs Zerox')
    ax.grid(alpha=0.3, axis='y')
    ax.legend(fontsize=9, loc='upper left')
    fig.tight_layout()
    ART_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(ART_DIR / 'plateau_decomposition.png', dpi=130)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--row', nargs=3, metavar=('MOTOR', 'MODEL', 'CAP'),
                         help="run a single row (for chunked execution), e.g. "
                              "'--row zerox ma2020 0.0'; appends to "
                              "RESULTS_JSONL and returns without plotting")
    parser.add_argument('--plots-only', action='store_true',
                         help='skip running sims; regenerate the 4 PNGs from '
                              'cached C:/tmp/zerox_axis/*.npy')
    args = parser.parse_args()

    if args.row:
        motor, model, cap_s = args.row
        cap = float(cap_s)
        fn = zerox_row if motor == 'zerox' else chunc_row
        r = fn(model, cap)
        r.update(motor=motor, erosive_model=model, cap=cap)
        _append_jsonl(r)
        print(json.dumps(r, indent=2))
        return

    if not args.plots_only:
        exp = load_zerox_experimental()
        te, Pe = exp
        m_exp = zerox_metrics(te, Pe)
        print(f"Zerox experimental (in-script, reversed-column recipe): "
              f"peak {m_exp['peak']:.3f} MPa @ {m_exp['t_peak']:.2f} s "
              f"(brief: {ZEROX_EXP_PEAK} @ 0.43), ref_1_2 {m_exp['ref_1_2']:.3f} "
              f"MPa (brief: {ZEROX_EXP_REF_1_2}), ratio {m_exp['ratio']:.3f} "
              f"(brief: {ZEROX_EXP_RATIO}), burn_duration "
              f"{m_exp['burn_duration']:.3f} s (brief: "
              f"{ZEROX_EXP_BURN_DURATION_BAND})")
        run_all_rows()

    print_verdict()
    plot_zerox_full()
    plot_zerox_startup()
    plot_chunc()
    plot_plateau_decomposition()
    print(f"plots written to {ART_DIR}")


if __name__ == '__main__':
    main()
