"""Baseline bit-identity gate (B0) for the solver-verification tasks.

Re-runs the three canonical ``ma2020`` / ``erosive_mp_cap=0`` rows of
``docs/v0_7_4/probes/mp_cap_crossmotor.py`` (Chunc, BALLSstick, Hasegawa A;
committed .ric files, canonical knobs, library BPNV igniter) and compares
them to the recorded values (see the BASELINE comment: re-recorded by T10) at FULL float precision plus exact step counts.

A task whose spec says "default behaviour unchanged" passes B0 only if every
float matches within ``--rtol`` (default 1e-5) and every step count within
``--steps-rtol`` (default 5e-4 relative). Why not bitwise: with
``fastmath=True`` the Numba machine code depends on the compile state (which
functions were compiled in-process vs loaded from the on-disk cache), so
IDENTICAL source gives slightly different results between a warm cache and a
clean recompile -- deterministically per state. History of the measured
noise floor:
  * before T10: 1.5e-16..1.3e-15 relative, exact step counts (2026-10-05);
  * after T10: up to 1.6e-10 (Chunc plateau), exact step counts (2026-10-06);
  * after T11 (energy-coupled iteration, default K=1): up to 2.4e-7 (Chunc
    peak) and step counts differing by up to 4.9e-5 relative (Hasegawa A,
    115 of 2.3M steps) (2026-10-11). The iteration's convergence test
    (max|dT|/T < 1e-6) decides whether the extra pass runs; last-digit
    differences flip that decision on marginal steps, and the adaptive time
    step carries the flip into the trajectory. A warm cache reproduces the
    recorded values exactly.
  * after IG3 (Peretz ignition h, 2026-10-12): a fresh recompile moved Chunc
    t_peak by 4.6e-5 relative (1.7 us, of the order of one time step) while every
    pressure stayed within 1e-5 and two warm runs reproduced the recorded
    values exactly. t_peak is the time of the largest SAMPLE (an argmax over
    steps), so it moves in whole steps like n_steps: it is compared with
    ``--steps-rtol``, not ``--rtol``.
1e-5 / 5e-4 keep ~40x / ~10x margin over that floor and stay far below the
smallest real perturbation seen (T10: -0.5% peak; T11: -11%). Anything beyond
them means the change perturbed the default trajectory -- report it, do not
loosen the tolerance to make it pass. Use ``--rtol 1e-9 --steps-rtol 0`` for a
same-compile-state A/B.

This script puts ITS OWN checkout first on ``sys.path`` so that, run from a
git worktree, it imports the worktree's ``srm_1d`` (not the editable install
of the main checkout). It prints the resolved package path as provenance.

Usage (from the checkout root, pyenv 3.10.5 interpreter):
    python docs/solver_verification/probes/baseline_gate.py
    python docs/solver_verification/probes/baseline_gate.py --only chunc
Exit code 0 = all identical, 1 = mismatch.
"""
import argparse
import hashlib
import importlib.util
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import srm_1d  # noqa: E402  (after the sys.path fix on purpose)

_spec = importlib.util.spec_from_file_location(
    "mp_cap_crossmotor", ROOT / "docs" / "v0_7_4" / "probes" / "mp_cap_crossmotor.py")
_probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_probe)

# RE-RECORDED 2026-10-11 by pyrogen-realism IG3 (stage 1): the default ignition
# heat-transfer closure became 'peretz' (Peretz et al. 1973, Eq. III-22, local
# port area) instead of the Ma/Gnielinski chain (ignition_heat_transfer='gnielinski').
# Only the first-light delay changes (Chunc 1.72 -> 5.81 ms, BALLSstick 1.75 -> 7.22,
# Hasegawa A 2.77 -> 7.69); the plateau is unchanged and the startup peak moves
# -1.3 % (Chunc) / -3.3 % (BALLSstick) / ~0 (Hasegawa A). Values below are the IG3
# results (docs/pyrogen_realism/results/ig3_B0_*.txt). The T11 values (the previous
# BASELINE; reproduce them with ignition_heat_transfer='gnielinski') were:
#   chunc: peak 10.08270103870215, t_peak 0.02791011521855539,
#          plateau 8.146643086031123, n_steps 531701
#   balls: peak 7.257654400712118, trough 4.343098504353659,
#          mid 4.868134347450266, n_steps 843077
#   haseg: peak 6.1072529740207635, t_peak 2.405585530808029,
#          n_steps 2336674
#
# (earlier) RE-RECORDED 2026-10-11 by solver-verification T11: the PISO second
# corrector now includes the storage the first corrector took
# (piso_corr2_storage=True; without it every transient pressure increment
# was applied ~twice) and the step iterates energy/pressure once
# (piso_energy_coupling_iters=1). Values below are the T11 results
# (C:/tmp/sv_t11/core_*.jsonl, worklog/T11.md). Chunc t_peak moving
# 16.3 -> 27.9 ms tripped T11's 15% startup-metric STOP rule; the reviewer
# ruled it expected (P5 measured the same shift before adoption: the
# doubled increments had made the startup pressurise too fast).
# Pre-T11 (T10, 2026-10-05) values:
#   chunc: peak 11.28243546264288, t_peak 0.016286760473606336,
#          plateau 8.150933877646011, n_steps 531883
#   balls: peak 8.331331016878805, trough 4.343222420159596,
#          mid 4.879629094113075, n_steps 843274
#   haseg: peak 6.108574454139173, t_peak 2.391606230911318, n_steps 2336140
# Pre-T10 (2026-08-08, erosive_model_registry WORKLOG T2) values:
#   chunc: peak 11.472831280605046, t_peak 0.016454109339050758,
#          plateau 8.222465813350324, n_steps 531651
#   balls: peak 8.42348927652058, trough 4.366702806070186,
#          mid 4.913378178011638, n_steps 842922
#   haseg: peak 6.143912554793205, t_peak 2.3816074267208034, n_steps 2335770
BASELINE = {
    'chunc': dict(peak=9.950855636512404, t_peak=0.036735041862973064,
                  plateau=8.141968067255977, n_steps=532381),
    'balls': dict(peak=7.020318627616664, trough=4.343390826697384,
                  mid=4.855316005870968, n_steps=843726),
    'haseg': dict(peak=6.107268310023869, t_peak=2.411042627819566,
                  n_steps=2336271),
}
RIC_FILE = {'chunc': 'machbusterNew.ric', 'balls': 'BALLSstick.ric',
            'haseg': 'hasegawa_a.ric'}
ROWS = {'chunc': _probe.chunc_row, 'balls': _probe.balls_row,
        'haseg': _probe.haseg_row}


def _print_provenance(names):
    """Output-only provenance (P6): git HEAD + dirty flag, library versions,
    and the sha256 of each .ric file this gate will actually run (resolved
    by the probe's own committed-copy logic). Never affects pass/fail."""
    import numpy
    import numba
    def _git(*a):
        try:
            return subprocess.run(['git', *a], cwd=str(ROOT), capture_output=True,
                                  text=True, check=True).stdout.strip()
        except Exception as exc:  # provenance must never break the gate
            return f"<unavailable: {exc}>"
    head = _git('rev-parse', 'HEAD')
    dirty = bool(_git('status', '--porcelain', '--untracked-files=no'))
    print(f"[provenance] git HEAD {head}  dirty={dirty}")
    print(f"[provenance] numpy {numpy.__version__}  numba {numba.__version__}  "
          f"python {sys.version.split()[0]}")
    for name in names:
        try:
            path = _probe._committed_ric_path(RIC_FILE[name])
            digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
            print(f"[provenance] {RIC_FILE[name]} sha256 {digest}  ({path})")
        except Exception as exc:
            print(f"[provenance] {RIC_FILE[name]} sha256 <unavailable: {exc}>")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default='chunc,balls,haseg',
                    help='comma-separated subset of chunc,balls,haseg')
    ap.add_argument('--rtol', type=float, default=1e-5,
                    help='float tolerance (relative)')
    ap.add_argument('--steps-rtol', type=float, default=5e-4,
                    help='relative tolerance on n_steps (0 = exact)')
    args = ap.parse_args()

    print(f"[provenance] srm_1d imported from {Path(srm_1d.__file__).parent}")
    _print_provenance([n.strip() for n in args.only.split(',')])
    ok = True
    for name in args.only.split(','):
        name = name.strip()
        t0 = time.time()
        row = ROWS[name]('ma2020', 0.0)
        wall = time.time() - t0
        print(f"\n[{name}] wall {wall:.1f} s, termination={row['termination']!r}")
        for key, ref in BASELINE[name].items():
            got = row[key]
            if key == 'n_steps':
                rel = abs(int(got) - int(ref)) / ref
                same = rel <= args.steps_rtol
            elif key == 't_peak':
                # argmax over time steps: discrete, so the step tolerance
                rel = abs(got - ref) / abs(ref)
                same = rel <= args.steps_rtol
            else:
                rel = abs(got - ref) / abs(ref)
                same = rel <= args.rtol
            ok &= same
            print(f"  {key:8s} got {got!r:>24}  ref {ref!r:>24}  "
                  f"rel {rel:.3e}  {'OK' if same else 'MISMATCH'}")
    print("\nB0 RESULT:",
          f"PASS (floats within rtol {args.rtol:g}, n_steps and t_peak within "
          f"{args.steps_rtol:g})" if ok
          else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
