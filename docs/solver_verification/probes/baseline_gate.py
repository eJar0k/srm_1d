"""Baseline bit-identity gate (B0) for the solver-verification tasks.

Re-runs the three canonical ``ma2020`` / ``erosive_mp_cap=0`` rows of
``docs/v0_7_4/probes/mp_cap_crossmotor.py`` (Chunc, BALLSstick, Hasegawa A;
committed .ric files, canonical knobs, library BPNV igniter) and compares
them to the values recorded in ``docs/erosive_model_registry/WORKLOG.md``
(T2, 2026-08-08) at FULL float precision plus exact step counts.

A task whose spec says "default behaviour unchanged" passes B0 only if every
compared value is identical (``--rtol 0``, the default). Any difference means
the change perturbed the default trajectory -- report it, do not loosen the
tolerance to make it pass.

This script puts ITS OWN checkout first on ``sys.path`` so that, run from a
git worktree, it imports the worktree's ``srm_1d`` (not the editable install
of the main checkout). It prints the resolved package path as provenance.

Usage (from the checkout root, pyenv 3.10.5 interpreter):
    python docs/solver_verification/probes/baseline_gate.py
    python docs/solver_verification/probes/baseline_gate.py --only chunc
Exit code 0 = all identical, 1 = mismatch.
"""
import argparse
import importlib.util
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

# Recorded 2026-08-08 (erosive_model_registry WORKLOG, T2 results.jsonl).
BASELINE = {
    'chunc': dict(peak=11.472831280605046, t_peak=0.016454109339050758,
                  plateau=8.222465813350324, n_steps=531651),
    'balls': dict(peak=8.42348927652058, trough=4.366702806070186,
                  mid=4.913378178011638, n_steps=842922),
    'haseg': dict(peak=6.143912554793205, t_peak=2.3816074267208034,
                  n_steps=2335770),
}
ROWS = {'chunc': _probe.chunc_row, 'balls': _probe.balls_row,
        'haseg': _probe.haseg_row}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default='chunc,balls,haseg',
                    help='comma-separated subset of chunc,balls,haseg')
    ap.add_argument('--rtol', type=float, default=0.0)
    args = ap.parse_args()

    print(f"[provenance] srm_1d imported from {Path(srm_1d.__file__).parent}")
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
                same = int(got) == int(ref)
                rel = 0.0 if same else abs(int(got) - int(ref)) / ref
            else:
                rel = abs(got - ref) / abs(ref)
                same = rel <= args.rtol
            ok &= same
            print(f"  {key:8s} got {got!r:>24}  ref {ref!r:>24}  "
                  f"rel {rel:.3e}  {'OK' if same else 'MISMATCH'}")
    print("\nB0 RESULT:", "PASS (identical)" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
