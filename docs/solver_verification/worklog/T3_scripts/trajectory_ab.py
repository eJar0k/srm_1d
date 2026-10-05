"""Bit-identity A/B tool for 'accounting only' changes (T3).

B0 only compares peak/plateau/n_steps of the 100-cell canonical runs, where the
1 kPa floor never fires. This tool dumps EVERY per-step history array (plus final
state) of a short Chunc run from any srm_1d tree and compares two dumps
byte-for-byte, so a change can be shown trajectory-identical in a regime where
the clamp is active (>= 400 cells).

    # dump (run once per tree, e.g. a `git archive HEAD srm_1d` copy and the worktree)
    PYTHONPATH=<pkg_root> python trajectory_ab.py dump <pkg_root> <out_dir> <tag> 100,400,800 [t_max]
    # compare
    python trajectory_ab.py cmp <a.npz> <b.npz> <label>

IMPORTANT (measured in T3): compare like with like -- two runs that both
JIT-compiled fresh, or two that both loaded the Numba cache. A fresh-compile run
and a cache-loaded run of the SAME source differ at the ~1e-16..1e-9 level
(fastmath compile-state noise), which this tool will (correctly) report.
"""
import sys
from pathlib import Path

import numpy as np

RIC = str(Path(__file__).resolve().parents[4] / "motors" / "machbusterNew.ric")


def dump(root, out, tag, cells, t_max):
    sys.path.insert(0, root)
    import srm_1d
    from srm_1d.openmotor_adapter import run_from_ric
    print("[provenance]", tag, srm_1d.__file__, flush=True)
    canon = dict(roughness=32e-6, kappa=0.44, T_ignition=756.0, pyrogen='bpnv',
                 pyrogen_mass=None, print_interval=0.0, verbose=False)
    for tpc in cells:
        res = run_from_ric(RIC, t_max=t_max, target_propellant_cells=tpc, **canon)[0]
        d = {k: res.raw[k] for k in res.raw.keys()
             if isinstance(res.raw[k], np.ndarray) and res.raw[k].dtype != object}
        d['__steps'] = np.array([res['summary']['steps']])
        np.savez(Path(out) / f"{tag}_{tpc}.npz", **d)
        print(tpc, "cells", res['summary']['cells'], "steps",
              res['summary']['steps'], flush=True)


def cmp(fa, fb, label):
    a, b = np.load(fa), np.load(fb)
    common = sorted(set(a.files) & set(b.files))
    bad = [k for k in common
           if a[k].shape != b[k].shape or a[k].tobytes() != b[k].tobytes()]
    mx = 0.0
    for k in bad:
        if a[k].shape == b[k].shape and a[k].dtype.kind == 'f':
            x, y = a[k], b[k]
            # denominator floored at 1e-300: 0-vs-1e-15 pairs report huge ratios
            mx = max(mx, float(np.max(np.abs(x - y) / np.maximum(np.abs(x), 1e-300))))
    only_a = sorted(set(a.files) - set(b.files))
    only_b = sorted(set(b.files) - set(a.files))
    print(f"{label}: {len(common)} common arrays, {len(bad)} byte-mismatched, "
          f"max rel diff {mx:.2e}" + ("  IDENTICAL" if not bad else "")
          + (f"  only-in-A={only_a}" if only_a else "")
          + (f"  only-in-B={only_b}" if only_b else ""))
    return bad


if __name__ == '__main__':
    if sys.argv[1] == 'dump':
        t = float(sys.argv[6]) if len(sys.argv) > 6 else 0.05
        dump(sys.argv[2], sys.argv[3], sys.argv[4],
             [int(c) for c in sys.argv[5].split(',')], t)
    else:
        cmp(sys.argv[2], sys.argv[3], sys.argv[4])
