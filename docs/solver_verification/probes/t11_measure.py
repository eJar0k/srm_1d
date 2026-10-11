"""T11 measurement driver: old (K=0,S=False) vs new (K=1,S=True) PISO defaults.

Job names (comma separated on the command line):
  motor:<chunc|balls|haseg>:<erosive_model>:<cap>:<old|new>
  grid:<target_cells>:<old|new>          (Chunc, CFL 0.3, t_max 1.0, lib BPNV)
Rows are appended to a JSONL accumulator (resumable: a job already present is
skipped). The piso flags are passed EXPLICITLY for both configs, so the result
does not depend on the library defaults.

Run from the worktree root, pyenv 3.10.5, SRM1D_OPENMOTOR_PATH exported, e.g.
  python docs/solver_verification/probes/t11_measure.py --out C:/tmp/sv_t11/results.jsonl \
      motor:chunc:ma2020:0.0:old grid:106:new
"""
import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np  # noqa: E402


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CFG = {'old': dict(piso_energy_coupling_iters=0, piso_corr2_storage=False),
       'new': dict(piso_energy_coupling_iters=1, piso_corr2_storage=True)}


def run_job(job):
    f = job.split(':')
    if f[0] == 'motor':
        _, motor, model, cap, cfg = f
        probe = _load("mp_cap_crossmotor", "docs/v0_7_4/probes/mp_cap_crossmotor.py")
        probe.CANON.update(CFG[cfg])
        fn = {'chunc': probe.chunc_row, 'balls': probe.balls_row,
              'haseg': probe.haseg_row}[motor]
        t0 = time.time()
        row = fn(model, float(cap))
        row['wall_s'] = time.time() - t0
        if row.get('n_steps'):
            row['steps_per_s'] = row['n_steps'] / row['wall_s']
    elif f[0] == 'grid':
        _, cells, cfg = f
        gc = _load("grid_convergence", "docs/solver_verification/probes/grid_convergence.py")
        gc.CANON.update(CFG[cfg])
        row = gc.run_row('A', 'lib_bpnv', int(cells), 0.3, 1.0)
    else:
        raise ValueError(job)
    row['job'] = job
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='C:/tmp/sv_t11/results.jsonl')
    ap.add_argument('jobs', nargs='+')
    a = ap.parse_args()
    out = Path(a.out)
    for job in [j for arg in a.jobs for j in arg.split(',')]:
        done = set()
        if out.exists():
            done = {json.loads(l)['job'] for l in out.read_text().splitlines() if l.strip()}
        if job in done:
            print('skip', job, flush=True)
            continue
        import srm_1d
        print('[run]', job, srm_1d.__file__, flush=True)
        row = run_job(job)
        with open(out, 'a') as fh:
            fh.write(json.dumps(row) + '\n')
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
