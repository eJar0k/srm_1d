"""GS1 driver: runs gs1_row.py jobs, each as its own guarded job, with a small worker pool.
Usage: python docs/pyrogen_realism/probes/gs1_run.py <phase A|B|C> <out.jsonl> <workers> <mem_gb_per_job>
Phases (one cold-compile warm-up job should have been run first so the numba cache is warm):
  A  N in {100,200,400}: clean rows (t_max 0.15) + threshold rows (t_max 0.3)
  B  N = 800: clean rows (t_max 0.15) + threshold rows (t_max 0.2)
  C  full burns at N in {100,200} for the threshold charges (plateau / late peak): chunc 1 s, balls 2 s, haseg 10 s
Job output goes to <out.jsonl>; the guard line (wall_s, peak_job_gb) is appended to <out>.guard.jsonl."""
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
phase, out, nw, mem = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3]), sys.argv[4]

CLEAN = [('chunc_fp10', 'chunc', 10.0, 'forward_plenum'), ('balls_fp10', 'balls', 10.0, 'forward_plenum'),
         ('haseg_fp30', 'haseg', 30.0, 'forward_plenum'), ('chunc_hb4p5', 'chunc', 4.5, 'head_basket')]
THRESH = [('chunc_th2p62', 'chunc', 2.62, 'forward_plenum'), ('balls_th4p83', 'balls', 4.83, 'forward_plenum'),
          ('haseg_th14p98', 'haseg', 14.98, 'forward_plenum')]
FULL_T = dict(chunc=1.0, balls=2.0, haseg=10.0)


def specs():
    out_ = []
    if phase == 'A':
        for N in (100, 200, 400):
            for nm, m, g, tp in CLEAN:
                out_.append(dict(tag=f'{nm}_N{N}', motor=m, g=g, topo=tp, N=N, t_max=0.15, plateau_win=[0.1, 0.15]))
            for nm, m, g, tp in THRESH:
                out_.append(dict(tag=f'{nm}_N{N}', motor=m, g=g, topo=tp, N=N, t_max=0.3, plateau_win=[0.2, 0.3]))
    elif phase == 'B':
        for nm, m, g, tp in CLEAN:
            out_.append(dict(tag=f'{nm}_N800', motor=m, g=g, topo=tp, N=800, t_max=0.15, plateau_win=[0.1, 0.15]))
        for nm, m, g, tp in THRESH:
            out_.append(dict(tag=f'{nm}_N800', motor=m, g=g, topo=tp, N=800, t_max=0.2, plateau_win=[0.15, 0.2]))
    elif phase == 'C':
        for N in (100, 200):
            for nm, m, g, tp in THRESH:
                t = FULL_T[m]
                out_.append(dict(tag=f'{nm}_N{N}_full', motor=m, g=g, topo=tp, N=N, t_max=t,
                                 plateau_win=[0.4 * t, 0.9 * t] if m == 'chunc' else [0.5 * t, 0.9 * t]))
    return out_


def job(spec):
    done = out.exists() and any(json.loads(l)['spec']['tag'] == spec['tag'] for l in out.read_text().splitlines() if l.strip())
    if done:
        return
    cmd = [sys.executable, '-m', 'srm_1d.tools.guarded_run', '--mem-gb', mem, '--timeout-min', '60', '--',
           sys.executable, str(ROOT / 'docs/pyrogen_realism/probes/gs1_row.py'), json.dumps(spec), str(out)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT))
    g = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith('[guarded_run]')]
    with open(str(out) + '.guard.jsonl', 'a') as f:
        f.write(json.dumps(dict(tag=spec['tag'], guard=g[-1] if g else None, rc=r.returncode, wall=round(time.time() - t0, 1),
                                err=(r.stderr[-300:] if r.returncode else None))) + '\n')
    print(spec['tag'], g[-1] if g else r.returncode, flush=True)


with ThreadPoolExecutor(nw) as ex:
    list(ex.map(job, specs()))
