"""MEM1: run mem1_probe.py once per committed motor, SEQUENTIALLY, each as its own guarded job
(4 GB / 30 min), appending a JSON line (probe output + guard summary) to the out file.
Usage: mem1_sweep.py <out.jsonl> [motor_substring ...]   (run from the worktree root, PYTHONPATH=$PWD)"""
import json, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
out = Path(sys.argv[1]); only = sys.argv[2:]
# (ric, extra probe args): B0 rows at their gate settings, pathfin54 on the GUI path, the rest at run_from_ric defaults
JOBS = [('motors/machbusterNew.ric', ['1.0', '10000']), ('motors/BALLSstick.ric', ['2.0', '10000']),
        ('motors/hasegawa_a.ric', ['10.0', '50000']),
        ('motors/2025.12.25 54-2800 ST2.0.ric', ['-', '-', 'gui'])]
done = {p.name for p in (ROOT / 'motors').glob('*.ric')}
for r in sorted((ROOT / 'motors').glob('*.ric')):
    rel = 'motors/' + r.name
    if rel not in [j[0] for j in JOBS]:
        JOBS.append((rel, ['-', '-']))
for ric, extra in JOBS:
    if only and not any(o in ric for o in only):
        continue
    cmd = [sys.executable, '-m', 'srm_1d.tools.guarded_run', '--mem-gb', '4', '--timeout-min', '30', '--',
           sys.executable, str(ROOT / 'docs/pyrogen_realism/probes/mem1_probe.py'), ric] + extra
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT))
    line = [l for l in r.stdout.splitlines() if l.startswith('{')]
    guard = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith('[guarded_run]')]
    row = dict(ric=ric, extra=extra, probe=json.loads(line[-1]) if line else None, guard=guard[-1] if guard else None,
               err=None if line else (r.stderr[-300:] or r.stdout[-300:]))
    with open(out, 'a') as f:
        f.write(json.dumps(row) + '\n')
    print(ric, row['guard'], flush=True)
