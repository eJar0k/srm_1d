"""IG3: run the canonical-motor validation rows (ig3_valid.py) for local/head/aft in parallel.
Usage: ig3_run_valid.py <out.jsonl> [workers=6] [demar=step]. Run from the worktree root with PYTHONPATH=$PWD."""
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
out = sys.argv[1]
nw = int(sys.argv[2]) if len(sys.argv) > 2 else 6
demar = sys.argv[3] if len(sys.argv) > 3 else 'step'
specs = [dict(motor=m, ap_ref=ap, demar=demar) for m in ('haseg', 'chunc', 'balls') for ap in ('local', 'head', 'aft')]


def go(spec):
    r = subprocess.run([sys.executable, str(HERE / 'ig3_valid.py'), json.dumps(spec), out],
                       capture_output=True, text=True)
    print(spec, 'rc', r.returncode, flush=True)
    if r.returncode:
        print(r.stderr[-2000:], flush=True)


with ThreadPoolExecutor(nw) as ex:
    list(ex.map(go, specs))
