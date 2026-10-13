"""SM2 sweeps. Usage: sm2_sweep.py <group> <out.jsonl> [workers<=8]
groups: A (base, 3 models x 4 cases, forcing recorded), H (gnielinski ign h), D (DeMar on), C (auto charge),
        K (calibration k_solid x T_ignition, short windows), G (200 cells).
Run from the repo root with PYTHONPATH=$PWD. Resumable (tag-keyed)."""
import sys, json, time
from pathlib import Path
from multiprocessing import get_context
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).parent))

MOTORS = ('chunc', 'chunc_hb', 'balls', 'haseg')
MODELS3 = ('goodman', 'peretz', 'conduction')
MODELS2 = ('peretz', 'conduction')
KS = (0.26, 0.271, 0.32)
TIGN = (756.0, 850.0, 950.0)          # v0.7.5 LHS band T_ignition in [750, 950] K (docs/v0_7_5/RESULT.md); 756 = optimum
SHORT = {'chunc': 0.3, 'chunc_hb': 0.2, 'balls': 0.3, 'haseg': 0.3}


def specs(group):
    out = []
    if group == 'A':
        for mo in MOTORS:
            for m in MODELS3:
                out.append(dict(tag=f'A_{mo}_{m}', motor=mo, solid=m, rec=120000))
    elif group == 'H':
        for mo in MOTORS:
            for m in MODELS2:
                out.append(dict(tag=f'H_{mo}_{m}', motor=mo, solid=m, ign_h='gnielinski', rec=60000))
    elif group == 'D':
        for mo in MOTORS:
            for m in MODELS2:
                out.append(dict(tag=f'D_{mo}_{m}', motor=mo, solid=m, demar='on', rec=60000))
    elif group == 'C':
        for mo in ('chunc', 'balls', 'haseg'):
            for m in MODELS2:
                out.append(dict(tag=f'C_{mo}_{m}', motor=mo, solid=m, charge='auto', rec=60000))
    elif group == 'K':
        for mo in MOTORS:
            for m in MODELS2:
                for k in KS:
                    for T in TIGN:
                        out.append(dict(tag=f'K_{mo}_{m}_k{k}_T{int(T)}', motor=mo, solid=m, k_solid=k, T_ign=T,
                                        t_max=SHORT[mo]))
    elif group == 'G':
        for mo in MOTORS:
            for m in MODELS3:
                out.append(dict(tag=f'G_{mo}_{m}', motor=mo, solid=m, cells=200))
    else:
        raise SystemExit('unknown group')
    return out


def _work(spec):
    import sm2_case
    try:
        return sm2_case.run_case(spec)
    except Exception as e:
        import traceback
        return dict(tag=spec['tag'], error=repr(e), tb=traceback.format_exc()[-600:], **{k: v for k, v in spec.items() if k != 'tag'})


if __name__ == '__main__':
    group, out = sys.argv[1], Path(sys.argv[2])
    workers = min(int(sys.argv[3]) if len(sys.argv) > 3 else 8, 8)
    done = set()
    if out.exists():
        done = {json.loads(l).get('tag') for l in out.read_text().splitlines() if l.strip()}
    todo = [s for s in specs(group) if s['tag'] not in done]
    # longest first (haseg, then balls) for load balance
    order = {'haseg': 0, 'balls': 1, 'chunc': 2, 'chunc_hb': 3}
    todo.sort(key=lambda s: (order[s['motor']], s['tag']))
    t0 = time.time()
    with get_context('spawn').Pool(workers) as pool, open(out, 'a') as f:
        for m in pool.imap_unordered(_work, todo):
            f.write(json.dumps(m) + '\n'); f.flush()
            print(m['tag'], 'ERR ' + m['error'] if 'error' in m else 'wall %.0f' % m['wall_s'], flush=True)
    print('group', group, 'n', len(todo), 'wall', time.time() - t0)
