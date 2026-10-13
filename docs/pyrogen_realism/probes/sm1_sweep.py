"""SM1 Part 3 sweeps. Usage: sm1_sweep.py <realistic|threshold> <out.jsonl> [workers<=8]
Run from the worktree root with PYTHONPATH=$PWD. Resumable (tag-keyed). Chunc (machbusterNew.ric), OB1 defaults
(Peretz ignition h, I-A pellets, exact-cylinder basket law)."""
import sys, json, time
from pathlib import Path
from multiprocessing import get_context
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).parent))

MODELS = ('goodman', 'peretz', 'conduction')


def _work(spec):
    import sm1_case
    try:
        return sm1_case.run_case(spec)
    except Exception as e:
        return dict(tag=spec['tag'], error=repr(e), **{k: v for k, v in spec.items() if k != 'tag'})


def specs_realistic():
    out = []
    for solid in MODELS:
        for pyro in ('bpnv', 'mtv'):
            for g in (4.5, 6.0):
                for d in ('on', 'off'):
                    out.append(dict(tag=f'rc_{solid}_{pyro}_{g}_{d}', topo='head_basket', g=g, pyrogen=pyro,
                                    demar=d, solid=solid, t_max=0.2, plateau_win=[0.1, 0.2]))
    return out


def _bisect(args):
    """Minimum charge (g) with >= crit of the grain lit within t_max, to 3 %. args = (topo, solid, start_g, crit)."""
    import sm1_case
    topo, solid, g0, crit = args
    rows = []

    def run(g):
        m = sm1_case.run_case(dict(tag=f'bi{crit}_{topo}_{solid}_{g:.5g}', topo=topo, g=g, demar='off', solid=solid,
                                   t_max=0.2, plateau_win=[0.1, 0.2], P_cutoff=0.05e6))
        rows.append(m)
        return m['frac_ignited'] >= crit - 1e-9
    hi = g0 * 4.0
    while not run(hi):
        hi *= 2.0
        if hi > 200:
            return rows
    lo = None
    g = hi
    while g > 2e-3:
        g /= 2.0
        if run(g):
            hi = g
        else:
            lo = g
            break
    if lo is None:
        return rows
    while hi / lo > 1.03:
        mid = (lo * hi) ** 0.5
        if run(mid):
            hi = mid
        else:
            lo = mid
    return rows


if __name__ == '__main__':
    mode, out = sys.argv[1], Path(sys.argv[2])
    workers = min(int(sys.argv[3]) if len(sys.argv) > 3 else 8, 8)
    t0 = time.time()
    done = set()
    if out.exists():
        done = {json.loads(l).get('tag') for l in out.read_text().splitlines() if l.strip()}
    with get_context('spawn').Pool(workers) as pool, open(out, 'a') as f:
        if mode == 'realistic':
            specs = [s for s in specs_realistic() if s['tag'] not in done]
            for m in pool.imap_unordered(_work, specs):
                f.write(json.dumps(m) + '\n'); f.flush(); print(m['tag'], m.get('wall_s'), flush=True)
        elif mode == 'threshold':
            jobs = [(t, s, g, crit) for crit in (0.99, 1.0) for t, g in (('head_basket', 3.0), ('forward_plenum', 2.27)) for s in MODELS]
            for rows in pool.imap_unordered(_bisect, jobs):
                for m in rows:
                    f.write(json.dumps(m) + '\n')
                f.flush(); print(rows[0]['topo'], rows[0]['solid'], len(rows), flush=True)
    print('wall', time.time() - t0)
