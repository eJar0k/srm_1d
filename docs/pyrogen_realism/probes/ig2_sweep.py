"""IG2 resumable sweeps. Usage: ig2_sweep.py <threshold|bisect|realistic|chunc_valid> <out.jsonl> [workers]
Appends one JSON line per finished case (tag-keyed, skipped if present). Run from the worktree root with
PYTHONPATH=$PWD (see docs/solver_verification/TASKS.md)."""
import sys, json, time
from pathlib import Path
from multiprocessing import get_context
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).parent))

HS = ['gnielinski', 'peretz', 'keller']
DM = ['step', 'off']
GRID_G = [0.05, 0.1, 0.2, 0.3, 0.36, 0.45, 0.6, 0.9, 1.36, 1.8, 2.27, 2.93, 3.6, 4.54, 6.3, 9.0, 13.6]


def _work(spec):
    import ig2_case
    try:
        return ig2_case.run_case(spec)
    except Exception as e:           # keep the sweep going; record the failure
        return dict(tag=spec['tag'], error=repr(e), **{k: v for k, v in spec.items() if k != 'tag'})


def specs_threshold():
    out = []
    for topo in ('head_basket', 'forward_plenum'):
        for h in HS:
            for d in DM:
                for g in GRID_G:
                    out.append(dict(tag=f'th_{topo}_{h}_{d}_{g}', topo=topo, g=g, ign_h=h, demar=d, t_max=0.2))
    return out


def specs_realistic():
    out = []
    for pyro in ('bpnv', 'mtv'):
        for g in (4.5, 6.0):
            for h in HS:
                for d in DM:
                    out.append(dict(tag=f'rc_{pyro}_{g}_{h}_{d}', topo='head_basket', g=g, pyrogen=pyro, cubes=True,
                                    ign_h=h, demar=d, t_max=0.2))
    return out


def _bisect(args):
    """Find the minimum charge (g) with frac_ignited >= 0.99 within ~2 % for one combo.
    args = (topo, h, d, lo_fail, hi_ok); lo_fail may be None (then extend downward)."""
    import ig2_case
    topo, h, d, lo, hi = args
    rows = []
    def run(g):
        sp = dict(tag=f'bi_{topo}_{h}_{d}_{g:.5g}', topo=topo, g=g, ign_h=h, demar=d, t_max=0.2)
        m = ig2_case.run_case(sp); rows.append(m); return m['frac_ignited'] >= 0.99
    if lo is None:                      # all grid points ignite: walk down by halving
        g = hi
        while g > 2e-4:
            g = g / 2.0
            if not run(g):
                lo = g; break
        else:
            return rows
    while hi / lo > 1.02:
        mid = (lo * hi) ** 0.5
        if run(mid): hi = mid
        else: lo = mid
    return rows


def load_done(path):
    done = set()
    if path.exists():
        for ln in path.read_text().splitlines():
            if ln.strip():
                done.add(json.loads(ln)['tag'])
    return done


def main():
    kind, out = sys.argv[1], Path(sys.argv[2]); nw = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    if kind == 'threshold':
        specs = specs_threshold()
    elif kind == 'realistic':
        specs = specs_realistic()
    elif kind == 'bisect':
        base = [json.loads(l) for l in open(sys.argv[4])]
        combos = {}
        for r in base:
            combos.setdefault((r['topo'], r['ign_h'], r['demar']), []).append(r)
        jobs = []
        for (topo, h, d), rs in combos.items():
            rs = sorted(rs, key=lambda r: r['g'])
            ok = [r['g'] for r in rs if r['frac_ignited'] >= 0.99]
            bad = [r['g'] for r in rs if r['frac_ignited'] < 0.99]
            hi = min(ok); lo = max([b for b in bad if b < hi], default=None)
            jobs.append((topo, h, d, lo, hi))
        t0 = time.time()
        with get_context('spawn').Pool(nw) as pool, open(out, 'a') as f:
            for rows in pool.imap_unordered(_bisect, jobs):
                for m in rows:
                    f.write(json.dumps(m) + chr(10))
                f.flush()
                print(rows[0]['topo'], rows[0]['ign_h'], rows[0]['demar'], 'runs', len(rows), flush=True)
        print(f'total wall {time.time()-t0:.0f} s'); return
    elif kind == 'file':          # extra specs from a json list file given as argv[4]
        specs = json.load(open(sys.argv[4]))
    else:
        raise SystemExit('unknown sweep kind')
    done = load_done(out)
    todo = [s for s in specs if s['tag'] not in done]
    print(f'{len(specs)} specs, {len(todo)} to run, {nw} workers', flush=True)
    t0 = time.time()
    ctx = get_context('spawn')
    with ctx.Pool(nw) as pool, open(out, 'a') as f:
        for i, m in enumerate(pool.imap_unordered(_work, todo), 1):
            f.write(json.dumps(m) + '\n'); f.flush()
            print(f"[{i}/{len(todo)}] {m['tag']} peak={m.get('peak_MPa')} frac={m.get('frac_ignited')} "
                  f"wall={m.get('wall_s')} err={m.get('error')}", flush=True)
    print(f'total wall {time.time()-t0:.0f} s')


if __name__ == '__main__':
    main()
