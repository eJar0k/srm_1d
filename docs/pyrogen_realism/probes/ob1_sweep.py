"""OB1 resumable sweeps (post-rebase). Usage: ob1_sweep.py <name> <out.jsonl> [workers]
names: realistic (5.1, t_max 0.2 s) | long (5.1, t_max 1.0 s) | thresh (5.2 bisection) | mass (peak vs charge) |
extra (attribution / ablation probes). Appends one JSON line per finished case (tag-keyed, skipped if present).
Run from the worktree root with PYTHONPATH=$PWD. "pre_ob1" = the library pellet of the previous library (3.2 mm,
L/D 1) with the pooled constant-area basket law = the state before OB1 (the forward-plenum law is identical there)."""
import sys, json, time
from pathlib import Path
from multiprocessing import get_context
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))

DEMAR = ['on', 'off']
# model name -> (geom, area, law)
MODELS = {
    'pre_ob1':      ('old32', 'constant', 'lib'),
    'lib_constant': ('lib', 'constant', 'lib'),
    'lib':          ('lib', 'shrinking', 'lib'),
    'lib_capZ':     ('lib', 'shrinking', 'capZ_lib'),
    'cube3_constant': ('cube3', 'constant', 'lib'),
    'cube3':        ('cube3', 'shrinking', 'lib'),
    'cube3_capZ':   ('cube3', 'shrinking', 'capZ_lib'),
}


def _work(spec):
    import ob1_case
    try:
        return ob1_case.run_case(spec)
    except Exception as e:
        return dict(tag=spec['tag'], error=repr(e), **{k: v for k, v in spec.items() if k != 'tag'})


def S(tag, model, **kw):
    geom, area, law = MODELS[model]
    return dict(tag=tag, model=model, geom=geom, area=area, law=law, **kw)


def specs_realistic():
    out = []
    for pyro in ('bpnv', 'mtv'):
        for g in (4.5, 6.0):
            for d in DEMAR:
                for m in MODELS:
                    if pyro == 'mtv' and m.endswith('capZ'):
                        continue
                    out.append(S(f're_{pyro}_{g}_{d}_{m}', m, pyrogen=pyro, g=g, demar=d, t_max=0.2))
    return out


def specs_long():
    out = []
    for pyro in ('bpnv', 'mtv'):
        for g in (4.5, 6.0):
            for d in DEMAR:
                for m in ('pre_ob1', 'lib', 'cube3'):
                    out.append(S(f'lg_{pyro}_{g}_{d}_{m}', m, pyrogen=pyro, g=g, demar=d, t_max=1.0,
                                 plateau_win=(0.4, 0.9), snap=5e-3))
    return out


def specs_mass():
    out = []
    for d in DEMAR:
        for m in ('pre_ob1', 'lib', 'cube3'):
            for g in (0.6, 0.9, 1.36, 1.8, 2.27, 2.93, 3.6, 4.54, 6.3, 9.0):
                out.append(S(f'ms_{d}_{m}_{g}', m, g=g, demar=d, t_max=0.2))
    return out


def specs_extra():
    out = []
    for m in ('pre_ob1', 'lib', 'lib_capZ'):
        for em in ('saint_robert', 'mukunda_paul'):
            out.append(S(f'ex_{m}_{em}', m, g=4.5, erosive_model=em, t_max=0.2))
        out.append(S(f'ex_{m}_pyro_only', m, g=4.5, T_ign=1.0e5, t_max=0.1))
    # ablations on the new default (lib, shrinking), BPNV 4.5 g
    for law in ('zhu_lo', 'zhu_hi', 'zhu_mean', 'cap7_lib'):
        out.append(dict(tag=f'ex_lib_{law}', model='lib', geom='lib', area='shrinking', law=law, g=4.5, t_max=0.2))
    for v in (3.0, 10.0, 30.0):
        out.append(dict(tag=f'ex_lib_bed{v:g}', model='lib', geom='lib', area='shrinking', law='lib', bed_v=v, g=4.5,
                        t_max=0.2))
    return out


GRID_G = [0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.9, 1.36, 1.8, 2.27, 2.93, 3.6, 4.54, 6.3, 9.0, 13.6, 20.0, 30.0]
# (topo, demar, model) series for the threshold search
THRESH_SERIES = (
    [('head_basket', d, m) for d in DEMAR for m in ('pre_ob1', 'lib_constant', 'lib', 'cube3')] +
    [('forward_plenum', d, m) for d in DEMAR for m in ('pre_ob1', 'lib', 'cube3')]
)


def _thresh_series(args):
    """Lowest charge on a grid scan with >= 99 % of the 100 grain cells lit within t_max, bisected to 2 %
    between the last failing and first lit grid point. The scan continues past the first lit point to flag
    non-monotonic (ragged) boundaries: 'relit' lists lit charges above any later failure."""
    import ob1_case
    topo, d, m = args
    rows = []
    grid = ([0.005, 0.01, 0.02, 0.03, 0.04] + GRID_G) if topo == 'forward_plenum' else GRID_G

    def run(g):
        sp = S(f'thg_{topo}_{d}_{m}_{g:.5g}', m, topo=topo, g=g, demar=d, t_max=0.2)
        r = ob1_case.run_case(sp)
        rows.append(r)
        return r['frac_ignited'] >= 0.99
    lo, hi, flags = None, None, []
    for g in grid:
        ok = run(g)
        flags.append((g, ok))
        if ok and hi is None:
            hi = g
        elif not ok and hi is None:
            lo = g
        if ok and g >= 1.0 and sum(1 for _, o in flags if o) >= 3:
            break
    thr = None
    if hi is not None:
        if lo is None:
            thr = ('<=', hi)
        else:
            while hi / lo > 1.02:
                mid = (hi * lo) ** 0.5
                if run(mid):
                    hi = mid
                else:
                    lo = mid
            thr = ('', hi)
    lit = [g for g, o in flags if o]
    notlit_after = [g for g, o in flags if not o and hi is not None and g > hi]
    summary = dict(tag=f'thr_{topo}_{d}_{m}', topo=topo, demar=d, model=m,
                   threshold_g=(None if thr is None else thr[1]), threshold_note=(thr[0] if thr else 'not lit up to grid top'),
                   grid_lit=lit, grid_fail_above_threshold=notlit_after, n_runs=len(rows))
    return rows + [summary]


def run_thresh(out, nw, force_plenum_on=False):
    done = set()
    if out.exists() and not force_plenum_on:
        done = {json.loads(l)['tag'] for l in out.read_text().splitlines() if l.strip()}
    jobs = [j for j in THRESH_SERIES if f'thr_{j[0]}_{j[1]}_{j[2]}' not in done]
    if force_plenum_on:
        jobs = [j for j in THRESH_SERIES if j[0] == 'forward_plenum' and j[1] == 'on']
    with get_context('spawn').Pool(nw) as pool:
        for rows in pool.imap_unordered(_thresh_series, jobs):
            with open(out, 'a') as f:
                for r in rows:
                    f.write(json.dumps(r) + chr(10))
            print(rows[-1], flush=True)


SPECS = {'realistic': specs_realistic, 'long': specs_long, 'mass': specs_mass, 'extra': specs_extra}

if __name__ == '__main__':
    name, out = sys.argv[1], Path(sys.argv[2])
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    if name == 'thresh':
        run_thresh(out, nw)
        sys.exit(0)
    if name == 'thresh_plenum_low':      # re-run the DeMar-on plenum series with a finer low grid (last summary row wins)
        run_thresh(out, nw, force_plenum_on=True)
        sys.exit(0)
    done = set()
    if out.exists():
        done = {json.loads(l)['tag'] for l in out.read_text().splitlines() if l.strip()}
    specs = [s for s in SPECS[name]() if s['tag'] not in done]
    print(f'{len(specs)} to run ({len(done)} done)', flush=True)
    t0 = time.time()
    with get_context('spawn').Pool(nw) as pool:
        for m in pool.imap_unordered(_work, specs):
            with open(out, 'a') as f:
                f.write(json.dumps(m) + chr(10))
            print(m['tag'], m.get('peak_MPa'), m.get('error', ''), f'{time.time()-t0:.0f}s', flush=True)
