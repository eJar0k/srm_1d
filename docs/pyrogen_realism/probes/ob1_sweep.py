"""OB1 resumable sweeps. Usage: ob1_sweep.py <name> <out.jsonl> [workers]
names: rank (Part-3 ranking probes, BPNV 4.5 g, both combos) | realistic (Part 5.1 matrix, t_max 0.2 s)
| long (Part 5.1 t_max 1.0 s) | extra (supplementary probes, see specs_extra).
Appends one JSON line per finished case (tag-keyed, skipped if present). Run from the worktree root with
PYTHONPATH=$PWD."""
import sys, json, time
from pathlib import Path
from multiprocessing import get_context
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))

COMBOS = ['current', 'forward']
# name -> (law, area, bed_v)
OPTS = {
    'base':        ('lib', 'constant', 0.0),
    'area23':      ('lib', 'shrinking', 0.0),
    'zhu':         ('zhu_mean', 'constant', 0.0),
    'zhu+area23':  ('zhu_mean', 'shrinking', 0.0),
    'cap7':        ('cap7_lib', 'constant', 0.0),
    'capZ':        ('capZ_lib', 'constant', 0.0),
    'capZ+area23': ('capZ_lib', 'shrinking', 0.0),
}


def _work(spec):
    import ob1_case
    try:
        return ob1_case.run_case(spec)
    except Exception as e:
        return dict(tag=spec['tag'], error=repr(e), **{k: v for k, v in spec.items() if k != 'tag'})


def specs_rank():
    out = []
    for combo in COMBOS:
        for name, (law, area, bv) in OPTS.items():
            out.append(dict(tag=f'rk_{combo}_{name}', combo=combo, law=law, area=area, bed_v=bv, g=4.5, t_max=0.2))
        # diagnostic: bed flame speed (no sourced value) on top of the area law
        for v in (3.0, 10.0, 30.0):
            out.append(dict(tag=f'rk_{combo}_area23+bed{v:g}', combo=combo, law='lib', area='shrinking', bed_v=v,
                            g=4.5, t_max=0.2))
    return out


def specs_realistic():
    out = []
    for pyro in ('bpnv', 'mtv'):
        for g in (4.5, 6.0):
            for combo in COMBOS:
                for name, (law, area, bv) in OPTS.items():
                    if pyro == 'mtv' and law != 'lib':
                        continue                  # no sourced plateau law for MTV
                    out.append(dict(tag=f're_{pyro}_{g}_{combo}_{name}', pyrogen=pyro, g=g, combo=combo,
                                    law=law, area=area, bed_v=bv, t_max=0.2))
    return out


def specs_long():
    out = []
    for pyro in ('bpnv', 'mtv'):
        for g in (4.5, 6.0):
            for combo in COMBOS:
                for name in ('base', 'area23', 'zhu+area23', 'capZ+area23'):
                    law, area, bv = OPTS[name]
                    if pyro == 'mtv' and law != 'lib':
                        continue
                    out.append(dict(tag=f'lg_{pyro}_{g}_{combo}_{name}', pyrogen=pyro, g=g, combo=combo,
                                    law=law, area=area, bed_v=bv, t_max=1.0, plateau_win=(0.4, 0.9), snap=5e-3))
    return out


def specs_extra():
    """Attribution probes on the best grounded model (capZ+area23) and the base model, BPNV 4.5 g, current combo."""
    out = []
    for name in ('base', 'capZ+area23'):
        law, area, bv = OPTS[name]
        for em in ('saint_robert', 'mukunda_paul'):
            out.append(dict(tag=f'ex_{name}_{em}', law=law, area=area, g=4.5, erosive_model=em, t_max=0.2))
        out.append(dict(tag=f'ex_{name}_pyro_only', law=law, area=area, g=4.5, T_ign=1.0e5, t_max=0.1))
    # sensitivity to the measured spread of Zhu's r0 (Table 2: 5.27e-3 .. 6.85e-3 m/s/MPa)
    for law in ('zhu_lo', 'zhu_hi'):
        out.append(dict(tag=f'ex_{law}+area23', law=law, area='shrinking', g=4.5, t_max=0.2))
    return out


GRID_G = [0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.9, 1.36, 1.8, 2.27, 2.93, 3.6, 4.54, 6.3, 9.0, 13.6]
THRESH_OPTS = ('base', 'area23', 'capZ+area23', 'zhu+area23')


def _thresh_series(args):
    """Minimum head-basket charge (g) with >= 99 % of the grain cells lit within t_max, bisected to 2 %.
    Returns the rows run plus a summary row (tag 'thr_<combo>_<opt>')."""
    import ob1_case
    combo, name = args
    law, area, bv = OPTS[name]
    rows = []

    def run(g):
        sp = dict(tag=f'thg_{combo}_{name}_{g:.5g}', g=g, combo=combo, law=law, area=area, bed_v=bv, t_max=0.2)
        m = ob1_case.run_case(sp)
        rows.append(m)
        return m['frac_ignited'] >= 0.99
    lo, hi = None, None
    for g in GRID_G:
        if run(g):
            hi = g
            break
        lo = g
    if hi is None:
        thr = None
    elif lo is None:
        thr = ('<=', hi)
    else:
        while hi / lo > 1.02:
            mid = (hi * lo) ** 0.5
            if run(mid):
                hi = mid
            else:
                lo = mid
        thr = ('', hi)
    summary = dict(tag=f'thr_{combo}_{name}', combo=combo, opt=name, threshold_g=(None if thr is None else thr[1]),
                   threshold_note=('not lit up to 13.6 g' if thr is None else thr[0]), n_runs=len(rows))
    return rows + [summary]


def run_thresh(out, nw):
    jobs = [(c, n) for c in COMBOS for n in THRESH_OPTS]
    with get_context('spawn').Pool(nw) as pool:
        for rows in pool.imap_unordered(_thresh_series, jobs):
            with open(out, 'a') as f:
                for r in rows:
                    f.write(json.dumps(r) + chr(10))
            print(rows[-1], flush=True)


def specs_mass():
    """Peak vs charge mass above threshold (BPNV head basket, 3 mm cubes)."""
    out = []
    for combo in COMBOS:
        for name in ('base', 'area23', 'capZ+area23'):
            law, area, bv = OPTS[name]
            for g in (0.6, 0.9, 1.36, 1.8, 2.27, 2.93, 3.6, 4.54, 6.3, 9.0):
                out.append(dict(tag=f'ms_{combo}_{name}_{g}', g=g, combo=combo, law=law, area=area, bed_v=bv, t_max=0.2))
    return out


SPECS = {'rank': specs_rank, 'realistic': specs_realistic, 'long': specs_long, 'extra': specs_extra, 'mass': specs_mass}

if __name__ == '__main__':
    name, out = sys.argv[1], Path(sys.argv[2])
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    if name == 'thresh':
        run_thresh(out, nw)
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
