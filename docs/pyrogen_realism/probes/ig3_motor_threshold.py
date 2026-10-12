"""IG3: minimum forward-plenum BPNV charge that lights each canonical motor (>=99 % of grain lit), per combo.
Usage: ig3_motor_threshold.py <out.jsonl>. Run from the worktree root with PYTHONPATH=$PWD."""
import sys, json, time
from pathlib import Path
from multiprocessing import get_context
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).parent))

def _job(a):
    import ig3_case
    ric, motor, h, d, sutton_g, ap = a
    rows = []
    def run(g):
        m = ig3_case.run_case(dict(tag=f'mt_{motor}_{h}_{ap}_{d}_{g:.5g}', ric=ric, topo='forward_plenum', g=g, ign_h=h, ap_ref=ap,
                                   demar=d, t_max=0.25, plateau_win=[0.1, 0.25], P_cutoff=0.05e6)); m['motor'] = motor
        rows.append(m); return m['frac_ignited'] >= 0.99
    hi = sutton_g * 16.0
    while not run(hi): hi *= 2.0
    g = hi
    lo = None
    while g > 1e-3:
        g /= 2.0
        if run(g): hi = g
        else: lo = g; break
    if lo is None: return rows
    while hi / lo > 1.03:
        mid = (lo * hi) ** 0.5
        if run(mid): hi = mid
        else: lo = mid
    return rows

if __name__ == '__main__':
    out = Path(sys.argv[1])
    motors = [('machbusterNew.ric', 'chunc', 0.904), ('BALLSstick.ric', 'balls', 1.386), ('hasegawa_a.ric', 'haseg', 4.184)]
    jobs = [('motors/' + r, m, h, d, g, ap) for r, m, g in motors for h, d, ap in (('peretz', 'step', 'local'), ('peretz', 'step', 'head'), ('peretz', 'step', 'aft'), ('peretz', 'off', 'local'), ('peretz', 'off', 'head'), ('peretz', 'off', 'aft'))]
    t0 = time.time()
    with get_context('spawn').Pool(6) as pool, open(out, 'a') as f:
        for rows in pool.imap_unordered(_job, jobs):
            for m in rows: f.write(json.dumps(m) + chr(10))
            f.flush(); print(rows[0]['motor'], rows[0]['ap_ref'], rows[0]['demar'], len(rows), flush=True)
    print('wall', time.time() - t0)
