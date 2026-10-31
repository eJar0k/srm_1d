"""HE1: minimum forward-plenum BPNV charge that lights >= 99 % of the grain within 0.25 s
(IG2 section 4.1 criterion), conduction solid, sensible/legacy/static, head_end_jet_h off vs max.
    python he1_threshold.py <motor> <out.json>      motor = chunc | balls | haseg
Doubling bracket from 1 g up to 64 g, then bisection to 3 %. Bounded: <= 20 runs per jet mode."""
import json, sys, time, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
warnings.simplefilter('ignore')
from srm_1d.openmotor_adapter import run_from_ric

RIC = {'chunc': 'machbusterNew.ric', 'balls': 'BALLSstick.ric', 'haseg': 'hasegawa_a.ric'}
motor, out = sys.argv[1], sys.argv[2]
rows = []


def lit(g, jet):
    r = run_from_ric(str(ROOT / 'motors' / RIC[motor]), t_max=0.25, P_cutoff=0.05e6, pyrogen='bpnv',
                     pyrogen_mass=g * 1e-3, roughness=32e-6, kappa=0.44, T_ignition=756.0,
                     print_interval=0.0, verbose=False, solid_heating_model='conduction',
                     head_end_jet_h=jet, snapshot_interval=0.05)[0]
    f = float(r['summary']['grain_fraction_lit'])
    rows.append(dict(motor=motor, jet=jet, g=g, frac=f, never=int(r['summary']['grain_cells_never_lit'])))
    print(motor, jet, round(g, 4), round(f, 3), flush=True)
    return f >= 0.99


res = {}
for jet in ('off', 'max'):
    g = 1.0
    hi = None
    for _ in range(7):
        if lit(g, jet):
            hi = g
            break
        g *= 2.0
    if hi is None:
        res[jet] = dict(threshold_g=None, note='not lit up to 64 g (>= 99 % of grain)')
        continue
    lo = None
    g = hi
    while g > 1e-3 and lo is None:
        g /= 2.0
        if lit(g, jet):
            hi = g
        else:
            lo = g
    if lo is None:
        res[jet] = dict(threshold_g=hi, note='lit down to 1 mg')
        continue
    n = 0
    while hi / lo > 1.03 and n < 8:
        mid = (lo * hi) ** 0.5
        if lit(mid, jet):
            hi = mid
        else:
            lo = mid
        n += 1
    res[jet] = dict(threshold_g=hi, bracket=[lo, hi])
Path(out).write_text(json.dumps(dict(result=res, rows=rows), indent=1))
print(res)
