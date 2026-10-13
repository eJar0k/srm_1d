# SM1 Part 3.4: wall-time cost of each solid model, serial, warm numba cache, best of 2 (steps/s).
import sys, time, dataclasses
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
nodemar = dataclasses.replace(load_pyrogen('bpnv'), heat_delivery_mode='none')
CASES = {'BALLSstick dead zone 4.83 g (2 s)': dict(ric='BALLSstick.ric', pyrogen=nodemar, pyrogen_mass=4.83e-3, t_max=2.0, P_cutoff=0.01e6),
         'Chunc head basket 4.5 g BPNV (0.2 s)': dict(ric='machbusterNew.ric', pyrogen='bpnv', pyrogen_mass=4.5e-3, t_max=0.2, injection_topology='head_basket'),
         'Chunc canonical (1 s)': dict(ric='machbusterNew.ric', pyrogen='bpnv', t_max=1.0, P_cutoff=0.01e6)}
for name, c in CASES.items():
    c = dict(c); ric = c.pop('ric')
    for m in ('goodman', 'peretz', 'conduction'):
        best = None
        for _ in range(2):
            t0 = time.time()
            r = run_from_ric(str(ROOT / 'motors' / ric), roughness=32e-6, kappa=0.44, T_ignition=756.0, print_interval=0.0,
                             verbose=False, solid_heating_model=m, **c)[0]
            w = time.time() - t0; st = r['summary']['steps']
            best = (w, st) if best is None or w < best[0] else best
        print(f'{name} | {m} | wall {best[0]:.2f} s | steps {best[1]} | {best[1]/best[0]/1e3:.1f} k steps/s', flush=True)
