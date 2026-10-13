"""Print the EN1 duct-ladder table from worklog/EN1_results/duct_*.json."""
import json, glob, sys
from pathlib import Path
d = Path(__file__).resolve().parents[1] / 'worklog' / 'EN1_results'
rows = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(d / 'duct_*.json')))]
def tag(r): return ('total' if r['total'] else 'sensible') + ('+cons' if r['cons_mom'] else '')
for kind in ('venturi', 'shock'):
    print(f"\n== {kind} ==")
    for r in sorted([r for r in rows if r['kind'] == kind], key=lambda r: (tag(r), r['N'])):
        if r['diverged']:
            print(tag(r), r['N'], 'DIVERGED'); continue
        m = r['metrics']
        s = f"{tag(r):15s} N={r['N']:4d} resid {r['resid']:.1e} L1 M/T/P {m['L1_M']:.4f} {m['L1_T']:.4f} {m['L1_P']:.4f}"
        if kind == 'shock':
            xs = m['x_shock']
            s += (f" | up-L1 M/T/P {m['L1up_M']:.4f} {m['L1up_T']:.4f} {m['L1up_P']:.4f}"
                  f" | xs {xs[0]:.4f} (ex {xs[1]:.4f}) postP {m['post_shock_P_ratio'][0]:.4f}/{m['post_shock_P_ratio'][1]:.4f}"
                  f" postT {m['post_shock_T_ratio'][0]:.4f}/{m['post_shock_T_ratio'][1]:.4f}"
                  f" postM {m['post_shock_M'][0]:.3f}/{m['post_shock_M'][1]:.3f} preM {m['pre_shock_M'][0]:.3f}/{m['pre_shock_M'][1]:.3f}"
                  f" p0ratio {m['p0_ratio_across_shock'][0]:.4f}/{m['p0_ratio_across_shock'][1]:.4f}")
        else:
            s += f" | Mthroat {m['M_throat'][0]:.4f}/{m['M_throat'][1]:.4f} T/T0 throat {m['T_over_T0_throat'][0]:.4f}/{m['T_over_T0_throat'][1]:.4f} exit {m['T_over_T0_exit'][0]:.4f}/{m['T_over_T0_exit'][1]:.4f}"
        print(s)
