"""Print the RT1 grid (ignition driver temperature) from worklog/RT1_results/*.json."""
import glob, json
from pathlib import Path
d = Path(__file__).resolve().parents[1] / 'worklog' / 'RT1_results'
rows = {}
for f in sorted(glob.glob(str(d / '*.json'))):
    r = json.loads(Path(f).read_text())
    rows[(r['motor'], r['solid'], 'auto' if r['mass'] == 'auto' else 'fixed', r['case'], r['driver'])] = r
def fmt(x, f='{:.2f}'):
    return 'n/a' if x is None or x != x else f.format(x)
for motor in ('chunc', 'balls', 'haseg'):
    for solid in ('goodman', 'conduction'):
        for mass in ('fixed', 'auto'):
            sub = {k: v for k, v in rows.items() if k[:3] == (motor, solid, mass)}
            if not sub: continue
            print(f"\n## {motor} | solid={solid} | charge={mass}")
            for case in ('sensible', 'total_cons'):
                for drv in ('static', 'recovery'):
                    r = sub.get((motor, solid, mass, case, drv))
                    if not r: continue
                    ig = r['ignition']
                    extra = ''
                    if motor == 'chunc': extra = f" plateau {fmt(r.get('plateau'),'{:.3f}')} ratio {fmt(r.get('ratio'),'{:.3f}')}"
                    if motor == 'balls': extra = f" trough {fmt(r.get('trough'),'{:.3f}')} mid {fmt(r.get('mid'),'{:.3f}')} pk/tr {fmt(r.get('pk_tr'),'{:.3f}')}"
                    if motor == 'haseg': extra = f" (t_max {r['t_max']} s)"
                    print(f"{case:10s} {drv:8s} charge {r['pyrogen_mass_final_g']:.2f} g att {r['ignition_attempts']} | first {fmt(ig.get('t_first_ms'))} t90 {fmt(ig.get('t90_ms'))} ms | never lit {ig['grain_cells_never_lit']} {ig['grain_never_lit_cell_range']} ignited {ig['ignited']} | peak {r['peak']:.3f} MPa @ {r['t_peak']*1e3:.1f} ms{extra}")
