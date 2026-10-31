"""Print the HE1 grid from worklog/HE1_results/*.json (fixed-charge and auto-size rows)."""
import glob, json, sys
from pathlib import Path
d = Path(__file__).resolve().parents[1] / 'worklog' / 'HE1_results'
rows = {}
for f in sorted(glob.glob(str(d / '*.json'))):
    r = json.loads(Path(f).read_text())
    key = (r['motor'], r['solid'], r['case'], r['driver'], r['jet'], 'auto' if r['mass'] == 'auto' else 'fixed')
    rows[key] = r


def fmt(x, f='{:.2f}'):
    return 'n/a' if x is None or x != x else f.format(x)


for mass in ('fixed', 'auto'):
    for motor in ('chunc', 'balls', 'haseg'):
        for solid in ('conduction', 'goodman'):
            sub = {k: v for k, v in rows.items() if k[0] == motor and k[1] == solid and k[5] == mass}
            if not sub:
                continue
            print(f"\n## {motor} | solid={solid} | charge={mass}")
            for case, drv in (('sensible', 'static'), ('total_cons', 'recovery')):
                for jet in ('off', 'max', 'replace'):
                    r = sub.get((motor, solid, case, drv, jet, mass))
                    if not r:
                        continue
                    ig = r['ignition']
                    extra = ''
                    if motor == 'chunc':
                        extra = f" plateau {fmt(r.get('plateau'),'{:.3f}')} ratio {fmt(r.get('ratio'),'{:.3f}')}"
                    if motor == 'balls':
                        extra = f" tr {fmt(r.get('trough'),'{:.3f}')} mid {fmt(r.get('mid'),'{:.3f}')} pk/tr {fmt(r.get('pk_tr'),'{:.3f}')} pk/mid {fmt(r.get('pk_mid'),'{:.3f}')}"
                    if motor == 'haseg':
                        extra = f" late {fmt(r.get('late_peak'),'{:.3f}')} @ {fmt(r.get('t_late_peak'))} s"
                    print(f"{case:10s}/{drv:8s} jet={jet:7s} charge {r['pyrogen_mass_final_g']:.2f} g att {r['ignition_attempts']} | first {fmt(ig.get('t_first_ms'))} t90 {fmt(ig.get('t90_ms'))} ms | never lit {ig['grain_cells_never_lit']} {ig['grain_never_lit_cell_range']} | peak {r['peak']:.3f} @ {r['t_peak']*1e3:.1f} ms{extra}")
