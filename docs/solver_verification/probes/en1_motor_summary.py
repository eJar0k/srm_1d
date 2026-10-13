"""Print the EN1 motor tables from worklog/EN1_results/motors/*.json."""
import glob, json
from pathlib import Path
d = Path(__file__).resolve().parents[1] / 'worklog' / 'EN1_results' / 'motors'
rows = {}
for f in sorted(glob.glob(str(d / '*.json'))):
    r = json.loads(Path(f).read_text())
    rows[(r['motor'], r['case'], r['mass'])] = r
order = ['sensible', 'total', 'total_cons', 'sens_cons']
def g(r, k, fmt='{:.4f}'):
    v = r.get(k)
    return 'n/a' if v is None else fmt.format(v)
for motor in ('chunc', 'balls', 'haseg'):
    for mass in sorted({k[2] for k in rows if k[0] == motor}):
        print(f"\n## {motor}  charge={mass}")
        for case in order:
            r = rows.get((motor, case, mass))
            if not r: continue
            ig = r['ignition']
            line = (f"{case:11s} charge {r['pyrogen_mass_final_g']:.3f} g (attempts {r['ignition_attempts']}) peak {r['peak']:.4f} MPa @ {r['t_peak']*1e3:.2f} ms")
            if motor == 'chunc': line += f" plateau {r['plateau']:.4f} ratio {r['ratio']:.4f}"
            if motor == 'balls': line += f" trough {r['trough']:.4f} mid {r['mid']:.4f} pk/tr {r['pk_tr']:.4f} pk/mid {r['pk_mid']:.4f}"
            line += (f" | steps {r['steps']} wall {r['wall_time_s']:.0f}s | fillMach {g(r,'fill_mach_max','{:.2f}')} floor {r['floor_activations']}"
                     f" | first {ig.get('t_first_ms',float('nan')):.2f} ms t90 {ig.get('t90_ms',float('nan')):.2f} ms lit {ig['n_lit']}/{ig['n_cells']} neverlit {ig['grain_cells_never_lit']} {ig['grain_never_lit_cell_range']} ignited {ig['ignited']}")
            print(line)
