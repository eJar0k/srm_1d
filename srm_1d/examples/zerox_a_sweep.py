"""
zerox_a_sweep.py — Sensitivity sweep of the propellant burn-rate `a`
=====================================================================

Sweeps a × {0.80, 0.85, ..., 1.05} on every PropellantTab in the loaded
zerox.ric and overlays the resulting head-end pressure traces against
experimental.

Parallelized via ProcessPoolExecutor (matches the pattern in
srm_1d/tools/sensitivity.py — Numba's cache=True is shared across
worker processes, so warmup is amortized).

Usage:
    "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" \\
        -m srm_1d.examples.zerox_a_sweep
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.openmotor_adapter import load_ric, load_transport, ric_to_sim_args
from srm_1d.simulation import run_simulation
from srm_1d.plotting import ZEROX_EXPERIMENTAL


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox.ric'
TRANSPORT_PATH = MOTOR_PATH.with_name('zerox.transport.yaml')

A_SCALES = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05]

SIM_KWARGS = dict(
    roughness=20e-6,
    kappa=0.45,
    igniter_mass=0.0024,
    igniter_tau=0.1269,
    ignition_ramp_tau=0.0136,
    P_ignition=0.042e6,
    cfl_target=0.5,
    dt_max=1e-4,
    t_max=8.0,
    P_cutoff=0.01e6,
    snapshot_interval=2.0,
    print_interval=10.0,   # suppress per-step prints — workers run in parallel
)


def _run_one_a_scale(args):
    """Top-level worker (must pickle into ProcessPoolExecutor children)."""
    a_scale, motor_path, transport_path, sim_kwargs = args
    motor = load_ric(motor_path)
    gas_props = load_transport(transport_path)
    sim_args = ric_to_sim_args(motor, gas_props=gas_props, **sim_kwargs)

    prop = sim_args['propellant']
    for tab in prop.tabs:
        tab.a = tab.a * a_scale

    geo = sim_args.pop('geo')
    propellant = sim_args.pop('propellant')
    result = run_simulation(geo, propellant, **sim_args)

    # Strip to what we plot — keeps return-value pickling cheap.
    return {
        'a_scale': a_scale,
        't': np.asarray(result['time']),
        'P_head': np.asarray(result['P_head']),
        'D_throat': np.asarray(result['D_throat']),
        'mass_produced': result['summary']['mass_produced'],
        'mass_propellant': result['summary']['propellant_mass'],
    }


def plateau_mean(t, P, t_lo=1.0, t_hi=2.5):
    mask = (t >= t_lo) & (t <= t_hi)
    return float(P[mask].mean()) if mask.any() else float('nan')


def main():
    n_workers = min(len(A_SCALES), os.cpu_count() or 1)
    print(f"Launching {len(A_SCALES)} parallel sims on {n_workers} workers...")

    work = [(scale, str(MOTOR_PATH), str(TRANSPORT_PATH), SIM_KWARGS)
            for scale in A_SCALES]
    runs = []
    with ProcessPoolExecutor(max_workers=n_workers) as exe:
        futures = {exe.submit(_run_one_a_scale, w): w[0] for w in work}
        for fut in as_completed(futures):
            scale = futures[fut]
            try:
                runs.append(fut.result())
                print(f"  done: a × {scale:.2f}")
            except Exception as exc:
                print(f"  FAILED: a × {scale:.2f} — {exc}")
    runs.sort(key=lambda r: r['a_scale'])

    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.viridis(np.linspace(0.10, 0.90, len(A_SCALES)))

    plateau_table = []
    for color, run in zip(colors, runs):
        t = run['t']
        P_mpa = run['P_head'] / 1e6
        ax.plot(t, P_mpa, '-', color=color, linewidth=1.8,
                label=f'a × {run["a_scale"]:.2f}')
        plateau_table.append((
            run['a_scale'],
            plateau_mean(t, P_mpa),
            float(P_mpa.max()),
            run['mass_produced'],
            run['mass_propellant'],
        ))

    t_exp = np.array(ZEROX_EXPERIMENTAL['time']) \
        + ZEROX_EXPERIMENTAL.get('time_offset', 0.0)
    p_exp = np.array(ZEROX_EXPERIMENTAL['pressure'])
    ax.plot(t_exp, p_exp, 'k-', linewidth=2.2, marker='o', markersize=3.5,
            markevery=3, label='Experimental')

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Head-End Pressure [MPa]', fontsize=12)
    ax.set_title('Zerox — Saint-Robert `a` sensitivity '
                 '(igniter/transport/roughness held fixed)',
                 fontsize=13)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 8.5)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    save_path = "zerox_a_sweep.png"
    fig.savefig(save_path, dpi=150)
    print(f"\nSaved {save_path}")
    plt.close(fig)

    # Compute integrated impulse-time as a conservation diagnostic.
    print("\nResults:")
    header = (f"  {'a × scale':<11}{'P_plateau [MPa]':<18}{'P_peak [MPa]':<14}"
              f"{'∫P dt [MPa·s]':<16}{'mass_prod [kg]':<16}{'mass_prop [kg]':<14}")
    print(header)
    for scale, P_plat, P_peak, m_prod, m_prop in plateau_table:
        run = next(r for r in runs if r['a_scale'] == scale)
        impulse_pt = float(np.trapezoid(run['P_head'] / 1e6, run['t']))
        print(f"  {scale:<11.2f}{P_plat:<18.2f}{P_peak:<14.2f}"
              f"{impulse_pt:<16.2f}{m_prod:<16.3f}{m_prop:<14.3f}")
    p_exp_plat = plateau_mean(t_exp, p_exp)
    impulse_exp = float(np.trapezoid(p_exp, t_exp))
    print(f"  {'experimental':<11}{p_exp_plat:<18.2f}{float(p_exp.max()):<14.2f}"
          f"{impulse_exp:<16.2f}{'(unknown)':<16}{'(unknown)':<14}")


if __name__ == '__main__':
    main()
