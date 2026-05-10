"""
zerox_erosion_sweep.py — Throat-erosion sensitivity sweep
==========================================================

The a-sweep (zerox_a_sweep.py) revealed that ∫P dt is nearly
conservation-invariant across `a × {0.80..1.05}` (sim 19.5-19.7
MPa·s) but experimental is 16.4 MPa·s — a 20% gap. Per 0D
equilibrium ∫P dt ≈ m_total · c* / Ā_t. With m_total fixed and
c* drop unphysical, the gap must come from throat-area history.

This script tests whether `nozzle.erosion_coeff × scale` for
scale ∈ {1.0, 1.5, 2.0, 2.5, 3.0} closes the gap, holding
`a × 1.00` and the borrowed Hasegawa A igniter fixed.

Parallelized via ProcessPoolExecutor.

Usage:
    "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" \\
        -m srm_1d.examples.zerox_erosion_sweep
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.openmotor_adapter import (
    build_pyrogen_chamber,
    load_pyrogen,
    load_ric,
    load_transport,
    ric_to_sim_args,
)
from srm_1d.simulation import run_simulation
from srm_1d.plotting import ZEROX_EXPERIMENTAL


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox.ric'
TRANSPORT_PATH = MOTOR_PATH.with_name('zerox.transport.yaml')

EROSION_SCALES = [1.0, 1.5, 2.0, 2.5, 3.0]

SIM_KWARGS = dict(
    roughness=20e-6,
    kappa=0.45,
    T_ignition=850.0,
    cfl_target=0.5,
    dt_max=1e-4,
    t_max=8.0,
    P_cutoff=0.01e6,
    snapshot_interval=2.0,
    print_interval=10.0,
)


def _run_one_erosion(args):
    """Top-level worker (must pickle into ProcessPoolExecutor children)."""
    erosion_scale, motor_path, transport_path, sim_kwargs = args
    motor = load_ric(motor_path)
    gas_props = load_transport(transport_path)
    sim_args = ric_to_sim_args(motor, gas_props=gas_props, **sim_kwargs)
    sim_args['pyrogen_chamber'] = build_pyrogen_chamber(
        load_pyrogen('bpnv'), sim_args['geo'], sim_args['nozzle']
    )

    nozzle = sim_args['nozzle']
    nozzle.erosion_coeff = nozzle.erosion_coeff * erosion_scale

    geo = sim_args.pop('geo')
    propellant = sim_args.pop('propellant')
    result = run_simulation(geo, propellant, **sim_args)

    return {
        'erosion_scale': erosion_scale,
        'erosion_coeff': nozzle.erosion_coeff,
        't': np.asarray(result['time']),
        'P_head': np.asarray(result['P_head']),
        'D_throat': np.asarray(result['D_throat']),
        'mass_produced': result['summary']['mass_produced'],
        'mass_propellant': result['summary']['propellant_mass'],
        'D_throat_final': result['summary']['D_throat_final'],
    }


def plateau_mean(t, P, t_lo=1.0, t_hi=2.5):
    mask = (t >= t_lo) & (t <= t_hi)
    return float(P[mask].mean()) if mask.any() else float('nan')


def main():
    n_workers = min(len(EROSION_SCALES), os.cpu_count() or 1)
    print(f"Launching {len(EROSION_SCALES)} parallel sims on {n_workers} workers...")

    work = [(scale, str(MOTOR_PATH), str(TRANSPORT_PATH), SIM_KWARGS)
            for scale in EROSION_SCALES]
    runs = []
    with ProcessPoolExecutor(max_workers=n_workers) as exe:
        futures = {exe.submit(_run_one_erosion, w): w[0] for w in work}
        for fut in as_completed(futures):
            scale = futures[fut]
            try:
                runs.append(fut.result())
                print(f"  done: erosion × {scale:.2f}")
            except Exception as exc:
                print(f"  FAILED: erosion × {scale:.2f} — {exc}")
    runs.sort(key=lambda r: r['erosion_scale'])

    fig, axes = plt.subplots(2, 1, figsize=(11, 11), sharex=True)
    ax_p, ax_t = axes
    colors = plt.cm.plasma(np.linspace(0.10, 0.85, len(EROSION_SCALES)))

    for color, run in zip(colors, runs):
        t = run['t']
        P_mpa = run['P_head'] / 1e6
        D_mm = run['D_throat'] * 1000
        ax_p.plot(t, P_mpa, '-', color=color, linewidth=1.8,
                  label=f'erosion × {run["erosion_scale"]:.1f}')
        ax_t.plot(t, D_mm, '-', color=color, linewidth=1.8,
                  label=f'erosion × {run["erosion_scale"]:.1f}')

    t_exp = np.array(ZEROX_EXPERIMENTAL['time']) \
        + ZEROX_EXPERIMENTAL.get('time_offset', 0.0)
    p_exp = np.array(ZEROX_EXPERIMENTAL['pressure'])
    ax_p.plot(t_exp, p_exp, 'k-', linewidth=2.2, marker='o', markersize=3.5,
              markevery=3, label='Experimental')

    ax_p.set_ylabel('Head-End Pressure [MPa]', fontsize=12)
    ax_p.set_title('Zerox — throat-erosion sensitivity '
                   '(a × 1.00, igniter/transport/roughness held fixed)',
                   fontsize=13)
    ax_p.legend(loc='upper right', fontsize=10)
    ax_p.grid(True, alpha=0.3)
    ax_p.set_ylim(bottom=0)

    ax_t.set_xlabel('Time [s]', fontsize=12)
    ax_t.set_ylabel('Throat diameter [mm]', fontsize=12)
    ax_t.set_title('Throat opening rate vs t', fontsize=12)
    ax_t.legend(loc='lower right', fontsize=9)
    ax_t.grid(True, alpha=0.3)
    ax_t.set_xlim(0, 8.5)

    plt.tight_layout()
    save_path = "zerox_erosion_sweep.png"
    fig.savefig(save_path, dpi=150)
    print(f"\nSaved {save_path}")
    plt.close(fig)

    # Conservation table.
    print("\nResults:")
    header = (f"  {'erosion×':<10}{'P_plateau [MPa]':<18}{'P_peak [MPa]':<14}"
              f"{'∫P dt [MPa·s]':<16}{'D_t_final [mm]':<16}{'mass_prod [kg]':<14}")
    print(header)
    for run in runs:
        t = run['t']
        P_mpa = run['P_head'] / 1e6
        impulse_pt = float(np.trapezoid(P_mpa, t))
        print(f"  {run['erosion_scale']:<10.2f}"
              f"{plateau_mean(t, P_mpa):<18.2f}"
              f"{float(P_mpa.max()):<14.2f}"
              f"{impulse_pt:<16.2f}"
              f"{run['D_throat_final']*1000:<16.2f}"
              f"{run['mass_produced']:<14.3f}")
    p_exp_plat = plateau_mean(t_exp, p_exp)
    impulse_exp = float(np.trapezoid(p_exp, t_exp))
    print(f"  {'experimental':<10}{p_exp_plat:<18.2f}"
          f"{float(p_exp.max()):<14.2f}"
          f"{impulse_exp:<16.2f}"
          f"{'(unknown)':<16}{'(unknown)':<14}")


if __name__ == '__main__':
    main()
