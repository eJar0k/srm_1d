"""
hasegawa_motor_a.py — Validate against Hasegawa et al. (2006) Motor A
=======================================================================

Runs the Hasegawa Motor A simulation and compares the head-end pressure
trace against digitized experimental data from Ma et al. (2020) Fig. 10.

Expected output:
    P_peak ≈ 6.2 MPa (experiment: 6.4 MPa, within 10%)
    Total impulse within 2% of experimental
    Burn time ≈ 4.2 s (experiment: ≈ 4.7 s)

Usage:
    python -m srm_1d.examples.hasegawa_motor_a

Outputs:
    hasegawa_a_pressure.png  — pressure trace vs experimental
    hasegawa_a_flow.png      — flow field snapshot at t ≈ 2s
    hasegawa_a_summary.png   — 4-panel summary
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for script use
import matplotlib.pyplot as plt

from srm_1d import run_simulation
from srm_1d.propellant import make_hasegawa_propellant_1
from srm_1d.grain_geometry import make_hasegawa_motor_A_geo
from srm_1d.nozzle import Nozzle, compute_motor_performance, print_performance_summary
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
    HASEGAWA_MOTOR_A_EXPERIMENTAL,
)


def main():
    # ============================================================
    # Run simulation
    # ============================================================
    geo = make_hasegawa_motor_A_geo()
    prop = make_hasegawa_propellant_1()
    prop.k_gas = 0.37      # between frozen (0.37) and effective (0.65)
    prop.Cp_gas = 2060.0    # between frozen (2060) and effective (1800)

    nozzle = Nozzle(
        D_throat=0.034,
        D_exit=0.050,
        div_angle=15.0,
        efficiency=0.95,
    )

    result = run_simulation(
        geo, prop, nozzle,
        roughness=15e-6,
        kappa=0.45,
        P_ignition=0.05e6,
        ignition_ramp_tau=0.010,
        P_cutoff=0.05e6,
        snapshot_interval=0.2,
        print_interval=0.2,
        igniter_mass = 0.010,
        # igniter_a,
        # igniter_n,
        # igniter_rho,
        igniter_A_burn = 10e-3,
    )

    # ============================================================
    # Nozzle performance
    # ============================================================
    perf = compute_motor_performance(result, nozzle, prop)
    print_performance_summary(perf, nozzle)

    # ============================================================
    # Plots
    # ============================================================
    # Pressure trace vs experimental
    plot_pressure(
        result,
        title="Motor A — 1D PISO vs Experimental (Hasegawa 2006)",
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        save_path="hasegawa_a_pressure.png",
    )

    # Flow snapshot at t ≈ 2s (mid-burn, before burnthrough)
    plot_flow_snapshot(
        result,
        t_target=2.0,
        save_path="hasegawa_a_flow.png",
    )

    # 4-panel summary
    plot_summary(
        result,
        performance=perf,
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        title="Hasegawa Motor A — Simulation Summary",
        save_path="hasegawa_a_summary.png",
    )

    plt.close('all')
    print("\nAll plots saved.")


if __name__ == '__main__':
    main()
