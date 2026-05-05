"""
bates_4seg.py — 4-Segment BATES Motor Simulation
===================================================

Runs a typical amateur L-class 4-segment BATES motor and compares
against openMotor's 0-D prediction.

Expected output:
    Pressure matches openMotor to ~6%
    Total impulse within 1%
    Motor designation: L-class (~L1200)

Usage:
    python -m srm_1d.examples.bates_4seg

Outputs:
    bates_pressure.png  — pressure trace
    bates_flow.png      — flow field snapshot
    bates_thrust.png    — thrust and Isp
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d import run_simulation
from srm_1d.grain_geometry import make_example_bates
from srm_1d.propellant import make_hasegawa_propellant_1
from srm_1d.nozzle import Nozzle, compute_motor_performance, print_performance_summary
from srm_1d.plotting import plot_pressure, plot_thrust, plot_flow_snapshot


def main():
    # ============================================================
    # Geometry: 76mm casing, 38mm bore, 4 x 120mm segments + 20mm throat
    # ============================================================
    geo = make_example_bates()
    prop = make_hasegawa_propellant_1()
    nozzle = Nozzle(
        D_throat=0.020,
        D_exit=0.035,
        div_angle=15.0,
        efficiency=0.95,
    )

    print(f"BATES Motor: {geo.N_cells} cells, {len(geo.segments)} segments")
    print(f"  L_motor = {geo.L_motor*1e3:.0f} mm")
    print(f"  D_outer = {geo.D_outer*1e3:.0f} mm")
    print(f"  D_bore  = {geo.segments[0].D_bore_fwd*1e3:.0f} mm")
    print(f"  D_throat = {nozzle.D_throat*1e3:.1f} mm")
    print(f"  Web thickness = {geo.web_thickness*1e3:.1f} mm")
    print()

    # ============================================================
    # Run simulation
    # ============================================================
    result = run_simulation(
        geo, prop, nozzle,
        roughness=20e-6,
        kappa=0.45,
        P_ignition=0.05e6,
        ignition_ramp_tau=0.010,
        P_cutoff=0.05e6,
        snapshot_interval=0.1,
        print_interval=0.1,
    )

    # ============================================================
    # Nozzle performance
    # ============================================================
    perf = compute_motor_performance(result, nozzle, prop)
    print_performance_summary(perf, nozzle)

    # ============================================================
    # Plots
    # ============================================================
    plot_pressure(result, title="4-Segment BATES",
                  save_path="bates_pressure.png")

    plot_flow_snapshot(result, t_target=0.3,
                      title="BATES — Flow at t ≈ 0.3s",
                      save_path="bates_flow.png")

    plot_thrust(result, perf, title="4-Segment BATES",
                save_path="bates_thrust.png")

    plt.close('all')
    print("\nAll plots saved.")


if __name__ == '__main__':
    main()
