"""
chunc.py — Validate against chunc data
=======================================================================

He's boosking it so hard

Expected output:
    Bomb

Usage:
    python -m srm_1d.examples.chunc

Outputs:
    chunc_pressure.png  — pressure trace vs experimental
    chunc_flow.png      — flow field snapshot at t ≈ 2s
    chunc_summary.png   — 4-panel summary
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for script use
import matplotlib.pyplot as plt
import numpy as np

from srm_1d import run_simulation
from srm_1d.propellant import Propellant
from srm_1d.grain_geometry import MotorGeometry, GrainSegment
from srm_1d.nozzle import Nozzle, compute_motor_performance, print_performance_summary
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
)

CHUNC_EXPERIMENTAL = {
    'label': 'Experimental (Chunc)',
    'time': np.array([
        0.0, 0.01, 0.044, 0.085, 0.126, 0.167, 0.207, 0.248, 0.289, 0.33,
        0.37, 0.411, 0.452, 0.493, 0.533, 0.574, 0.615, 0.655, 0.696,
        0.737, 0.778, 0.818, 0.859, 0.9, 0.941, 0.981, 1.022, 1.063,
        1.104, 1.144, 1.185, 1.226, 1.267, 1.307, 1.348, 1.389, 1.429,
        1.47, 1.511, 1.552, 1.592, 1.633, 1.674, 1.715, 1.755, 1.796,
        1.837, 1.878, 1.918, 1.959, 2, 2.041, 2.081, 2.122, 2.163, 2.203,
        2.244, 2.285, 2.326, 0, 
    ]),
    'pressure': np.array([
        0.0, 8.466000756, 8.735058654, 8.807471236, 8.829845452, 8.863448346,
        8.815005926, 8.807753257, 8.800600028, 8.83391275, 8.833599756,
        8.882640451, 8.840643762, 8.807647295, 8.788535068, 8.76042588,
        8.709474614, 8.670301396, 8.616598386, 8.569848093, 8.609606545,
        8.602280518, 8.591972891, 8.520044473, 8.411780979, 8.106008174,
        7.254074637, 6.251002249, 5.644352069, 5.241037462, 4.852299922,
        4.272797125, 3.858902652, 3.498468291, 3.111628281, 2.699448757,
        2.439607895, 2.199316148, 1.955314112, 1.690312101, 1.491945315,
        1.307850753, 1.149263612, 0.979278255, 0.802036967, 0.668733914,
        0.553478579, 0.449849686, 0.354674905, 0.288253234, 0.236817804,
        0.192674168, 0.150211248, 0.124883146, 0.10305993, 0.08138017,
        0.06160283, 0.048840154, 0.041499454, 0, 
    
    ]),
    'time_offset': 0.04,  # Align ignition events
}

def main():
    # ============================================================
    # Run simulation
    # ============================================================
    geo = make_chunc_geo()
    prop = make_boosk_prop()
    prop.k_gas = 0.6584      # between frozen (0.3654) and effective (0.6584)
    prop.Cp_gas = 2826.0    # between frozen (2050) and effective (2826)

    nozzle = Nozzle(
        D_throat=0.017017,
        D_exit=0.047625,
        div_angle=15.0,
        efficiency=0.95,
        erosion_coeff=184.2,
    )

    result = run_simulation(
        geo, prop, nozzle,
        roughness=30e-6,
        kappa=0.45,
        igniter_mass=0.001,  # 1g — chunc-specific calibration for the static fire trace
        P_ignition=0.05e6,
        ignition_ramp_tau=0.010,
        P_cutoff=0.05e6,
        snapshot_interval=0.2,
        print_interval=0.2,
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
        title="Chunc Sim vs Firing",
        experimental=CHUNC_EXPERIMENTAL,
        save_path="chunc_pressure.png",
    )

    # Flow snapshot at t ≈ 2s (mid-burn, before burnthrough)
    plot_flow_snapshot(
        result,
        t_target=0.5,
        save_path="chunc_flow.png",
    )

    # 4-panel summary
    plot_summary(
        result,
        performance=perf,
        experimental=CHUNC_EXPERIMENTAL,
        title="Chunc — Simulation Summary",
        save_path="chunc_summary.png",
    )

    plt.close('all')
    print("\nAll plots saved.")

def make_boosk_prop():
    """
    8448802384023475908345
    Data sources:
        openMotor gooning

    Gas transport from RPA at 5 MPa (frozen values at nozzle inlet):
        μ  = 9.125e-5 Pa·s,  k = 0.3654 W/(m·K),  Cp = 2050 J/(kg·K)

    Surface properties (from AP/HTPB/Al characterization literature):
        T_surface ≈ 1000 K (thermocouple measurement, typical for AP/HTPB)
        Cps ≈ 1500 J/(kg·K) (DSC measurement)

    These surface properties are NOT arbitrary, but the erosive burning
    prediction is moderately sensitive to (T_surface - T_initial).
    This pair is the least-constrained parameter in the model and a
    candidate for sensitivity analysis.
    """

    return Propellant(
        name="Boosk Prop 84",
        a=0.01569483e-3,
        n=0.4,
        rho_propellant=1700.0,
        Cps=1500.0,
        T_surface=1000.0,
        T_initial=293.0,
        T_flame=3105.0,
        gamma=1.19,
        molecular_weight=0.02498,
        mu_gas=9.125e-5,
        k_gas=0.3654,
        Cp_gas=2050.0,
    )

def make_chunc_geo():

    spacing = 0.001

    segments = []

    segments.append(GrainSegment(
            x_start=0, length=127e-3,
            D_bore_fwd=12.7e-3, D_outer=47.625e-3,
            D_bore_aft=12.e-3,
            inhibit_fwd=True, inhibit_aft=True,
        )),
    
    segments.append(GrainSegment(
            x_start=127e-3, length=685.8e-3,
            D_bore_fwd=14.478e-3, D_outer=47.625e-3,
            D_bore_aft=22.860e-3,
            inhibit_fwd=True, inhibit_aft=True,
        )),

    return MotorGeometry(
        L_motor=812.80e-3,
        D_outer=47.625e-3,
        segments=segments,
        N_cells=180,
    )


if __name__ == '__main__':
    main()