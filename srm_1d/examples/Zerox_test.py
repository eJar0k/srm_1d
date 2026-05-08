# basic test run with BATES motor

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.openmotor_adapter import run_from_ric, save_csv
from srm_1d.plotting import plot_pressure, plot_thrust, plot_flow_snapshot

ZEROX_EXPERIMENTAL = {
    'label': 'Experimental (Hasegawa)',
    'time': np.array([
        0, 0.155, 0.292, 0.43, 0.574, 0.719, 0.874, 1.025, 1.179, 1.316, 1.478, 1.61,
        1.754, 1.898, 2.032, 2.17, 2.31, 2.446, 2.588, 2.73, 2.867, 3.002, 3.147, 3.298,
        3.438, 3.572, 3.711, 3.852, 4.002, 4.143, 4.278, 4.417, 4.558, 4.7, 4.842, 4.974,
        5.115, 5.266, 5.41, 5.559, 5.694, 5.832, 5.985, 6.135, 6.276, 6.426, 6.56, 6.717,
        6.849, 6.991, 7.132, 7.284, 7.425, 7.559, 7.701, 7.834, 7.973, 8.113

    ]),
    'pressure': np.array([
        0, 0.0239, 0.4591, 3.9938, 3.863, 3.6684, 3.5838, 3.489, 3.4492, 3.3936, 3.3551, 3.2957,
        3.2805, 3.2319, 3.2073, 3.1725, 3.0986, 3.0645, 3.0449, 2.9944, 2.964, 2.9394, 2.868,
        2.8017, 2.6892, 2.6406, 2.5509, 2.4397, 2.3771, 2.2855, 2.2476, 2.2047, 2.1472, 2.1042,
        2.0562, 1.9899, 1.9779, 1.9507, 1.8686, 1.8029, 1.7043, 1.5546, 1.3101, 1.1705, 0.9961, 0.806,
        0.6777, 0.4794, 0.4017, 0.365, 0.2905, 0.2197, 0.1679, 0.1313, 0.1073, 0.051, 0.0599, 0.0302

    
    ]),
    'time_offset': -0.3,  # Align ignition events
}

result, perf, nozzle, geo, prop = run_from_ric(
    "Zerox V1.ric", 
    gas_props={'mu': 8.604e-5, 'k': 0.3455, 'Cp': 2052.0}, #N_cells = 177,
    # --- Igniter ---
    # 1. The Physical Erosive Friction
        roughness=20e-6,        # 37.1 microns
        kappa=0.45,               # Standard gas entrance effect
        
        # 2. The Numerical Ignition Smoothing
        igniter_mass=0.0024,      # 2.4 g
        igniter_tau=0.1269,       # 126.9 ms
        ignition_ramp_tau=0.0136, # 13.6 ms
        P_ignition=0.042e6,       # 0.042 MPa
        
        # General Solver Limits
        cfl_target=0.5,
        dt_max=1e-4,
        t_max=8.0,                # Adjust based on your expected BATES burn time
        P_cutoff=0.01e6,
        snapshot_interval=0.5,    # Take snapshots to watch the gap flows!
        print_interval=0.2,
)

# ============================================================
# Plots
# ============================================================
plot_pressure(result, title="Zerox Test",
                experimental=ZEROX_EXPERIMENTAL,
                time_offset=ZEROX_EXPERIMENTAL.get('time_offset', 0.0),
                save_path="Zerox_pressure.png")

plot_flow_snapshot(result, t_target=0.3,
                    title="Zerox — Flow at t ≈ 0.3s",
                    save_path="Zerox_flow.png")

plot_thrust(result, perf, title="Zerox Test",
            save_path="Zerox_thrust.png")

plt.close('all')
print("\nAll plots saved.")

save_csv("output.csv", result, perf, geo=geo, propellant=prop)
print("\nCSV saved.")