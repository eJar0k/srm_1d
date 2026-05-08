import numpy as np
import concurrent.futures
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import qmc
from scipy.interpolate import interp1d

from srm_1d import run_simulation
from srm_1d.grain_geometry import make_hasegawa_motor_A_geo, make_hasegawa_motor_A_nozzle
from srm_1d.propellant import make_hasegawa_propellant_1
from srm_1d.openmotor_adapter import load_ric, ric_to_sim_args

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

# =====================================================================
# 1. Experimental Benchmark Setup (for Error Calculation)
# =====================================================================
t_exp = ZEROX_EXPERIMENTAL['time'] + ZEROX_EXPERIMENTAL['time_offset']
p_exp = ZEROX_EXPERIMENTAL['pressure']
exp_interp = interp1d(t_exp, p_exp, bounds_error=False, fill_value=0.0)

# =====================================================================
# 2. Worker Function
# =====================================================================
def run_single_lhs_case(config):
    idx, params = config
    
    # Unpack all 7 variables strictly matching the bounds dictionary order
    (roughness, ign_mass, ign_ramp, p_ign, 
     ign_tau) = params
    
    motor = load_ric("Zerox V1.ric")
    args = ric_to_sim_args(
    motor, 
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
    
    geo = args.pop('geo')
    prop = args.pop('propellant')
    nozzle = args['nozzle']
    
    
    result = run_simulation(
        geo=geo, propellant=prop, nozzle=nozzle,
        
        # The 7 Sampled Variables
        roughness=roughness,
        igniter_mass=ign_mass,
        ignition_ramp_tau=ign_ramp,
        P_ignition=p_ign,
        igniter_tau=ign_tau,
        kappa=0.45,
        
        # Locked constants
        cfl_target=0.5,
        dt_max=1e-4,
        t_max=8.0,
        P_cutoff=0.05e6,
        snapshot_interval=2.0, 
        print_interval=20.0, # Silence output
    )
    
    # --- Fitness Evaluation (Mean Squared Error) ---
    t_sim = result['time']
    p_sim = result['P_head'] / 1e6  
    
    if len(t_sim) < 100 or t_sim[-1] < 1.0:
        error = 1e6 
    else:
        p_sim_at_exp_times = np.interp(t_exp, t_sim, p_sim)
        valid_idx = t_exp > 0.01
        sq_errors = (p_sim_at_exp_times[valid_idx] - p_exp[valid_idx])**2
        error = np.mean(sq_errors)
        
    return idx, params, error, result

# =====================================================================
# 3. Main LHS Driver
# =====================================================================
def main():
    N_SAMPLES = 50  # Increased for the 7-dimensional space
    
    # Define Parameter Bounds: [Lower Bound, Upper Bound]
    bounds = {
        'roughness':         [5e-6, 50e-6],
        'igniter_mass':      [0.001, 0.050],
        'ignition_ramp_tau': [0.001, 0.10],
        'P_ignition':        [0.005e6, 0.1e6],
        'igniter_tau':       [0.001, 0.20],
        # 'slag_coeff':        [0.0, 0.005],   
        # 'kappa':             [0.35, 0.55]    
    }
    
    print(f"Generating {N_SAMPLES} Latin Hypercube samples across {len(bounds)} dimensions...")
    
    sampler = qmc.LatinHypercube(d=len(bounds))
    sample = sampler.random(n=N_SAMPLES)
    
    l_bounds = [b[0] for b in bounds.values()]
    u_bounds = [b[1] for b in bounds.values()]
    scaled_samples = qmc.scale(sample, l_bounds, u_bounds)
    
    configs = [(i, tuple(row)) for i, row in enumerate(scaled_samples)]
    results_list = []

    print(f"Executing {N_SAMPLES} parallel simulations...")
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = executor.map(run_single_lhs_case, configs)
        
        for i, res in enumerate(futures):
            idx, params, error, sim_data = res
            results_list.append((error, params, sim_data))
            if (i+1) % 25 == 0:
                print(f"  Completed {i+1}/{N_SAMPLES}")
                
    # =====================================================================
    # 4. Sorting and Plotting
    # =====================================================================
    results_list.sort(key=lambda x: x[0])
    
    print("\n" + "="*50)
    print("--- TOP 5 BEST FITS ---")
    print("="*50)
    
    # Comprehensive Console Output
    for rank in range(5):
        err, p, _ = results_list[rank]
        print(f"Rank {rank+1} (MSE: {err:.4f}):")
        print(f"  Roughness    = {p[0]*1e6:.1f} μm")
        print(f"  Ign Mass     = {p[1]*1000:.1f} g")
        print(f"  Ign Ramp Tau = {p[2]*1000:.1f} ms")
        print(f"  P_ignition   = {p[3]/1e6:.3f} MPa")
        print(f"  Ign Tau      = {p[4]*1000:.1f} ms")
        # print(f"  Slag Coeff   = {p[5]:.5f}")
        # print(f"  Kappa        = {p[6]:.3f}")
        print("-" * 30)

    # Plot the Top 5
    plt.figure(figsize=(14, 9))
    plt.plot(t_exp, p_exp, 'k.-', linewidth=2.5, zorder=10, label='Experimental (Zerox)')
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for rank in range(5):
        err, params, sim_data = results_list[rank]
        t_sim = sim_data['time']
        p_sim = sim_data['P_head'] / 1e6
        
        # Multi-line plot legend to handle all 7 variables cleanly
        lbl = (f"Rank {rank+1} | R={params[0]*1e6:.0f}μm, m={params[1]*1000:.0f}g, "
               f"τ_r={params[2]*1000:.0f}ms\n"
               f"             P_i={params[3]/1e6:.2f}MPa, τ_i={params[4]*1000:.0f}ms")
               
        plt.plot(t_sim, p_sim, color=colors[rank], linewidth=1.5, alpha=0.9, label=lbl)

    plt.tight_layout()
    plt.title("Zerox — 5-Variable LHS Optimization (Top 5 Results)", fontsize=16)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("Head-End Pressure [MPa]", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 5.5)
    
    # Adjust legend size and spacing so it doesn't overlap the curves too badly
    plt.legend(loc='upper right', fontsize=9, labelspacing=1.2)
    
    plt.tight_layout()
    save_path = "zerox_LHS_5var.png"
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"\nOptimization complete! Plot saved to {save_path}.")

if __name__ == "__main__":
    main()