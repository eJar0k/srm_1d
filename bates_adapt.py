# basic test run with BATES motor

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric, save_csv
from srm_1d.plotting import plot_pressure, plot_thrust, plot_flow_snapshot

result, perf, nozzle, geo, prop = run_from_ric(
    "BATES Test.ric",
    gas_props={'mu': 8.842e-5, 'k': 0.3685, 'Cp': 2060.0},
)

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

save_csv("output.csv", result, perf, geo=geo, propellant=prop)
print("\nCSV saved.")