"""
srm_1d — 1D SRM Internal Ballistics Simulator
================================================

A transient 1D finite-volume solver for solid rocket motor internal
ballistics with the Ma et al. (2020) erosive burning model.

Quick start:
    from srm_1d import run_simulation
    from srm_1d.propellant import make_hasegawa_propellant_1
    from srm_1d.grain_geometry import make_hasegawa_motor_A_geo
    from srm_1d.nozzle import Nozzle, compute_motor_performance

    geo = make_hasegawa_motor_A_geo()
    prop = make_hasegawa_propellant_1()
    result = run_simulation(geo, prop, roughness=20e-6)
"""

__version__ = "0.5.0"

# Main entry point
from .simulation import run_simulation

# Geometry (grain side; nozzle is separate — see srm_1d.nozzle)
from .grain_geometry import (
    MotorGeometry,
    GrainSegment,
    make_bates_motor,
    make_single_cylinder,
    make_conical_grain,
    make_stepped_motor,
    make_hasegawa_motor_A_geo,
    make_hasegawa_motor_B_geo,
    make_hasegawa_motor_C_geo,
    make_hasegawa_motor_A_nozzle,
    make_hasegawa_motor_B_nozzle,
    make_hasegawa_motor_C_nozzle,
    make_example_bates,
)

# Propellant
from .propellant import (
    Propellant,
    PropellantTab,
    GasProperties,
    make_hasegawa_propellant_1,
    make_king_propellant_4525,
    create_gas_properties,
    critical_flow_function,
    characteristic_velocity,
)

# Nozzle
from .nozzle import (
    Nozzle,
    compute_motor_performance,
    print_performance_summary,
)

# Plotting (optional — requires matplotlib)
try:
    from .plotting import (
        plot_pressure,
        plot_thrust,
        plot_flow_snapshot,
        plot_summary,
        plot_comparison,
        plot_grain_regression,
        load_experimental_csv,
        HASEGAWA_MOTOR_A_EXPERIMENTAL,
    )
except ImportError:
    pass

# openMotor adapter (optional — requires pyyaml)
try:
    from .openmotor_adapter import (
        load_ric,
        run_from_ric,
        ric_to_sim_args,
        result_to_csv,
        save_csv,
        print_ric_summary,
        compute_grain_metrics,
        load_openmotor_csv,
    )
except ImportError:
    pass
