"""
srm_1d — 1D SRM Internal Ballistics Simulator
================================================

A transient 1D finite-volume solver for solid rocket motor internal
ballistics with the Ma et al. (2020) erosive burning model.

Quick start (v0.6.0+):
    from srm_1d.openmotor_adapter import run_from_ric

    result, perf, nozzle, geo, prop = run_from_ric(
        'srm_1d/motors/hasegawa_a.ric',
        roughness=37.1e-6, igniter_tau=0.1269, igniter_mass=0.0024,
    )

The ``srm_1d/motors/`` directory ships canonical motor specifications
as ``.ric`` files (openMotor schema) with sibling ``.transport.yaml``
files supplying combustion gas transport (mu, k, Cp).

To define a parametric geometry programmatically, build it directly via
``build_snapped_geometry`` from ``grain_geometry``.
"""

__version__ = "0.6.0"

# Main entry point
from .simulation import run_simulation

# Geometry (grain side; nozzle is separate — see srm_1d.nozzle)
from .grain_geometry import (
    MotorGeometry,
    GrainSegment,
    build_snapped_geometry,
)

# Propellant
from .propellant import (
    Propellant,
    PropellantTab,
    GasProperties,
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
        load_transport,
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
