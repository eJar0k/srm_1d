"""
srm_1d — 1D SRM Internal Ballistics Simulator
================================================

A transient 1D finite-volume solver for solid rocket motor internal
ballistics with the Ma et al. (2020) erosive burning model.

Quick start (v0.6.0+):
    from srm_1d.openmotor_adapter import run_from_ric

    result, perf, nozzle, geo, prop = run_from_ric(
        'srm_1d/motors/hasegawa_a.ric',
        roughness=37.1e-6, pyrogen='bpnv', T_ignition=850.0,
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

# v0.8.0 Phase 1 — channel result model (openMotor-aligned).
# Additive: run_simulation() still returns the results dict; build_channels
# re-shapes it into channels for the frontend-integration path.
from .channels import (
    Channel,
    AxialChannel,
    SimulationChannels,
    build_channels,
    as_channels,
)

# v0.8.x station-viz — headless backend for the per-station axial panel
# (payload extraction + default fore/mid/aft station model). Qt-free; the
# GUI panel on the openMotor-fork side consumes these.
from .station_viz import (
    AxialPayload,
    Station,
    build_axial_payload,
    default_stations,
    make_station,
    grain_cell_spans,
    gap_cell_indices,
    cell_categories,
    grain_role,
    classify_cell,
    station_full_label,
)

# Igniter
from .igniter_plenum import PyrogenChamber

# Geometry (grain side; nozzle is separate — see srm_1d.nozzle)
from .grain_geometry import (
    MotorGeometry,
    GrainSegment,
    build_snapped_geometry,
)

# FMM grains + parametric axial tapers (optional — taper resolution needs
# the openMotor checkout + scikit-fmm, but importing the module is light).
try:
    from .fmm_grain import (
        FmmTable,
        from_openmotor,
        from_ric_grain,
        TaperSpec,
        linear_taper,
        taper_profile,
        resolve_taper,
    )
except ImportError:
    pass

# Propellant
from .propellant import (
    Propellant,
    PropellantTab,
    Pyrogen,
    GasSpecies,
    GasProperties,
    create_gas_properties,
    species_array,
    ambient_air_species,
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
        plot_channels,
        plot_pressure,
        plot_thrust,
        plot_flow_snapshot,
        plot_flow_snapshots,
        plot_field_heatmap,
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
        load_pyrogen,
        run_from_ric,
        ric_to_sim_args,
        result_to_csv,
        save_csv,
        print_ric_summary,
        compute_grain_metrics,
        load_openmotor_csv,
        migrate_ric_transport,
        migrate_all_motors,
        build_transport_library,
        load_igniter,
        default_igniter_block,
        build_pyrogen_chamber,
    )
except ImportError:
    pass
