"""
srm1d_plugin.py — v0.8.0 Phase 5: srm_1d as an openMotor solver plugin (D1).

Registers srm_1d's 1-D PISO transient solver against openMotor's solver
registry (``motorlib.solvers``) so the GUI / simulation manager can select
it alongside the built-in quasi-steady solver. Importing this module
performs the registration.

The solver consumes an openMotor ``Motor`` (via ``Motor.getDict()`` — the
same ``{nozzle, propellant, grains, config}`` shape the ``.ric`` adapter
already converts, now carrying per-tab transport from Phase 3), runs
``run_simulation``, and maps the resulting channels back into an openMotor
``SimulationResult`` so the GUI renders it natively.

Scope: the scalar per-step channels (time, pressure, force, exit pressure,
throat change, Kn) and the per-grain multi-value channels (mass, massFlow,
massFlux, regression, web, machNumber) are populated — the latter aggregated
from srm_1d's per-cell axial snapshots via the cell→segment map and
interpolated onto the per-step time grid (Phase 6 task 3). The igniter is
read from the motor's own ``data.igniter`` block (Phase 6 task 4), falling
back to the built-in BPNV pyrogen for motors that predate the block. The GUI
solver picker (D6) is Phase 6 task 1, verified live in the openMotor app.
"""

import math as _math
import threading
import time as _time

import numpy as np

from .fmm_grain import _setup_openmotor_path
from .igniter_plenum import _VALID_TOPOLOGIES, _VALID_AREA_LAWS

SOLVER_NAME = 'srm_1d-transient'
ROUGHNESS_KEY = 'roughness'   # edited in um in the GUI; run_simulation wants m
IGN_H_KEY = 'ignition_heat_transfer'   # choice field; default read from run_simulation
# Allowed values of ``ignition_heat_transfer`` (srm_1d.simulation._IGN_H_MODES;
# pinned by tests/test_solver_config_upgrade.py). Default first = the GUI default.
_IGNITION_H_CHOICES = ('peretz', 'gnielinski', 'keller')
SOLID_KEY = 'solid_heating_model'      # choice field; default read from run_simulation
# D1 (2026-10-31): EN1/RT1/HE1 options, choice fields whose defaults are read from
# run_simulation ('total' / 'conservative' / 'recovery' / 'max'); each tuple's first entry is
# the current default and the rest are the alternatives (the pre-D1 behaviour included).
ENERGY_KEY = 'energy_formulation'
MOMENTUM_KEY = 'momentum_formulation'
DRIVER_KEY = 'ignition_driver_temperature'
JET_KEY = 'head_end_jet_h'
NOZZLE_BC_KEY = 'nozzle_bc_stagnation'   # NB1 choice field; default read from run_simulation
_ENERGY_CHOICES = ('total', 'sensible')
_MOMENTUM_CHOICES = ('conservative', 'legacy')
_DRIVER_CHOICES = ('recovery', 'static')
_JET_CHOICES = ('max', 'replace', 'off')
_NOZZLE_BC_CHOICES = ('isentropic', 'static')   # NB1; first = default
# srm_1d.simulation._SOLID_MODES (pinned by tests/test_solver_config_upgrade.py); default first.
_SOLID_HEATING_CHOICES = ('conduction', 'peretz', 'goodman')

# The transient solver emits one sample per solver step (10^5+ for a typical
# burn). openMotor's graph widget re-converts (per-point, in Python) and
# re-plots the whole series on every channel switch, so handing it the raw
# resolution causes multi-second hangs. We decimate the per-step channels to
# this many points for the GUI SimulationResult (a display artifact), while
# always preserving the peak-pressure / peak-thrust / endpoint samples so
# peak P, burn time, and impulse stay faithful. srm_1d's own full-resolution
# channels are untouched (analysis / CSV use those).
GUI_MAX_POINTS = 5000


def _om():
    """Lazily import the openMotor modules this plugin builds on."""
    _setup_openmotor_path()
    from motorlib.solvers import SolverPlugin, register_solver  # type: ignore
    from motorlib.simResult import SimulationResult  # type: ignore
    return SolverPlugin, register_solver, SimulationResult


def _transient_config_schema():
    """Build the srm_1d transient solver's run-parameter schema as an
    openMotor ``PropertyCollection``. Property keys match ``run_simulation``
    keyword arguments so the GUI-collected dict can be passed straight through
    as the solver ``config`` (→ ``simulate_motor`` overrides). Defaults mirror
    ``run_simulation``'s. ``ambPressure`` is intentionally absent — it stays in
    the shared global config (the adapter reads it from the motor)."""
    _setup_openmotor_path()
    from motorlib.properties import (  # type: ignore
        PropertyCollection, FloatProperty, EnumProperty,
    )
    schema = PropertyCollection()
    schema.props['t_max'] = FloatProperty('Max Simulation Time', 's', 0.01, 600.0)
    schema.props['P_cutoff'] = FloatProperty('Pressure Cutoff', 'Pa', 1.0e3, 5.0e6)
    schema.props['cfl_target'] = FloatProperty('CFL Target', '', 0.01, 1.0)
    # Roughness is edited in micrometers (its own unit category, default um),
    # converted to metres at the run boundary (see ROUGHNESS_KEY handling).
    schema.props['roughness'] = FloatProperty('Surface Roughness', 'um', 0.0, 1000.0)
    # kappa is the Gnielinski (Ma Eq. 9) temperature-ratio exponent in the
    # turbulent Nusselt correction (T_gas/T_surface)^kappa; recommended ~0.45
    # for heated gas. It is NOT an erosive-burning coefficient.
    schema.props['kappa'] = FloatProperty('Gnielinski Temp-Ratio Exponent', '', 0.0, 2.0)
    schema.props['T_ignition'] = FloatProperty('Ignition Temperature', 'K', 300.0, 2000.0)
    # IG3: closure for the convective h that heats the UNIGNITED propellant
    # surface (choice field; run_simulation validates the value). The Ma 2020
    # erosive burn-rate path is unaffected by this choice.
    schema.props[IGN_H_KEY] = EnumProperty(
        'Ignition Heat Transfer', list(_IGNITION_H_CHOICES))
    # How the unignited propellant surface temperature is advanced toward T_ignition
    # (conduction = exact 1-D transient conduction, the default since 2026-10-13).
    schema.props[SOLID_KEY] = EnumProperty(
        'Solid Heating Model', list(_SOLID_HEATING_CHOICES))
    # D1: total-energy PISO step, conservative momentum predictor, recovery-temperature ignition
    # driver and head-end igniter-jet heat transfer (docs/solver_verification/worklog/{EN1,RT1}.md,
    # docs/pyrogen_realism/HE1_REPORT.md). Each alternative reproduces the pre-D1 behaviour.
    schema.props[ENERGY_KEY] = EnumProperty('Energy Equation', list(_ENERGY_CHOICES))
    schema.props[MOMENTUM_KEY] = EnumProperty('Momentum Predictor', list(_MOMENTUM_CHOICES))
    schema.props[DRIVER_KEY] = EnumProperty('Ignition Driver Temperature', list(_DRIVER_CHOICES))
    schema.props[JET_KEY] = EnumProperty('Head-End Igniter-Jet Heat Transfer', list(_JET_CHOICES))
    # NB1: nozzle-boundary stagnation state ('static' = the pre-NB1 static-as-stagnation boundary).
    schema.props[NOZZLE_BC_KEY] = EnumProperty('Nozzle Boundary Stagnation', list(_NOZZLE_BC_CHOICES))
    # Seed defaults: DERIVED from run_simulation's signature (so they cannot
    # drift from the solver again -- they once did, 50 um / 0.45 / 850 K were
    # the pre-v0.7.5 values), roughness converted m -> um at this boundary.
    # NOTE: transport variant is a per-propellant property (propellant editor /
    # .ric), not a run-control param.
    schema.setProperties(_solver_config_defaults())
    return schema


# Pre-v0.7.5 run defaults (GUI units: roughness in um), superseded by the
# v0.7.5 re-LHS (roughness 32 um / kappa 0.44 / T_ignition 756 K). A saved GUI
# solver config (preferences or per-motor) that still holds exactly one of these
# is an UNTOUCHED stale default and is upgraded to the current default; any
# other value is a deliberate user choice and is never touched. Append, never
# edit, when a default changes again.
SUPERSEDED_SOLVER_CONFIG = {
    'roughness': [50.0],
    'kappa': [0.45],
    'T_ignition': [850.0],
}


def _solver_config_defaults():
    """``{key: default}`` in GUI units, read from ``run_simulation``'s
    signature for every key of the config schema."""
    import inspect
    from .simulation import run_simulation
    params = inspect.signature(run_simulation).parameters
    keys = ('t_max', 'P_cutoff', 'cfl_target', 'roughness', 'kappa', 'T_ignition')
    out = {k: float(params[k].default) for k in keys}
    out[ROUGHNESS_KEY] = round(out[ROUGHNESS_KEY] * 1.0e6, 9)   # m -> um
    out[IGN_H_KEY] = str(params[IGN_H_KEY].default)             # IG3 choice field
    out[SOLID_KEY] = str(params[SOLID_KEY].default)             # choice field
    for key in (ENERGY_KEY, MOMENTUM_KEY, DRIVER_KEY, JET_KEY, NOZZLE_BC_KEY):  # D1/NB1 choice fields
        d = params[key].default
        out[key] = 'off' if d is False else str(d)
    return out


def upgrade_solver_config(values):
    """Upgrade stale GUI solver-config defaults.

    ``values`` is a saved run-config dict (GUI units). Returns
    ``(new_values, upgraded_keys)``: each key whose value equals (rel tol 1e-9)
    a :data:`SUPERSEDED_SOLVER_CONFIG` default becomes the current default;
    everything else is untouched. Pure (never mutates ``values``) and
    idempotent."""
    import math
    new = dict(values or {})
    current = _solver_config_defaults()
    upgraded = []
    for key, old_values in SUPERSEDED_SOLVER_CONFIG.items():
        v = new.get(key)
        if isinstance(v, (int, float)) and not isinstance(v, bool) and any(
                math.isclose(float(v), o, rel_tol=1e-9, abs_tol=0.0)
                for o in old_values):
            new[key] = current[key]
            upgraded.append(key)
    return new, upgraded


# ----------------------------------------------------------------
# Plugin-provided igniter (docs/igniter_provider/DESIGN.md, D1/D2)
# ----------------------------------------------------------------
# srm_1d owns the igniter model: the material fields (``Pyrogen``), the chamber
# sizing (``_IGNITER_SIZING_DEFAULTS``), the built-in material library and the
# GUI metadata below. The openMotor fork hardcodes none of it; it asks the
# registered solver through the optional ``get_igniter_*`` hooks.

# Allowed enum values. ``injection_topology`` / ``burn_law`` /
# ``heat_delivery_mode`` are VALIDATED by srm_1d itself, so they are sourced
# from (or pinned by tests/test_igniter_provider.py to) srm_1d's validators:
# ``igniter_plenum._VALID_TOPOLOGIES``, ``igniter_plenum._burn_law_code`` and
# ``simulation._heat_delivery_code``. ``form`` is informational only (not
# validated anywhere in srm_1d), so its list is the fork's former one.
_PYROGEN_FORMS = ('pellets', 'powder', 'chunks')
_HEAT_DELIVERY_MODES = ('demar', 'radiation', 'none')
_BURN_LAWS = ('0d', 'end_burning')
_BASKET_AREA_LAWS = _VALID_AREA_LAWS

# Per ``Pyrogen`` init field -> GUI metadata. Kinds are implied by the keys:
# ``choices`` -> EnumProperty, ``min``/``max`` -> FloatProperty, else a
# StringProperty. Labels / units / ranges are ported verbatim from the fork's
# former ``motorlib/igniter.py`` ``Pyrogen`` class; ``gas_mass_fraction`` is new.
# Order = ``Pyrogen`` dataclass order (minus the passthrough fields).
_PYROGEN_FIELD_META = {
    'name': {'label': 'Name'},
    'a': {'label': 'Burn rate Coefficient', 'unit': 'm/(s*Pa^n)',
          'min': 1e-9, 'max': 2},
    'n': {'label': 'Burn rate Exponent', 'unit': '', 'min': -0.99, 'max': 0.99},
    'rho': {'label': 'Density', 'unit': 'kg/m^3', 'min': 1, 'max': 10000},
    'T_flame': {'label': 'Flame Temperature', 'unit': 'K',
                'min': 1, 'max': 10000},
    'M': {'label': 'Molar Mass', 'unit': 'kg/mol', 'min': 1e-6, 'max': 1},
    'gamma': {'label': 'Specific Heat Ratio', 'unit': '',
              'min': 1 + 1e-6, 'max': 10},
    'impetus_W': {'label': 'Impetus', 'unit': 'psi*in^3/g',
                  'min': 0, 'max': 1e6},
    'heat_flux_cal_cm2_s': {'label': 'Heat Flux', 'unit': 'cal/(cm^2*s)',
                            'min': 0, 'max': 1e5},
    'gas_mass_fraction': {'label': 'Gas Mass Fraction', 'unit': '',
                          'min': 0, 'max': 1},
    'kappa_jet': {'label': 'Jet Decay Coefficient', 'unit': '',
                  'min': 0, 'max': 100},
    'form': {'label': 'Form', 'choices': _PYROGEN_FORMS},
    'particle_diameter_m': {'label': 'Particle Diameter', 'unit': 'm',
                            'min': 0, 'max': 1},
    'particle_LD_ratio': {'label': 'Particle L/D', 'unit': '',
                          'min': 0, 'max': 1000},
    'heat_delivery_mode': {'label': 'Heat Delivery',
                           'choices': _HEAT_DELIVERY_MODES},
    'pellet_emissivity': {'label': 'Pellet Emissivity', 'unit': '',
                          'min': 0, 'max': 1},
    'radiation_absorption_length_m': {
        'label': 'Radiation Absorption Length', 'unit': 'm',
        'min': 0, 'max': 1000},
    'P_plateau': {'label': 'Burn Rate Plateau Pressure', 'unit': 'Pa',
                  'min': 0, 'max': 1e10},
}

# ``Pyrogen`` fields the GUI schema deliberately does NOT carry. They still
# round-trip losslessly through the raw igniter dict (DESIGN D3/D4): the GUI
# never rewrites a key the user did not edit.
_PYROGEN_GUI_PASSTHROUGH = {
    'Cp_gas': 'optional product-gas Cp; None = derive from gamma/M, and a '
              'FloatProperty cannot hold None',
}

# Per ``_IGNITER_SIZING_DEFAULTS`` key -> GUI metadata for the chamber sizing
# (ported verbatim from the fork's former ``Igniter`` class). -1 is the "auto"
# sentinel on the sizing fields.
_IGNITER_CHAMBER_META = {
    'mass': {'label': 'Pyrogen Mass', 'unit': 'kg', 'min': -1, 'max': 100},
    'throat_area': {'label': 'Throat Area', 'unit': 'm^2',
                    'min': -1, 'max': 1},
    'volume': {'label': 'Volume', 'unit': 'm^3', 'min': -1, 'max': 1},
    'burn_area': {'label': 'Burn Area', 'unit': 'm^2', 'min': -1, 'max': 100},
    'burn_law': {'label': 'Burn Law', 'choices': _BURN_LAWS},
    'injection_topology': {'label': 'Injection Topology',
                           'choices': _VALID_TOPOLOGIES},
    'cartridge_length_m': {'label': 'Cartridge Length', 'unit': 'm',
                           'min': -1, 'max': 100},
    'basket_fill_fraction': {'label': 'Basket Fill Fraction', 'unit': '',
                             'min': 0, 'max': 1},
    'pellet_packing_fraction': {'label': 'Pellet Packing Fraction', 'unit': '',
                                'min': 0, 'max': 1},
    'basket_area_law': {'label': 'Basket Burn Area Law',
                        'choices': _BASKET_AREA_LAWS},
    # G1: -1 = by the first grain's forward-face inhibition
    'head_gap_m': {'label': 'Head Gap (closure to grain, -1 = auto)',
                   'unit': 'm', 'min': -1, 'max': 1},
}


def _igniter_schema(meta, defaults):
    """A FRESH openMotor ``PropertyCollection`` built from a metadata table,
    then seeded with ``defaults`` (a dict; keys outside ``meta`` are ignored)."""
    _setup_openmotor_path()
    from motorlib.properties import (  # type: ignore
        PropertyCollection, FloatProperty, EnumProperty, StringProperty,
    )
    schema = PropertyCollection()
    for key, m in meta.items():
        if 'choices' in m:
            schema.props[key] = EnumProperty(m['label'], list(m['choices']))
        elif 'min' in m:
            schema.props[key] = FloatProperty(
                m['label'], m['unit'], m['min'], m['max'])
        else:
            schema.props[key] = StringProperty(m['label'])
    schema.setProperties({k: v for k, v in defaults.items()
                          if k in meta and v is not None})
    return schema


def _pyrogen_default_values():
    """The ``Pyrogen`` dataclass's own field defaults (``kappa_jet`` 8.0, ...),
    so a material missing those keys displays what srm_1d would use."""
    import dataclasses
    from .propellant import Pyrogen
    return {f.name: f.default for f in dataclasses.fields(Pyrogen)
            if f.init and f.default is not dataclasses.MISSING}


def _igniter_material_schema():
    """Schema for ONE igniter material (``data.igniter.pyrogen``)."""
    return _igniter_schema(_PYROGEN_FIELD_META, _pyrogen_default_values())


def _igniter_chamber_schema():
    """Schema for the per-motor chamber sizing / topology."""
    from .openmotor_adapter import _IGNITER_SIZING_DEFAULTS
    return _igniter_schema(_IGNITER_CHAMBER_META, _IGNITER_SIZING_DEFAULTS)


def simulate_motor(motor, igniter_pyrogen=None, callback=None, **sim_overrides):
    """Run srm_1d on an openMotor ``Motor`` and return an openMotor
    ``SimulationResult``. Reusable headlessly without the registry.

    The igniter is read from the motor's own ``data.igniter`` block (Phase 4
    self-describing motors, now carried by openMotor's ``Motor`` — Phase 6
    task 4). Passing ``igniter_pyrogen='<name>'`` forces a named library
    pyrogen with auto chamber sizing instead, overriding the motor's block.
    Motors that predate the block fall back to the BPNV default.

    ``callback`` (the openMotor progress callback) is driven live: the @njit
    time loop publishes a 0..1 progress metric into a shared array while it
    runs (GIL released via ``nogil``), a poller thread here forwards it to
    ``callback``, and a truthy callback return requests a cooperative cancel.
    """
    from .openmotor_adapter import (
        ric_to_sim_args, build_pyrogen_chamber, load_pyrogen, load_igniter,
    )
    from .simulation import run_simulation
    from .nozzle import compute_motor_performance

    motor_dict = motor.getDict()
    # Transport is read from the propellant's per-tab schema (gas_props=None);
    # a sentinel/missing transport hard-faults (D7), consistent with the .ric path.
    args = ric_to_sim_args(motor_dict, gas_props=None, **sim_overrides)
    geo = args.pop('geo')
    prop = args.pop('propellant')
    nozzle = args['nozzle']

    # Igniter resolution: explicit override > motor's own block > BPNV default.
    if igniter_pyrogen is not None:
        pyro, sizing = load_pyrogen(igniter_pyrogen), {}
    else:
        loaded = load_igniter(motor_dict)
        if loaded is not None:
            pyro, sizing = loaded
        else:
            pyro, sizing = load_pyrogen('bpnv'), {}
    args['pyrogen_chamber'] = build_pyrogen_chamber(pyro, geo, nozzle, **sizing)
    args.setdefault('verbose', False)

    # v0.8.x station-viz: snapshots are the source of BOTH the per-cell axial
    # payload and the per-grain channels, so the run_simulation default 0.2 s
    # interval is far too coarse for the GUI — it heavily interpolates the
    # per-grain channels and under-samples the axial fields, making the station
    # traces jagged and eating the ignition transient. 1 ms frames resolve the
    # ms-scale fill/ignition transient in the slice viewer (whose slider steps
    # frame by frame). The GUI passes t_max = the "Max Simulation Time" setting
    # (default 10 s), not the burn time, so the frame budget is set against it:
    # 1 ms up to t_max = 20 s, coarser beyond. The snapshot buffer is allocated
    # with np.empty, so frames after burnout are never written.
    t_max_eff = float(args.get('t_max', 10.0))
    args.setdefault('snapshot_interval', max(0.001, t_max_eff / 20000.0))

    from .simulation import PyrogenAutoSizeError
    try:
        if callback is None:
            result = run_simulation(geo, prop, **args)
        else:
            result = _run_with_progress(run_simulation, geo, prop, args, callback)
    except PyrogenAutoSizeError as exc:
        # The auto-sized igniter charge never lit the motor: return a failed result with an
        # ERROR alert instead of letting the exception reach the GUI.
        return _failed_igniter_result(motor, exc)

    perf = compute_motor_performance(result, nozzle, prop)

    # openMotor's post-run consumers — motor-stats ``getPortRatio`` and the
    # grain burnback cross-section in ``resultsWidget`` — call grain-geometry
    # methods (``getPortArea`` → ``getFaceArea``) that need each grain's FMM
    # regression map / ``faceArea``, which are populated ONLY by
    # ``simulationSetup``. srm_1d runs its own FMM bridge and never calls it,
    # so the motor's own grain objects keep ``faceArea = None`` and crash
    # ``getFaceArea`` for FMM-class grains (Finocyl/Star/Custom/...). Mirror
    # openMotor's solver setup loop. (Base ``Grain.simulationSetup`` is a
    # no-op, so this is safe for every grain type.)
    for grain in motor.grains:
        grain.simulationSetup(motor.config)

    return _result_to_om_simresult(motor, result, perf, geo, prop)


def _run_with_progress(run_simulation, geo, prop, args, callback,
                       poll_interval=0.05):
    """Run ``run_simulation`` in a worker thread while polling its shared
    ``progress_state`` array and forwarding progress to ``callback``. A truthy
    callback return sets the cancel flag the @njit loop reads each step.

    The loop releases the GIL (``nogil=True``), so this thread runs
    concurrently and sees progress update live. Returns the worker's result;
    re-raises any exception the worker hit."""
    progress_state = np.zeros(2, dtype=np.float64)
    args['progress_state'] = progress_state
    holder = {}

    def _worker():
        try:
            holder['result'] = run_simulation(geo, prop, **args)
        except BaseException as exc:  # propagate to the caller's thread
            holder['error'] = exc

    # Progress display has three regimes, because no single signal predicts the
    # wall-clock cost of the whole run:
    #   • WARMUP — before the @njit loop publishes anything (JIT compile + the
    #     first geometry/burn setup), creep slowly to WARMUP_CAP so the bar
    #     shows life instead of sitting frozen at 0.
    #   • BURN (phys < TAIL_GATE) — follow the physics fraction. The in-loop
    #     metric is now web-depth-weighted (total web burned / total web), so it
    #     tracks the BULK burn and no longer snaps to the gate the instant a
    #     thin-web OD-taper tip cell burns through.
    #   • TAIL (phys >= TAIL_GATE) — the end-of-burn / depressurization. Keep
    #     following the (still-moving) physics tail, but add an ASYMPTOTIC
    #     wall-clock floor toward TAIL_CAP so the bar keeps creeping even if the
    #     physics value flattens — always moving, never a hard stall at a fixed
    #     %, decelerating over a horizon scaled by the burn-to-gate time.
    #     Everything is capped at TAIL_CAP while the worker runs; completion
    #     snaps to 1.0.
    WARMUP_CAP = 0.05
    WARMUP_DUR = 8.0          # ~JIT-compile time; creep 0 -> WARMUP_CAP over it
    TAIL_GATE = 0.9
    TAIL_CAP = 0.995          # asymptote; only true completion fills the bar
    displayed = 0.0
    tail_start = None
    tail_tau = None
    start = _time.monotonic()

    worker = threading.Thread(target=_worker, daemon=True)
    worker.start()
    while worker.is_alive():
        phys = float(progress_state[0])
        now = _time.monotonic()
        if phys <= 0.0:
            # Setup / JIT-compile: indeterminate slow creep (bounded).
            target = min(WARMUP_CAP, (now - start) / WARMUP_DUR * WARMUP_CAP)
        else:
            target = phys                  # follow the (representative) burn
            if phys >= TAIL_GATE:
                if tail_start is None:     # entering the tail — calibrate
                    tail_start = now
                    tail_tau = max(1.5, 0.5 * (now - start))
                # Asymptotic anti-stall floor — decelerates but never stops.
                floor = TAIL_GATE + (TAIL_CAP - TAIL_GATE) * (
                    1.0 - _math.exp(-(now - tail_start) / tail_tau))
                if floor > target:
                    target = floor
            if target > TAIL_CAP:
                target = TAIL_CAP          # only true completion fills the bar
        if target > displayed:
            displayed = target             # monotonic
        if callback(displayed):
            progress_state[1] = 1.0         # request cooperative cancel
        _time.sleep(poll_interval)
    worker.join()
    if progress_state[1] <= 0.5:
        callback(1.0)                       # completed → fill the bar

    if 'error' in holder:
        raise holder['error']
    return holder['result']


def _per_grain_series(result, geo, prop):
    """Aggregate srm_1d's per-cell axial snapshot fields to per-grain series,
    sampled on the snapshot time base (Phase 6 task 3).

    Returns ``(snap_times, series, vol_loading)`` where ``series`` maps each
    openMotor multi-value channel name to an ``(n_frames, n_grains)`` array
    and ``vol_loading`` is the ``(n_frames,)`` scalar volume-loading series:

    - ``mass``       — solid propellant mass per grain, integrated from the
      per-cell port diameter: ρ·dx·Σ (A_outer − A_port).
    - ``regression`` / ``web`` — per-grain radial regression / remaining web
      (already computed per-grain in ``result['grains']``).
    - ``machNumber`` — peak core Mach over each grain's cells.
    - ``massFlow``   — cumulative downstream mass-loss rate (forward→aft sum
      of each grain's −d(mass)/dt), matching openMotor's convention.
    - ``massFlux``   — ``massFlow`` divided by the grain's aft port area.
    - ``vol_loading`` — 100 · (solid propellant volume / grain bounding
      volume), the per-step scalar openMotor reports as Volume Loading.

    Returns ``None`` if the result carries no snapshots (nothing to map).
    """
    snapshots = result['snapshots'] if 'snapshots' in result else None
    grains = result['grains'] if 'grains' in result else None
    if not snapshots or not grains:
        return None

    ga = geo.compile_geometry_arrays()
    cell_seg = np.asarray(ga['cell_segment_id'])
    bscale = np.asarray(ga['cell_burn_scale'])
    n_grains = len(grains)
    n_frames = len(snapshots)
    dx = geo.dx
    A_outer = np.pi / 4.0 * geo.D_outer ** 2
    rho = float(prop.rho_propellant)

    # Per-grain cell masks + aft-most cell index (for the aft port area).
    seg_cells = [np.flatnonzero(cell_seg == k) for k in range(n_grains)]
    seg_aft = [int(cells[-1]) if cells.size else -1 for cells in seg_cells]

    snap_times = np.array([snap['t'] for snap in snapshots])
    mass = np.zeros((n_frames, n_grains))
    mach = np.zeros((n_frames, n_grains))
    aft_port_area = np.zeros((n_frames, n_grains))
    for s, snap in enumerate(snapshots):
        d_port = np.asarray(snap['D_port'])
        mach_field = np.asarray(snap['Mach'])
        a_port = np.pi / 4.0 * d_port ** 2
        solid = np.maximum(A_outer - a_port, 0.0)
        for k in range(n_grains):
            cells = seg_cells[k]
            if cells.size:
                # G2: x burn_scale (L_true / L_snapped) so the grain's mass is the .ric one
                mass[s, k] = rho * dx * float(np.sum(solid[cells] * bscale[cells]))
                mach[s, k] = float(np.max(mach_field[cells]))
            if seg_aft[k] >= 0:
                aft_port_area[s, k] = float(a_port[seg_aft[k]])

    regression = np.stack(
        [np.asarray(grains[k]['regression']) for k in range(n_grains)], axis=1)
    web = np.stack(
        [np.asarray(grains[k]['web']) for k in range(n_grains)], axis=1)

    # massFlow: each grain's mass-loss rate, accumulated forward→aft so a
    # grain reports the total flow passing through it (openMotor convention).
    dmass = np.zeros_like(mass)
    if n_frames > 1:
        dt = np.diff(snap_times)
        dt[dt <= 0.0] = np.nan  # guard against duplicate snapshot times
        rate = np.maximum(-(np.diff(mass, axis=0)) / dt[:, None], 0.0)
        dmass[1:] = np.nan_to_num(rate)
    mass_flow = np.cumsum(dmass, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        mass_flux = np.where(aft_port_area > 0.0,
                             mass_flow / aft_port_area, 0.0)

    # Volume loading: solid propellant volume / grain bounding volume. The
    # bounding volume is the grain cells' cylinder (A_outer · dx · n_cells),
    # matching openMotor's "fraction of the chamber occupied by propellant".
    n_grain_cells = float(np.sum(bscale[cell_seg >= 0]))      # G2: true length / dx
    bounding_volume = A_outer * dx * max(n_grain_cells, 1.0)
    solid_volume = mass.sum(axis=1) / rho
    vol_loading = 100.0 * solid_volume / bounding_volume

    return snap_times, {
        'mass': mass, 'regression': regression, 'web': web,
        'machNumber': mach, 'massFlow': mass_flow, 'massFlux': mass_flux,
    }, vol_loading


def _decimate_indices(n, p_head, thrust, max_points):
    """Return a sorted, unique index array selecting at most ~``max_points``
    samples from a length-``n`` series, always including the first, last, and
    peak-pressure / peak-thrust samples so the decimated trace preserves peak
    P, burn time, and (closely) impulse. Returns all indices when ``n`` is
    already under the cap."""
    if n <= max_points:
        return np.arange(n)
    stride = int(np.ceil(n / max_points))
    keep = set(range(0, n, stride))
    keep.update((0, n - 1))
    if p_head is not None and p_head.size:
        keep.add(int(np.argmax(p_head)))
    if thrust is not None and thrust.size:
        keep.add(int(np.argmax(thrust)))
    return np.array(sorted(keep))


def _axial_payload_for_gui(result):
    """Build the per-station axial payload the GUI station panel consumes,
    as **plain, GUI-friendly structures** (numpy arrays + dicts; no srm_1d
    types) attached to the SimulationResult as ``sr.srm1d_axial``.

    This is the v0.8.x station-viz data contract (design phase 3): the
    capability-gated station panel slices these per-cell field matrices on
    demand instead of routing through openMotor's fixed per-grain channels.
    The default fore/mid/aft station model is embedded so the GUI needs no
    srm_1d import to populate its selector.

    Returns ``None`` when the result lacks the axial contract (e.g. results
    produced before v0.8.x), so the GUI cleanly falls back to per-grain mode.
    """
    if 'cell_segment_id' not in result or not result.get('snapshots'):
        return None
    from .station_viz import build_axial_payload, default_stations
    # Carry the snapshots at (nearly) full resolution — the snapshot interval
    # (set in simulate_motor) already bounds the frame count (~2000), so a high
    # cap avoids a second decimation stage that would re-introduce jaggedness.
    payload = build_axial_payload(result, max_frames=4000)
    if payload is None:
        return None
    # G2: station positions are PHYSICAL (cell axis + the lumped inter-grain gaps upstream)
    stations = default_stations(payload.cell_segment_id, payload.x_cell + payload.x_offset)
    return {
        'snap_times': payload.snap_times,
        'x_cell': payload.x_cell,
        'x_offset': payload.x_offset,   # G2: physical axis = x_cell + x_offset
        'cell_segment_id': payload.cell_segment_id,
        'fields': payload.fields,
        # Roadmap #2 longitudinal motor-slice geometry (constant per run).
        'dx': payload.dx,
        'D_outer': payload.D_outer,
        'cell_D_outer': payload.cell_D_outer,   # per-cell OD-tapered casing
        'cell_wall_web': payload.cell_wall_web,
        'seg_geom': payload.seg_geom,   # per-segment end-face burnback
        'stations': [
            {
                'grain': s.grain, 'cell_index': s.cell_index,
                'position_m': s.position_m, 'active': s.active,
                'role': s.role, 'label': s.label,
            }
            for s in stations
        ],
    }


def _failed_igniter_result(motor, exc):
    """A failed openMotor SimulationResult (``success=False``) with an ERROR alert naming
    the last (largest) auto charge tried; used when the igniter auto-sizing gave up."""
    _, _, SimulationResult = _om()
    from motorlib.simResult import SimAlert, SimAlertLevel, SimAlertType  # type: ignore
    sr = SimulationResult(motor)
    sr.success = False
    masses = getattr(exc, 'masses_kg', ())
    last = ("{:.3g} g".format(masses[-1] * 1e3) if masses else "an unknown charge")
    first = ("{:.3g} g".format(masses[0] * 1e3) if masses else "?")
    sr.addAlert(SimAlert(
        SimAlertLevel.ERROR, SimAlertType.VALUE,
        "srm_1d: the igniter never lit the motor -- the auto-sized charge was raised "
        "from {} to {} ({} attempts) without success. Set the pyrogen charge explicitly "
        "in the igniter settings (the motor may need a much larger charge, or an igniter "
        "topology that reaches the head end) and re-run.".format(first, last, len(masses)),
        'srm_1d'))
    return sr


def _result_to_om_simresult(motor, result, perf, geo=None, prop=None):
    """Map srm_1d results (SimulationChannels) + performance into an
    openMotor SimulationResult. Populates the scalar per-step channels and —
    when ``geo``/``prop`` are supplied — the per-grain multi-value channels
    (interpolated from the snapshot time base onto the per-step time grid).

    Also attaches ``sr.srm1d_axial`` (the per-station axial payload) when the
    result carries the v0.8.x axial contract, so the GUI's station panel can
    slice per-cell fields without re-deriving them."""
    _, _, SimulationResult = _om()
    sr = SimulationResult(motor)

    time = np.asarray(result['time'])
    p_head = np.asarray(result['P_head'])
    p_exit = np.asarray(result['P_exit'])
    d_throat = np.asarray(result['D_throat'])
    thrust = np.asarray(perf['thrust'])
    kn = np.asarray(result['Kn']) if 'Kn' in result else None
    d0 = float(d_throat[0]) if d_throat.size else 0.0

    # Per-grain channels live on the sparser snapshot time base; interpolate
    # each grain's series onto the per-step time grid so every channel shares
    # the 'time' channel's length (an openMotor SimulationResult invariant).
    grain_step = None
    vol_loading_step = None
    if geo is not None and prop is not None:
        pg = _per_grain_series(result, geo, prop)
        if pg is not None:
            snap_times, series, vol_loading = pg
            grain_step = {
                name: np.column_stack([
                    np.interp(time, snap_times, arr[:, k])
                    for k in range(arr.shape[1])
                ]) for name, arr in series.items()
            }
            vol_loading_step = np.interp(time, snap_times, vol_loading)

    # Decimate to a GUI-friendly sample count, keeping ballistics-critical
    # samples (peak pressure, peak thrust, first/last).
    idx = _decimate_indices(time.size, p_head, thrust, GUI_MAX_POINTS)

    for i in idx:
        sr.channels['time'].addData(float(time[i]))
        sr.channels['pressure'].addData(float(p_head[i]))
        sr.channels['exitPressure'].addData(float(p_exit[i]))
        sr.channels['force'].addData(float(thrust[i]))
        sr.channels['dThroat'].addData(float(d_throat[i] - d0))
        if kn is not None:
            sr.channels['kn'].addData(float(kn[i]))
        if vol_loading_step is not None:
            sr.channels['volumeLoading'].addData(float(vol_loading_step[i]))
        if grain_step is not None:
            for name, arr in grain_step.items():
                sr.channels[name].addData(tuple(float(v) for v in arr[i]))

    summary = result['summary'] if 'summary' in result else {}
    # Run-health watchdog: surface failed runs instead of silently plotting
    # nothing. A run is UNHEALTHY when it (a) numerically collapsed (code 4),
    # (b) was canceled (code 5), or (c) reached t_max (code 0) having never
    # pressurized — i.e. it never ignited (peak pressure near ambient). Note
    # ``t_burn`` is just the run duration, so it can't detect non-ignition;
    # peak pressure is the right signal. The ~0.3 MPa floor (≈3× ambient) is
    # far below any real solid-motor chamber pressure.
    term_code = summary.get('termination_code', 0)
    p_peak_mpa = float(summary.get('P_peak', 0.0) or 0.0) / 1.0e6
    canceled = term_code == 5
    collapsed = term_code == 4
    # IG3: code 6 = stopped early because ignition is impossible.
    no_ignition = (term_code == 6) or (
        (term_code == 0) and (p_peak_mpa < 0.3))
    # IG4: the run's own success flag (the grain, not the igniter, drives the
    # chamber after the igniter is spent). A hang-fire -- a lit cell or two
    # holding the chamber at the igniter's dribble pressure, ended by the
    # pressure cutoff -- has P_peak above the check above, so it is flagged here.
    dud = (summary.get('ignited') is False) and not (canceled or collapsed or no_ignition)
    sr.success = not (canceled or collapsed or no_ignition or dud)
    if not sr.success and not canceled:
        from motorlib.simResult import SimAlert, SimAlertLevel, SimAlertType
        if dud:
            n_lit = summary.get('grain_cells_ignited')
            n_all = summary.get('grain_cells')
            why = ("the motor did not reach self-sustaining combustion "
                   "(only {} of {} grain cells lit and the grain never out-"
                   "produced the igniter; a hang-fire)".format(n_lit, n_all))
        elif no_ignition:
            why = "the motor never ignited"
        else:
            why = "the solution collapsed"
        sr.addAlert(SimAlert(
            SimAlertLevel.ERROR, SimAlertType.VALUE,
            "srm_1d: no usable trace -- {} (P_peak {:.2f} MPa). Check the "
            "igniter / propellant / geometry.".format(why, p_peak_mpa),
            'srm_1d'))

    # AS1: the grain lit and sustains combustion but not PROMPTLY (fewer than the required
    # fraction of the grain ignited within the prompt window after igniter firing): the trace
    # is usable, so this is a WARNING, not an ERROR (a never-lit / hang-fire motor was already
    # reported above as an ERROR).
    if (summary.get('prompt_start') is False
            and not (canceled or collapsed or no_ignition or dud)):
        from motorlib.simResult import SimAlert, SimAlertLevel, SimAlertType
        sr.addAlert(SimAlert(
            SimAlertLevel.WARNING, SimAlertType.VALUE,
            "srm_1d: slow ignition -- only {:.0f} % of the grain lit within {:g} s of igniter "
            "firing (required {:.0f} %; first light {:.3g} ms, t90 {:.3g} ms, final igniter "
            "charge {:.3g} g). The motor does ignite, but the startup is not prompt: a larger "
            "igniter charge may be needed.".format(
                100.0 * float(summary.get('prompt_lit_fraction', 0.0)),
                float(summary.get('prompt_window_effective_s', 0.1)),
                100.0 * float(summary.get('prompt_fraction_required', 0.95)),
                1.0e3 * float(summary.get('first_ignition_time_s', float('nan'))),
                1.0e3 * float(summary.get('ignition_t90_s', float('nan'))),
                1.0e3 * float(summary.get('pyrogen_mass_final_kg', 0.0))),
            'srm_1d'))

    # v0.8.x station-viz: attach the per-station axial payload (capability-
    # gated side attribute; the GUI station panel consumes it, QS ignores it).
    axial = _axial_payload_for_gui(result)
    if axial is not None:
        sr.srm1d_axial = axial
    return sr


def register():
    """Register the srm_1d transient solver with openMotor's registry."""
    SolverPlugin, register_solver, _ = _om()

    class Srm1dTransientSolver(SolverPlugin):
        name = SOLVER_NAME
        capabilities = {
            'transient': True,
            'axial_fields': True,
            'needs_transport': True,
            'igniter': True,  # uses the pyrogen igniter (drives the GUI Igniter row)
        }

        def simulate(self, motor, config=None, callback=None):
            overrides = dict(config) if isinstance(config, dict) else {}
            # The config screen edits roughness in micrometers (its own unit
            # category); run_simulation wants metres. Convert at this boundary
            # so simulate_motor / run_simulation stay SI.
            if 'roughness' in overrides:
                overrides['roughness'] = overrides['roughness'] * 1.0e-6
            return simulate_motor(motor, callback=callback, **overrides)

        def get_config_schema(self):
            return _transient_config_schema()

        def upgrade_solver_config(self, values):
            return upgrade_solver_config(values)

        # ---- plugin-provided igniter hooks (DESIGN D2) ----
        def get_igniter_material_schema(self):
            return _igniter_material_schema()

        def get_igniter_chamber_schema(self):
            return _igniter_chamber_schema()

        def get_igniter_materials(self):
            from .openmotor_adapter import library_materials
            return library_materials()

        def get_default_igniter(self):
            from .openmotor_adapter import default_igniter_block
            return default_igniter_block('bpnv')

        def upgrade_igniter_block(self, block):
            from .openmotor_adapter import upgrade_igniter_block
            return upgrade_igniter_block(block)

    return register_solver(Srm1dTransientSolver())


# Register on import so `import srm_1d.srm1d_plugin` makes the solver
# selectable via motorlib.solvers.get_solver('srm_1d-transient').
try:
    register()
except Exception:  # pragma: no cover - openMotor checkout may be absent
    pass
