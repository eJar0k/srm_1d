"""
openmotor_adapter.py — openMotor .ric File Adapter
=====================================================

Reads openMotor .ric files (YAML format), converts motor definitions
into srm_1d data structures, runs simulations, and exports results
as openMotor-compatible CSV.

Standalone — does not require motorlib to be installed.

Usage:
    from srm_1d.openmotor_adapter import run_from_ric

    result, perf = run_from_ric(
        "my_motor.ric",
        gas_props={'mu': 8.842e-5, 'k': 0.3685, 'Cp': 2060.0},
    )

    # Or step-by-step:
    from srm_1d.openmotor_adapter import load_ric, ric_to_sim_args, result_to_csv

    motor = load_ric("my_motor.ric")
    args = ric_to_sim_args(motor, gas_props={...})
    nozzle = args['nozzle']
    geo = args.pop('geo')
    prop = args.pop('propellant')
    result = run_simulation(geo, prop, **args)
    csv_str = result_to_csv(result, perf)

Unit conversions from openMotor internal units:
    erosionCoeff: m/(s·Pa) → μm/(s·MPa)  (multiply by 1e12)
    slagCoeff:    kept as-is, (m·MPa)/s — verify against openMotor source
    propellant.m: g/mol → kg/mol  (divide by 1000)
    propellant.a: m/s per Pa^n  (same as ours, no conversion)
"""

import warnings
import numpy as np

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

from .propellant import Propellant, PropellantTab
from .grain_geometry import MotorGeometry, GrainSegment, make_bates_motor
from .nozzle import Nozzle, compute_motor_performance, print_performance_summary
from .simulation import run_simulation


# ================================================================
# .ric file reader
# ================================================================

def load_ric(filepath):
    """
    Load an openMotor .ric file.

    Parameters
    ----------
    filepath : str
        Path to the .ric file.

    Returns
    -------
    dict with keys 'grains', 'nozzle', 'propellant', 'config',
    'version'. Raw openMotor format, not yet converted.
    """
    if not HAS_YAML:
        raise ImportError(
            "PyYAML is required to read .ric files. "
            "Install with: pip install pyyaml"
        )

    with open(filepath, 'r') as f:
        # openMotor .ric files contain Python-specific YAML tags
        # (!!python/object/apply, !!python/tuple) that safe_load
        # rejects. Use a custom loader that converts these to plain
        # Python objects — we only need the data dict.
        class _RicLoader(yaml.SafeLoader):
            pass
        # Handle any unknown tag by returning its constructed value
        _RicLoader.add_multi_constructor(
            'tag:yaml.org,2002:python/',
            lambda loader, suffix, node: loader.construct_mapping(node)
            if isinstance(node, yaml.MappingNode)
            else loader.construct_sequence(node)
            if isinstance(node, yaml.SequenceNode)
            else loader.construct_scalar(node),
        )
        raw = yaml.load(f, Loader=_RicLoader)

    if 'data' not in raw:
        raise ValueError(f"Invalid .ric file: missing 'data' key in {filepath}")

    data = raw['data']
    required = ['grains', 'nozzle', 'propellant', 'config']
    for key in required:
        if key not in data:
            raise ValueError(f"Invalid .ric file: missing '{key}' in {filepath}")

    return {
        'grains': data['grains'],
        'nozzle': data['nozzle'],
        'propellant': data['propellant'],
        'config': data['config'],
        'version': raw.get('version', (0, 0, 0)),
    }


# ================================================================
# Conversion: openMotor → srm_1d data structures
# ================================================================

# Maps openMotor inhibitedEnds strings to (inhibit_fwd, inhibit_aft)
_INHIBIT_MAP = {
    'Neither': (False, False),
    'Top': (True, False),       # Top = forward (head end)
    'Bottom': (False, True),    # Bottom = aft (nozzle end)
    'Both': (True, True),
}


def _ric_tab_to_srm(ric_tab):
    """Convert one openMotor tab dict to srm_1d PropellantTab."""
    return PropellantTab(
        min_pressure=ric_tab.get('minPressure', 0.0),
        max_pressure=ric_tab.get('maxPressure', 2.0e7),
        a=ric_tab['a'],                 # m/s per Pa^n (same convention)
        n=ric_tab['n'],
        gamma=ric_tab['k'],             # openMotor: 'k' = γ
        T_flame=ric_tab['t'],           # openMotor: 't' = T_flame
        molecular_weight=ric_tab['m'] / 1000.0,  # g/mol → kg/mol
    )


def _representative_ric_tab(tabs):
    """Pick the tab spanning the widest pressure range — used to source
    transport-property estimation when gas_props isn't supplied."""
    if len(tabs) == 1:
        return tabs[0]
    return max(tabs,
               key=lambda t: t.get('maxPressure', 0) - t.get('minPressure', 0))


def convert_propellant(ric_propellant, gas_props=None):
    """
    Convert an openMotor propellant dict to srm_1d Propellant. All
    `tabs` in the .ric file are preserved one-to-one as PropellantTab
    entries (full multi-tab support since v0.4.0).

    Parameters
    ----------
    ric_propellant : dict
        From the .ric file's 'propellant' key.
    gas_props : dict or None
        Transport properties: {'mu': Pa·s, 'k': W/(m·K), 'Cp': J/(kg·K)}.
        If None, estimated properties are derived from the
        widest-range tab (less accurate).

    Returns
    -------
    Propellant
    """
    ric_tabs = ric_propellant['tabs']
    tabs = [_ric_tab_to_srm(t) for t in ric_tabs]

    name = ric_propellant.get('name', 'openMotor propellant')
    density = ric_propellant['density']

    if gas_props is not None:
        mu = gas_props['mu']
        k_gas = gas_props['k']
        Cp = gas_props['Cp']
    else:
        # Estimate transport from the widest-range tab
        from .propellant import create_gas_properties_estimated
        rep = _representative_ric_tab(ric_tabs)
        gamma = rep['k']
        mw_kgmol = rep['m'] / 1000.0
        est = create_gas_properties_estimated(gamma, mw_kgmol, rep['t'])
        mu = est.mu
        k_gas = est.k_thermal
        Cp = est.Cp

    return Propellant(
        name=name,
        tabs=tabs,
        rho_propellant=density,
        Cps=1500.0,        # Default solid Cp — not in .ric
        T_surface=1000.0,  # Default surface temp — not in .ric
        T_initial=293.0,
        mu_gas=mu,
        k_gas=k_gas,
        Cp_gas=Cp,
    )


def convert_geometry(ric_grains, N_cells=None, spacing=None,
                     fmm_map_dim=1001):
    """
    Convert openMotor grain list to srm_1d MotorGeometry.
    The throat lives on a separate Nozzle object — see convert_nozzle().

    Supports BATES + Conical (analytic) and all 7 of openMotor's FMM
    grain types (Finocyl, Star, Moonburner, X, C, D, Custom). FMM grains
    have their regression maps built via openMotor (see srm_1d.fmm_grain)
    and attached as `GrainSegment.fmm_table`.

    Parameters
    ----------
    ric_grains : list of dict
        From the .ric file's 'grains' key.
    N_cells : int or None
        Number of cells. If None, auto-computed.
    spacing : float or None
        Inter-segment gap [m]. If None, defaults to 5% of grain
        outer diameter.
    fmm_map_dim : int
        FMM regression-map resolution for FMM grain types. openMotor
        default is 1001. Higher = more accurate perimeter/port-area
        sampling but quadratically slower setup (skfmm.distance is
        O(mapDim²)). Sim-level config knob.

    Returns
    -------
    MotorGeometry
    """
    # BATES and Conical are analytic — handled directly. All FMM types
    # are dispatched to srm_1d.fmm_grain.from_ric_grain.
    segments_data = []

    for i, grain in enumerate(ric_grains):
        gtype = grain['type']
        props = grain['properties']
        inhibit_str = props.get('inhibitedEnds', 'Neither')
        if inhibit_str not in _INHIBIT_MAP:
            raise ValueError(
                f"Grain {i} has unknown inhibitedEnds '{inhibit_str}'. "
                f"Valid values: {list(_INHIBIT_MAP.keys())}"
            )
        inh_fwd, inh_aft = _INHIBIT_MAP[inhibit_str]

        if gtype == 'BATES':
            segments_data.append({
                'kind': 'analytic',
                'D_bore_fwd': props['coreDiameter'],
                'D_bore_aft': props['coreDiameter'],
                'D_outer': props['diameter'],
                'length': props['length'],
                'inhibit_fwd': inh_fwd,
                'inhibit_aft': inh_aft,
                'fmm_table': None,
            })
        elif gtype == 'Conical':
            segments_data.append({
                'kind': 'analytic',
                'D_bore_fwd': props['forwardCoreDiameter'],
                'D_bore_aft': props['aftCoreDiameter'],
                'D_outer': props['diameter'],
                'length': props['length'],
                'inhibit_fwd': inh_fwd,
                'inhibit_aft': inh_aft,
                'fmm_table': None,
            })
        else:
            # Try FMM dispatch. from_ric_grain raises with a clear
            # message if `gtype` isn't a registered FMM type.
            from .fmm_grain import from_ric_grain
            try:
                fmm_table = from_ric_grain(grain, map_dim=fmm_map_dim)
            except ValueError as e:
                raise ValueError(
                    f"Grain {i} has unsupported type '{gtype}'. "
                    f"BATES uses analytic; FMM types must be registered. "
                    f"Inner error: {e}"
                ) from e
            segments_data.append({
                'kind': 'fmm',
                # FMM has no circular bore; D_outer placeholder gets
                # overwritten by FmmTable in compile_geometry_arrays.
                'D_bore_fwd': props['diameter'],
                'D_bore_aft': props['diameter'],
                'D_outer': props['diameter'],
                'length': props['length'],
                'inhibit_fwd': inh_fwd,
                'inhibit_aft': inh_aft,
                'fmm_table': fmm_table,
            })

    if not segments_data:
        raise ValueError("No supported grain segments found in .ric file.")

    D_outer = segments_data[0]['D_outer']

    # Default gap: 5% of outer diameter, minimum 3mm
    if spacing is None:
        spacing = max(0.003, D_outer * 0.05)

    # Build geometry
    N_segments = len(segments_data)
    L_motor = sum(s['length'] for s in segments_data) + (N_segments + 1) * spacing

    segments = []
    x_cursor = spacing
    for sd in segments_data:
        segments.append(GrainSegment(
            x_start=x_cursor,
            length=sd['length'],
            D_bore_fwd=sd['D_bore_fwd'],
            D_bore_aft=sd['D_bore_aft'],
            D_outer=sd['D_outer'],
            inhibit_fwd=sd['inhibit_fwd'],
            inhibit_aft=sd['inhibit_aft'],
            fmm_table=sd['fmm_table'],
        ))
        x_cursor += sd['length'] + spacing

    # Auto cell count: ~25 cells per segment, ~3 cells per gap
    if N_cells is None:
        cells_per_seg = 25
        cells_per_gap = 3
        N_cells = N_segments * cells_per_seg + (N_segments + 1) * cells_per_gap
        N_cells = max(N_cells, 50)

    return MotorGeometry(
        L_motor=L_motor,
        D_outer=D_outer,
        segments=segments,
        N_cells=N_cells,
    )


def convert_nozzle(ric_nozzle):
    """
    Convert an openMotor nozzle dict to srm_1d Nozzle.

    Field mapping (openMotor → srm_1d):
        throat        → D_throat              (m, same units)
        exit          → D_exit                (m, same units)
        efficiency    → efficiency            (—)
        divAngle      → div_angle             (deg, same)
        convAngle     → conv_angle            (deg, same)
        throatLength  → throat_length         (m, same)
        erosionCoeff  → erosion_coeff         m/(s·Pa) → μm/(s·MPa) (×1e12)
        slagCoeff     → slag_coeff            (m·MPa)/s, same convention
    """
    # openMotor stores erosionCoeff in m/(s·Pa); our convention: μm/(s·MPa).
    # 1 m/(s·Pa) = 1e12 μm/(s·MPa).
    erosion_ours = ric_nozzle.get('erosionCoeff', 0.0) * 1e12
    # Slag matches the openMotor convention.
    slag_ours = ric_nozzle.get('slagCoeff', 0.0)

    return Nozzle(
        D_throat=ric_nozzle['throat'],
        D_exit=ric_nozzle['exit'],
        efficiency=ric_nozzle.get('efficiency', 0.95),
        div_angle=ric_nozzle.get('divAngle', 15.0),
        conv_angle=ric_nozzle.get('convAngle', 30.0),
        throat_length=ric_nozzle.get('throatLength', 0.0),
        erosion_coeff=erosion_ours,
        slag_coeff=slag_ours,
    )


# ================================================================
# High-level convenience functions
# ================================================================

def ric_to_sim_args(motor, gas_props=None, N_cells=None, spacing=None,
                    **sim_overrides):
    """
    Convert a loaded .ric motor dict to run_simulation keyword arguments.

    Parameters
    ----------
    motor : dict
        Output from load_ric().
    gas_props : dict or None
        Transport properties: {'mu': Pa·s, 'k': W/(m·K), 'Cp': J/(kg·K)}.
    N_cells : int or None
        Override cell count.
    spacing : float or None
        Override inter-segment gap [m].
    **sim_overrides
        Additional keyword arguments passed to run_simulation
        (e.g. roughness, kappa, igniter params).

    Returns
    -------
    dict of keyword arguments for run_simulation. Includes 'geo',
    'propellant', 'nozzle', 'P_ambient', and 'P_cutoff' by default.
    """
    prop = convert_propellant(motor['propellant'], gas_props)
    geo = convert_geometry(motor['grains'], N_cells=N_cells, spacing=spacing)
    nozzle = convert_nozzle(motor['nozzle'])

    P_amb = motor['config'].get('ambPressure', 101325.0)

    args = {
        'geo': geo,
        'propellant': prop,
        'nozzle': nozzle,
        'P_ambient': P_amb,
        'P_cutoff': P_amb * 5,
    }
    args.update(sim_overrides)

    return args


def run_from_ric(filepath, gas_props=None, N_cells=None, spacing=None,
                 **sim_overrides):
    """
    Load a .ric file, run the 1D simulation, compute performance.

    Returns
    -------
    result, perf, nozzle, geo, prop
    """
    motor = load_ric(filepath)
    args = ric_to_sim_args(
        motor, gas_props=gas_props,
        N_cells=N_cells, spacing=spacing, **sim_overrides,
    )

    geo = args.pop('geo')
    prop = args.pop('propellant')
    nozzle = args['nozzle']
    P_amb = args.get('P_ambient', 101325.0)
    result = run_simulation(geo, prop, **args)

    perf = compute_motor_performance(result, nozzle, prop, P_ambient=P_amb)
    print_performance_summary(perf, nozzle)

    return result, perf, nozzle, geo, prop


# ================================================================
# Per-grain metrics
# ================================================================

def compute_grain_metrics(result, geo, propellant):
    """
    Compute per-grain time histories from snapshot data.

    Returns regression depth, web remaining, and mass remaining for
    each grain at each snapshot time. Also computes Kn from the
    pressure trace using the equilibrium relation.

    Parameters
    ----------
    result : dict
        Output from run_simulation (must have 'snapshots').
    geo : MotorGeometry
        The geometry used for the simulation.
    propellant : Propellant
        The propellant used.

    Returns
    -------
    dict with:
        'snap_times': ndarray (n_snaps,)
        'regression': ndarray (n_snaps, n_grains) — radial regression [m]
        'web': ndarray (n_snaps, n_grains) — web remaining [m]
        'grain_mass': ndarray (n_snaps, n_grains) — propellant mass [kg]
        'kn': ndarray (n_time,) — Kn from equilibrium at every time step
        'kn_times': ndarray (n_time,) — corresponding times
    """
    from .propellant import critical_flow_function, R_UNIVERSAL

    snapshots = result.get('snapshots', [])
    ga = geo.compile_geometry_arrays()
    N_seg = ga['N_seg']
    dx = geo.dx

    n_snaps = len(snapshots)
    regression = np.zeros((n_snaps, N_seg))
    web = np.zeros((n_snaps, N_seg))
    grain_mass = np.zeros((n_snaps, N_seg))

    for s, snap in enumerate(snapshots):
        D_port = snap['D_port']
        for k in range(N_seg):
            # Cells belonging to this grain
            mask = ga['cell_segment_id'] == k
            if not np.any(mask):
                continue

            D_bore_init = np.mean(ga['cell_D_bore_init'][mask])
            D_cells = D_port[mask]
            avg_D = np.mean(D_cells)

            regression[s, k] = (avg_D - D_bore_init) / 2.0
            web[s, k] = max(0.0, (geo.D_outer - avg_D) / 2.0)

            # Mass remaining: annular volume × density
            # Use per-cell D_port for accuracy
            cell_volume = np.sum(
                np.pi / 4.0 * (geo.D_outer**2 - D_cells**2) * dx
            )
            grain_mass[s, k] = cell_volume * propellant.rho_propellant

    # Kn from equilibrium: P^(1-n) = rho_p * a * Kn * c*
    # → Kn = P^(1-n) / (rho_p * a * c*)
    rep_tab = propellant.representative_tab()
    R_spec = R_UNIVERSAL / rep_tab.molecular_weight
    Gamma = critical_flow_function(rep_tab.gamma)
    c_star = np.sqrt(R_spec * rep_tab.T_flame) / Gamma

    t_arr = result['time']
    P_arr = result['P_head']
    # NOTE: Kn equilibrium uses representative-tab a/n. For multi-tab
    # propellants, pressure-dependent Kn would require per-sample tab
    # lookup — left as a follow-up since most exports use a single tab.
    kn_denom = propellant.rho_propellant * rep_tab.a * c_star
    kn = np.where(
        P_arr > 1e4,
        np.power(P_arr, 1.0 - rep_tab.n) / kn_denom,
        0.0,
    )

    snap_times = np.array([s['t'] for s in snapshots])

    return {
        'snap_times': snap_times,
        'regression': regression,
        'web': web,
        'grain_mass': grain_mass,
        'kn': kn,
        'kn_times': t_arr,
    }


# ================================================================
# CSV export (openMotor-compatible)
# ================================================================

def result_to_csv(result, perf=None, geo=None, propellant=None,
                  dt_sample=None, separator=','):
    """
    Export simulation results as openMotor-compatible CSV.

    Parameters
    ----------
    result : dict
        Output from run_simulation.
    perf : dict or None
        Output from compute_motor_performance.
    geo : MotorGeometry or None
        If provided (with propellant), per-grain columns are included.
    propellant : Propellant or None
        Required alongside geo for grain metrics.
    dt_sample : float or None
        Resample to this timestep [s] for manageable file size.
        If None, uses 0.001s (1 kHz). openMotor default is 0.025s.
    separator : str
        Column separator.

    Returns
    -------
    str : CSV content.
    """
    t = result['time']
    P = result['P_head']
    P_exit = result['P_exit']
    D_throat = result['D_throat']

    # Compute Kn and grain metrics if geometry is available
    grain_met = None
    if geo is not None and propellant is not None:
        grain_met = compute_grain_metrics(result, geo, propellant)

    # Downsample
    if dt_sample is None:
        dt_sample = 0.001
    if dt_sample > 0 and len(t) > 1:
        dt_actual = t[1] - t[0]
        stride = max(1, int(dt_sample / dt_actual))
    else:
        stride = 1
    indices = np.arange(0, len(t), stride)

    # Build headers
    headers = ['Time (s)', 'Kn', 'Pressure (Pa)', 'Exit Pressure (Pa)',
               'dThroat (m)']
    if perf is not None:
        headers.extend(['Force (N)', 'Isp (s)'])
    if grain_met is not None:
        n_grains = grain_met['regression'].shape[1]
        for k in range(n_grains):
            headers.append(f'Regression G{k} (m)')
        for k in range(n_grains):
            headers.append(f'Web G{k} (m)')

    lines = [separator.join(headers)]

    # Interpolate per-grain data to time history if available
    reg_interp = None
    web_interp = None
    if grain_met is not None and len(grain_met['snap_times']) > 1:
        n_grains = grain_met['regression'].shape[1]
        reg_interp = np.zeros((len(t), n_grains))
        web_interp = np.zeros((len(t), n_grains))
        for k in range(n_grains):
            reg_interp[:, k] = np.interp(
                t, grain_met['snap_times'], grain_met['regression'][:, k]
            )
            web_interp[:, k] = np.interp(
                t, grain_met['snap_times'], grain_met['web'][:, k]
            )

    for i in indices:
        # Kn from equilibrium
        kn_val = grain_met['kn'][i] if grain_met is not None else 0.0

        row = [
            f"{t[i]:.6f}",
            f"{kn_val:.2f}",
            f"{P[i]:.2f}",
            f"{P_exit[i]:.2f}",
            f"{D_throat[i] - D_throat[0]:.8f}",
        ]
        if perf is not None:
            row.append(f"{perf['thrust'][i]:.2f}")
            row.append(f"{perf['Isp'][i]:.2f}")
        if reg_interp is not None:
            for k in range(n_grains):
                row.append(f"{reg_interp[i, k]:.8f}")
            for k in range(n_grains):
                row.append(f"{web_interp[i, k]:.8f}")
        lines.append(separator.join(row))

    return '\n'.join(lines)


def save_csv(filepath, result, perf=None, geo=None, propellant=None,
             dt_sample=None, separator=','):
    """Write result_to_csv output to a file."""
    csv = result_to_csv(result, perf, geo, propellant, dt_sample, separator)
    with open(filepath, 'w') as f:
        f.write(csv)
    n_lines = csv.count('\n')
    print(f"Saved {filepath} ({n_lines} rows, dt={dt_sample or 0.001:.3f}s)")


# ================================================================
# openMotor CSV loader (for comparison)
# ================================================================

def load_openmotor_csv(filepath, separator=','):
    """
    Load an openMotor CSV export for comparison plotting.

    Parameters
    ----------
    filepath : str
        Path to the CSV file exported from openMotor.

    Returns
    -------
    dict with numpy arrays keyed by column name (lowercase, no units).
    Always includes 'time' and 'pressure'. Other keys depend on what
    openMotor exported.
    """
    with open(filepath, 'r') as f:
        header_line = f.readline().strip()

    headers_raw = [h.strip() for h in header_line.split(separator)]

    # Normalize header names: "Time (s)" → "time", "Pressure (Pa)" → "pressure"
    def _normalize(h):
        h = h.split('(')[0].strip().lower()
        h = h.replace(' ', '_')
        return h

    col_names = [_normalize(h) for h in headers_raw]
    data = np.loadtxt(filepath, delimiter=separator, skiprows=1)

    result = {}
    for i, name in enumerate(col_names):
        if i < data.shape[1]:
            result[name] = data[:, i]

    return result


# ================================================================
# Summary / comparison helper
# ================================================================

def print_ric_summary(filepath):
    """Print a human-readable summary of a .ric file."""
    motor = load_ric(filepath)

    grains = motor['grains']
    noz = motor['nozzle']
    prop = motor['propellant']

    print(f"openMotor file: {filepath}")
    print(f"  Propellant: {prop.get('name', 'unnamed')}")
    print(f"    density = {prop['density']:.0f} kg/m³")
    for i, tab in enumerate(prop['tabs']):
        print(f"    tab {i}: a={tab['a']:.4e}  n={tab['n']:.2f}  "
              f"γ={tab['k']:.2f}  MW={tab['m']:.1f}g/mol  "
              f"T_flame={tab['t']:.0f}K  "
              f"P=[{tab.get('minPressure',0)/1e6:.1f}-{tab.get('maxPressure',0)/1e6:.1f}]MPa")

    print(f"  Grains: {len(grains)}")
    for i, g in enumerate(grains):
        p = g['properties']
        if g['type'] == 'Conical':
            bore_str = (f"D_bore_fwd={p['forwardCoreDiameter']*1e3:.1f}mm  "
                        f"D_bore_aft={p['aftCoreDiameter']*1e3:.1f}mm")
        elif 'coreDiameter' in p:
            bore_str = f"D_bore={p['coreDiameter']*1e3:.1f}mm"
        else:
            bore_str = "(FMM core)"
        print(f"    [{i}] {g['type']}: {bore_str}  "
              f"D_outer={p['diameter']*1e3:.1f}mm  L={p['length']*1e3:.1f}mm  "
              f"inhibited={p.get('inhibitedEnds', 'Neither')}")

    print(f"  Nozzle: throat={noz['throat']*1e3:.1f}mm  "
          f"exit={noz['exit']*1e3:.1f}mm  "
          f"div={noz.get('divAngle',15):.0f}°  "
          f"efficiency={noz.get('efficiency',0.95):.3f}")

    erosion = noz.get('erosionCoeff', 0)
    slag = noz.get('slagCoeff', 0)
    if erosion > 0 or slag > 0:
        print(f"    erosion={erosion:.2e} m/(s·Pa) "
              f"[= {erosion*1e12:.1f} μm/(s·MPa)]  "
              f"slag={slag:.2e} (m·MPa)/s")
    print()
