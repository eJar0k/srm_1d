"""
plotting.py — Visualization Utilities
=======================================

Quick plotting for simulation results. This is a convenience module,
not part of the solver — it has no effect on numerics.

Usage:

    from srm_1d import run_simulation
    from srm_1d.propellant import make_hasegawa_propellant_1
    from srm_1d.grain_geometry import make_hasegawa_motor_A_geo
    from srm_1d.nozzle import Nozzle, compute_motor_performance
    from srm_1d.plotting import plot_pressure, plot_thrust, plot_flow_snapshot

    result = run_simulation(make_hasegawa_motor_A_geo(),
                            make_hasegawa_propellant_1(), roughness=20e-6)

    # Pressure trace with experimental overlay
    plot_pressure(result, title="Motor A",
                  experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL)

    # Thrust (requires nozzle post-processing)
    nozzle = Nozzle(D_throat=0.034, D_exit=0.050)
    perf = compute_motor_performance(result, nozzle, prop)
    plot_thrust(result, perf, title="Motor A")

    # Flow field snapshot at a specific time
    plot_flow_snapshot(result, t_target=2.0)
"""

import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ================================================================
# Experimental Data
# ================================================================
# Paste experimental data here as dict with 'time' and 'pressure'
# arrays. Pressure in MPa, time in seconds.
#
# To add your own data: either define a new dict following the same
# format, or use load_experimental_csv() to read from a file.

HASEGAWA_MOTOR_A_EXPERIMENTAL = {
    'label': 'Experimental (Hasegawa)',
    'time': np.array([
        0.000, 0.045, 0.053, 0.064, 0.068, 0.092, 0.116, 0.147, 0.206,
        0.358, 0.541, 0.953, 1.074, 1.214, 1.545, 1.822, 2.133,
        2.258, 2.348, 2.422, 2.555, 2.688, 2.754, 2.891, 3.036,
        3.168, 3.328, 3.481, 3.633, 3.852, 4.113, 4.297, 4.430,
        4.609, 4.675, 4.738,
    ]),
    'pressure': np.array([
        0.143, 2.299, 5.634, 6.219, 6.436, 5.922, 5.489, 5.297, 5.088,
        4.983, 5.014, 5.284, 5.388, 5.411, 5.578, 5.705, 5.887,
        5.951, 5.990, 5.854, 5.348, 4.762, 4.441, 3.935, 3.550,
        3.253, 2.923, 2.602, 2.256, 1.894, 1.219, 0.778, 0.464,
        0.191, 0.110, 0.110,
    
    ]),
    'time_offset': 0.07,  # Align ignition events
}


def load_experimental_csv(filepath, time_col=0, pressure_col=1,
                          delimiter=',', skip_header=1,
                          pressure_unit='MPa', label=None):
    """
    Load experimental data from a CSV file.

    Parameters
    ----------
    filepath : str
        Path to CSV file.
    time_col : int
        Column index for time data (0-based). Default: 0.
    pressure_col : int
        Column index for pressure data (0-based). Default: 1.
    delimiter : str
        Column delimiter. Default: ','.
    skip_header : int
        Number of header rows to skip. Default: 1.
    pressure_unit : str
        Unit of pressure in the file. One of 'Pa', 'kPa', 'MPa', 'psi'.
        Data will be converted to MPa for plotting.
    label : str or None
        Legend label. If None, uses the filename.

    Returns
    -------
    dict with 'time', 'pressure' (in MPa), and 'label'.
    """
    data = np.loadtxt(filepath, delimiter=delimiter, skiprows=skip_header)
    time = data[:, time_col]
    pressure = data[:, pressure_col]

    # Convert to MPa
    conversion = {
        'Pa': 1e-6,
        'kPa': 1e-3,
        'MPa': 1.0,
        'psi': 0.00689476,
    }
    if pressure_unit not in conversion:
        raise ValueError(
            f"Unknown pressure unit '{pressure_unit}'. "
            f"Use one of: {list(conversion.keys())}"
        )
    pressure = pressure * conversion[pressure_unit]

    if label is None:
        label = filepath.split('/')[-1].split('\\')[-1]

    return {
        'time': time,
        'pressure': pressure,
        'label': label,
    }


# ================================================================
# Pressure Plot
# ================================================================

def plot_pressure(result, title="Head-End Pressure",
                  experimental=None, save_path=None, ax=None,
                  n_head_cells=1):
    """
    Plot head-end pressure vs time.

    Parameters
    ----------
    result : dict
        Output from run_simulation. Must have 'time' and 'P_head'.
    title : str
        Plot title.
    experimental : dict or list of dicts, optional
        Experimental data to overlay. Each dict must have 'time'
        (seconds), 'pressure' (MPa), and optionally 'label' and
        'time_offset' (seconds, added to experimental time axis
        to align ignition events).
    save_path : str or None
        If provided, save figure to this path.
    ax : matplotlib Axes or None
        If provided, plot on this axes. Otherwise create a new figure.
    n_head_cells : int
        Number of head-end cells to plot. 1 = just P[0] (default).
        Values > 1 show pressure at cells 0..n-1 to reveal axial
        variation during the ignition transient.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes (None, None if no matplotlib).
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.get_figure()

    t = result['time']

    if n_head_cells <= 1:
        # Single head-end trace
        P = result['P_head'] / 1e6
        ax.plot(t, P, 'b-', linewidth=2, label='1D PISO (this work)')
    else:
        # Multi-cell traces from snapshots
        snapshots = result.get('snapshots', [])
        if snapshots and 'P' in snapshots[0]:
            snap_t = np.array([s['t'] for s in snapshots])
            cell_colors = plt.cm.viridis(np.linspace(0, 0.7, n_head_cells))
            for ci in range(n_head_cells):
                P_cell = np.array([s['P'][ci] / 1e6 for s in snapshots])
                x_mm = snapshots[0]['x'][ci] * 1000
                ax.plot(snap_t, P_cell, '-', color=cell_colors[ci],
                        linewidth=1.5, label=f'Cell {ci} (x={x_mm:.0f}mm)')
        # Always show P_head as the primary trace
        P = result['P_head'] / 1e6
        ax.plot(t, P, 'b-', linewidth=2, alpha=0.4, label='P[0] (full res)')

    # Experimental overlay(s)
    if experimental is not None:
        if isinstance(experimental, dict):
            experimental = [experimental]
        colors = ['k', 'r', 'g', 'orange', 'purple']
        markers = ['o', 's', '^', 'D', 'v']
        for i, expt in enumerate(experimental):
            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]
            label = expt.get('label', f'Experimental {i+1}')
            t_offset = expt.get('time_offset', 0.0)
            ax.plot(expt['time'] + t_offset, expt['pressure'],
                    color=color, linewidth=2, marker=marker,
                    markersize=3, markevery=max(1, len(expt['time'])//20),
                    label=label)

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Head-End Pressure [MPa]', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)

    if own_fig:
        plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved {save_path}")

    return fig, ax


# ================================================================
# Thrust Plot
# ================================================================

def plot_thrust(result, performance, title="Thrust",
                save_path=None, ax=None):
    """
    Plot thrust and Isp vs time.

    Parameters
    ----------
    result : dict
        Output from run_simulation.
    performance : dict
        Output from compute_motor_performance. Must have 'thrust' and 'Isp'.
    title : str
        Plot title.
    save_path : str or None
        If provided, save figure to this path.
    ax : matplotlib Axes or None
        If provided, plot on this axes (thrust only). Otherwise create
        a two-panel figure with thrust and Isp.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes.
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    t = result['time']
    thrust = performance['thrust']
    Isp = performance['Isp']

    own_fig = ax is None
    if own_fig:
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        ax_thrust = axes[0]
        ax_isp = axes[1]
    else:
        fig = ax.get_figure()
        ax_thrust = ax
        ax_isp = None
        axes = [ax]

    # Thrust
    ax_thrust.plot(t, thrust, 'b-', linewidth=2)
    ax_thrust.set_ylabel('Thrust [N]', fontsize=12)
    ax_thrust.set_title(title, fontsize=14)
    ax_thrust.grid(True, alpha=0.3)
    ax_thrust.set_ylim(bottom=0)

    # Add designation annotation
    desig = performance.get('motor_designation', '')
    impulse = performance.get('total_impulse', 0)
    if desig:
        ax_thrust.annotate(
            f"{desig}\n{impulse:.0f} N·s",
            xy=(0.98, 0.95), xycoords='axes fraction',
            ha='right', va='top', fontsize=12,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8),
        )

    # Isp
    if ax_isp is not None:
        ax_isp.plot(t, Isp, 'r-', linewidth=2)
        ax_isp.set_xlabel('Time [s]', fontsize=12)
        ax_isp.set_ylabel('Isp [s]', fontsize=12)
        ax_isp.grid(True, alpha=0.3)
        ax_isp.set_ylim(bottom=0)

    if own_fig:
        plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved {save_path}")

    return fig, axes


# ================================================================
# Flow Snapshot Plot
# ================================================================

def plot_flow_snapshot(result, t_target=None, snap_index=None,
                      title=None, save_path=None):
    """
    Plot a flow field snapshot (P, Mach, burn rate, port diameter).

    Parameters
    ----------
    result : dict
        Output from run_simulation. Must have 'snapshots'.
    t_target : float or None
        Target simulation time [s]. Uses the nearest snapshot.
    snap_index : int or None
        Snapshot index directly. Overrides t_target.
    title : str or None
        Plot title. Auto-generated if None.
    save_path : str or None
        If provided, save figure to this path.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes (2x2 grid).
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    snapshots = result.get('snapshots', [])
    if not snapshots:
        print("No snapshots in result dict.")
        return None, None

    # Select snapshot
    if snap_index is not None:
        snap_index = max(0, min(snap_index, len(snapshots) - 1))
    elif t_target is not None:
        snap_times = [s['t'] for s in snapshots]
        snap_index = min(range(len(snap_times)),
                         key=lambda i: abs(snap_times[i] - t_target))
    else:
        # Default: ~30% through the burn
        snap_index = len(snapshots) // 3

    snap = snapshots[snap_index]
    x_mm = snap['x'] * 1000

    if title is None:
        title = f"Flow Snapshot at t = {snap['t']:.3f} s"

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Pressure
    axes[0, 0].plot(x_mm, snap['P'] / 1e6, 'b-', linewidth=1.5)
    axes[0, 0].set_ylabel('Pressure [MPa]')
    axes[0, 0].set_title('Pressure')
    axes[0, 0].grid(True, alpha=0.3)

    # Mach
    axes[0, 1].plot(x_mm, snap['Mach'], 'r-', linewidth=1.5)
    axes[0, 1].set_ylabel('Mach Number')
    axes[0, 1].set_title('Mach Number')
    axes[0, 1].grid(True, alpha=0.3)

    # Burn rate
    axes[1, 0].plot(x_mm, snap['r_total'] * 1000, 'b-', linewidth=1.5,
                    label='Total')
    axes[1, 0].plot(x_mm, snap['r_erosive'] * 1000, 'r--', linewidth=1.5,
                    label='Erosive')
    r_normal = (snap['r_total'] - snap['r_erosive']) * 1000
    axes[1, 0].plot(x_mm, r_normal, 'g:', linewidth=1.5, label='Normal')
    axes[1, 0].set_xlabel('Position [mm]')
    axes[1, 0].set_ylabel('Burn Rate [mm/s]')
    axes[1, 0].set_title('Burn Rate')
    axes[1, 0].legend(fontsize=9)
    axes[1, 0].grid(True, alpha=0.3)

    # Show end-face mass injection locations if available
    if 'endface_msource' in snap:
        ef = snap['endface_msource']
        ef_mask = ef > 0
        if np.any(ef_mask):
            ax2 = axes[1, 0].twinx()
            ax2.bar(x_mm[ef_mask], ef[ef_mask], width=x_mm[1]-x_mm[0],
                    alpha=0.25, color='orange', label='End-face source')
            ax2.set_ylabel('End-face [kg/(m·s)]', fontsize=8, color='orange')
            ax2.tick_params(axis='y', labelcolor='orange', labelsize=7)

    # Port diameter
    axes[1, 1].plot(x_mm, snap['D_port'] * 1000, 'k-', linewidth=1.5)
    axes[1, 1].set_xlabel('Position [mm]')
    axes[1, 1].set_ylabel('Port Diameter [mm]')
    axes[1, 1].set_title('Port Diameter')
    axes[1, 1].grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved {save_path}")

    return fig, axes


# ================================================================
# Combined Summary Plot
# ================================================================

def plot_summary(result, performance=None, experimental=None,
                 title="Simulation Summary", save_path=None):
    """
    Combined plot: pressure, thrust (if available), and a flow snapshot.

    Parameters
    ----------
    result : dict
        Output from run_simulation.
    performance : dict or None
        Output from compute_motor_performance. If None, thrust panel
        is replaced with exit pressure.
    experimental : dict or list of dicts, optional
        Experimental pressure data for overlay.
    title : str
        Overall figure title.
    save_path : str or None
        If provided, save figure to this path.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes.
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    t = result['time']

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (0,0) Pressure trace
    plot_pressure(result, title="Head-End Pressure",
                  experimental=experimental, ax=axes[0, 0])

    # (0,1) Thrust or exit pressure
    if performance is not None:
        axes[0, 1].plot(t, performance['thrust'], 'b-', linewidth=2)
        axes[0, 1].set_ylabel('Thrust [N]', fontsize=11)
        axes[0, 1].set_title('Thrust', fontsize=12)
        desig = performance.get('motor_designation', '')
        if desig:
            axes[0, 1].annotate(
                desig, xy=(0.98, 0.95), xycoords='axes fraction',
                ha='right', va='top', fontsize=12,
                bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='wheat', alpha=0.8),
            )
    else:
        axes[0, 1].plot(t, result['P_exit'] / 1e6, 'r-', linewidth=2)
        axes[0, 1].set_ylabel('Exit Pressure [MPa]', fontsize=11)
        axes[0, 1].set_title('Nozzle-End Pressure', fontsize=12)
    axes[0, 1].set_xlabel('Time [s]', fontsize=11)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_ylim(bottom=0)
    axes[0, 1].set_xlim(left=0)

    # (1,0) and (1,1) Flow snapshot at ~30% burn
    snapshots = result.get('snapshots', [])
    if snapshots:
        snap = snapshots[len(snapshots) // 3]
        x_mm = snap['x'] * 1000

        axes[1, 0].plot(x_mm, snap['r_total'] * 1000, 'b-', lw=1.5,
                        label='Total')
        axes[1, 0].plot(x_mm, snap['r_erosive'] * 1000, 'r--', lw=1.5,
                        label='Erosive')
        axes[1, 0].set_xlabel('Position [mm]', fontsize=11)
        axes[1, 0].set_ylabel('Burn Rate [mm/s]', fontsize=11)
        axes[1, 0].set_title(f'Burn Rate at t={snap["t"]:.2f}s', fontsize=12)
        axes[1, 0].legend(fontsize=9)
        axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].plot(x_mm, snap['Mach'], 'r-', lw=1.5)
        axes[1, 1].set_xlabel('Position [mm]', fontsize=11)
        axes[1, 1].set_ylabel('Mach Number', fontsize=11)
        axes[1, 1].set_title(f'Mach at t={snap["t"]:.2f}s', fontsize=12)
        axes[1, 1].grid(True, alpha=0.3)
    else:
        for a in [axes[1, 0], axes[1, 1]]:
            a.text(0.5, 0.5, 'No snapshots', ha='center', va='center',
                   transform=a.transAxes, fontsize=14, color='gray')

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved {save_path}")

    return fig, axes


# ================================================================
# Comparison Plot (1D PISO vs openMotor or experimental)
# ================================================================

def plot_comparison(result, perf=None, reference=None,
                    ref_label="openMotor", title="1D PISO vs Reference",
                    save_path=None):
    """
    Plot pressure and thrust comparison between our sim and a reference.

    Parameters
    ----------
    result : dict
        Output from run_simulation.
    perf : dict or None
        Output from compute_motor_performance.
    reference : dict or None
        Reference data with at least 'time' and 'pressure' arrays.
        Can be from load_openmotor_csv() or load_experimental_csv().
        If it has 'force', thrust is also compared.
        Pressure should be in Pa (or if < 100, assumed MPa).
    ref_label : str
        Legend label for the reference data.
    title : str
        Plot title.
    save_path : str or None
        If provided, save figure.

    Returns
    -------
    fig, axes
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available.")
        return None, None

    has_thrust = perf is not None and reference is not None and 'force' in reference
    n_rows = 2 if has_thrust else 1
    fig, axes_raw = plt.subplots(n_rows, 1, figsize=(10, 5 * n_rows), sharex=True)
    if n_rows == 1:
        axes = [axes_raw]
    else:
        axes = list(axes_raw)

    t_sim = result['time']
    P_sim = result['P_head'] / 1e6

    # Pressure panel
    axes[0].plot(t_sim, P_sim, 'b-', linewidth=2, label='1D PISO')

    if reference is not None:
        t_ref = reference['time']
        P_ref = reference['pressure']
        if np.max(P_ref) > 100:  # Assume Pa, convert to MPa
            P_ref = P_ref / 1e6
        axes[0].plot(t_ref, P_ref, 'r--', linewidth=2, label=ref_label)

    axes[0].set_ylabel('Pressure [MPa]', fontsize=12)
    axes[0].set_title(title, fontsize=14)
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(bottom=0)

    # Thrust panel
    if has_thrust:
        axes[1].plot(t_sim, perf['thrust'], 'b-', linewidth=2, label='1D PISO')
        t_ref = reference['time']
        F_ref = reference['force']
        axes[1].plot(t_ref, F_ref, 'r--', linewidth=2, label=ref_label)
        axes[1].set_ylabel('Thrust [N]', fontsize=12)
        axes[1].legend(fontsize=11)
        axes[1].grid(True, alpha=0.3)
        axes[1].set_ylim(bottom=0)

    axes[-1].set_xlabel('Time [s]', fontsize=12)
    axes[-1].set_xlim(left=0)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved {save_path}")

    return fig, axes


# ================================================================
# Per-Grain Regression Plot
# ================================================================

def plot_grain_regression(grain_metrics, geo, title="Per-Grain Regression",
                          save_path=None):
    """
    Plot per-grain regression depth and web remaining vs time.

    Parameters
    ----------
    grain_metrics : dict
        Output from compute_grain_metrics().
    geo : MotorGeometry
        The geometry (for segment labels).
    title : str
        Plot title.
    save_path : str or None
        If provided, save figure.

    Returns
    -------
    fig, axes
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available.")
        return None, None

    t = grain_metrics['snap_times']
    reg = grain_metrics['regression'] * 1000  # to mm
    web = grain_metrics['web'] * 1000
    n_grains = reg.shape[1]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    colors = plt.cm.tab10(np.linspace(0, 1, max(n_grains, 1)))

    for k in range(n_grains):
        axes[0].plot(t, reg[:, k], color=colors[k], linewidth=1.5,
                     label=f'Grain {k}')
        axes[1].plot(t, web[:, k], color=colors[k], linewidth=1.5,
                     label=f'Grain {k}')

    axes[0].set_xlabel('Time [s]', fontsize=11)
    axes[0].set_ylabel('Regression [mm]', fontsize=11)
    axes[0].set_title('Radial Regression', fontsize=12)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel('Time [s]', fontsize=11)
    axes[1].set_ylabel('Web Remaining [mm]', fontsize=11)
    axes[1].set_title('Web Remaining', fontsize=12)
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_ylim(bottom=0)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved {save_path}")

    return fig, axes
