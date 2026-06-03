"""
zerox_diagnostic.py — Localize the Zerox over-pressure
========================================================

Phases A + B of the Zerox tuning workflow. Runs the current borrowed
Hasegawa-A calibration on Zerox with dense snapshots and produces
diagnostic figures that classify the dominant error source:

  - Igniter-driven spike    → P_head peaks before r_erosive develops
  - Erosive-driven spike    → r_erosive(x,t) ramps up in finocyl region
                              (x < 25mm) during the 4-6 MPa window
  - FMM-boundary step       → D_port(x,t) shows discontinuous jump in
                              the finocyl cells; head-end P kinks
  - Throat erosion stuck    → D_throat(t) flat despite erosion_coeff > 0

Outputs:
    zerox_diagnostic_pressure.png  — head-end P + experimental overlay
    zerox_diagnostic_xt.png        — x-t pcolormesh of r_erosive, D_port
    zerox_diagnostic_throat.png    — D_throat(t) and A_throat ratio
    zerox_diagnostic_snapshots.png — flow snapshots at 6 times

Usage:
    "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" \\
        -m examples.zerox_diagnostic
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import plot_pressure, plot_flow_snapshot, ZEROX_EXPERIMENTAL


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox.ric'

# Forward Finocyl is 254mm of the ~911mm grain stack (~28% of the motor)
# — annotate the finocyl/BATES boundary on x-t plots.
FINOCYL_LENGTH_MM = 254.0


def run_baseline():
    """Run Zerox with current Hasegawa-A calibration + dense snapshots."""
    return run_from_ric(
        str(MOTOR_PATH),
        # Borrowed Hasegawa A LHS rank-1 — known-wrong baseline.
        roughness=20e-6,
        kappa=0.45,
        pyrogen='bpnv',
        pyrogen_mass=None,
        T_ignition=850.0,
        cfl_target=0.5,
        dt_max=1e-4,
        t_max=8.0,
        P_cutoff=0.01e6,
        snapshot_interval=0.05,   # dense — ~160 frames over 8s burn
        print_interval=1.0,
    )


def stack_snapshots(snapshots, key):
    """Stack snapshot[key] arrays into a (n_snaps, N) 2D array."""
    return np.array([s[key] for s in snapshots])


def plot_xt_diagnostic(result, save_path):
    """
    x-t pcolormesh of r_erosive (mm/s) and D_port (mm).
    Annotates the finocyl/BATES segment boundary.
    """
    snapshots = result['snapshots']
    if not snapshots:
        print("No snapshots — cannot build x-t diagnostic.")
        return

    t_snap = np.array([s['t'] for s in snapshots])
    x_mm = snapshots[0]['x'] * 1000

    r_erosive = stack_snapshots(snapshots, 'r_erosive') * 1000  # m/s -> mm/s
    D_port = stack_snapshots(snapshots, 'D_port') * 1000        # m -> mm

    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

    # r_erosive (mm/s) — should localize spike to finocyl region if erosive
    pcm0 = axes[0].pcolormesh(
        t_snap, x_mm, r_erosive.T,
        shading='nearest', cmap='inferno',
    )
    cb0 = fig.colorbar(pcm0, ax=axes[0])
    cb0.set_label('r_erosive [mm/s]')
    axes[0].axhline(FINOCYL_LENGTH_MM, color='cyan', linestyle='--',
                    linewidth=1.2, label=f'Finocyl/BATES @ {FINOCYL_LENGTH_MM:.1f}mm')
    axes[0].set_ylabel('Axial position x [mm]')
    axes[0].set_title('Erosive burn-rate increment vs (x, t)')
    axes[0].legend(loc='upper right', fontsize=9)

    # D_port (mm) — should show the burnout step in finocyl cells
    pcm1 = axes[1].pcolormesh(
        t_snap, x_mm, D_port.T,
        shading='nearest', cmap='viridis',
    )
    cb1 = fig.colorbar(pcm1, ax=axes[1])
    cb1.set_label('D_port [mm]')
    axes[1].axhline(FINOCYL_LENGTH_MM, color='red', linestyle='--',
                    linewidth=1.2, label=f'Finocyl/BATES @ {FINOCYL_LENGTH_MM:.1f}mm')
    axes[1].set_ylabel('Axial position x [mm]')
    axes[1].set_xlabel('Time [s]')
    axes[1].set_title('Effective port diameter vs (x, t)')
    axes[1].legend(loc='upper right', fontsize=9)

    plt.tight_layout()
    fig.savefig(save_path, dpi=150)
    print(f"Saved {save_path}")
    plt.close(fig)


def plot_throat_diagnostic(result, save_path):
    """D_throat(t) and A_throat(t)/A_throat(0) — verify nozzle erosion."""
    t = result['time']
    D_throat = result.get('D_throat')
    if D_throat is None or len(D_throat) == 0:
        print("No D_throat history — cannot build throat diagnostic.")
        return

    D_throat_mm = D_throat * 1000
    A_ratio = (D_throat / D_throat[0]) ** 2

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(t, D_throat_mm, 'b-', linewidth=2)
    axes[0].set_ylabel('D_throat [mm]')
    axes[0].set_title(
        f'Throat diameter vs t  (initial {D_throat_mm[0]:.2f} mm, '
        f'final {D_throat_mm[-1]:.2f} mm, '
        f'Δ {D_throat_mm[-1] - D_throat_mm[0]:+.2f} mm)'
    )
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, A_ratio, 'r-', linewidth=2)
    axes[1].axhline(1.0, color='k', linestyle=':', linewidth=0.8)
    axes[1].set_ylabel('A_throat(t) / A_throat(0)')
    axes[1].set_xlabel('Time [s]')
    axes[1].set_title(f'Throat-area ratio  (final {A_ratio[-1]:.3f})')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(save_path, dpi=150)
    print(f"Saved {save_path}")
    plt.close(fig)


def main():
    print(f"Running baseline Zerox sim ({MOTOR_PATH.name})...")
    result, perf, _, _, _ = run_baseline()

    print("\nGenerating diagnostic plots...")
    plot_pressure(
        result, title="Zerox — Baseline (Hasegawa-A calibration)",
        experimental=ZEROX_EXPERIMENTAL,
        time_offset=ZEROX_EXPERIMENTAL.get('time_offset', 0.0),
        save_path="zerox_diagnostic_pressure.png",
    )
    plot_xt_diagnostic(result, "zerox_diagnostic_xt.png")
    plot_throat_diagnostic(result, "zerox_diagnostic_throat.png")
    # Single mid-spike snapshot — x-t plot above already covers full axial-temporal evolution.
    plot_flow_snapshot(
        result, t_target=0.30,
        title='Zerox flow at t ≈ 0.30s (mid-spike)',
        save_path="zerox_diagnostic_flow.png",
    )

    plt.close('all')
    print("\nDone. Inspect:")
    print("  zerox_diagnostic_pressure.png — head-end P vs experimental")
    print("  zerox_diagnostic_xt.png       — r_erosive(x,t), D_port(x,t)")
    print("  zerox_diagnostic_throat.png   — D_throat(t), A_throat ratio")
    print("  zerox_diagnostic_flow.png     — flow snapshot at t ≈ 0.30s")


if __name__ == '__main__':
    main()
