"""
zerox_fmm_dump.py — Visualize the Zerox finocyl FMM table
==========================================================

Phase C of the Zerox tuning workflow. Loads the Zerox .ric, extracts
the forward Finocyl grain's FmmTable, and plots:

  1. perimeter(reg)        — should reveal the snap-to-zero step at wall_web
  2. port_area(reg)         — jump to casting_area at burnout
  3. dPerimeter/dReg        — quantifies sharpness of the burnout step

The aim is to judge whether the "sharp step" symptom in the simulated
pressure trace is consistent with this discontinuity (per the burnout
clamp at fmm_grain.py:386-394 — when reg >= wall_web, perimeter is
forced to 0 and port_area to casting_area).

Usage:
    "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" \\
        -m srm_1d.examples.zerox_fmm_dump
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.fmm_grain import from_ric_grain
from srm_1d.openmotor_adapter import load_ric


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox.ric'


def main():
    motor = load_ric(str(MOTOR_PATH))
    grains = motor['grains']

    # Find the first FMM-type grain (Finocyl in Zerox); BATES is analytic.
    fmm_grain = next((g for g in grains if g['type'] != 'BATES'), None)
    if fmm_grain is None:
        raise SystemExit("No FMM-type grain found in zerox.ric")

    print(f"Building FmmTable for grain type {fmm_grain['type']!r}...")
    table = from_ric_grain(fmm_grain, map_dim=750)  # match Zerox config mapDim

    reg_mm = table.reg_depth * 1000
    perim_mm = table.perimeter * 1000
    port_cm2 = table.port_area * 1e4
    wall_web_mm = table.wall_web * 1000

    casting_area_cm2 = (np.pi / 4 * table.grain_outer_diameter ** 2) * 1e4

    # Numerical derivative of perimeter wrt reg, for sharpness panel.
    dPerim_dReg = np.gradient(perim_mm, reg_mm)

    # Quantify the snap: last finite vs final (zero) sample.
    perim_jump_mm = perim_mm[-2] if len(perim_mm) >= 2 else 0.0
    port_jump_cm2 = casting_area_cm2 - port_cm2[-2] if len(port_cm2) >= 2 else 0.0

    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)

    axes[0].plot(reg_mm, perim_mm, 'b-', linewidth=2)
    axes[0].axvline(wall_web_mm, color='red', linestyle='--', linewidth=1.2,
                    label=f'wall_web = {wall_web_mm:.3f} mm')
    axes[0].annotate(
        f'snap-to-zero  Δ ≈ {perim_jump_mm:.1f} mm',
        xy=(wall_web_mm, perim_jump_mm),
        xytext=(wall_web_mm * 0.6, perim_jump_mm * 0.8),
        arrowprops=dict(arrowstyle='->', color='black'),
        fontsize=10,
    )
    axes[0].set_ylabel('Perimeter [mm]')
    axes[0].set_title(f'Zerox forward {table.geom_name} FMM table'
                      f'  (n_samples={len(table.reg_depth)})')
    axes[0].legend(loc='upper right')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(reg_mm, port_cm2, 'g-', linewidth=2)
    axes[1].axhline(casting_area_cm2, color='gray', linestyle=':',
                    linewidth=1.0, label=f'casting_area = {casting_area_cm2:.2f} cm²')
    axes[1].axvline(wall_web_mm, color='red', linestyle='--', linewidth=1.2)
    axes[1].annotate(
        f'jump to casting_area  Δ ≈ {port_jump_cm2:.2f} cm²',
        xy=(wall_web_mm, casting_area_cm2),
        xytext=(wall_web_mm * 0.55, casting_area_cm2 * 0.7),
        arrowprops=dict(arrowstyle='->', color='black'),
        fontsize=10,
    )
    axes[1].set_ylabel('Port area [cm²]')
    axes[1].legend(loc='lower right')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(reg_mm, dPerim_dReg, 'm-', linewidth=2)
    axes[2].axvline(wall_web_mm, color='red', linestyle='--', linewidth=1.2)
    axes[2].axhline(0.0, color='k', linestyle=':', linewidth=0.7)
    axes[2].set_ylabel('dPerimeter / dReg [mm / mm]')
    axes[2].set_xlabel('Regression depth [mm]')
    axes[2].set_title('Sharpness of the burnout step '
                      '(numerical gradient — large negative spike at wall_web → discontinuity)')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    save_path = "zerox_fmm_dump.png"
    fig.savefig(save_path, dpi=150)
    print(f"Saved {save_path}")
    plt.close(fig)

    print(f"\nSummary:")
    print(f"  geom = {table.geom_name}")
    print(f"  grain OD = {table.grain_outer_diameter*1e3:.2f} mm")
    print(f"  grain length = {table.grain_length*1e3:.2f} mm")
    print(f"  wall_web = {wall_web_mm:.3f} mm")
    print(f"  initial perimeter = {perim_mm[0]:.1f} mm")
    print(f"  perimeter just before burnout = {perim_jump_mm:.1f} mm  "
          f"(then snaps to 0)")
    print(f"  initial port_area = {port_cm2[0]:.2f} cm²")
    print(f"  port_area just before burnout = {port_cm2[-2]:.2f} cm²  "
          f"(then jumps to {casting_area_cm2:.2f})")


if __name__ == '__main__':
    main()
