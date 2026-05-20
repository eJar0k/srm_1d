"""Dump the snapped grid layout for Hasegawa A at several cell counts.

Reports dx, cell-segment mapping, head/aft grain cells, and the
trailing-gas cell width relative to the throat diameter. Used by
docs/v0_7_0/audits/2026-05-20_radiation_collapse_localT.md to check
whether the snapped discretization contributes to the PISO/throat
numerical-front sensitivity.
"""

from pathlib import Path

import numpy as np

from srm_1d.openmotor_adapter import (
    convert_geometry,
    convert_nozzle,
    load_ric,
)


MOTOR_PATH = Path(__file__).resolve().parent.parent / "motors" / "hasegawa_a.ric"


def _summarize(cells: int) -> dict:
    ric = load_ric(str(MOTOR_PATH))
    geo = convert_geometry(ric["grains"], target_propellant_cells=cells)
    nozzle = convert_nozzle(ric["nozzle"])

    dx = geo.L_motor / geo.N_cells
    x_centers = (np.arange(geo.N_cells) + 0.5) * dx

    n_seg = len(geo.segments)
    seg_x_start = np.array([s.x_start for s in geo.segments])
    seg_length = np.array([s.length for s in geo.segments])

    cell_seg = np.full(geo.N_cells, -1, dtype=int)
    for i, x in enumerate(x_centers):
        x_lo = x - 0.5 * dx
        x_hi = x + 0.5 * dx
        best = 0.0
        for k in range(n_seg):
            ovl = max(0.0, min(x_hi, seg_x_start[k] + seg_length[k]) - max(x_lo, seg_x_start[k]))
            if ovl > best:
                best = ovl
                cell_seg[i] = k
        if best <= 1e-12 * dx:
            cell_seg[i] = -1

    grain_cells = np.flatnonzero(cell_seg >= 0)
    head_grain = int(grain_cells[0]) if grain_cells.size else -1
    aft_grain = int(grain_cells[-1]) if grain_cells.size else -1

    leading_gas = int(grain_cells[0]) if grain_cells.size else 0
    trailing_gas = int(geo.N_cells - 1 - grain_cells[-1]) if grain_cells.size else 0

    return {
        "cells_target": cells,
        "N_cells": geo.N_cells,
        "L_motor_m": geo.L_motor,
        "dx_mm": dx * 1000.0,
        "n_segments": n_seg,
        "first_seg_xstart_mm": seg_x_start[0] * 1000.0,
        "first_seg_length_mm": seg_length[0] * 1000.0,
        "head_grain_cell": head_grain,
        "aft_grain_cell": aft_grain,
        "leading_gas_cells": leading_gas,
        "trailing_gas_cells": trailing_gas,
        "trailing_gas_width_mm": trailing_gas * dx * 1000.0,
        "D_throat_mm": nozzle.D_throat * 1000.0,
        "trailing_gas_width_over_throat": (trailing_gas * dx) / nozzle.D_throat if nozzle.D_throat > 0 else None,
    }


def main():
    print("Hasegawa A snapped grid layout vs target cell count:")
    print()
    header = (
        "target", "N_cells", "L_motor_m", "dx_mm",
        "lead_gas", "head_grain", "aft_grain", "trail_gas",
        "trail_gas_mm", "trail/throat",
    )
    print("  ".join(f"{h:>12}" for h in header))
    for cells in (50, 100, 200):
        s = _summarize(cells)
        row = (
            f"{s['cells_target']}", f"{s['N_cells']}", f"{s['L_motor_m']:.4f}",
            f"{s['dx_mm']:.3f}",
            f"{s['leading_gas_cells']}", f"{s['head_grain_cell']}",
            f"{s['aft_grain_cell']}", f"{s['trailing_gas_cells']}",
            f"{s['trailing_gas_width_mm']:.2f}",
            f"{s['trailing_gas_width_over_throat']:.3f}",
        )
        print("  ".join(f"{v:>12}" for v in row))


if __name__ == "__main__":
    main()
