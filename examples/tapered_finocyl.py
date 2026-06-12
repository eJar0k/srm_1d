"""
tapered_finocyl.py — Parametric axial-taper FMM grain demo
===========================================================

Demonstrates srm_1d's parametric FMM tapering (viz/geometry roadmap #3):
a single finocyl grain whose fins GROW along the axis, built from a start
and end cross-section without hand-authoring stepped segments.

The taper is authored with ``linear_taper`` (forward → aft cross-section)
and resolved into a stack of REAL per-station FMM regression tables by
``build_snapped_geometry`` — one genuine cross-section per axial station,
with the station count automatically matched to the snapped mesh.

The forward grain ALSO carries an **OD / end taper** (round-4 transient
support): a forward elliptical dome that shrinks the CASING diameter over
the head end (a domed closure). The per-cell ``cell_D_outer`` the solver
builds from it tapers the outer wall — visible in the GUI longitudinal
slice viewer's casing, and reflected here in the printed casing range.

This example reuses the Zerox motor's propellant / nozzle / igniter /
transport (Zerox is itself a 4-fin finocyl + BATES) and swaps in a
tapered finocyl for the forward grain so the run is physically sensible.

Usage:
    "C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe" \\
        -m examples.tapered_finocyl

Outputs (under artifacts/tapered_finocyl/<timestamp>_<sha>[-dirty]/):
    pressure.png  — head-end pressure trace
    flow.png      — flow-field snapshot (port diameter tapers along x)
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from srm_1d.openmotor_adapter import (
    load_ric, load_transport, ric_to_sim_args,
    build_pyrogen_chamber, load_pyrogen,
)
from srm_1d.grain_geometry import build_snapped_geometry
from srm_1d.fmm_grain import linear_taper
from srm_1d.simulation import run_simulation
from srm_1d.plotting import plot_pressure, plot_flow_snapshot
from srm_1d.run_artifacts import artifact_dir


MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'zerox.ric'
TRANSPORT_PATH = MOTOR_PATH.with_name('zerox.transport.yaml')

# FMM cross-sectional resolution per station (radial; independent of the
# axial mesh). 401 is plenty for a demo and keeps the ~dozen station
# solves fast; production runs use openMotor's 1001 default.
MAP_DIM = 401


def build_tapered_geometry(d_outer):
    """A forward 4-fin finocyl whose fins grow 0.6in -> 1.2in tip length,
    followed by Zerox's aft BATES segment. Core diameter, fin width, and
    fin count are held constant (only finLength tapers). The forward grain
    also has a forward elliptical OD dome (a head-end closure) that shrinks
    the casing over the first 50 mm."""
    IN = 0.0254
    base = dict(
        diameter=d_outer,
        length=0.254,            # snapped below; nominal forward-grain length
        inhibitedEnds='Bottom',  # constant across the taper (non-interpolable)
        coreDiameter=0.0330,
        numFins=4,
        finWidth=0.0038,
        invertedFins=False,
    )
    props_fwd = {**base, 'finLength': 0.6 * IN}
    props_aft = {**base, 'finLength': 1.2 * IN}

    # OD / end taper: a forward elliptical dome shrinking the casing from the
    # full OD to 60% over the head end. Carried on BOTH the TaperSpec (so each
    # per-station FMM table is clipped to the local casing) and the segment
    # spec (so compile builds the per-cell cell_D_outer); the domed end is
    # auto-inhibited (bonded closure).
    od_fwd_dome = [{'end': 'fwd', 'length': 0.050,
                    'endDiameter': 0.60 * d_outer, 'profile': 'elliptical'}]

    taper = linear_taper('Finocyl', props_fwd, props_aft, map_dim=MAP_DIM,
                         od_ends=od_fwd_dome, grain_length=0.254)

    segments_spec = [
        # Forward grain: the tapered finocyl + forward OD dome. D_bore_fwd is a
        # placeholder for FMM cells (port comes from the regression table).
        {'length': 0.254, 'taper': taper, 'od_ends': od_fwd_dome,
         'inhibit_fwd': True},
        # Aft grain: Zerox's cylindrical BATES (touches the finocyl, so the
        # shared faces auto-inhibit).
        {'length': 0.657, 'D_bore_fwd': 0.03175},
    ]
    return build_snapped_geometry(segments_spec, D_outer=d_outer,
                                  target_propellant_cells=80)


def main():
    # Borrow Zerox's propellant / nozzle / transport / igniter sizing.
    motor = load_ric(str(MOTOR_PATH))
    gas_props = load_transport(str(TRANSPORT_PATH))
    sim_args = ric_to_sim_args(
        motor, gas_props=gas_props,
        roughness=20e-6, kappa=0.45, T_ignition=850.0,
        cfl_target=0.5, dt_max=1e-4, t_max=8.0, P_cutoff=0.05e6,
        snapshot_interval=0.25, print_interval=0.5,
    )

    d_outer = sim_args['geo'].D_outer

    # Swap in the tapered geometry; rebuild the pyrogen chamber against it.
    print(f"Building tapered finocyl geometry (map_dim={MAP_DIM})...")
    geo = build_tapered_geometry(d_outer)
    tabs = geo.segments[0].fmm_tables
    n_stations = len(tabs)
    print(f"  forward grain resolved to {n_stations} FMM stations; "
          f"initial burn perimeter "
          f"{tabs[0].initial_perimeter*1e3:.1f} -> "
          f"{tabs[-1].initial_perimeter*1e3:.1f} mm, "
          f"initial port area "
          f"{tabs[0].initial_port_area*1e6:.0f} -> "
          f"{tabs[-1].initial_port_area*1e6:.0f} mm^2 (fwd -> aft, "
          f"fins grow 0.6 -> 1.2 in)")
    print(f"  total propellant volume: "
          f"{geo.total_propellant_volume()*1e6:.1f} cm^3")
    ga = geo.compile_geometry_arrays()
    cD = ga['cell_D_outer']
    print(f"  casing OD (cell_D_outer): {cD.min()*1e3:.1f} -> {cD.max()*1e3:.1f} "
          f"mm (forward dome shrinks the head-end casing)")

    sim_args.pop('geo')
    propellant = sim_args.pop('propellant')
    sim_args['pyrogen_chamber'] = build_pyrogen_chamber(
        load_pyrogen('bpnv'), geo, sim_args['nozzle'],
    )

    result = run_simulation(geo, propellant, **sim_args)

    out = artifact_dir('tapered_finocyl')
    plot_pressure(result, title="Tapered Finocyl — Head-End Pressure",
                  save_path=str(out / 'pressure.png'))
    plot_flow_snapshot(result, t_target=1.0,
                       title="Tapered Finocyl — Flow at t ≈ 1.0s "
                             "(port diameter tapers along x)",
                       save_path=str(out / 'flow.png'))
    plt.close('all')
    print(f"\nPlots saved to {out}")


if __name__ == '__main__':
    main()
