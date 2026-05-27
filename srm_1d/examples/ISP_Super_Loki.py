"""
ISP_Super_Loki.py -- run the 'ISP_Super_Loki.ric' through head_basket.

**Provenance note (v0.7.4 Phase C.1)**: this example models the
**RCS Rocket Motor Components surplus Super Loki recreation** with
the **8522 HTPB analog propellant** (developed by Industrial Solid
Propulsion to replace the original 1968-1993 polysulfide formula).
The "ISP" in the filename refers to Industrial Solid Propulsion,
NOT a single "ISP Corporation". See lit dive at
``srm_1d/docs/v0_7_3/references/super_loki_igniter_lit_dive.md``.

**Igniter assumption**: the RCS Super Loki does NOT ship with an
igniter; the user assembles their own. This example models an
**amateur head-end BKNO3 pellet pack glued to the forward bulkhead**
(=  head_basket topology). This is one of several plausible amateur
configurations — see PyrogenChamber docstring L82-100 for the full
discussion. The factory 1968-era Super Loki igniter was nozzle-
inserted; that configuration would map to the deferred
``aft_fore_firing`` topology, not head_basket.

**Validation status**: no verified Super Loki experimental data is
in the repo. The previous overlay (commented-out array in this
file's history) was mis-labeled Chunc data. NASA CR-61238 (1968)
contains a polysulfide-era pressure trace that could be digitized
as a partial reference (15-20% peak deviation expected from the
HTPB analog).

Usage:
    python -m srm_1d.examples.ISP_Super_Loki
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
from srm_1d.plotting import (
    plot_pressure, plot_flow_snapshot, plot_summary,
    plot_flow_snapshots, plot_field_heatmap,
)
from srm_1d.run_artifacts import artifact_dir


CASE_NAME = 'ISP_Super_Loki'
MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'ISP_Super_Loki.ric'


def _run_one(mode, output_dir, particle_d=5.0e-3, particle_LD=3.0):
    """Run Super Loki with one heat-delivery mode and save artifacts.

    v0.7.4 Phase C.1: particle_d and particle_LD are exposed as run-
    script knobs (default to MTV-pellet 5 mm × L/D=3). Overriding them
    here is the cleanest way to A/B sweep particle geometry without
    editing the YAML.
    """
    # Build the pyrogen Pyrogen object explicitly so we can override
    # heat_delivery_mode + particle geometry per A/B run.
    pyrogen_obj = load_pyrogen('mtv')
    pyrogen_obj.heat_delivery_mode = mode
    pyrogen_obj.particle_diameter_m = particle_d
    pyrogen_obj.particle_LD_ratio = particle_LD

    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        roughness=30e-6,
        kappa=0.45,
        pyrogen=pyrogen_obj,
        pyrogen_mass=None,
        T_ignition=850.0,
        P_cutoff=0.01e6,
        injection_topology='head_basket',
        cartridge_length_m=-1.0,
        snapshot_interval=0.01,
        print_interval=0.2,
        verbose=False,
    )

    summary = result['summary']
    label = f"ISP_Super_Loki [head_basket / {mode}]"
    print(
        f"{label}: P_peak={summary['P_peak'] / 1e6:.2f} MPa, "
        f"t_peak={summary['t_peak']:.3f} s, "
        f"t_burn={summary['t_burn']:.3f} s, "
        f"impulse={perf['total_impulse']:.1f} N*s, "
        f"designation={perf['motor_designation']}"
    )

    # v0.7.3 Phase B (2026-05-25): NO verified Super Loki experimental
    # dataset in the repo. The pressure trace prior to today was
    # comparing against a mis-labeled Chunc dataset that had been
    # copy-pasted into this example years ago. Plot without overlay
    # until verified Super Loki static-fire data is sourced.
    plot_pressure(
        result,
        title=f"ISP_Super_Loki — head_basket / {mode}",
        save_path=str(output_dir / f"pressure_{mode}.png"),
    )
    plot_flow_snapshots(
        result,
        t_targets=[0.01, 0.05, 0.20, 0.50, 1.00],
        fields=('P', 'u', 'T', 'is_burning'),
        title=f"ISP_Super_Loki — flow evolution ({mode})",
        save_path=str(output_dir / f"flow_multi_{mode}.png"),
    )
    plot_field_heatmap(
        result,
        fields=('P', 'u', 'T', 'T_surf', 'is_burning'),
        t_max=0.5,
        title=f"ISP_Super_Loki — x-t heatmap (ignition, {mode})",
        save_path=str(output_dir / f"heatmap_ignition_{mode}.png"),
    )
    plot_summary(
        result,
        performance=perf,
        title=f"ISP_Super_Loki Summary ({mode})",
        save_path=str(output_dir / f"summary_{mode}.png"),
    )
    plt.close('all')
    return summary


def main():
    OUTPUT_DIR = artifact_dir(CASE_NAME)
    print(f"v0.7.3 Phase B.6 — ISP Super Loki A/B (head_basket)")
    print(f"  NOTE: no verified Super Loki experimental data in repo")
    print(f"  (the previous overlay was mis-labeled Chunc data;")
    print(f"   removed 2026-05-25 when provenance was traced).")
    print()

    # A/B independent modes per docs/v0_7_3/PHASE_B_SCOPE.md §B.4
    # (don't stack — DeMar already includes ~83% of radiation
    # contribution; stacking would double-count). 'none' is the
    # control case (no pyrogen surface heat flux at all) — if this
    # also reaches ignition, the structural problem is solved by
    # B.0 + B.2 alone and B.4 is just trimming.
    summary_none      = _run_one('none', OUTPUT_DIR)
    summary_demar     = _run_one('demar', OUTPUT_DIR)
    summary_radiation = _run_one('radiation', OUTPUT_DIR)

    print()
    print(f"Plots saved under {OUTPUT_DIR}")
    print()
    print("A/B summary:")
    print(f"  none:      P_peak = {summary_none['P_peak']/1e6:7.3f} MPa "
          f"@ t={summary_none['t_peak']:.3f} s")
    print(f"  demar:     P_peak = {summary_demar['P_peak']/1e6:7.3f} MPa "
          f"@ t={summary_demar['t_peak']:.3f} s")
    print(f"  radiation: P_peak = {summary_radiation['P_peak']/1e6:7.3f} MPa "
          f"@ t={summary_radiation['t_peak']:.3f} s")
    print(f"  reference: no Super Loki experimental data in repo")


if __name__ == '__main__':
    main()
