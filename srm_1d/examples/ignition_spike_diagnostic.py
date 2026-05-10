"""
Run ignition-spike diagnostics for a .ric motor.

This example is diagnostic-only. It compares the current default startup
path with opt-in isolation runs and writes plots/CSVs under
``artifacts/ignition_diagnostics/<case>/``.

Usage:
    python -m srm_1d.examples.ignition_spike_diagnostic --case hasegawa_a
    python -m srm_1d.examples.ignition_spike_diagnostic --case BALLSstick --t-max 3.0
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.tools.ignition_diagnostics import (
    analyze_ignition_spike,
    classification_report,
    plot_diagnostic_figures,
    write_diagnostic_outputs,
)


ROOT = Path(__file__).resolve().parents[2]
MOTORS_DIR = Path(__file__).resolve().parents[1] / "motors"
ARTIFACT_ROOT = Path("artifacts") / "ignition_diagnostics"


VARIANTS = {
    "baseline": {},
    "ambient_initial_gas": {"initial_gas_temperature": 293.0},
    "no_erosive": {"diagnostic_disable_erosive": True},
    "no_endfaces": {"diagnostic_disable_endfaces": True},
}


def _motor_path(case: str, explicit: str | None) -> Path:
    if explicit:
        return Path(explicit)
    return MOTORS_DIR / f"{case}.ric"


def _run_variant(args, motor_path: Path, variant: str, overrides: dict):
    output_dir = ARTIFACT_ROOT / args.case / variant
    output_dir.mkdir(parents=True, exist_ok=True)

    sim_kwargs = {
        "roughness": args.roughness,
        "kappa": args.kappa,
        "pyrogen": args.pyrogen,
        "pyrogen_mass": args.pyrogen_mass,
        "pyrogen_throat_area": args.pyrogen_throat_area,
        "T_ignition": args.T_ignition,
        "P_cutoff": args.P_cutoff,
        "t_max": args.t_max,
        "dt_max": args.dt_max,
        "snapshot_interval": args.snapshot_interval,
        "print_interval": args.print_interval,
        "verbose": args.verbose,
        **overrides,
    }

    result, perf, _nozzle, geo, prop = run_from_ric(str(motor_path), **sim_kwargs)
    diagnostics = analyze_ignition_spike(result, geo=geo, propellant=prop)
    write_diagnostic_outputs(diagnostics, output_dir, f"{args.case}_{variant}")
    plot_diagnostic_figures(result, diagnostics, output_dir, f"{args.case}_{variant}")
    plt.close("all")
    return diagnostics, result, perf


def main():
    parser = argparse.ArgumentParser(description="Run ignition spike diagnostics.")
    parser.add_argument("--case", default="hasegawa_a",
                        help="Case name. Defaults to srm_1d/motors/<case>.ric")
    parser.add_argument("--motor-path", default=None,
                        help="Explicit .ric path. Overrides --case lookup.")
    parser.add_argument("--variants", nargs="+",
                        default=["baseline", "ambient_initial_gas", "no_erosive", "no_endfaces"],
                        choices=sorted(VARIANTS),
                        help="Diagnostic variants to run.")
    parser.add_argument("--pyrogen", default="bpnv")
    parser.add_argument("--pyrogen-mass", type=float, default=None)
    parser.add_argument("--pyrogen-throat-area", type=float, default=None)
    parser.add_argument("--roughness", type=float, default=30e-6)
    parser.add_argument("--kappa", type=float, default=0.45)
    parser.add_argument("--T-ignition", dest="T_ignition", type=float, default=850.0)
    parser.add_argument("--t-max", type=float, default=3.0)
    parser.add_argument("--dt-max", type=float, default=1.0e-4)
    parser.add_argument("--P-cutoff", dest="P_cutoff", type=float, default=0.01e6)
    parser.add_argument("--snapshot-interval", type=float, default=0.005)
    parser.add_argument("--print-interval", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    motor_path = _motor_path(args.case, args.motor_path)
    if not motor_path.exists():
        raise FileNotFoundError(f"Motor file not found: {motor_path}")

    case_dir = ARTIFACT_ROOT / args.case
    case_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for variant in args.variants:
        print(f"Running {args.case}/{variant}...")
        diagnostics, _result, _perf = _run_variant(
            args, motor_path, variant, VARIANTS[variant]
        )
        print(classification_report(diagnostics))
        cls = diagnostics["classification"]
        pressure = diagnostics["pressure"]
        spread = diagnostics["ignition_spread"]
        rows.append({
            "variant": variant,
            "primary_driver": cls["primary_driver"],
            "peak_time_s": pressure["peak_time_s"],
            "peak_pressure_mpa": pressure["peak_pressure_mpa"],
            "spread_10_90_s": spread["spread_10_90_s"],
            "instant_ignition_collapse": spread["instant_ignition_collapse"],
        })

    with open(case_dir / f"{args.case}_variant_summary.csv", "w", newline="") as f:
        fieldnames = [
            "variant", "primary_driver", "peak_time_s", "peak_pressure_mpa",
            "spread_10_90_s", "instant_ignition_collapse",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote diagnostics under {case_dir}")


if __name__ == "__main__":
    main()
