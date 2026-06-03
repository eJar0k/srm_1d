"""
Run ignition-spike diagnostics for a .ric motor.

This example is diagnostic-only. It compares the current default startup
path with opt-in isolation runs and writes plots/CSVs under
``artifacts/ignition_diagnostics/<case>/``.

Usage:
    python -m examples.ignition_spike_diagnostic --case hasegawa_a
    python -m examples.ignition_spike_diagnostic --case BALLSstick --t-max 3.0
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from srm_1d.nozzle import compute_motor_performance
from srm_1d.openmotor_adapter import (
    build_pyrogen_chamber,
    load_pyrogen,
    load_ric,
    load_transport,
    ric_to_sim_args,
    run_from_ric,
)
from srm_1d.propellant import Pyrogen
from srm_1d.simulation import run_simulation
from srm_1d.tools.ignition_diagnostics import (
    analyze_ignition_spike,
    classification_report,
    literature_evaluation_report,
    plot_diagnostic_figures,
    write_diagnostic_outputs,
)


ROOT = Path(__file__).resolve().parents[2]
MOTORS_DIR = Path(__file__).resolve().parents[1] / "motors"
ARTIFACT_ROOT = Path("artifacts") / "ignition_diagnostics"


VARIANTS = {
    "baseline": {},
    "ambient_initial_gas": {"initial_gas_temperature": 293.0},
    "ambient_no_surface_heating": {
        "initial_gas_temperature": 293.0,
        "diagnostic_disable_pyrogen_surface_heating": True,
    },
    "ambient_no_radiation": {
        "initial_gas_temperature": 293.0,
        "diagnostic_disable_adjacent_radiation": True,
    },
    "no_erosive": {"diagnostic_disable_erosive": True},
    "no_endfaces": {"diagnostic_disable_endfaces": True},
    "no_momentum": {"diagnostic_disable_momentum": True},
}


DEFAULT_VARIANTS = [
    "baseline", "ambient_initial_gas",
    "ambient_no_surface_heating",
    "ambient_no_radiation", "no_erosive",
    "no_endfaces", "no_momentum",
]


RADIATION_PROBE_VARIANTS = {
    "ambient_nominal_T850_dt1e-4": {
        "initial_gas_temperature": 293.0,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "ambient_no_radiation_T850_dt1e-4": {
        "initial_gas_temperature": 293.0,
        "diagnostic_disable_adjacent_radiation": True,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "ambient_no_surface_heating_T850_dt1e-4": {
        "initial_gas_temperature": 293.0,
        "diagnostic_disable_pyrogen_surface_heating": True,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "baseline_hotfill_T850_dt1e-4": {
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "no_erosive_hotfill_T850_dt1e-4": {
        "diagnostic_disable_erosive": True,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "ambient_emissivity_0_T850_dt1e-4": {
        "initial_gas_temperature": 293.0,
        "propellant_radiation_emissivity": 0.0,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "ambient_emissivity_high_T850_dt1e-4": {
        "initial_gas_temperature": 293.0,
        "propellant_radiation_emissivity": 0.9,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "ambient_nominal_T650_dt1e-4": {
        "initial_gas_temperature": 293.0,
        "T_ignition": 650.0,
        "dt_max": 1.0e-4,
    },
    "ambient_nominal_T750_dt1e-4": {
        "initial_gas_temperature": 293.0,
        "T_ignition": 750.0,
        "dt_max": 1.0e-4,
    },
    "ambient_nominal_T850_dt2e-5": {
        "initial_gas_temperature": 293.0,
        "T_ignition": 850.0,
        "dt_max": 2.0e-5,
    },
}


RADIATION_PROBE_DEFAULT_VARIANTS = list(RADIATION_PROBE_VARIANTS)


def _ambient_emissivity_variant(emissivity: float, **overrides) -> dict:
    variant = {
        "initial_gas_temperature": 293.0,
        "propellant_radiation_emissivity": float(emissivity),
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    }
    variant.update(overrides)
    return variant


RADIATION_COLLAPSE_VARIANTS = {
    "ambient_emissivity_0": _ambient_emissivity_variant(0.0),
    "ambient_no_radiation": {
        "initial_gas_temperature": 293.0,
        "diagnostic_disable_adjacent_radiation": True,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "ambient_emissivity_0p05": _ambient_emissivity_variant(0.05),
    "ambient_emissivity_0p10": _ambient_emissivity_variant(0.10),
    "ambient_emissivity_0p20": _ambient_emissivity_variant(0.20),
    "ambient_emissivity_0p30": _ambient_emissivity_variant(0.30),
    "ambient_emissivity_0p40": _ambient_emissivity_variant(0.40),
    "ambient_emissivity_0p45": _ambient_emissivity_variant(0.45),
    "ambient_emissivity_0p50": _ambient_emissivity_variant(0.50),
    "ambient_emissivity_0p60": _ambient_emissivity_variant(0.60),
    "ambient_emissivity_0p75": _ambient_emissivity_variant(0.75),
    "ambient_emissivity_0p90": _ambient_emissivity_variant(0.90),
    "baseline_hotfill": {
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "no_erosive_hotfill": {
        "diagnostic_disable_erosive": True,
        "T_ignition": 850.0,
        "dt_max": 1.0e-4,
    },
    "ambient_rad045_no_erosive": _ambient_emissivity_variant(
        0.45, diagnostic_disable_erosive=True
    ),
    "ambient_rad045_receiver_heat_no_sink": _ambient_emissivity_variant(
        0.45, diagnostic_disable_radiation_gas_sink=True
    ),
    "ambient_rad045_no_surface_heating": _ambient_emissivity_variant(
        0.45, diagnostic_disable_pyrogen_surface_heating=True
    ),
    "ambient_rad045_cfl025": _ambient_emissivity_variant(0.45, cfl_target=0.25),
    "ambient_rad045_cfl010": _ambient_emissivity_variant(0.45, cfl_target=0.10),
    "ambient_rad090_cfl025": _ambient_emissivity_variant(0.90, cfl_target=0.25),
    "ambient_rad090_cfl010": _ambient_emissivity_variant(0.90, cfl_target=0.10),
    "ambient_rad045_dt2e-5": _ambient_emissivity_variant(0.45, dt_max=2.0e-5),
    "ambient_rad090_dt2e-5": _ambient_emissivity_variant(0.90, dt_max=2.0e-5),
    "ambient_rad045_cells50": _ambient_emissivity_variant(0.45, target_propellant_cells=50),
    "ambient_rad045_cells200": _ambient_emissivity_variant(0.45, target_propellant_cells=200),
    "ambient_rad090_cells50": _ambient_emissivity_variant(0.90, target_propellant_cells=50),
    "ambient_rad090_cells200": _ambient_emissivity_variant(0.90, target_propellant_cells=200),
}


RADIATION_COLLAPSE_DEFAULT_VARIANTS = list(RADIATION_COLLAPSE_VARIANTS)


# ---------------------------------------------------------------------------
# Ignition-tuning Cartesian sweep (Step 4 of
# continue-with-the-numerical-zippy-dawn.md). Sweep T_ignition and k_solid
# at two radiation_emissivity values to determine whether existing
# ignition-timing knobs can match the Hasegawa A experimental trace
# without introducing a burn-establishment ramp model.
# ---------------------------------------------------------------------------

def _ignition_tuning_variant(T_ignition: float, k_solid: float,
                             emissivity: float, **overrides) -> dict:
    """Cartesian ignition-tuning variant. T_ignition [K], k_solid [W/m/K]."""
    variant = {
        "T_ignition": float(T_ignition),
        "propellant_k_solid": float(k_solid),
        "propellant_radiation_emissivity": float(emissivity),
        "dt_max": 1.0e-4,
    }
    variant.update(overrides)
    return variant


def _ignition_tuning_name(T_ignition: float, k_solid: float,
                          emissivity: float) -> str:
    k_str = f"{k_solid:.1f}".replace(".", "p")
    eps_str = f"{emissivity:.2f}".replace(".", "p")
    return f"T{int(T_ignition)}_k{k_str}_eps{eps_str}"


IGNITION_TUNING_T_IGNITION = (650.0, 750.0, 850.0)
IGNITION_TUNING_K_SOLID = (0.2, 0.3, 0.5)
IGNITION_TUNING_EMISSIVITY = (0.0, 0.45)

IGNITION_TUNING_VARIANTS = {
    _ignition_tuning_name(T, k, e): _ignition_tuning_variant(T, k, e)
    for T in IGNITION_TUNING_T_IGNITION
    for k in IGNITION_TUNING_K_SOLID
    for e in IGNITION_TUNING_EMISSIVITY
}

IGNITION_TUNING_DEFAULT_VARIANTS = list(IGNITION_TUNING_VARIANTS)


def _motor_path(case: str, explicit: str | None) -> Path:
    if explicit:
        return Path(explicit)
    return MOTORS_DIR / f"{case}.ric"


def _run_from_ric_with_propellant_overrides(
    filepath, *, propellant_overrides=None, gas_props=None,
    transport_path=None, pyrogen=None, pyrogen_mass=None,
    pyrogen_throat_area=None, pyrogen_volume=None,
    pyrogen_burn_area=None, pyrogen_burn_law="0d",
    T_ignition=850.0, verbose=True, **sim_overrides,
):
    """Local diagnostic runner for probe-only propellant property overrides."""
    filepath = str(filepath)
    motor = load_ric(filepath)
    stem = str(Path(filepath).with_suffix(""))

    if gas_props is None:
        if transport_path is None:
            candidate = stem + ".transport.yaml"
            if Path(candidate).exists():
                transport_path = candidate
        if transport_path is not None:
            gas_props = load_transport(transport_path)

    args = ric_to_sim_args(motor, gas_props=gas_props, **sim_overrides)
    geo = args.pop("geo")
    prop = args.pop("propellant")
    nozzle = args["nozzle"]
    P_amb = args.get("P_ambient", 101325.0)

    for key, value in (propellant_overrides or {}).items():
        if not hasattr(prop, key):
            raise AttributeError(f"Unknown propellant override: {key}")
        setattr(prop, key, value)

    if pyrogen is None:
        candidate = stem + ".pyrogen.yaml"
        if Path(candidate).exists():
            pyrogen_obj = load_pyrogen(candidate)
        else:
            raise ValueError(
                f"No pyrogen specified for {filepath}. Pass pyrogen='bpnv', "
                "pyrogen=<Pyrogen>, or add a sibling <motor>.pyrogen.yaml."
            )
    elif isinstance(pyrogen, Pyrogen):
        pyrogen_obj = pyrogen
    else:
        pyrogen_obj = load_pyrogen(pyrogen)

    args["pyrogen_chamber"] = build_pyrogen_chamber(
        pyrogen_obj, geo, nozzle,
        pyrogen_mass=pyrogen_mass,
        pyrogen_throat_area=pyrogen_throat_area,
        pyrogen_volume=pyrogen_volume,
        pyrogen_burn_area=pyrogen_burn_area,
        pyrogen_burn_law=pyrogen_burn_law,
    )
    args["T_ignition"] = T_ignition
    args["verbose"] = verbose

    result = run_simulation(geo, prop, **args)
    perf = compute_motor_performance(result, nozzle, prop, P_ambient=P_amb)
    return result, perf, nozzle, geo, prop


def _run_variant(args, motor_path: Path, variant: str, overrides: dict):
    output_case = getattr(args, "output_case", args.case)
    output_dir = ARTIFACT_ROOT / output_case / variant
    output_dir.mkdir(parents=True, exist_ok=True)

    overrides = dict(overrides)
    propellant_overrides = {}
    if "propellant_radiation_emissivity" in overrides:
        propellant_overrides["radiation_emissivity"] = float(
            overrides.pop("propellant_radiation_emissivity")
        )
    if "propellant_k_solid" in overrides:
        propellant_overrides["k_solid"] = float(
            overrides.pop("propellant_k_solid")
        )

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
        "cfl_target": args.cfl_target,
        "diagnostic_history_capacity": args.diagnostic_history_capacity,
        "snapshot_interval": args.snapshot_interval,
        "print_interval": args.print_interval,
        "verbose": args.verbose,
        **overrides,
    }
    if args.target_propellant_cells is not None and "target_propellant_cells" not in sim_kwargs:
        sim_kwargs["target_propellant_cells"] = args.target_propellant_cells

    runner = run_from_ric
    if propellant_overrides:
        runner = _run_from_ric_with_propellant_overrides
        sim_kwargs["propellant_overrides"] = propellant_overrides

    result, perf, _nozzle, geo, prop = runner(str(motor_path), **sim_kwargs)
    diagnostics = analyze_ignition_spike(result, geo=geo, propellant=prop)
    write_diagnostic_outputs(diagnostics, output_dir, f"{args.case}_{variant}")
    plot_diagnostic_figures(result, diagnostics, output_dir, f"{args.case}_{variant}")
    plt.close("all")
    return diagnostics, result, perf


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run ignition spike diagnostics.")
    parser.add_argument("--case", default="hasegawa_a",
                        help="Case name. Defaults to srm_1d/motors/<case>.ric")
    parser.add_argument("--motor-path", default=None,
                        help="Explicit .ric path. Overrides --case lookup.")
    all_variants = (
        set(VARIANTS)
        | set(RADIATION_PROBE_VARIANTS)
        | set(RADIATION_COLLAPSE_VARIANTS)
        | set(IGNITION_TUNING_VARIANTS)
    )
    parser.add_argument(
        "--mode",
        choices=["standard", "radiation-probe", "radiation-collapse",
                 "ignition-tuning"],
        default="standard",
        help="Diagnostic matrix to run.",
    )
    parser.add_argument("--variants", nargs="+",
                        default=DEFAULT_VARIANTS,
                        choices=sorted(all_variants),
                        help="Diagnostic variants to run.")
    parser.add_argument("--pyrogen", default="bpnv")
    parser.add_argument("--pyrogen-mass", type=float, default=0.03)
    parser.add_argument("--pyrogen-throat-area", type=float, default=None)
    parser.add_argument("--roughness", type=float, default=30e-6)
    parser.add_argument("--kappa", type=float, default=0.45)
    parser.add_argument("--cfl-target", type=float, default=0.5)
    parser.add_argument("--target-propellant-cells", type=int, default=None)
    parser.add_argument("--T-ignition", dest="T_ignition", type=float, default=850.0)
    parser.add_argument("--t-max", type=float, default=3.0)
    parser.add_argument("--dt-max", type=float, default=1.0e-4)
    parser.add_argument("--diagnostic-history-cap", type=int, default=None,
                        dest="diagnostic_history_capacity",
                        help="Optional diagnostic-only max history rows.")
    parser.add_argument("--P-cutoff", dest="P_cutoff", type=float, default=0.01e6)
    parser.add_argument("--snapshot-interval", type=float, default=0.005)
    parser.add_argument("--print-interval", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    motor_path = _motor_path(args.case, args.motor_path)
    if not motor_path.exists():
        raise FileNotFoundError(f"Motor file not found: {motor_path}")

    variant_map = VARIANTS
    variants = list(args.variants)
    args.output_case = args.case
    if args.mode == "radiation-probe":
        variant_map = RADIATION_PROBE_VARIANTS
        args.output_case = f"{args.case}_radiation_probe"
        if variants == DEFAULT_VARIANTS:
            variants = list(RADIATION_PROBE_DEFAULT_VARIANTS)
    elif args.mode == "radiation-collapse":
        variant_map = RADIATION_COLLAPSE_VARIANTS
        args.output_case = f"{args.case}_radiation_collapse"
        if variants == DEFAULT_VARIANTS:
            variants = list(RADIATION_COLLAPSE_DEFAULT_VARIANTS)
    elif args.mode == "ignition-tuning":
        variant_map = IGNITION_TUNING_VARIANTS
        args.output_case = f"{args.case}_ignition_tuning"
        if variants == DEFAULT_VARIANTS:
            variants = list(IGNITION_TUNING_DEFAULT_VARIANTS)

    unknown = [variant for variant in variants if variant not in variant_map]
    if unknown:
        raise ValueError(
            f"Variants {unknown} are not valid for --mode {args.mode}."
        )

    case_dir = ARTIFACT_ROOT / args.output_case
    case_dir.mkdir(parents=True, exist_ok=True)

    # Optional MSE-vs-experimental for ignition-tuning mode. The plan's
    # Step 4 pass criterion is MSE < 0.15 MPa^2 vs the Hasegawa A
    # experimental trace; pre-compute the target so per-variant MSE is
    # cheap.
    experimental_target = None
    if args.case == "hasegawa_a" and args.mode == "ignition-tuning":
        try:
            from srm_1d.plotting import HASEGAWA_MOTOR_A_EXPERIMENTAL as _hex
            experimental_target = (_hex["time"], _hex["pressure"])
        except Exception:  # noqa: BLE001 -- plotting module is optional
            experimental_target = None

    rows = []
    for variant in variants:
        print(f"Running {args.output_case}/{variant}...")
        diagnostics, _result, _perf = _run_variant(
            args, motor_path, variant, variant_map[variant]
        )
        print(classification_report(diagnostics))
        cls = diagnostics["classification"]
        pressure = diagnostics["pressure"]
        spread = diagnostics["ignition_spread"]
        early = diagnostics["early"]
        mse_vs_exp = float("nan")
        if experimental_target is not None:
            import numpy as _np
            t_sim = _np.asarray(_result["time"])
            p_sim_mpa = _np.asarray(_result["P_head"]) / 1.0e6
            t_exp, p_exp_mpa = experimental_target
            # Compare only over the experimental support window
            mask = (t_exp >= max(t_sim.min(), 0.01)) & (t_exp <= t_sim.max())
            if mask.any():
                p_at_exp = _np.interp(t_exp[mask], t_sim, p_sim_mpa)
                mse_vs_exp = float(
                    _np.mean((p_at_exp - p_exp_mpa[mask]) ** 2)
                )
        rows.append({
            "variant": variant,
            "mse_vs_experimental_mpa2": mse_vs_exp,
            "primary_driver": cls["primary_driver"],
            "startup_window_peak_time_s": pressure["startup_window_peak_time_s"],
            "startup_window_peak_pressure_mpa": pressure["startup_window_peak_pressure_mpa"],
            "global_peak_time_s": pressure["global_peak_time_s"],
            "global_peak_pressure_mpa": pressure["global_peak_pressure_mpa"],
            "spread_10_90_s": spread["spread_10_90_s"],
            "spread_metric_source": spread["spread_metric_source"],
            "exact_spread_metrics": spread["exact_spread_metrics"],
            "instant_ignition_collapse": spread["instant_ignition_collapse"],
            "diagnostic_failure_mode": early["diagnostic_failure_mode"],
            "collapse_class": early["collapse_class"],
            "collapse_branch_suspect": early["collapse_branch_suspect"],
            "collapse_detected": early["collapse_detected"],
            "termination": early["termination"],
            "history_cap_reached": early["history_cap_reached"],
            "steps": early["steps"],
            "final_time_s": early["final_time_s"],
            "dt_min_s": early["dt_min_s"],
            "dt_median_s": early["dt_median_s"],
            "first_dt_collapse_time_s": early["first_dt_collapse_time_s"],
            "first_pressure_collapse_time_s": early["first_pressure_collapse_time_s"],
            "first_mach_collapse_time_s": early["first_mach_collapse_time_s"],
            "first_clipping_dominated_time_s": early["first_clipping_dominated_time_s"],
            "first_collapse_time_s": early["first_collapse_time_s"],
            "pressure_time_1_mpa_s": early["pressure_time_1_mpa_s"],
            "pressure_time_2_mpa_s": early["pressure_time_2_mpa_s"],
            "pressure_time_5_mpa_s": early["pressure_time_5_mpa_s"],
            "pressure_time_10_mpa_s": early["pressure_time_10_mpa_s"],
            "first_ignition_time_s": early["first_ignition_time_s"],
            "first_ignition_cell": early["first_ignition_cell"],
            "first_sustained_erosive_time_s": early["first_sustained_erosive_time_s"],
            "first_full_grain_ignition_time_s": early["first_full_grain_ignition_time_s"],
            "active_radiation_emissivity": early["active_radiation_emissivity"],
            "max_radiation_emitters": early["max_radiation_emitters"],
            "max_radiation_receivers": early["max_radiation_receivers"],
            "first_radiation_time_s": early["first_radiation_time_s"],
            "radiation_enabled_zero_activity": early["radiation_enabled_zero_activity"],
            "radiation_energy_mismatch": early["radiation_energy_mismatch"],
            "clipping_dominated_energy": early["clipping_dominated_energy"],
            "max_clipping_power_abs_w": early["max_clipping_power_abs_w"],
            "max_convective_scalar_flux_power_abs_w": early["max_convective_scalar_flux_power_abs_w"],
            "max_thermal_source_power_abs_w": early["max_thermal_source_power_abs_w"],
            "max_energy_residual_abs_w": early["max_energy_residual_abs_w"],
            "energy_residual_relative": early["energy_residual_relative"],
            "max_momentum_residual_abs_n": early["max_momentum_residual_abs_n"],
            "residuals_closed": early["residuals_closed"],
            "max_nozzle_massflow_kg_s": early["max_nozzle_massflow_kg_s"],
            "max_pressure_pa": early["max_pressure_pa"],
            "max_gas_temperature_k": early["max_gas_temperature_k"],
            "max_surface_temperature_k": early["max_surface_temperature_k"],
            "max_mach": early["max_mach"],
        })

    with open(case_dir / f"{args.output_case}_variant_summary.csv", "w", newline="") as f:
        fieldnames = [
            "variant", "mse_vs_experimental_mpa2", "primary_driver",
            "startup_window_peak_time_s", "startup_window_peak_pressure_mpa",
            "global_peak_time_s", "global_peak_pressure_mpa",
            "spread_10_90_s", "spread_metric_source", "exact_spread_metrics",
            "instant_ignition_collapse",
            "diagnostic_failure_mode", "collapse_class",
            "collapse_branch_suspect", "collapse_detected",
            "termination", "history_cap_reached",
            "steps", "final_time_s", "dt_min_s", "dt_median_s",
            "first_dt_collapse_time_s", "first_pressure_collapse_time_s",
            "first_mach_collapse_time_s", "first_clipping_dominated_time_s",
            "first_collapse_time_s",
            "pressure_time_1_mpa_s", "pressure_time_2_mpa_s",
            "pressure_time_5_mpa_s", "pressure_time_10_mpa_s",
            "first_ignition_time_s", "first_ignition_cell",
            "first_sustained_erosive_time_s", "first_full_grain_ignition_time_s",
            "active_radiation_emissivity", "max_radiation_emitters",
            "max_radiation_receivers", "first_radiation_time_s",
            "radiation_enabled_zero_activity", "radiation_energy_mismatch",
            "clipping_dominated_energy", "max_clipping_power_abs_w",
            "max_convective_scalar_flux_power_abs_w",
            "max_thermal_source_power_abs_w", "max_energy_residual_abs_w",
            "energy_residual_relative", "max_momentum_residual_abs_n",
            "residuals_closed", "max_nozzle_massflow_kg_s",
            "max_pressure_pa", "max_gas_temperature_k",
            "max_surface_temperature_k", "max_mach",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if args.mode in ("radiation-probe", "radiation-collapse"):
        with open(case_dir / f"{args.output_case}_literature_check.txt", "w") as f:
            f.write(literature_evaluation_report(args.output_case, rows))
    print(f"Wrote diagnostics under {case_dir}")


if __name__ == "__main__":
    main()
