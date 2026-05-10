"""
Ignition-spike diagnostics for v0.7.0 startup investigations.

The functions here are post-processing only. They consume a simulation
``result`` plus optional geometry/propellant objects and estimate which
source family is active when the head-end pressure spike occurs.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:  # pragma: no cover - exercised only without matplotlib
    HAS_MATPLOTLIB = False


def stack_snapshots(result: dict, key: str) -> np.ndarray:
    """Return a ``(n_snapshots, n_cells)`` array for a snapshot key."""
    snapshots = result.get("snapshots", [])
    if not snapshots:
        return np.empty((0, 0))
    if key not in snapshots[0]:
        return np.empty((len(snapshots), 0))
    return np.array([np.asarray(s[key]) for s in snapshots])


def _snapshot_times(result: dict) -> np.ndarray:
    return np.array([float(s["t"]) for s in result.get("snapshots", [])])


def _infer_dx(result: dict, geo=None) -> float:
    if geo is not None and getattr(geo, "dx", None) is not None:
        return float(geo.dx)
    snapshots = result.get("snapshots", [])
    if not snapshots:
        return 1.0
    x = np.asarray(snapshots[0].get("x", []), dtype=float)
    if len(x) > 1:
        return float(np.median(np.diff(x)))
    return 1.0


def _rho_propellant(propellant=None) -> float:
    if propellant is None:
        return 1.0
    return float(getattr(propellant, "rho_propellant", 1.0))


def _first_time_at_fraction(times: np.ndarray, frac: np.ndarray,
                            threshold: float) -> float:
    idx = np.flatnonzero(frac >= threshold)
    if idx.size == 0:
        return float("nan")
    return float(times[int(idx[0])])


def pressure_landmarks(result: dict) -> dict[str, float]:
    """Pressure peak, takeoff, and post-peak trough timing."""
    t = np.asarray(result.get("time", []), dtype=float)
    p = np.asarray(result.get("P_head", []), dtype=float)
    if t.size == 0 or p.size == 0:
        return {
            "takeoff_time_s": float("nan"),
            "peak_time_s": float("nan"),
            "peak_pressure_pa": float("nan"),
            "peak_pressure_mpa": float("nan"),
            "post_peak_trough_time_s": float("nan"),
            "post_peak_trough_pressure_pa": float("nan"),
        }

    peak_idx = int(np.argmax(p))
    ambient = float(result.get("P_ambient", p[0]))
    threshold = ambient + 0.05 * max(float(p[peak_idx]) - ambient, 0.0)
    takeoff_candidates = np.flatnonzero(p >= threshold)
    takeoff_time = float(t[int(takeoff_candidates[0])]) if takeoff_candidates.size else float("nan")

    post = p[peak_idx:]
    trough_rel = int(np.argmin(post)) if post.size else 0
    trough_idx = peak_idx + trough_rel
    return {
        "takeoff_time_s": takeoff_time,
        "peak_time_s": float(t[peak_idx]),
        "peak_pressure_pa": float(p[peak_idx]),
        "peak_pressure_mpa": float(p[peak_idx] / 1.0e6),
        "post_peak_trough_time_s": float(t[trough_idx]),
        "post_peak_trough_pressure_pa": float(p[trough_idx]),
        "post_peak_trough_pressure_mpa": float(p[trough_idx] / 1.0e6),
    }


def pyrogen_landmarks(result: dict, peak_time_s: float | None = None) -> dict[str, float | bool]:
    """Pyrogen pressure/mass-flow timing from time-history arrays."""
    t = np.asarray(result.get("time", []), dtype=float)
    mdot = np.asarray(result.get("mdot_ig", []), dtype=float)
    pig = np.asarray(result.get("P_ig", []), dtype=float)
    if t.size == 0 or mdot.size == 0:
        return {
            "mdot_peak_time_s": float("nan"),
            "mdot_peak_kg_s": float("nan"),
            "pyrogen_active_end_s": float("nan"),
            "P_ig_peak_time_s": float("nan"),
            "P_ig_peak_pa": float("nan"),
            "pyrogen_active_at_peak": False,
            "mdot_at_pressure_peak_kg_s": float("nan"),
        }

    mdot_peak = float(np.max(mdot))
    mdot_threshold = max(1.0e-9, 0.01 * mdot_peak)
    active = mdot > mdot_threshold
    mdot_peak_idx = int(np.argmax(mdot))
    active_end = float(t[int(np.flatnonzero(active)[-1])]) if np.any(active) else float("nan")
    pig_peak_idx = int(np.argmax(pig)) if pig.size else 0

    mdot_at_peak = float("nan")
    active_at_peak = False
    if peak_time_s is not None and np.isfinite(peak_time_s):
        mdot_at_peak = float(np.interp(peak_time_s, t, mdot))
        active_at_peak = bool(mdot_at_peak > mdot_threshold)

    return {
        "mdot_peak_time_s": float(t[mdot_peak_idx]),
        "mdot_peak_kg_s": mdot_peak,
        "mdot_peak_g_s": mdot_peak * 1000.0,
        "pyrogen_active_end_s": active_end,
        "P_ig_peak_time_s": float(t[pig_peak_idx]) if pig.size else float("nan"),
        "P_ig_peak_pa": float(pig[pig_peak_idx]) if pig.size else float("nan"),
        "P_ig_peak_mpa": float(pig[pig_peak_idx] / 1.0e6) if pig.size else float("nan"),
        "pyrogen_active_at_peak": active_at_peak,
        "mdot_at_pressure_peak_kg_s": mdot_at_peak,
        "mdot_at_pressure_peak_g_s": mdot_at_peak * 1000.0 if np.isfinite(mdot_at_peak) else float("nan"),
    }


def ignition_spread_metrics(result: dict) -> dict[str, Any]:
    """Burning-fraction time history and inferred first-burning times."""
    times = _snapshot_times(result)
    if times.size == 0:
        return {
            "times_s": times,
            "burning_fraction": np.array([]),
            "first_ignition_time_s": float("nan"),
            "t10_s": float("nan"),
            "t50_s": float("nan"),
            "t90_s": float("nan"),
            "t100_s": float("nan"),
            "spread_10_90_s": float("nan"),
            "instant_ignition_collapse": False,
            "first_burning_time_by_cell_s": np.array([]),
            "axial_ignition_order": np.array([], dtype=int),
        }

    burning = stack_snapshots(result, "is_burning").astype(bool)
    is_grain = stack_snapshots(result, "is_grain").astype(bool)
    grain_mask = np.any(is_grain, axis=0) if is_grain.size else np.ones(burning.shape[1], dtype=bool)
    n_grain = max(int(np.sum(grain_mask)), 1)
    burning_fraction = np.sum(burning[:, grain_mask], axis=1) / n_grain

    first_by_cell = np.full(burning.shape[1], np.nan)
    for i in range(burning.shape[1]):
        idx = np.flatnonzero(burning[:, i])
        if idx.size:
            first_by_cell[i] = times[int(idx[0])]
    finite_cells = np.flatnonzero(np.isfinite(first_by_cell))
    order = finite_cells[np.argsort(first_by_cell[finite_cells])]

    t10 = _first_time_at_fraction(times, burning_fraction, 0.10)
    t50 = _first_time_at_fraction(times, burning_fraction, 0.50)
    t90 = _first_time_at_fraction(times, burning_fraction, 0.90)
    t100 = _first_time_at_fraction(times, burning_fraction, 0.999)
    if times.size > 1:
        snapshot_dt = float(np.median(np.diff(times)))
    else:
        snapshot_dt = float("inf")
    spread_10_90 = t90 - t10 if np.isfinite(t10) and np.isfinite(t90) else float("nan")
    instant = bool(np.isfinite(spread_10_90) and spread_10_90 <= 2.0 * snapshot_dt)

    ignited = np.flatnonzero(burning_fraction > 0.0)
    first_ignition = float(times[int(ignited[0])]) if ignited.size else float("nan")
    return {
        "times_s": times,
        "burning_fraction": burning_fraction,
        "first_ignition_time_s": first_ignition,
        "t10_s": t10,
        "t50_s": t50,
        "t90_s": t90,
        "t100_s": t100,
        "spread_10_90_s": spread_10_90,
        "instant_ignition_collapse": instant,
        "first_burning_time_by_cell_s": first_by_cell,
        "axial_ignition_order": order,
    }


def source_timeseries(result: dict, geo=None, propellant=None) -> dict[str, np.ndarray]:
    """Estimate source-family mass-flow rates at snapshot times."""
    times = _snapshot_times(result)
    if times.size == 0:
        empty = np.array([])
        return {
            "times_s": empty,
            "normal_sidewall_kg_s": empty,
            "erosive_sidewall_kg_s": empty,
            "total_sidewall_kg_s": empty,
            "endface_kg_s": empty,
            "pyrogen_kg_s": empty,
            "total_estimated_kg_s": empty,
            "normal_fraction": empty,
            "erosive_fraction": empty,
            "endface_fraction": empty,
            "pyrogen_fraction": empty,
        }

    dx = _infer_dx(result, geo)
    rho_p = _rho_propellant(propellant)
    c_burn = stack_snapshots(result, "C_burn")
    r_total = stack_snapshots(result, "r_total")
    r_erosive = stack_snapshots(result, "r_erosive")
    endface = stack_snapshots(result, "endface_msource")

    normal_rate = np.maximum(r_total - r_erosive, 0.0)
    normal = rho_p * np.sum(normal_rate * c_burn, axis=1) * dx
    erosive = rho_p * np.sum(np.maximum(r_erosive, 0.0) * c_burn, axis=1) * dx
    sidewall = rho_p * np.sum(np.maximum(r_total, 0.0) * c_burn, axis=1) * dx
    endface_total = np.sum(np.maximum(endface, 0.0), axis=1) * dx if endface.size else np.zeros_like(times)

    hist_t = np.asarray(result.get("time", []), dtype=float)
    mdot = np.asarray(result.get("mdot_ig", []), dtype=float)
    if hist_t.size and mdot.size:
        pyrogen = np.interp(times, hist_t, mdot)
    else:
        pyrogen = np.zeros_like(times)

    total = normal + erosive + endface_total + pyrogen
    denom = np.where(total > 0.0, total, 1.0)
    return {
        "times_s": times,
        "normal_sidewall_kg_s": normal,
        "erosive_sidewall_kg_s": erosive,
        "total_sidewall_kg_s": sidewall,
        "endface_kg_s": endface_total,
        "pyrogen_kg_s": pyrogen,
        "total_estimated_kg_s": total,
        "normal_fraction": normal / denom,
        "erosive_fraction": erosive / denom,
        "endface_fraction": endface_total / denom,
        "pyrogen_fraction": pyrogen / denom,
    }


def _nearest_source_at_peak(sources: dict[str, np.ndarray], peak_time_s: float) -> dict[str, float]:
    times = sources.get("times_s", np.array([]))
    if times.size == 0 or not np.isfinite(peak_time_s):
        return {
            "normal_fraction_at_peak": float("nan"),
            "erosive_fraction_at_peak": float("nan"),
            "endface_fraction_at_peak": float("nan"),
            "pyrogen_fraction_at_peak": float("nan"),
            "source_snapshot_time_s": float("nan"),
        }
    idx = int(np.argmin(np.abs(times - peak_time_s)))
    return {
        "normal_fraction_at_peak": float(sources["normal_fraction"][idx]),
        "erosive_fraction_at_peak": float(sources["erosive_fraction"][idx]),
        "endface_fraction_at_peak": float(sources["endface_fraction"][idx]),
        "pyrogen_fraction_at_peak": float(sources["pyrogen_fraction"][idx]),
        "source_snapshot_time_s": float(times[idx]),
    }


def classify_driver(pressure: dict[str, float], pyrogen: dict[str, Any],
                    spread: dict[str, Any], sources_at_peak: dict[str, float]) -> dict[str, Any]:
    """Classify the most likely pressure-spike driver."""
    pyrogen_active = bool(pyrogen.get("pyrogen_active_at_peak", False))
    pyro_frac = float(sources_at_peak.get("pyrogen_fraction_at_peak", np.nan))
    erosive_frac = float(sources_at_peak.get("erosive_fraction_at_peak", np.nan))
    endface_frac = float(sources_at_peak.get("endface_fraction_at_peak", np.nan))
    instant = bool(spread.get("instant_ignition_collapse", False))

    flags = {
        "pyrogen_combustion": pyrogen_active and (not np.isfinite(pyro_frac) or pyro_frac >= 0.25),
        "erosive_snap_on": np.isfinite(erosive_frac) and erosive_frac >= 0.35 and (
            not pyrogen_active or erosive_frac >= pyro_frac
        ),
        "endface_sources": np.isfinite(endface_frac) and endface_frac >= 0.35 and (
            not pyrogen_active or endface_frac >= pyro_frac
        ),
        "immediate_goodman_spread": instant,
    }

    if flags["pyrogen_combustion"]:
        primary = "pyrogen_combustion"
    elif flags["erosive_snap_on"]:
        primary = "erosive_snap_on"
    elif flags["endface_sources"]:
        primary = "endface_sources"
    elif flags["immediate_goodman_spread"]:
        primary = "immediate_goodman_spread"
    else:
        primary = "mixed_or_unclassified"

    return {
        "primary_driver": primary,
        **flags,
        "peak_time_s": pressure.get("peak_time_s", float("nan")),
    }


def analyze_ignition_spike(result: dict, geo=None, propellant=None) -> dict[str, Any]:
    """Run the full diagnostic reduction for one simulation result."""
    pressure = pressure_landmarks(result)
    pyrogen = pyrogen_landmarks(result, pressure["peak_time_s"])
    spread = ignition_spread_metrics(result)
    sources = source_timeseries(result, geo=geo, propellant=propellant)
    sources_at_peak = _nearest_source_at_peak(sources, pressure["peak_time_s"])
    classification = classify_driver(pressure, pyrogen, spread, sources_at_peak)
    return {
        "pressure": pressure,
        "pyrogen": pyrogen,
        "ignition_spread": spread,
        "sources": sources,
        "sources_at_peak": sources_at_peak,
        "classification": classification,
    }


def _flatten_scalars(prefix: str, data: dict[str, Any], rows: list[tuple[str, Any]]) -> None:
    for key, value in data.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            _flatten_scalars(name, value, rows)
        elif np.isscalar(value) or isinstance(value, (str, bool)):
            rows.append((name, value))


def classification_report(diagnostics: dict[str, Any]) -> str:
    """Compact human-readable report for one diagnostic run."""
    pressure = diagnostics["pressure"]
    pyrogen = diagnostics["pyrogen"]
    spread = diagnostics["ignition_spread"]
    sources = diagnostics["sources_at_peak"]
    cls = diagnostics["classification"]
    lines = [
        f"primary_driver: {cls['primary_driver']}",
        f"pressure_peak: {pressure['peak_pressure_mpa']:.3g} MPa at {pressure['peak_time_s']:.6g} s",
        f"pyrogen_active_at_peak: {pyrogen['pyrogen_active_at_peak']} "
        f"(mdot={pyrogen['mdot_at_pressure_peak_g_s']:.3g} g/s)",
        f"ignition_spread_10_90: {spread['spread_10_90_s']:.6g} s "
        f"(instant={spread['instant_ignition_collapse']})",
        "source_fractions_at_peak: "
        f"normal={sources['normal_fraction_at_peak']:.3g}, "
        f"erosive={sources['erosive_fraction_at_peak']:.3g}, "
        f"endface={sources['endface_fraction_at_peak']:.3g}, "
        f"pyrogen={sources['pyrogen_fraction_at_peak']:.3g}",
    ]
    return "\n".join(lines) + "\n"


def write_diagnostic_outputs(diagnostics: dict[str, Any], output_dir,
                             case_name: str) -> None:
    """Write summary and time-series CSVs plus a text report."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, Any]] = []
    _flatten_scalars("", {
        "pressure": diagnostics["pressure"],
        "pyrogen": diagnostics["pyrogen"],
        "sources_at_peak": diagnostics["sources_at_peak"],
        "classification": diagnostics["classification"],
    }, rows)
    spread = diagnostics["ignition_spread"]
    for key in ("first_ignition_time_s", "t10_s", "t50_s", "t90_s",
                "t100_s", "spread_10_90_s", "instant_ignition_collapse"):
        rows.append((f"ignition_spread.{key}", spread[key]))

    with open(output_dir / f"{case_name}_summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)

    sources = diagnostics["sources"]
    source_keys = [k for k in sources if k != "times_s"]
    with open(output_dir / f"{case_name}_sources.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s"] + source_keys)
        for i, t in enumerate(sources["times_s"]):
            writer.writerow([t] + [sources[k][i] for k in source_keys])

    first_by_cell = spread["first_burning_time_by_cell_s"]
    with open(output_dir / f"{case_name}_ignition_times.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell", "first_burning_time_s"])
        for i, value in enumerate(first_by_cell):
            writer.writerow([i, value])

    with open(output_dir / f"{case_name}_classification.txt", "w") as f:
        f.write(classification_report(diagnostics))


def plot_diagnostic_figures(result: dict, diagnostics: dict[str, Any],
                            output_dir, case_name: str) -> None:
    """Write overview and x-t diagnostic plots."""
    if not HAS_MATPLOTLIB:
        return
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t = np.asarray(result.get("time", []), dtype=float)
    p = np.asarray(result.get("P_head", []), dtype=float) / 1.0e6
    sources = diagnostics["sources"]
    spread = diagnostics["ignition_spread"]

    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=False)
    axes[0].plot(t, p, "b-", linewidth=1.8)
    axes[0].axvline(diagnostics["pressure"]["peak_time_s"], color="k", linestyle=":", linewidth=1.0)
    axes[0].set_ylabel("P_head [MPa]")
    axes[0].set_title(f"{case_name}: pressure and startup drivers")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, np.asarray(result.get("mdot_ig", [])) * 1000.0,
                 color="tab:orange", linewidth=1.4, label="mdot_ig")
    axes[1].set_ylabel("mdot_ig [g/s]")
    axes[1].grid(True, alpha=0.3)

    if sources["times_s"].size:
        st = sources["times_s"]
        axes[2].plot(st, sources["normal_fraction"], label="normal")
        axes[2].plot(st, sources["erosive_fraction"], label="erosive")
        axes[2].plot(st, sources["endface_fraction"], label="endface")
        axes[2].plot(st, sources["pyrogen_fraction"], label="pyrogen")
    axes[2].set_ylabel("source fraction")
    axes[2].legend(loc="best", fontsize=8)
    axes[2].grid(True, alpha=0.3)

    if spread["times_s"].size:
        axes[3].plot(spread["times_s"], spread["burning_fraction"],
                     "g-", linewidth=1.6)
    axes[3].set_ylabel("burning fraction")
    axes[3].set_xlabel("time [s]")
    axes[3].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / f"{case_name}_overview.png", dpi=160)
    plt.close(fig)

    snapshots = result.get("snapshots", [])
    if not snapshots:
        return
    snap_t = _snapshot_times(result)
    x_mm = np.asarray(snapshots[0]["x"], dtype=float) * 1000.0
    maps = [
        ("T", "T_gas [K]", stack_snapshots(result, "T")),
        ("T_surf", "T_surf [K]", stack_snapshots(result, "T_surf")),
        ("is_burning", "is_burning [-]", stack_snapshots(result, "is_burning").astype(float)),
        ("r_erosive", "r_erosive [mm/s]", stack_snapshots(result, "r_erosive") * 1000.0),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    for ax, (_key, title, data) in zip(axes.flat, maps):
        pcm = ax.pcolormesh(snap_t, x_mm, data.T, shading="nearest")
        fig.colorbar(pcm, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("time [s]")
        ax.set_ylabel("x [mm]")
    fig.tight_layout()
    fig.savefig(output_dir / f"{case_name}_xt.png", dpi=160)
    plt.close(fig)
