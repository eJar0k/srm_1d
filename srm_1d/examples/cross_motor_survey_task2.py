"""
cross_motor_survey_task2.py — v0.7.1 Phase 5 Task 2 cross-motor spike survey
=============================================================================

Runs the four fired motors (Hasegawa A, Zerox, BALLSstick, Chunc) at the
SAME default-knob configuration to test whether the v0.7.1 build
systematically over-predicts ignition-spike magnitudes across motors.
All runs use:

    roughness = 35 µm        (Ma 2020 mid-band literature)
    kappa     = 0.45         (Gnielinski mid-band)
    T_ignition= 900 K        (Goodman mid-band)
    k_solid   = 0.4 W/(m·K)  (upper edge of AP/HTPB+Al 0.20-0.40 literature)
    pyrogen sizing = Sutton default (m=0.12·V_F^0.7, throat & volume scaled)
    transport YAML = CURRENT (frozen) per-motor YAMLs

Pyrogen species are motor-specific to honor the user's existing scripts:
    Hasegawa A    → bpnv  (canonical Hasegawa pyrogen)
    Zerox         → bpnv
    BALLSstick    → bpnv
    Chunc         → mtv   (user's pre-staged choice for machbusterNew)

Outputs:
    artifacts/cross_motor_survey_task2/<stamp>/
        <motor>_pressure.png   per-motor head-end pressure (with exp overlay if available)
        summary.txt            P_peak / t_peak / steady plateau / spike ratio per motor
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import (
    plot_pressure,
    HASEGAWA_MOTOR_A_EXPERIMENTAL,
    ZEROX_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir


MOTORS_DIR = Path(__file__).resolve().parents[1] / 'motors'


# Same default-knob configuration applied to every motor.
DEFAULT_KNOBS = dict(
    roughness=35.0e-6,
    kappa=0.45,
    T_ignition=900.0,
    k_solid=0.4,
    pyrogen_mass=None,         # Sutton default
    pyrogen_throat_area=None,  # Sutton default
    pyrogen_volume=None,       # Sutton default
    P_cutoff=0.05e6,
    cfl_target=0.3,
    snapshot_interval=0.5,
    print_interval=0.5,
    verbose=False,
)


# Per-motor differences: ric path, pyrogen, time cutoff, experimental overlay.
MOTORS = [
    dict(
        label='Hasegawa A',
        slug='hasegawa_a',
        ric=MOTORS_DIR / 'hasegawa_a.ric',
        pyrogen='bpnv',
        t_max=3.0,
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=1.1,
    ),
    dict(
        label='Zerox',
        slug='zerox',
        ric=MOTORS_DIR / 'zerox.ric',
        pyrogen='bpnv',
        t_max=8.0,
        experimental=ZEROX_EXPERIMENTAL,
        time_offset=ZEROX_EXPERIMENTAL.get('time_offset', 0.0),
    ),
    dict(
        label='BALLSstick',
        slug='BALLSstick',
        ric=MOTORS_DIR / 'BALLSstick.ric',
        pyrogen='bpnv',
        t_max=5.0,
        experimental=None,
        time_offset=0.0,
    ),
    dict(
        label='Chunc (machbusterNew)',
        slug='machbusterNew',
        ric=MOTORS_DIR / 'machbusterNew.ric',
        pyrogen='mtv',
        t_max=3.0,
        experimental=None,
        time_offset=0.0,
    ),
]


def run_motor(motor):
    print(f"\n=== {motor['label']} ({motor['ric'].name}) ===")
    result, perf, nozzle, geo, prop = run_from_ric(
        str(motor['ric']),
        pyrogen=motor['pyrogen'],
        t_max=motor['t_max'],
        **DEFAULT_KNOBS,
    )
    return result


def metrics(result):
    """Pull spike-magnitude diagnostics from a run."""
    t = np.asarray(result['time'])
    P = np.asarray(result['P_head']) / 1e6  # MPa

    P_peak = float(P.max())
    i_peak = int(P.argmax())
    t_peak = float(t[i_peak])

    # Steady plateau: median P over the middle third of the trace
    # (skips ignition transient and tail-off).
    n = len(P)
    if n > 6:
        lo, hi = n // 3, 2 * n // 3
        P_plateau = float(np.median(P[lo:hi]))
    else:
        P_plateau = float(P.mean())

    spike_ratio = P_peak / P_plateau if P_plateau > 0.0 else float('nan')

    return dict(
        P_peak=P_peak, t_peak=t_peak,
        P_plateau=P_plateau, spike_ratio=spike_ratio,
        n_samples=n,
    )


def main():
    out = artifact_dir('cross_motor_survey_task2')
    print(f"Outputs -> {out}\n")

    rows = []

    for motor in MOTORS:
        result = run_motor(motor)
        m = metrics(result)
        rows.append((motor, m, result))

        # Per-motor pressure plot.
        png = out / f"{motor['slug']}_pressure.png"
        plot_pressure(
            result,
            title=(f"{motor['label']} — default knobs (k_solid=0.4, roughness=35µm, "
                   f"kappa=0.45, T_ign=900K, {motor['pyrogen']})"),
            experimental=motor['experimental'],
            time_offset=motor['time_offset'],
            save_path=str(png),
        )
        plt.close('all')
        print(f"  P_peak={m['P_peak']:.2f} MPa at t={m['t_peak']:.3f}s, "
              f"plateau={m['P_plateau']:.2f} MPa, spike/plateau={m['spike_ratio']:.2f}")
        print(f"  Plot: {png.name}")

    # ----- Summary table ---------------------------------------------------
    lines = [
        "v0.7.1 Phase 5 Task 2 — Cross-Motor Spike Survey",
        "=" * 78,
        "",
        "Default knobs (identical across all motors):",
        f"  roughness  = 35 um",
        f"  kappa      = 0.45",
        f"  T_ignition = 900 K",
        f"  k_solid    = 0.4 W/(m.K)   (upper edge of 0.20-0.40 AP/HTPB+Al band)",
        f"  pyrogen sizing = Sutton default",
        f"  transport YAML = current frozen YAMLs",
        "",
        f"{'Motor':<22} {'Pyrogen':<7} {'P_peak [MPa]':>14} {'t_peak [s]':>12} "
        f"{'P_plateau':>12} {'P_peak/P_plateau':>18}",
        "-" * 92,
    ]
    for motor, m, _ in rows:
        lines.append(
            f"{motor['label']:<22} {motor['pyrogen']:<7} "
            f"{m['P_peak']:>14.3f} {m['t_peak']:>12.4f} "
            f"{m['P_plateau']:>12.3f} {m['spike_ratio']:>18.2f}"
        )

    lines += [
        "",
        "Reading guide:",
        "  - P_peak >> P_plateau (ratio > ~1.3) suggests an ignition transient",
        "    that exceeds the steady-state burn, which the handoff memory flags",
        "    as the v0.7.1 cross-motor systematic.",
        "  - If 4/4 show spike_ratio > 1.3, the systematic is confirmed and",
        "    motivates re-LHS with effective k_gas (Task 3 path A).",
        "  - If only some over-predict, motor-size/L-D pattern narrows the",
        "    diagnosis.",
    ]

    summary_path = out / 'summary.txt'
    summary_path.write_text('\n'.join(lines))
    print(f"\nSummary: {summary_path}")
    print()
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
