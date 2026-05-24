"""
cross_motor_frozen_vs_effective.py — v0.7.2 cleanup
====================================================

Runs the four fired motors (Hasegawa A, Zerox, BALLSstick, Chunc) at
the SAME default-knob configuration TWICE per motor: once with the
FROZEN transport YAML, once with the EFFECTIVE YAML. Overlays the two
sim traces (plus experimental, where available) per motor so we can
see whether the Phase 5 structural ignition-kernel diagnosis (cross-
motor spike pattern is gas-transport-independent) generalizes beyond
Hasegawa A.

This is the cleanup for v0.7.1 Phase 5 Task 2's blind spot: Task 2 ran
each non-Hasegawa motor only once (at frozen defaults), so we never
established whether effective transport would shift their behavior the
way it shifted Hasegawa A's.

All four motors share knobs:
    roughness = 35 um
    kappa     = 0.45
    T_ignition= 900 K
    k_solid   = 0.4 W/(m.K)
    pyrogen sizing = Sutton default

Pyrogen species per motor: bpnv for Hasegawa A, Zerox, BALLSstick;
mtv for Chunc/machbusterNew (mirrors cross_motor_survey_task2.py).

Outputs:
    artifacts/cross_motor_frozen_vs_effective/<stamp>/
        <motor>_compare.png   per-motor frozen+effective+exp overlay
        summary.txt           tabulated P_peak / t_peak / spike_ratio per motor x transport
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import (
    HASEGAWA_MOTOR_A_EXPERIMENTAL,
    ZEROX_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir


MOTORS_DIR = Path(__file__).resolve().parents[1] / 'motors'


DEFAULT_KNOBS = dict(
    roughness=35.0e-6,
    kappa=0.45,
    T_ignition=900.0,
    k_solid=0.4,
    pyrogen_mass=None,
    pyrogen_throat_area=None,
    pyrogen_volume=None,
    P_cutoff=0.05e6,
    cfl_target=0.3,
    snapshot_interval=0.5,
    print_interval=0.5,
    verbose=False,
)


MOTORS = [
    dict(
        label='Hasegawa A',
        slug='hasegawa_a',
        ric=MOTORS_DIR / 'hasegawa_a.ric',
        frozen_yaml=MOTORS_DIR / 'hasegawa_a.frozen.transport.yaml',
        effective_yaml=MOTORS_DIR / 'hasegawa_a.transport.yaml',
        pyrogen='bpnv',
        t_max=3.0,
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=1.1,
    ),
    dict(
        label='Zerox',
        slug='zerox',
        ric=MOTORS_DIR / 'zerox.ric',
        frozen_yaml=MOTORS_DIR / 'zerox.frozen.transport.yaml',
        effective_yaml=MOTORS_DIR / 'zerox.transport.yaml',
        pyrogen='bpnv',
        t_max=8.0,
        experimental=ZEROX_EXPERIMENTAL,
        time_offset=ZEROX_EXPERIMENTAL.get('time_offset', 0.0),
    ),
    dict(
        label='BALLSstick',
        slug='BALLSstick',
        ric=MOTORS_DIR / 'BALLSstick.ric',
        frozen_yaml=MOTORS_DIR / 'BALLSstick.frozen.transport.yaml',
        effective_yaml=MOTORS_DIR / 'BALLSstick.transport.yaml',
        pyrogen='bpnv',
        t_max=5.0,
        experimental=None,
        time_offset=0.0,
    ),
    dict(
        label='Chunc (machbusterNew)',
        slug='machbusterNew',
        ric=MOTORS_DIR / 'machbusterNew.ric',
        frozen_yaml=MOTORS_DIR / 'machbusterNew.frozen.transport.yaml',
        effective_yaml=MOTORS_DIR / 'machbusterNew.transport.yaml',
        pyrogen='mtv',
        t_max=3.0,
        experimental=None,
        time_offset=0.0,
    ),
]


def run_motor_with_transport(motor, transport_path):
    result, perf, nozzle, geo, prop = run_from_ric(
        str(motor['ric']),
        transport_path=str(transport_path),
        pyrogen=motor['pyrogen'],
        t_max=motor['t_max'],
        **DEFAULT_KNOBS,
    )
    return result


def metrics(result):
    t = np.asarray(result['time'])
    P = np.asarray(result['P_head']) / 1e6
    P_peak = float(P.max())
    i_peak = int(P.argmax())
    t_peak = float(t[i_peak])
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
    )


def main():
    out = artifact_dir('cross_motor_frozen_vs_effective')
    print(f"Outputs -> {out}\n")

    rows = []

    for motor in MOTORS:
        print(f"\n=== {motor['label']} ===")
        print(f"  Running FROZEN...")
        r_frozen = run_motor_with_transport(motor, motor['frozen_yaml'])
        m_frozen = metrics(r_frozen)
        print(f"    P_peak={m_frozen['P_peak']:.2f} MPa at t={m_frozen['t_peak']:.3f}s, "
              f"plateau={m_frozen['P_plateau']:.2f} MPa, "
              f"spike/plateau={m_frozen['spike_ratio']:.2f}")

        print(f"  Running EFFECTIVE...")
        r_eff = run_motor_with_transport(motor, motor['effective_yaml'])
        m_eff = metrics(r_eff)
        print(f"    P_peak={m_eff['P_peak']:.2f} MPa at t={m_eff['t_peak']:.3f}s, "
              f"plateau={m_eff['P_plateau']:.2f} MPa, "
              f"spike/plateau={m_eff['spike_ratio']:.2f}")

        rows.append((motor, m_frozen, m_eff, r_frozen, r_eff))

        # Per-motor comparison plot
        fig, ax = plt.subplots(figsize=(12, 7))
        if motor['experimental'] is not None:
            exp = motor['experimental']
            ax.plot(exp['time'] + motor['time_offset'], exp['pressure'],
                    'k-', linewidth=1.5, marker='o', markersize=3,
                    markevery=max(1, len(exp['time']) // 25),
                    label=exp.get('label', 'Experimental'))
        ax.plot(r_frozen['time'], r_frozen['P_head'] / 1e6, 'b-', linewidth=2,
                label=f"Frozen (P_peak={m_frozen['P_peak']:.2f} MPa, "
                      f"ratio={m_frozen['spike_ratio']:.2f})")
        ax.plot(r_eff['time'], r_eff['P_head'] / 1e6, 'r-', linewidth=2,
                label=f"Effective (P_peak={m_eff['P_peak']:.2f} MPa, "
                      f"ratio={m_eff['spike_ratio']:.2f})")
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Head-End Pressure [MPa]', fontsize=12)
        ax.set_title(
            f"{motor['label']} — Frozen vs Effective Transport at Default Knobs\n"
            f"(roughness=35um, kappa=0.45, T_ign=900K, k_solid=0.4, "
            f"{motor['pyrogen']} Sutton-default)",
            fontsize=10,
        )
        ax.grid(alpha=0.3)
        ax.legend(loc='best', fontsize=10)
        png = out / f"{motor['slug']}_compare.png"
        fig.savefig(png, dpi=130, bbox_inches='tight')
        plt.close(fig)
        print(f"  Plot: {png.name}")

    # Summary table
    lines = [
        "v0.7.2 cleanup — Cross-Motor Frozen vs Effective Survey",
        "=" * 88,
        "",
        "Same default knobs across all motors (k_solid=0.4, roughness=35um,",
        "kappa=0.45, T_ign=900K, pyrogen sizing = Sutton default).",
        "Each motor run twice: once with frozen YAML, once with effective.",
        "",
        f"{'Motor':<22} {'YAML':<10} {'P_peak [MPa]':>14} {'t_peak [s]':>12} "
        f"{'plateau':>10} {'spike/plat':>12}",
        "-" * 90,
    ]
    for motor, m_frozen, m_eff, _, _ in rows:
        lines.append(
            f"{motor['label']:<22} {'frozen':<10} "
            f"{m_frozen['P_peak']:>14.3f} {m_frozen['t_peak']:>12.4f} "
            f"{m_frozen['P_plateau']:>10.3f} {m_frozen['spike_ratio']:>12.2f}"
        )
        lines.append(
            f"{'':<22} {'effective':<10} "
            f"{m_eff['P_peak']:>14.3f} {m_eff['t_peak']:>12.4f} "
            f"{m_eff['P_plateau']:>10.3f} {m_eff['spike_ratio']:>12.2f}"
        )
        # Delta row
        d_peak = m_eff['P_peak'] - m_frozen['P_peak']
        d_ratio = m_eff['spike_ratio'] - m_frozen['spike_ratio']
        lines.append(
            f"{'':<22} {'  delta':<10} "
            f"{d_peak:>+14.3f} {'':>12} "
            f"{'':>10} {d_ratio:>+12.2f}"
        )

    lines += [
        "",
        "Diagnostic reading guide:",
        "  - If effective YAML meaningfully REDUCES P_peak or spike_ratio for",
        "    Zerox/Chunc/BALLSstick (the way the v0.7.0 frozen->effective LHS",
        "    shifted Hasegawa A's k_solid), then the cross-motor systematic IS",
        "    partly gas-transport-driven and v0.7.2 should include per-motor",
        "    effective LHS recalibration before structural kernel work.",
        "  - If effective YAML AMPLIFIES the over-prediction (as the Zerox YAML",
        "    historical comment claims), the structural ignition-kernel diagnosis",
        "    locks: no transport choice resolves it, and v0.7.2 should go",
        "    straight to Z-N / spatial ignition front / pyrogen distribution.",
        "  - If effective leaves the spike essentially unchanged, that also",
        "    points to the structural diagnosis (the kernel artifact dominates",
        "    the gas-transport contribution).",
    ]
    summary_path = out / 'summary.txt'
    summary_path.write_text('\n'.join(lines))
    print(f"\nSummary: {summary_path}")
    print()
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
