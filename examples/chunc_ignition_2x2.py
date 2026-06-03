"""
chunc_ignition_2x2.py — v0.7.4 Phase FZ validation harness
==========================================================

Runs BOTH the Chunc (machbusterNew) and Hasegawa A motors under four
ignition-transient configurations and overlays them against static-fire
data:

    1. baseline          — flame-front OFF, Z-N OFF (the 16.85 MPa Chunc reference)
    2. flame_front        — Phase F only (bottom-up: gate ignition to a front)
    3. zn                 — Phase Z only (top-down: Z-N burn-rate relaxation)
    4. flame_front + zn   — combined

Primary target: Chunc's ignition spike should drop toward its flat
~8.5 MPa static-fire plateau WITHOUT depressing the plateau or taildown.
Hasegawa A runs alongside as a robustness check (must not regress).

Knobs are Hasegawa-canon-aligned (roughness 37.1 um, kappa 0.45,
T_ignition 850 K, k_solid default), cfl_target 0.3. NO pyrogen-throat
override — for head_basket the pyrogen orifice is vestigial. The intent
is to let Phase F / Z suppress the spike PHYSICALLY, not via knob tuning.

Outputs:
    artifacts/chunc_ignition_2x2/<stamp>/
        machbusterNew_2x2.png   4 configs + experimental
        hasegawa_a_2x2.png      4 configs + experimental
        summary.txt             P_peak / plateau / spike_ratio / ignition spread
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
from srm_1d.plotting import (
    HASEGAWA_MOTOR_A_EXPERIMENTAL,
    CHUNC_EXPERIMENTAL,
)
from srm_1d.run_artifacts import artifact_dir
from srm_1d.tools.ignition_diagnostics import ignition_spread_metrics


MOTORS_DIR = Path(__file__).resolve().parents[1] / 'motors'

# Hasegawa-canon-aligned knobs, shared by both motors (let F/Z do the
# spike suppression, not the knobs).
COMMON_KNOBS = dict(
    roughness=37.1e-6,
    kappa=0.45,
    T_ignition=850.0,
    P_cutoff=0.05e6,
    cfl_target=0.3,
    pyrogen_mass=None,
    pyrogen_volume=None,
    snapshot_interval=0.05,   # fine enough that spread metrics are meaningful
    print_interval=1.0,
    verbose=False,
)

# (label, flame_front_enabled, zn_enabled)
CONFIGS = [
    ('baseline', False, False),
    ('flame_front', True, False),
    ('zn', False, True),
    ('flame_front+zn', True, True),
]
CONFIG_COLORS = {
    'baseline': 'tab:red',
    'flame_front': 'tab:blue',
    'zn': 'tab:green',
    'flame_front+zn': 'tab:purple',
}


def _chunc_pyrogen():
    pyro = load_pyrogen('mtv')
    pyro.heat_delivery_mode = 'radiation'
    return pyro


MOTORS = [
    dict(
        label='Chunc (machbusterNew)',
        slug='machbusterNew',
        ric=MOTORS_DIR / 'machbusterNew.ric',
        transport=MOTORS_DIR / 'machbusterNew.frozen.transport.yaml',
        pyrogen=_chunc_pyrogen,          # callable → fresh object per run
        topology='head_basket',
        t_max=3.0,
        experimental=CHUNC_EXPERIMENTAL,
        time_offset=CHUNC_EXPERIMENTAL.get('time_offset', 0.0),
    ),
    dict(
        label='Hasegawa A',
        slug='hasegawa_a',
        ric=MOTORS_DIR / 'hasegawa_a.ric',
        transport=MOTORS_DIR / 'hasegawa_a.frozen.transport.yaml',
        pyrogen='bpnv',
        topology='forward_plenum',
        t_max=3.0,
        experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL,
        time_offset=1.1,
    ),
]


def run_config(motor, flame_front, zn):
    pyro = motor['pyrogen']() if callable(motor['pyrogen']) else motor['pyrogen']
    transport = motor['transport']
    transport_path = str(transport) if transport.exists() else None
    result, _perf, _nz, _geo, _prop = run_from_ric(
        str(motor['ric']),
        transport_path=transport_path,
        pyrogen=pyro,
        injection_topology=motor['topology'],
        t_max=motor['t_max'],
        flame_front_enabled=flame_front,
        zn_enabled=zn,
        kappa_zn=1.0,
        **COMMON_KNOBS,
    )
    return result


def metrics(result):
    t = np.asarray(result['time'])
    P = np.asarray(result['P_head']) / 1e6
    P_peak = float(P.max())
    t_peak = float(t[int(P.argmax())])
    n = len(P)
    if n > 6:
        lo, hi = n // 3, 2 * n // 3
        P_plateau = float(np.median(P[lo:hi]))
    else:
        P_plateau = float(P.mean())
    spike_ratio = P_peak / P_plateau if P_plateau > 0.0 else float('nan')
    spread = ignition_spread_metrics(result)
    return dict(
        P_peak=P_peak, t_peak=t_peak, P_plateau=P_plateau,
        spike_ratio=spike_ratio,
        spread_10_90_ms=1e3 * spread.get('spread_10_90_s', float('nan')),
    )


def main():
    out = artifact_dir('chunc_ignition_2x2')
    print(f"Outputs -> {out}\n")
    rows = []

    for motor in MOTORS:
        print(f"\n=== {motor['label']} ===")
        fig, ax = plt.subplots(figsize=(12, 7))
        if motor['experimental'] is not None:
            exp = motor['experimental']
            ax.plot(np.asarray(exp['time']) + motor['time_offset'],
                    exp['pressure'], 'k-', linewidth=1.6, marker='o',
                    markersize=3, markevery=max(1, len(exp['time']) // 25),
                    label=exp.get('label', 'Experimental'), zorder=10)

        for label, ff, zn in CONFIGS:
            print(f"  {label} (flame_front={ff}, zn={zn}) ...")
            result = run_config(motor, ff, zn)
            m = metrics(result)
            print(f"    P_peak={m['P_peak']:.2f} MPa @ t={m['t_peak']:.3f}s, "
                  f"plateau={m['P_plateau']:.2f} MPa, "
                  f"spike/plateau={m['spike_ratio']:.2f}, "
                  f"ign spread={m['spread_10_90_ms']:.1f} ms")
            rows.append((motor['label'], label, m))
            ax.plot(result['time'], np.asarray(result['P_head']) / 1e6,
                    '-', linewidth=1.8, color=CONFIG_COLORS[label],
                    label=f"{label} (peak {m['P_peak']:.2f}, "
                          f"ratio {m['spike_ratio']:.2f})")

        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Head-End Pressure [MPa]', fontsize=12)
        ax.set_title(
            f"{motor['label']} — v0.7.4 ignition-transient configs\n"
            f"(baseline / Phase F flame-front / Phase Z Z-N / F+Z; "
            f"canon knobs, frozen transport)",
            fontsize=10,
        )
        ax.grid(alpha=0.3)
        ax.legend(loc='best', fontsize=9)
        png = out / f"{motor['slug']}_2x2.png"
        fig.savefig(png, dpi=130, bbox_inches='tight')
        plt.close(fig)
        print(f"  Plot: {png.name}")

    # Summary table
    lines = [
        "v0.7.4 Phase FZ — Ignition-Transient Config Sweep (Chunc + Hasegawa A)",
        "=" * 92,
        "",
        "Knobs: roughness 37.1um, kappa 0.45, T_ign 850K, k_solid default, "
        "cfl 0.3, frozen transport.",
        "Chunc: head_basket + mtv(radiation). Hasegawa A: forward_plenum + bpnv.",
        "",
        f"{'Motor':<22} {'Config':<16} {'P_peak [MPa]':>13} {'plateau':>9} "
        f"{'spike/plat':>11} {'ign spread [ms]':>16}",
        "-" * 92,
    ]
    for motor_label, cfg_label, m in rows:
        lines.append(
            f"{motor_label:<22} {cfg_label:<16} "
            f"{m['P_peak']:>13.2f} {m['P_plateau']:>9.2f} "
            f"{m['spike_ratio']:>11.2f} {m['spread_10_90_ms']:>16.1f}"
        )

    lines += [
        "",
        "Validation criterion: Chunc spike -> ~8.5 MPa (no overshoot) while the",
        "plateau and taildown stay put; Hasegawa A must not regress. Phase F should",
        "widen Chunc's ignition spread from <5 ms toward tens of ms; if the spike",
        "persists despite that, Root B (memoryless burn rate) dominates and Phase Z",
        "is the load-bearing fix.",
    ]
    summary_path = out / 'summary.txt'
    summary_path.write_text('\n'.join(lines))
    print(f"\nSummary: {summary_path}")
    print()
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
