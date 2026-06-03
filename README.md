# srm_1d

A **1D transient finite-volume solid rocket motor internal-ballistics
simulator** with the Ma et al. (2020) erosive-burning model. A Numba-JIT
compiled PISO time loop (≈45–90k steps/s) resolves the axial pressure,
velocity, temperature, and per-cell grain regression through the full
firing — ignition transient, plateau, erosive lift, and tail-off.

## Features

- **PISO transient core** — pressure–velocity coupling, TDMA, adaptive CFL.
- **Ma 2020 erosive burning** — Haaland → Gnielinski → bisection closure.
- **Pyrogen ignition** — Goodman integral solid-heating subsolver, choked
  igniter plenum, and uncontained (`head_basket`/`aft_basket`) topologies.
- **N-species bore gas** with per-cell transport (frozen / effective).
- **openMotor integration** — loads `.ric` motors and registers as an
  openMotor solver plugin (`srm_1d.srm1d_plugin`).
- **FMM grain regression** via a bridge to a local openMotor checkout.

## Install

```bash
pip install -e .            # core
pip install -e ".[fmm,dev]" # + FMM grains (scikit-fmm) and pytest
```

Requires Python ≥ 3.10 (developed on 3.10.5). Core deps: numpy, scipy,
numba, pyyaml, matplotlib.

## Quick start

Run from the repo root:

```bash
python -m examples.hasegawa_motor_a   # runs motors/hasegawa_a.ric
python -m pytest tests/               # test suite
```

```python
from srm_1d import run_from_ric

result, perf, nozzle, geo, prop = run_from_ric(
    "motors/hasegawa_a.ric", pyrogen="bpnv"
)
print(perf["P_peak"], perf["total_impulse"])
```

## Layout

| Path | Role |
|------|------|
| `srm_1d/` | the importable package (ships): solver core, `tools/`, `pyrogens/` |
| `motors/` | motor data — `<motor>.ric` (transport embedded per-tab) |
| `examples/` | runnable studies (`python -m examples.<name>`) |
| `tests/` | pytest suite |
| `docs/` | design packages and development narrative |

See [`srm_1d/README.md`](srm_1d/README.md) for the full public API and
validated parameters, and [`CLAUDE.md`](CLAUDE.md) for an orientation map.

## License

MIT.
