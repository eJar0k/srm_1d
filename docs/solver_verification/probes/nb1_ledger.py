"""NB1 step 3: mass and energy ledgers with the new boundary (EN1 harness quantities).

Runs Chunc (committed .ric, canonical knobs, defaults of D1: total energy + conservative
momentum) for 0.1 s with each ``nozzle_bc_stagnation`` and reports the closure measures the
EN1 / T6 tests pin:

  * ``energy_residual``       per-step PISO ledger residual, relative to the peak thermal-source power
  * ``source_ledger_mismatch_W``  source ledger, relative to the peak thermal-source power
  * ``nozzle_ledger_mismatch_W``  PISO nozzle enthalpy flux vs the post-PISO reconstruction
  * ``nozzle_mass_mismatch_kg_s`` same for mass flow, relative to the peak nozzle mass flow
  * cumulative mass: produced - exited - stored  (summary ``mass_balance_error``)

    python docs/solver_verification/probes/nb1_ledger.py MODE OUT.json
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from srm_1d.openmotor_adapter import run_from_ric  # noqa: E402


def main():
    mode, out = sys.argv[1], sys.argv[2]
    res = run_from_ric(str(ROOT / 'motors' / 'machbusterNew.ric'), roughness=32e-6, kappa=0.44,
                       pyrogen='bpnv', pyrogen_mass=None, T_ignition=756.0, t_max=0.1,
                       P_cutoff=0.01e6, print_interval=0.0, verbose=False,
                       nozzle_bc_stagnation=mode)[0]
    a = lambda k: np.abs(np.asarray(res[k]))
    peak_src = float(a('thermal_source_power').max())
    peak_noz = float(a('nozzle_enthalpy_power').max())
    peak_md = float(a('nozzle_massflow_piso').max())
    row = dict(
        mode=mode, steps=int(res['summary']['steps']),
        energy_residual_over_peak_source=float(a('energy_residual').max() / peak_src),
        source_ledger_mismatch_over_peak_source=float(a('source_ledger_mismatch_W').max() / peak_src),
        nozzle_ledger_mismatch_over_peak_nozzle=float(a('nozzle_ledger_mismatch_W').max() / peak_noz),
        nozzle_mass_mismatch_over_peak_mdot=float(a('nozzle_mass_mismatch_kg_s').max() / peak_md),
        mass_balance_error=float(res['summary']['mass_balance_error']),
        peak_P_MPa=float(np.max(res['P_head'])) / 1e6)
    Path(out).write_text(json.dumps(row, indent=1))
    print(json.dumps(row, indent=1))


if __name__ == '__main__':
    main()
