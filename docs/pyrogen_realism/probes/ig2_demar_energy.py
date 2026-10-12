"""IG2 Part 3 support: energy the DeMar surface flux delivers vs the pyrogen gas enthalpy, in the sim.
Run from the worktree root with PYTHONPATH=$PWD."""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
import numpy as np
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen
out = {}
for topo, g in (('forward_plenum', 4.5), ('forward_plenum', 0.5), ('head_basket', 4.5)):
    pyro = load_pyrogen('bpnv')
    res = run_from_ric(str(ROOT / 'motors/machbusterNew.ric'), pyrogen=pyro, pyrogen_mass=g * 1e-3,
                       injection_topology=topo, roughness=32e-6, kappa=0.44, T_ignition=756.0, t_max=0.1,
                       print_interval=0.0, verbose=False)[0]
    t = np.asarray(res['time']); q = np.asarray(res['pyrogen_surface_heat_power'])
    E_demar = float(np.trapz(q, t))
    s = res['summary']
    m = float(s['pyrogen_mass_burned'])
    gmf = pyro.gas_mass_fraction or 1.0
    Cp_p = 1385.0  # approx; real value is per-species in the sim
    E_gas = m * gmf * pyro.T_flame * 0 + m * gmf * (pyro.gamma / (pyro.gamma - 1.0)) * (8314.462618 / (pyro.M * 1000.0)) * (pyro.T_flame - 293.0)
    out[f'{topo}_{g}g'] = dict(E_demar_J=E_demar, E_pyro_gas_enthalpy_J=E_gas, ratio=E_demar / E_gas,
                               impetus_energy_J=g * 565.0)
    print(topo, g, out[f'{topo}_{g}g'])
json.dump(out, open(ROOT / 'docs/pyrogen_realism/results/ig2_demar_energy.json', 'w'), indent=1)
