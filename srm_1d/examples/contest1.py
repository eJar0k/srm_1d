from srm_1d import run_simulation
from srm_1d.grain_geometry import make_conical_grain
from srm_1d.propellant import make_hasegawa_propellant_1
from srm_1d.nozzle import Nozzle

geo = make_conical_grain(0.030, 0.050, 0.080, 0.500, N_cells=100)
nozzle = Nozzle(D_throat=0.020, D_exit=0.035)
result = run_simulation(geo, make_hasegawa_propellant_1(), nozzle, roughness=20e-6)