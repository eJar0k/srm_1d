from srm_1d.openmotor_adapter import run_from_ric


result, perf, nozzle, geo, prop = run_from_ric(
    "srm_1d/motors/hasegawa_a.ric",
    roughness=20e-6,
    pyrogen="bpnv",
    T_ignition=850.0,
)
