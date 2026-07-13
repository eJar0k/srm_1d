"""Bisect why LHS calls collapse where the example doesn't.

Reference call (from hasegawa_motor_a.py) — known to run clean:
    run_from_ric('hasegawa_a.ric', roughness=37.1e-6, kappa=0.45,
                 pyrogen='bpnv', pyrogen_mass=None,
                 T_ignition=850.0, P_cutoff=0.05e6,
                 snapshot_interval=0.2, print_interval=0.2)

LHS call (collapsed):
    run_from_ric('hasegawa_a.ric', roughness=37.1e-6, kappa=0.45,
                 pyrogen=<Pyrogen obj with mode='demar'>,
                 transport_path='hasegawa_a.transport.yaml',
                 injection_topology='forward_plenum',
                 t_max=6.0, P_cutoff=0.05e6,
                 snapshot_interval=2.0, print_interval=20.0,
                 T_ignition=850, k_solid=0.30)

Bisecting: add LHS kwargs one at a time on top of the reference.
"""
from srm_1d.openmotor_adapter import run_from_ric, load_pyrogen


def report(label, result):
    s = result['summary']
    p_peak = s.get('P_peak', float('nan'))
    if p_peak is None: p_peak = float('nan')
    print(f"[{label}] "
          f"P_peak={p_peak/1e6:6.2f} MPa, "
          f"t_burn={s.get('t_burn', 0):6.4f}s, "
          f"term_code={s.get('termination_code', '?')}")


# Baseline: canonical knobs as hasegawa_motor_a.py runs them
print("\n=== Baseline canonical (hasegawa_motor_a.py style) ===")
r, *_ = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    roughness=37.1e-6, kappa=0.45,
    pyrogen='bpnv', pyrogen_mass=None,
    T_ignition=850.0, P_cutoff=0.05e6,
    snapshot_interval=0.2, print_interval=0.2,
    verbose=False,
)
report("baseline", r)

# +1: add explicit pyrogen object instead of string
print("\n=== +1: pyrogen as constructed object ===")
pyro = load_pyrogen('bpnv')
r, *_ = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    roughness=37.1e-6, kappa=0.45,
    pyrogen=pyro, pyrogen_mass=None,
    T_ignition=850.0, P_cutoff=0.05e6,
    snapshot_interval=0.2, print_interval=0.2,
    verbose=False,
)
report("pyrogen=obj", r)

# +2: add explicit injection_topology
print("\n=== +2: explicit injection_topology='forward_plenum' ===")
r, *_ = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    roughness=37.1e-6, kappa=0.45,
    pyrogen=pyro, pyrogen_mass=None,
    injection_topology='forward_plenum',
    T_ignition=850.0, P_cutoff=0.05e6,
    snapshot_interval=0.2, print_interval=0.2,
    verbose=False,
)
report("+topology", r)

# +3: add explicit transport_path
print("\n=== +3: explicit transport_path=effective ===")
r, *_ = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    transport_path='srm_1d/motors/hasegawa_a.transport.yaml',
    roughness=37.1e-6, kappa=0.45,
    pyrogen=pyro, pyrogen_mass=None,
    injection_topology='forward_plenum',
    T_ignition=850.0, P_cutoff=0.05e6,
    snapshot_interval=0.2, print_interval=0.2,
    verbose=False,
)
report("+transport_eff", r)

# +4: add explicit k_solid
print("\n=== +4: explicit k_solid=0.30 ===")
r, *_ = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    transport_path='srm_1d/motors/hasegawa_a.transport.yaml',
    roughness=37.1e-6, kappa=0.45,
    pyrogen=pyro, pyrogen_mass=None,
    injection_topology='forward_plenum',
    T_ignition=850.0, P_cutoff=0.05e6,
    snapshot_interval=0.2, print_interval=0.2,
    k_solid=0.30,
    verbose=False,
)
report("+k_solid=0.30", r)

# +5: switch to LHS's snapshot/print intervals
print("\n=== +5: snapshot_interval=2.0, print_interval=20.0 ===")
r, *_ = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    transport_path='srm_1d/motors/hasegawa_a.transport.yaml',
    roughness=37.1e-6, kappa=0.45,
    pyrogen=pyro, pyrogen_mass=None,
    injection_topology='forward_plenum',
    T_ignition=850.0, P_cutoff=0.05e6,
    snapshot_interval=2.0, print_interval=20.0,
    k_solid=0.30,
    verbose=False,
)
report("+slow_snap", r)

# +6: add explicit t_max=6.0
print("\n=== +6: t_max=6.0 (LHS uses this; example uses default) ===")
r, *_ = run_from_ric(
    'srm_1d/motors/hasegawa_a.ric',
    transport_path='srm_1d/motors/hasegawa_a.transport.yaml',
    roughness=37.1e-6, kappa=0.45,
    pyrogen=pyro, pyrogen_mass=None,
    injection_topology='forward_plenum',
    T_ignition=850.0, P_cutoff=0.05e6,
    snapshot_interval=2.0, print_interval=20.0,
    k_solid=0.30, t_max=6.0,
    verbose=False,
)
report("+tmax=6", r)
