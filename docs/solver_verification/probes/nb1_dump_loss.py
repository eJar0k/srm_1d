"""NB1 step 6 (report only): stagnation-pressure loss of a sudden expansion at the aft end.

Port exit (area A_e, Mach M_e, static P_e) dumps into the aft plenum / convergent entrance (area
A_c ~ the case ID).  1-D control volume from the exit plane to the plane just downstream of the
step (wall pressure on the annular step face = P_e, no friction, adiabatic):

    mass     rho_e u_e A_e = rho_c u_c A_c
    momentum P_e A_e + rho_e u_e^2 A_e + P_e (A_c - A_e) = P_c A_c + rho_c u_c^2 A_c
             -> P_c (1 + gamma M_c^2) = P_e (1 + gamma sigma M_e^2),     sigma = A_e/A_c
    energy   T0 conserved

Solved (subsonic root) for M_c; the loss is 1 - P0_c/P0_e.  The incompressible Borda-Carnot limit is
dP0 = (rho_e u_e^2 / 2)(1 - sigma)^2 -> dP0/P_e = (gamma/2) M_e^2 (1 - sigma)^2.  No free constants.

    python docs/solver_verification/probes/nb1_dump_loss.py chunc_isentropic.json ...
"""
import json
import math
import sys


def g(M, gam):
    return (1.0 + 0.5 * (gam - 1.0) * M * M) ** (gam / (gam - 1.0))


def dump_loss(M_e, sigma, gam):
    """Return (P0_c/P0_e, M_c, P_c/P_e) for a sudden expansion with area ratio sigma = A_e/A_c."""
    x = 0.5 * (gam - 1.0)
    Pi = 1.0 + gam * sigma * M_e * M_e                    # P_c (1 + gam M_c^2) / P_e
    rhs = sigma * M_e * math.sqrt(1.0 + x * M_e * M_e)    # P_c M_c sqrt(1 + x M_c^2) / P_e

    def f(M):
        return Pi * M * math.sqrt(1.0 + x * M * M) / (1.0 + gam * M * M) - rhs
    lo, hi = 0.0, M_e * sigma * 1.0 + 1e-12            # subsonic branch lies below M_e
    hi = 1.0
    # the function increases from 0 up to its sonic maximum: first root from the left
    M_prev, f_prev = 0.0, f(0.0)
    n = 4000
    for k in range(1, n + 1):
        M = k / n
        fm = f(M)
        if f_prev < 0.0 <= fm:
            lo, hi = M_prev, M
            break
        M_prev, f_prev = M, fm
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    M_c = 0.5 * (lo + hi)
    P_c_over_P_e = Pi / (1.0 + gam * M_c * M_c)
    return g(M_c, gam) * P_c_over_P_e / g(M_e, gam), M_c, P_c_over_P_e


def main():
    print(f"{'file':28s} {'M_e':>7s} {'A_e/A_c':>8s} {'A_e/A_t':>8s} {'P0/P iso':>9s} "
          f"{'dump dP0/P0':>12s} {'BC incompr':>11s}")
    for path in sys.argv[1:]:
        r = json.load(open(path))
        e = r['exit']
        gam, M_e = e['gamma'], e['M_exit']
        A_e = 0.25 * math.pi * e['D_port_exit_m'] ** 2
        A_c = 0.25 * math.pi * e['D_case_exit_m'] ** 2
        sigma = A_e / A_c
        ratio, M_c, _ = dump_loss(M_e, sigma, gam)
        bc_inc = 0.5 * gam * M_e ** 2 * (1.0 - sigma) ** 2 / g(M_e, gam)
        print(f"{Path_name(path):28s} {M_e:7.4f} {sigma:8.4f} {e['A_exit_over_A_throat']:8.3f} "
              f"{e['P0_over_P'] - 1.0:9.4%} {1.0 - ratio:12.4%} {bc_inc:11.4%}")


def Path_name(p):
    return p.replace('\\', '/').split('/')[-1]


if __name__ == '__main__':
    main()
