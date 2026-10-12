"""IG2 Part 1: decompose the ignition h at the first-lit cells (IG1 cases) into the
Gnielinski-chain terms, relative to a plain developed-flow smooth-pipe value, and compare with
Peretz 1974 Eq. III-22/23 and Keller 1966 Eq. 11.  Reads IG1's snapshot npz (C:/tmp/sv_ig1/npz).
Run from the worktree root with PYTHONPATH=$PWD."""
import sys, json
import numpy as np
from srm_1d.burn_rate import haaland_friction, gnielinski_nusselt

# Chunc (machbusterNew.ric) propellant gas, constants used by the sim (scalar mu/k/Cp)
MU, K, CP = 9.125e-5, 0.3654, 2050.0
PR = MU * CP / K
KAPPA, ROUGH = 0.44, 32e-6
CAL = 41840.0
D_THROAT = None  # filled below from the ric

def smooth_f(Re):                      # Petukhov smooth-tube friction factor
    return (0.790 * np.log(Re) - 1.64) ** -2

def nu_gn(Re, Pr, f, entrance, K):    # Gnielinski turbulent form (Ma Eq. 9), pieces switchable
    f8 = f / 8.0
    return f8 * (Re - 1000.0) * Pr / (1 + 12.7 * np.sqrt(f8) * (Pr ** (2 / 3) - 1)) * entrance * K

def peretz_h(P, u, Tg, Ts, R, d, x, Ap_At, gamma=1.2, Cp=CP, a6=2.56e-2):
    Taf = 0.5 * (Tg + Ts)
    W = 8314.462 / R                     # g/mol
    mu = 0.87e-7 * W ** 0.5 * Taf ** 0.65   # Pa s  (Eq. III-24, 0.87e-6 poise)
    rho = P / (R * Taf)
    Re = rho * abs(u) * d / mu
    Pr = gamma / (1.77 * gamma - 0.45)
    k = Cp * mu / Pr
    xr = max(x, d)
    Nu = a6 * (xr / d) ** -0.1 * Re ** 0.8 * Pr ** 0.4 * Ap_At ** 0.4
    return max(Nu, 3.66) * k / d, Re, Taf, mu

def keller_h(Tg, G):
    return 1.435e-4 * Tg ** 0.3 * (G * 0.1) ** 0.905 * CAL

def main():
    import srm_1d.openmotor_adapter as a, warnings
    warnings.simplefilter('ignore')
    m = a.load_ric('motors/machbusterNew.ric')
    nz = a.ric_to_sim_args(m, pyrogen='bpnv')['nozzle']
    At = np.pi / 4 * nz.D_throat ** 2
    out = {}
    for tag in ['th_hb_0.01', 'th_hb_0.001']:
        z = np.load(f'C:/tmp/sv_ig1/npz/{tag}.npz'); ign = z['ignition_time']; isg = z['snap_is_grain'][0] > .5
        idx = np.flatnonzero(isg); order = idx[np.argsort(ign[idx])][:8]
        print(f'\n== {tag}: first 8 lit cells (last pre-ignition snapshot); At={At*1e6:.1f} mm2')
        rows = []
        for c in order:
            ks = np.flatnonzero(z['snap_probe_h_c'][:, c] > 0); k = ks[-1]
            Tg = z['snap_T'][k, c]; Ts = z['snap_T_surf'][k, c]; u = z['snap_u'][k, c]
            rho = z['snap_rho'][k, c]; P = z['snap_P'][k, c]; d = z['snap_D_port'][k, c]
            x = z['snap_x'][c]; h_sim = z['snap_probe_h_c'][k, c]
            G = rho * abs(u); R = P / (rho * Tg)
            Re = rho * abs(u) * d / MU
            L = max(x, d)
            f_r = haaland_friction(Re, ROUGH, d); f_s = smooth_f(Re)
            ent = 1 + (d / L) ** (2 / 3)
            res = {}
            for lab, Tsx in (('last', Ts), ('cold293', 293.0)):
                K_ = (Tg / Tsx) ** KAPPA
                nu_ref = nu_gn(Re, PR, f_s, 1.0, 1.0)
                nu_db = 0.023 * Re ** 0.8 * PR ** 0.4
                res[lab] = dict(K=K_, ratio_Tg_Ts=Tg / Tsx)
            nu_ref = nu_gn(Re, PR, f_s, 1.0, 1.0)
            nu_db = 0.023 * Re ** 0.8 * PR ** 0.4
            fr = (f_r / f_s)   # roughness enters via f/8 (denominator effect small): compute exact factor
            nu_rough = nu_gn(Re, PR, f_r, 1.0, 1.0) / nu_ref
            Nu_sim = h_sim * d / K
            hp, Rep, Taf, mup = peretz_h(P, u, Tg, Ts, R, d, x, (np.pi / 4 * d * d) / At)
            hk = keller_h(Tg, G)
            rows.append(dict(cell=int(c), Tg=Tg, Ts=Ts, u=u, G=G, d=d, x=x, Re=Re, f_haal=f_r, f_smooth=f_s,
                             nu_db=nu_db, nu_ref=nu_ref, rough=nu_rough, ent=ent, K_last=(Tg / Ts) ** KAPPA,
                             K_cold=(Tg / 293.) ** KAPPA, Nu_sim=Nu_sim, tot_vs_ref=Nu_sim / nu_ref,
                             h_sim=h_sim, h_peretz=hp, Re_peretz=Rep, h_keller=hk,
                             Ap_At=(np.pi / 4 * d * d) / At))
            r = rows[-1]
            print(f"cell {c}: Tg {Tg:.0f} Ts {Ts:.0f} G {G:.0f} d {d*1e3:.1f}mm x {x*1e3:.0f}mm Re {Re:.0f} | "
                  f"Nu_DB {nu_db:.0f} Nu_ref(smooth Gn) {nu_ref:.0f} | rough x{nu_rough:.3f} ent x{ent:.3f} "
                  f"K(last) x{r['K_last']:.3f} K(cold) x{r['K_cold']:.3f} | Nu_sim {Nu_sim:.0f} = {Nu_sim/nu_ref:.3f} x ref "
                  f"(model product {nu_rough*ent*r['K_last']:.3f}) | h sim {h_sim:.0f} peretz {hp:.0f} keller {hk:.0f} "
                  f"(sim/peretz {h_sim/hp:.2f}, sim/keller {h_sim/hk:.2f}, Ap/At {r['Ap_At']:.2f})")
        out[tag] = rows
    json.dump(out, open('docs/pyrogen_realism/results/ig2_decomp.json', 'w'), indent=1, default=float)

if __name__ == '__main__':
    main()
