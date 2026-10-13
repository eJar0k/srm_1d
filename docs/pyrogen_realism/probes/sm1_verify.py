# SM1 Part 2: verification of the three solid-heating models against exact solutions.
# Usage (repo root):  python docs/pyrogen_realism/probes/sm1_verify.py [--convergence]
# Prints markdown tables; redirect to docs/pyrogen_realism/results/sm1_verify.md.
import sys
import numpy as np
from srm_1d import solid_thermal as st
from srm_1d.tools import solid_verify as sv

k = 0.271
rho, cps = 1800.0, 1500.0           # sim: Cps = 1500 fixed; rho of the order of the BPNV/APCP tabs
al = k / (rho * cps)
Ti, Tg = 293.0, 2800.0
Y = st.build_solid_grid(al, 2.0, st.CONDUCTION_DY0, st.CONDUCTION_RATIO)
AF = st.CONDUCTION_AGE_FRACTION
MODELS = ('goodman', 'peretz', 'conduction')


def run(model, hh, TT, dt, y=Y, f=AF):
    if model == 'goodman':
        return sv.run_goodman(hh, TT, dt, Ti, al, k)
    if model == 'peretz':
        return sv.run_peretz(hh, TT, dt, Ti, al, k)
    return sv.run_conduction(hh, TT, dt, Ti, al, k, y, f)


print(f"k = {k} W/(m K), alpha = {al:.4e} m^2/s, T_i = {Ti} K, conduction grid: "
      f"dy0 = {st.CONDUCTION_DY0:.1e} m, ratio {st.CONDUCTION_RATIO}, M = {len(Y)} nodes "
      f"to {Y[-1]*1e3:.2f} mm, BE sub-step <= {AF} x heating age, active depth 7 sqrt(alpha t)\n")

# ---- 1. constant h, T_g: exact Robin
print("### 1. Constant h, T_g = 2800 K: exact Robin (erfcx) solution\n")
print("Error = max |T_s - exact| over the window T_exact <= 900 K [K]; dt90% = time-to-756 K error "
      "[% of exact t756].\n")
print("| h [W/m2K] | t756 exact [ms] | dt [us] | " + " | ".join(f"{m} err K / t756 %" for m in MODELS) + " |")
print("|---|---|---|---|---|---|")
for h in (1000., 5000., 20000.):
    for dt in (2e-6, 2e-5):
        tend = 0.05 if h > 500 else 0.5
        n = int(tend / dt)
        t = (np.arange(n) + 1) * dt
        hh = np.full(n, h)
        TT = np.full(n, Tg)
        ex = sv.exact_robin(t, h, Tg, Ti, al, k)
        t756 = sv.time_to(ex, dt, 756., Ti)
        cells = []
        for m in MODELS:
            r = run(m, hh, TT, dt)
            msk = ex <= 900.
            e = np.abs(r - ex)[msk].max()
            et = (sv.time_to(r, dt, 756., Ti) - t756) / t756 * 100
            cells.append(f"{e:.2f} / {et:+.2f}%")
        print(f"| {h:.0f} | {t756*1e3:.3f} | {dt*1e6:.0f} | " + " | ".join(cells) + " |")

# ---- 2. constant flux
print("\n### 2. Constant surface flux: exact T_s = T_i + 2 q sqrt(alpha t / pi) / k\n")
print("Realised as h = 1 W/m2K, T_drv = T_i + q/h (flux droop = (T_s - T_i) / (q/h) < 5e-5). "
      "q = 1e7 W/m2.\n")
q = 1.0e7
h1 = 1.0
dt = 2e-6
n = int(4e-3 / dt)
t = (np.arange(n) + 1) * dt
hh = np.full(n, h1)
TT = np.full(n, Ti + q / h1)
ex = sv.exact_flux(t, q, Ti, al, k)
t756 = sv.time_to(ex, dt, 756., Ti)
print("| t [ms] | exact T_s [K] | " + " | ".join(f"{m} err [K]" for m in MODELS) + " |")
print("|---|---|---|---|---|")
res = {m: run(m, hh, TT, dt) for m in MODELS}
for tq in (0.1e-3, 0.3e-3, 0.6e-3, 1.0e-3, 1.5e-3):
    i = int(round(tq / dt)) - 1
    print(f"| {tq*1e3:.1f} | {ex[i]:.1f} | " + " | ".join(f"{res[m][i]-ex[i]:+.2f}" for m in MODELS) + " |")
print(f"\nTime to 756 K: exact {t756*1e3:.4f} ms; error " +
      ", ".join(f"{m} {(sv.time_to(res[m], dt, 756., Ti)-t756)/t756*100:+.2f}%" for m in MODELS))

# ---- 3. step changes
print("\n### 3. Step changes in h (T_g = 2800 K const.), dt = 2 us\n")
print("Reference = Duhamel/Volterra flux-superposition solve (independent of the conduction grid, "
      "Richardson-extrapolated from dt = 2 us and 1 us); the converged conduction solve "
      "(dy0 = 1.25e-8, ratio 1.03, sub-step 0.003 x age) is quoted as a cross-check.\n")
YF = st.build_solid_grid(al, 2.0, 1.25e-8, 1.03)
cases = (("step DOWN  h 5400 -> 5400/14 at 0.8 ms (cell-9 like)", 5400., 1 / 14., 0.8e-3),
         ("step UP    h 500 -> 5000 (x10) at 1.0 ms", 500., 10., 1.0e-3))
dt = 2e-6
te = 0.01
n = int(te / dt)
tt = (np.arange(n) + 1) * dt
for name, h1, fac, t1 in cases:
    hf = lambda t, h1=h1, fac=fac, t1=t1: np.where(t <= t1, h1, h1 * fac)
    tv, v = sv.volterra_richardson(hf, lambda t: np.full_like(t, Tg), te, dt, Ti, al, k)
    hh = hf(tt)
    TT = np.full(n, Tg)
    fine = sv.run_conduction(hh, TT, dt, Ti, al, k, YF, 0.003125)
    res = {m: run(m, hh, TT, dt) for m in MODELS}
    print(f"**{name}**\n")
    print("| t [ms] | Volterra ref [K] | converged cond. [K] | goodman err [K] | peretz err [K] | conduction err [K] |")
    print("|---|---|---|---|---|---|")
    for tq in (0.5e-3, t1, t1 + 0.1e-3, t1 + 0.5e-3, t1 + 1e-3, t1 + 2e-3, t1 + 4e-3, 9e-3):
        i = int(round(tq / dt)) - 1
        print(f"| {tq*1e3:.2f} | {v[i]:.1f} | {fine[i]:.1f} | " +
              " | ".join(f"{res[m][i]-v[i]:+.1f}" for m in MODELS) + " |")
    after = tt > t1
    summary = []
    for m in MODELS:
        e = np.abs(res[m] - v)
        summary.append(f"{m} {e[after].max():.1f}")
    print(f"\nmax |T_s error| after the step [K]: " + ", ".join(summary) +
          f"; converged-conduction vs Volterra {np.abs(fine - v)[after].max():.2f}")
    tr = sv.time_to(v, dt, 756., Ti)
    if np.isfinite(tr):
        print("time to 756 K: Volterra %.3f ms; " % (tr * 1e3) +
              ", ".join(f"{m} {sv.time_to(res[m], dt, 756., Ti)*1e3:.3f} ms" for m in MODELS))
    else:
        mx = v.max()
        print(f"T_s never reaches 756 K (Volterra max {mx:.0f} K); "
              + ", ".join(f"{m} max {res[m].max():.0f} K" for m in MODELS))
    print()

if '--convergence' in sys.argv:
    print("### 4. Conduction-grid convergence study (T_s <= 900 K; dt = 2 us / 20 us)\n")
    print("max |T_s - exact Robin| [K] / time-to-756 error [%] for h = 1000, 5000, 20000; "
          "'step' = max |conduction - converged conduction| [K] on the step-down and step-up cases "
          "(T_s <= 900 K window).\n")
    cands = [(1e-7, 1.15, 0.05), (1e-7, 1.15, 0.0125), (5e-8, 1.1, 0.05), (5e-8, 1.1, 0.025),
             (5e-8, 1.1, 0.0125), (5e-8, 1.1, 0.00625), (2.5e-8, 1.05, 0.0125),
             (2.5e-8, 1.05, 0.00625), (1.25e-8, 1.03, 0.003125)]
    for dt in (2e-6, 2e-5):
        n = int(te / dt)
        tt2 = (np.arange(n) + 1) * dt
        refs = []
        for name, h1, fac, t1 in cases:
            hh = np.where(tt2 <= t1, h1, h1 * fac)
            TT = np.full(n, Tg)
            refs.append((hh, TT, sv.run_conduction(hh, TT, dt, Ti, al, k, YF, 0.003125)))
        print(f"dt = {dt*1e6:.0f} us\n")
        print("| dy0 [m] | ratio | sub-step frac | nodes | h=1000 | h=5000 | h=20000 | step-down | step-up |")
        print("|---|---|---|---|---|---|---|---|---|")
        for dy0, r, f in cands:
            y = st.build_solid_grid(al, 2.0, dy0, r)
            cells = []
            for h in (1000., 5000., 20000.):
                tend = 0.05
                nn = int(tend / dt)
                t = (np.arange(nn) + 1) * dt
                ex = sv.exact_robin(t, h, Tg, Ti, al, k)
                rr = sv.run_conduction(np.full(nn, h), np.full(nn, Tg), dt, Ti, al, k, y, f)
                e = np.abs(rr - ex)[ex <= 900.].max()
                t756 = sv.time_to(ex, dt, 756., Ti)
                et = (sv.time_to(rr, dt, 756., Ti) - t756) / t756 * 100
                cells.append(f"{e:.2f} K / {et:+.2f}%")
            for hh, TT, ref in refs:
                c = sv.run_conduction(hh, TT, dt, Ti, al, k, y, f)
                msk = ref <= 900.
                cells.append(f"{np.abs(c - ref)[msk].max():.2f} K")
            mark = " (chosen)" if (dy0, r, f) == (5e-8, 1.1, 0.0125) else ""
            print(f"| {dy0:.2e} | {r} | {f} | {len(y)}{mark} | " + " | ".join(cells) + " |")
        print()
