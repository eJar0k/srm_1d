"""Baseline bit-identity gate (B0) for the solver-verification tasks.

Re-runs the three canonical ``ma2020`` / ``erosive_mp_cap=0`` rows of
``docs/v0_7_4/probes/mp_cap_crossmotor.py`` (Chunc, BALLSstick, Hasegawa A;
committed .ric files, canonical knobs, library BPNV igniter) and compares
them to the recorded values (see the BASELINE comment: re-recorded by T10) at FULL float precision plus exact step counts.

A task whose spec says "default behaviour unchanged" passes B0 only if every
float matches within ``--rtol`` (default 1e-5) and every step count within
``--steps-rtol`` (default 5e-4 relative). Why not bitwise: with
``fastmath=True`` the Numba machine code depends on the compile state (which
functions were compiled in-process vs loaded from the on-disk cache), so
IDENTICAL source gives slightly different results between a warm cache and a
clean recompile -- deterministically per state. History of the measured
noise floor:
  * before T10: 1.5e-16..1.3e-15 relative, exact step counts (2026-10-05);
  * after T10: up to 1.6e-10 (Chunc plateau), exact step counts (2026-10-06);
  * after T11 (energy-coupled iteration, default K=1): up to 2.4e-7 (Chunc
    peak) and step counts differing by up to 4.9e-5 relative (Hasegawa A,
    115 of 2.3M steps) (2026-10-11). The iteration's convergence test
    (max|dT|/T < 1e-6) decides whether the extra pass runs; last-digit
    differences flip that decision on marginal steps, and the adaptive time
    step carries the flip into the trajectory. A warm cache reproduces the
    recorded values exactly.
  * after IG3 (Peretz ignition h, 2026-10-12): a fresh recompile moved Chunc
    t_peak by 4.6e-5 relative (1.7 us, of the order of one time step) while every
    pressure stayed within 1e-5 and two warm runs reproduced the recorded
    values exactly. t_peak is the time of the largest SAMPLE (an argmax over
    steps), so it moves in whole steps like n_steps: it is compared with
    ``--steps-rtol``, not ``--rtol``.
1e-5 / 5e-4 keep ~40x / ~10x margin over that floor and stay far below the
smallest real perturbation seen (T10: -0.5% peak; T11: -11%). Anything beyond
them means the change perturbed the default trajectory -- report it, do not
loosen the tolerance to make it pass. Use ``--rtol 1e-9 --steps-rtol 0`` for a
same-compile-state A/B.

This script puts ITS OWN checkout first on ``sys.path`` so that, run from a
git worktree, it imports the worktree's ``srm_1d`` (not the editable install
of the main checkout). It prints the resolved package path as provenance.

Usage (from the checkout root, pyenv 3.10.5 interpreter):
    python docs/solver_verification/probes/baseline_gate.py
    python docs/solver_verification/probes/baseline_gate.py --only chunc
Exit code 0 = all identical, 1 = mismatch.
"""
import argparse
import hashlib
import importlib.util
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import srm_1d  # noqa: E402  (after the sys.path fix on purpose)

_spec = importlib.util.spec_from_file_location(
    "mp_cap_crossmotor", ROOT / "docs" / "v0_7_4" / "probes" / "mp_cap_crossmotor.py")
_probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_probe)

# RE-RECORDED 2026-11-11 by NB1 (user decision 2026-11-10 "Prefer a physical fix, so taking the isentropic relations to
# find actual stagnation properties at that choked boundary"): run_simulation(nozzle_bc_stagnation='isentropic') is the
# default. The last cell's (P, T) are the STATIC port-exit state and the throat flow is computed from the true stagnation
# state P0 = P g(M_e) (area-Mach root for choked flow), so for a given gas generation the chamber static pressure settles
# LOWER by about g(M_e) - 1 at the plateau (Chunc M_e 0.157: +1.5 % P0; BALLSstick 0.128: +1.0 %; Hasegawa A 0.113: +0.8 %)
# and by more where the port Mach is higher (the startup). Moves: Chunc peak -6.5 % (9.826), plateau -1.5 % (7.917), peak/plateau
# 1.308 -> 1.241; BALLSstick peak -4.2 %, trough -2.3 %, mid -1.8 %; Hasegawa A late peak -0.9 % (5.961 at 2.53 s). The auto
# charges are unchanged (Chunc 0.815 g, BALLSstick 1.982 g, Hasegawa A 7.861 g). The 'static' boundary reproduces the AS1 values
# below bit for bit (nozzle_bc_stagnation='static'). Recorded from a warm-cache run after a clean __pycache__; reproduced by a
# second run of the A/B probe (docs/solver_verification/probes/nb1_motor_ab.py). The AS1 values (the previous BASELINE) were:
#   chunc: peak 10.510377787816175, t_peak 0.25905381442784203,
#          plateau 8.034819676399866, n_steps 543997
#   balls: peak 8.43291738963611, trough 4.43144847227655,
#          mid 4.925227950671004, n_steps 844368
#   haseg: peak 6.017191292703041, t_peak 2.499738880086473,
#          n_steps 2334673
#
# RE-RECORDED 2026-11-10 by AS1 (prompt-start auto charge; user decision): the auto-sized charge now has to meet
# "at least 95 % of the grain cells ignited within 0.5 s of igniter firing" (attempts that miss it are failed at the
# window end and retried at x1.25; the IG4 self-sustaining test alone had accepted hang-fires). The DeMar 300 psi
# initial guess already meets it for Chunc (0.815 g, 99 % lit by 0.235 s) and Hasegawa A (7.86 g, 99 % by 0.11 s), so
# those two rows are UNCHANGED. BALLSstick's 1.268 g (G1) lit its grain over ~1 s (t90 971 ms, a hang-fire); it now
# needs 3 attempts (1.268 -> 1.585 -> 1.982 g, DeMar-equivalent 469 psi): first light 1.7 ms, t90 17 ms, a real startup
# (peak 8.433 MPa at 29 ms, trough 4.431, mid 4.925; pk/trough 1.90, pk/mid 1.71). The previous (G1) BALLSstick row was:
#   balls: peak 0.2945771944198506, trough 0.07628399187690121,
#          mid 6.7042510894841545, n_steps 855216
#
# (earlier) RE-RECORDED 2026-11-10 by G1 (user decisions 2026-10-31 "Let's attempt A ..."): the head and aft
# gas-buffer cells (MIN_BUFFER_CELLS = 3, full case diameter, 3 dx long) are gone. Cell 0 is the first
# grain cell against the closure; a head gap (none for these three motors: all grains are inhibited on the
# head face) and the nozzle convergent are lumped into the end cells' gas storage (aft_end='lump'). Every
# startup shifts: the 24 mm (N = 100) head pocket is gone, and the DeMar auto charge uses the geometry's own
# free volume (no fictitious buffer volume), so the initial auto charge falls 25 % (Chunc 1.07 -> 0.81 g).
# Chunc: auto charge 1.34 -> 0.81 g, peak 9.955 -> 10.511 MPa at 0.259 s (the progressive late hump of D1 is
# now earlier), cell 0 never lit. BALLSstick: auto charge 3.09 -> 1.27 g and the grain lights over ~1 s (t90
# 971 ms), so this row's "peak" (largest P in the first 0.3 s, 0.295 MPa) and "trough" (0.076 MPa) are the
# slow-ignition pressure, NOT a startup spike/trough: the auto-size criterion (grain out-produces the igniter)
# accepts a charge that is too small for a prompt start once the free volume is physical. Flagged to the
# user in worklog/G1.md. Hasegawa A: 9.58 -> 7.86 g, late peak 6.109 -> 6.017 MPa (t 2.64 -> 2.50 s).
# Recorded on top of QF1 (the first G1 recording, before QF1's KE term, was Chunc 10.510546 / BALLS mid
# 6.704349 / Hasegawa A 6.017199: a 1.5e-5 shift, as QF1 predicted) after a clean __pycache__.
# The QF1 values (the previous BASELINE) were:
#   chunc: peak 9.954942701677405, t_peak 0.4882626717511268,
#          plateau 8.10882300869407, n_steps 541753
#   balls: peak 7.56221990287067, trough 4.341708721229631,
#          mid 4.8594489133399525, n_steps 844548
#   haseg: peak 6.108548290011714, t_peak 2.6393358955391624,
#          n_steps 2343878
# and the D1 values before them:
#   chunc: peak 9.95509537258403, t_peak 0.48826256590025513,
#          plateau 8.108865458753694, n_steps 541753
#   balls: peak 7.56231250893722, trough 4.3417299619930585,
#          mid 4.8594671206858635, n_steps 844549
#   haseg: peak 6.10855599581484, t_peak 2.6393309440372743,
#          n_steps 2343883
#
# (earlier) RE-RECORDED 2026-10-31 by D1 (user decision "keep total energy opt-in, build the next piece",
# then "approved as DEFAULTS"): run_simulation now defaults to energy_formulation='total',
# momentum_formulation='conservative', ignition_driver_temperature='recovery',
# head_end_jet_h='max' (docs/solver_verification/worklog/{EN1,RT1}.md,
# docs/pyrogen_realism/HE1_REPORT.md). The jet lights the head-end dead zones, so the IG3 auto-size
# fallback now stops at smaller charges (Chunc 2.10 -> 1.34 g, BALLSstick 3.86 -> 3.09 g,
# Hasegawa A 14.98 -> 9.58 g); Chunc's peak is a late progressive hump (t_peak 0.488 s, t90 455 ms).
# Recorded from a fresh-compile run (reproduced by a warm run). The conduction-default values
# (the previous BASELINE; reproduce them with energy_formulation='sensible',
# momentum_formulation='legacy', ignition_driver_temperature='static', head_end_jet_h=False) were:
#   chunc: peak 9.744849573359163, t_peak 0.035353187954651435,
#          plateau 7.99119980103874, n_steps 534236
#   balls: peak 5.8149312202096075, trough 3.6073917578343226,
#          mid 3.97139852072466, n_steps 852772
#   haseg: peak 5.363255636176837, t_peak 2.6881129738957537, n_steps 2347930
#
# (earlier) RE-RECORDED 2026-10-21 by the conduction default flip (user decision 2026-10-13):
# run_simulation(solid_heating_model='conduction') is the default (was 'goodman'): the
# unignited surface is advanced by exact 1-D transient conduction instead of Goodman's
# algebraic cubic-profile integral (T_ignition / k_solid were fit with Goodman). Ignition is
# later and the fallback needs smaller charges: auto charges Chunc 2.623 -> 2.098 g (4
# attempts; 1 head-end cell never lit), BALLSstick 4.829 -> 3.863 g (4; 8 never lit),
# Hasegawa A 14.98 g (3; 8 never lit, unchanged). Chunc peak 8.698 -> 9.745 MPa (pk/pl 1.22),
# BALLSstick peak 7.193 -> 5.815 (pk/trough 2.04 -> 1.61), Hasegawa A late peak 5.362 -> 5.363.
# Values recorded from a cold-compile run (docs/pyrogen_realism/results/cond_B0_*). The previous
# (Goodman, IG4+MEM1) BASELINE was:
#   chunc: peak 8.698445703737404, t_peak 0.024529374887138196,
#          plateau 8.06619314574538, n_steps 533513
#   balls: peak 7.193312409822636, trough 3.522240937005205,
#          mid 3.8878973484404313, n_steps 850245
#   haseg: peak 5.362271099044967, t_peak 2.6818221964560043, n_steps 2347427
#
# (earlier) RE-RECORDED 2026-10-12 by pyrogen-realism IG3 (stage 2, user decision "self-sustaining +
# warn"): the fallback's success test is no longer >= 99 % of the grain lit; it is "the grain
# lit and keeps burning after the igniter is spent", and grain cells that never light are only
# WARNED about (convective-only ignition model, not an igniter-size problem). Final auto
# charges: Chunc 2.623 g (unchanged, all cells lit), BALLSstick 4.829 g (was 7.546; 9 head-end
# cells never lit), Hasegawa A 14.98 g (was 18.72; 8 head-end cells never lit). Chunc rows are
# unchanged; BALLSstick peak 11.18 -> 7.193 (pk/trough 2.64 -> 2.04); Hasegawa A's B0 global
# peak is now the LATE one, 5.362 MPa at t = 2.68 s (the 23 ms ignition spike of the 18.72 g
# charge is gone). Values recorded from a clean-compile run (see results/ig3s2b_*). The
# previous (99 %-lit) BASELINE was:
#   chunc: peak 8.698445703737404, t_peak 0.024529374887138196,
#          plateau 8.06619314574538, n_steps 533513
#   balls: peak 11.182269177449877, trough 4.235361601649791,
#          mid 4.770663098447565, n_steps 842363
#   haseg: peak 6.195206225864981, t_peak 0.02305180400111348, n_steps 2335529
#
# (earlier) RE-RECORDED 2026-10-12 by pyrogen-realism IG3 (stage 2), on top of OB1: the DeMar surface
# flux is OFF by default (library pyrogens heat_delivery_mode 'none') and the auto igniter
# charge is DeMar-sized (V_free * 300 psi / impetus_W) with the auto-upsize fallback
# (x1.25 until >= 99 % of the grain lights). The canonical rows use that auto path, so each
# motor now runs a different, larger charge than the Sutton sizes: Chunc 0.904 -> 2.623 g
# (5 attempts), BALLSstick 1.386 -> 7.546 g (7 attempts), Hasegawa A 4.184 -> 18.72 g
# (4 attempts). Chunc peak 9.977 -> 8.698 MPa (peak/plateau 1.078, exp 1.027), plateau
# -0.9 %; BALLSstick peak 7.061 -> 11.18 (pk/trough 2.64, exp ~1.5), trough -2.5 %; Hasegawa A
# global peak 6.107 (t = 2.41 s) -> 6.195 (an ignition peak at t = 23 ms; its late 2.4 s peak
# is ~6.1). Values recorded from a clean-compile run and reproduced EXACTLY by a warm run
# (docs/pyrogen_realism/results/ig3s2_B0_rebased_*.txt). The OB1 values (the previous BASELINE) were:
#   chunc: peak 9.977396119842126, t_peak 0.03519436869536591,
#          plateau 8.141928366040833, n_steps 532348
#   balls: peak 7.060894562817688, trough 4.343382326291922,
#          mid 4.855564928020005, n_steps 843680
#   haseg: peak 6.107263748109281, t_peak 2.4109782579175234, n_steps 2336318
#
# (earlier) RE-RECORDED 2026-10-12 by pyrogen-realism OB1: (i) the library pellet pyrogens
# (bpnv, mtv) now use the MIL-P-46994B Type I-A pellet geometry (3.175 mm dia,
# L/D 1.504; was 3.2 mm, L/D 1.0) and (ii) the forward plenum's '0d' burn-area
# law is the exact right-cylinder form function (``_cyl_area_fraction``) instead of
# A0 (m/m0)^(2/3) (identical for L/D = 1), and the sphere initial-area branch was
# removed (a single cylinder formula). The canonical motors use the Sutton-sized
# forward-plenum BPNV igniter, so they move: Chunc peak +0.27 %, t_peak -4.2 %
# (36.7 -> 35.2 ms), plateau unchanged; BALLSstick peak +0.58 %; Hasegawa A unchanged
# (7e-7). Values below were recorded from a WARM-cache run (compile-state noise:
# a clean-compile run differs at 1e-7..4e-7 in the peaks and a few steps; two warm
# runs reproduce these EXACTLY). The IG3 (stage 1) values (the previous BASELINE)
# were:
#   chunc: peak 9.950855636512404, t_peak 0.036735041862973064,
#          plateau 8.141968067255977, n_steps 532381
#   balls: peak 7.020318627616664, trough 4.343390826697384,
#          mid 4.855316005870968, n_steps 843726
#   haseg: peak 6.107268310023869, t_peak 2.411042627819566, n_steps 2336271
#
# (earlier) RE-RECORDED 2026-10-11 by pyrogen-realism IG3 (stage 1): the default ignition
# heat-transfer closure became 'peretz' (Peretz et al. 1973, Eq. III-22, local
# port area) instead of the Ma/Gnielinski chain (ignition_heat_transfer='gnielinski').
# Only the first-light delay changes (Chunc 1.72 -> 5.81 ms, BALLSstick 1.75 -> 7.22,
# Hasegawa A 2.77 -> 7.69); the plateau is unchanged and the startup peak moves
# -1.3 % (Chunc) / -3.3 % (BALLSstick) / ~0 (Hasegawa A). Values below are the IG3
# results (docs/pyrogen_realism/results/ig3_B0_*.txt). The T11 values (the previous
# BASELINE; reproduce them with ignition_heat_transfer='gnielinski') were:
#   chunc: peak 10.08270103870215, t_peak 0.02791011521855539,
#          plateau 8.146643086031123, n_steps 531701
#   balls: peak 7.257654400712118, trough 4.343098504353659,
#          mid 4.868134347450266, n_steps 843077
#   haseg: peak 6.1072529740207635, t_peak 2.405585530808029,
#          n_steps 2336674
#
# (earlier) RE-RECORDED 2026-10-11 by solver-verification T11: the PISO second
# corrector now includes the storage the first corrector took
# (piso_corr2_storage=True; without it every transient pressure increment
# was applied ~twice) and the step iterates energy/pressure once
# (piso_energy_coupling_iters=1). Values below are the T11 results
# (C:/tmp/sv_t11/core_*.jsonl, worklog/T11.md). Chunc t_peak moving
# 16.3 -> 27.9 ms tripped T11's 15% startup-metric STOP rule; the reviewer
# ruled it expected (P5 measured the same shift before adoption: the
# doubled increments had made the startup pressurise too fast).
# Pre-T11 (T10, 2026-10-05) values:
#   chunc: peak 11.28243546264288, t_peak 0.016286760473606336,
#          plateau 8.150933877646011, n_steps 531883
#   balls: peak 8.331331016878805, trough 4.343222420159596,
#          mid 4.879629094113075, n_steps 843274
#   haseg: peak 6.108574454139173, t_peak 2.391606230911318, n_steps 2336140
# Pre-T10 (2026-08-08, erosive_model_registry WORKLOG T2) values:
#   chunc: peak 11.472831280605046, t_peak 0.016454109339050758,
#          plateau 8.222465813350324, n_steps 531651
#   balls: peak 8.42348927652058, trough 4.366702806070186,
#          mid 4.913378178011638, n_steps 842922
#   haseg: peak 6.143912554793205, t_peak 2.3816074267208034, n_steps 2335770
# RE-RECORDED 2026-11-10 by QF1: KE term in the total-energy regression-fill debit; moved peaks ~1.5e-5.
# Recorded from a cold-compile run (reproduced by a warm run). The D1 values (the previous BASELINE) were:
#   chunc: peak 9.95509537258403, t_peak 0.48826256590025513,
#          plateau 8.108865458753694, n_steps 541753
#   balls: peak 7.56231250893722, trough 4.3417299619930585,
#          mid 4.8594671206858635, n_steps 844549
#   haseg: peak 6.10855599581484, t_peak 2.6393309440372743, n_steps 2343883
BASELINE = {
    'chunc': dict(peak=9.826276661918815, t_peak=0.25915098317205126,
                  plateau=7.917330440691093, n_steps=545969),
    'balls': dict(peak=8.081732502821485, trough=4.331017370696415,
                  mid=4.838013889209436, n_steps=847603),
    'haseg': dict(peak=5.960864370215009, t_peak=2.528567952802282,
                  n_steps=2341103),
}
RIC_FILE = {'chunc': 'machbusterNew.ric', 'balls': 'BALLSstick.ric',
            'haseg': 'hasegawa_a.ric'}
ROWS = {'chunc': _probe.chunc_row, 'balls': _probe.balls_row,
        'haseg': _probe.haseg_row}


def _print_provenance(names):
    """Output-only provenance (P6): git HEAD + dirty flag, library versions,
    and the sha256 of each .ric file this gate will actually run (resolved
    by the probe's own committed-copy logic). Never affects pass/fail."""
    import numpy
    import numba
    def _git(*a):
        try:
            return subprocess.run(['git', *a], cwd=str(ROOT), capture_output=True,
                                  text=True, check=True).stdout.strip()
        except Exception as exc:  # provenance must never break the gate
            return f"<unavailable: {exc}>"
    head = _git('rev-parse', 'HEAD')
    dirty = bool(_git('status', '--porcelain', '--untracked-files=no'))
    print(f"[provenance] git HEAD {head}  dirty={dirty}")
    print(f"[provenance] numpy {numpy.__version__}  numba {numba.__version__}  "
          f"python {sys.version.split()[0]}")
    for name in names:
        try:
            path = _probe._committed_ric_path(RIC_FILE[name])
            digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
            print(f"[provenance] {RIC_FILE[name]} sha256 {digest}  ({path})")
        except Exception as exc:
            print(f"[provenance] {RIC_FILE[name]} sha256 <unavailable: {exc}>")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default='chunc,balls,haseg',
                    help='comma-separated subset of chunc,balls,haseg')
    ap.add_argument('--rtol', type=float, default=1e-5,
                    help='float tolerance (relative)')
    ap.add_argument('--steps-rtol', type=float, default=5e-4,
                    help='relative tolerance on n_steps (0 = exact)')
    args = ap.parse_args()

    print(f"[provenance] srm_1d imported from {Path(srm_1d.__file__).parent}")
    _print_provenance([n.strip() for n in args.only.split(',')])
    ok = True
    for name in args.only.split(','):
        name = name.strip()
        t0 = time.time()
        row = ROWS[name]('ma2020', 0.0)
        wall = time.time() - t0
        print(f"\n[{name}] wall {wall:.1f} s, termination={row['termination']!r}")
        for key, ref in BASELINE[name].items():
            got = row[key]
            if key == 'n_steps':
                rel = abs(int(got) - int(ref)) / ref
                same = rel <= args.steps_rtol
            elif key == 't_peak':
                # argmax over time steps: discrete, so the step tolerance
                rel = abs(got - ref) / abs(ref)
                same = rel <= args.steps_rtol
            else:
                rel = abs(got - ref) / abs(ref)
                same = rel <= args.rtol
            ok &= same
            print(f"  {key:8s} got {got!r:>24}  ref {ref!r:>24}  "
                  f"rel {rel:.3e}  {'OK' if same else 'MISMATCH'}")
    print("\nB0 RESULT:",
          f"PASS (floats within rtol {args.rtol:g}, n_steps and t_peak within "
          f"{args.steps_rtol:g})" if ok
          else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
