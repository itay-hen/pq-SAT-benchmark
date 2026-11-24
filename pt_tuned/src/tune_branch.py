#!/usr/bin/env python3
"""Optimize the PT temperature ladder separately for every d, for one branch.

usage: tune_branch.py factor|cnf dmin dmax

Per d (training instances only, from train.tsv):
  1. single-replica scan on a geometric T grid: <E>, var(E), tau_int.
     T_max = smallest grid T with median tau_int <= 1 sweep.
  2. provisional M0 = 1 + (integral of sigma_E d beta over [beta_max, beta_min]) / 1.47
     (spacing that gives ~30% swap acceptance for Gaussian energy histograms).
  3. T_min candidates:
       factor: c * T*, c in {0.6, 0.8, 1.0, 1.25}, T* = smallest T(P0 = 0.5) over
               the training instances, from the exact density of states.
       cnf:    a fixed grid; adequacy judged from the measured ground-state
               occupation p0 at T_min (meaningful when round trips occur).
     For each: feedback-optimized ladder, M grown from M0 by factors of 1.4
     until every adjacent swap acceptance is >= 0.2, then training
     time to solution (4 instances x 3 seeds).
  4. T_min = candidate with the smallest median training TTS among adequate
     ones (factor: c <= 1; cnf: measured p0(T_min) >= 0.5 with >= 50 round
     trips; if none qualifies, TTS alone and the choice is flagged).
  5. M scan {0.7, 1, 1.4, 2} x M0 at that T_min, feedback-optimized; M = smallest
     cost per round trip with min swap acceptance >= 0.15.
Writes ladders/<branch>_d<d>.json and appends to tune_<branch>.log.
"""
import csv, glob, json, math, os, re, statistics as st, subprocess, sys
import numpy as np
from scipy.optimize import brentq
import tune_ladder as tl

branch, dmin, dmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
os.makedirs("ladders", exist_ok=True)
logf = open(f"tune_{branch}.log", "a")
def log(s):
    print(s, flush=True); logf.write(s + "\n"); logf.flush()

train = [r for r in csv.DictReader(open("train.tsv"), delimiter="\t")]
CFG = dict(
    factor=dict(grid=list(np.geomspace(0.15, 8, 22)), single_sweeps=60000, sweeps0=20000, iters=4,
                tts_limit=60, cands=[0.6, 0.8, 1.0, 1.25]),
    cnf=dict(grid=list(np.geomspace(0.08, 4, 20)), single_sweeps=20000, sweeps0=4000, iters=3,
             tts_limit=60, cands=[0.10, 0.13, 0.17, 0.22, 0.28]),
)[branch]


def insts_for(d):
    out = []
    for r in train:
        if int(r["d"]) != d: continue
        p, q = int(r["p"]), int(r["q"])
        out.append(dict(branch=branch, d=d, p=p, q=q, N=int(r["N"]), np=p.bit_length(), nq=q.bit_length(), cnf=r["cnf"]))
    return out


def tstar(d):
    Ts = []
    for f in glob.glob(f"dos/d{d}_*.txt"):
        g = {int(a): int(b) for a, b in (l.split() for l in open(f))}
        E = np.array(sorted(g)); w = np.array([g[e] for e in E], float)
        P0 = lambda T: (w[E == 0] * 1.0).sum() / (w * np.exp(-E / T)).sum()
        Ts.append(brentq(lambda T: P0(T) - 0.5, 0.01, 50))
    return min(Ts), Ts


def tts(insts, temps, limit, seeds=(1, 2, 3)):
    jobs = [(i, s) for i in insts for s in seeds]
    def one(a):
        i, s = a
        o = tl.run(tl.model_args(i) + ["--mode", "solve", "--temps", ",".join(f"{t:.6g}" for t in temps),
                                       "--time", str(limit), "--seed", str(100 + s)])
        x = o.split()
        ok = x[1] == "SAT"
        return (float(x[x.index("time") + 1]) if ok else math.inf), (int(x[x.index("work") + 1]) if ok else math.inf)
    res = list(tl.POOL.map(one, jobs))
    t = sorted(r[0] for r in res); w = sorted(r[1] for r in res)
    return st.median(w), sum(1 for x in t if x < math.inf), len(t), st.median(t)


for d in range(dmin, dmax + 1):
    insts = insts_for(d)
    if not insts or (branch == "cnf" and not insts[0]["cnf"]):
        continue
    log(f"=== {branch} d={d}: {len(insts)} training instances")
    # 1. single-replica scan
    jobs = [(i, T) for T in CFG["grid"] for i in insts[:2]]
    sres = list(tl.POOL.map(lambda a: tl.single(a[0], a[1], CFG["single_sweeps"]), jobs))
    # need var(E): rerun parse via meanE only is not enough -> compute from a tune-like single call
    byT = {}
    for (inst, T), r in zip(jobs, sres):
        byT.setdefault(T, []).append(r)
    tau_med = {T: st.median(x["tau"] for x in v) for T, v in byT.items()}
    # lowest T such that every hotter grid temperature also decorrelates within
    # one sweep (scanning down from the hot end; frozen chains at very low T
    # have zero variance and a meaningless tau_int of 0.5)
    tmax = CFG["grid"][-1]
    for T in sorted(CFG["grid"], reverse=True):
        if tau_med[T] <= 1.0: tmax = T
        else: break
    log(f"  T_max = {tmax:.3f}  (tau_int at grid: " + ", ".join(f"{T:.2f}:{tau_med[T]:.2g}" for T in CFG["grid"]) + ")")
    # 2. provisional M0 from energy fluctuations (single-replica variances on the grid)
    var = {}
    for T in CFG["grid"]:
        vs = []
        for inst in insts[:2]:
            o = tl.run(tl.model_args(inst) + ["--mode", "single", "--tmin", f"{float(T):.10g}", "--sweeps", "20000", "--burn", "2000"])
            x = o.split(); vs.append(float(x[x.index("varE") + 1]))
        var[T] = st.median(vs)
    def M0_for(tmin):
        Tg = [T for T in CFG["grid"] if tmin * 0.999 <= T <= tmax * 1.001]
        if len(Tg) < 2: return 4
        b = 1 / np.array(Tg); s = np.sqrt([max(var[T], 1e-9) for T in Tg])
        I = float(np.trapezoid(s[::-1], b[::-1]))
        # add the stretch between tmin and the first grid point using the edge sigma
        I += s[0] * max(0.0, 1 / tmin - 1 / Tg[0])
        return int(max(4, round(1 + I / 1.47)))
    # 3. T_min candidates
    if branch == "factor":
        Ts, allTs = tstar(d)
        cands = [c * Ts for c in CFG["cands"]]
        log(f"  T* (P0 = 0.5, exact) = {Ts:.3f}; per instance {[round(x, 3) for x in allTs]}")
    else:
        cands = list(CFG["cands"])
    results = []
    for tmin in cands:
        M0 = M0_for(tmin)
        # grow M until every adjacent pair of the feedback-optimized ladder
        # swaps at least 20% of the time (single-replica variances underestimate
        # the fluctuations at low T, where an isolated replica freezes)
        for _ in range(6):
            final, hist = tl.optimize_ladder(insts, tmin, tmax, M0, CFG["sweeps0"], CFG["iters"], log=None)
            if min(final["acc"]) >= 0.2: break
            M0 = int(math.ceil(M0 * 1.4))
        med, nsol, n, tmed = tts(insts, final["T"], CFG["tts_limit"])
        results.append(dict(tmin=tmin, M=M0, ladder=final, tts_median=med, tts_time_median=tmed, solved=nsol, runs=n))
        log(f"  T_min {tmin:.3f}: M0 {M0}, round trips {final['round_trips']}, min acc {min(final['acc']):.3f}, "
            f"p0(T_min) {final['p0_tmin']:.3f}, training median work {med:.4g} replica-sweeps, time {tmed:.3g} s ({nsol}/{n} solved)")
    if branch == "factor":
        adequate = [r for r in results if r["tmin"] <= Ts * 1.0001]
        note = "exact P0(T_min) >= 0.5 on every training instance"
    else:
        adequate = [r for r in results if r["ladder"]["p0_tmin"] >= 0.5 and r["ladder"]["round_trips"] >= 50]
        note = "measured p0(T_min) >= 0.5 with >= 50 round trips"
        if not adequate:
            adequate = results
            note = "NOT VERIFIABLE: no candidate reached p0 >= 0.5 with >= 50 round trips; chosen by TTS alone"
    best = min(adequate, key=lambda r: (r["tts_median"], -r["solved"]))
    tmin = best["tmin"]; M0 = best["M"]
    log(f"  chosen T_min = {tmin:.3f}  (adequacy: {note})")
    # 5. M scan
    Ms = sorted({max(3, int(round(c * M0))) for c in (0.7, 1.0, 1.4, 2.0)})
    bestM, allM = tl.choose_M(insts, tmin, tmax, Ms, CFG["sweeps0"], CFG["iters"], log=log)
    if bestM is None:
        bestM = max(allM, key=lambda m: (min(allM[m]["acc"]) >= 0.15, -allM[m]["cost"]))
        log("  no ladder met both criteria; taking the best available")
    fin = allM[bestM]
    log(f"  chosen M = {bestM}; temperatures " + ", ".join(f"{t:.4g}" for t in fin["T"]))
    log(f"  swap acceptance " + ", ".join(f"{a:.2f}" for a in fin["acc"]))
    log(f"  flow f(T)       " + ", ".join(f"{x:.2f}" for x in fin["f"]))
    json.dump(dict(branch=branch, d=d, tmax=tmax, tau_scan={f"{T:.4g}": tau_med[T] for T in CFG["grid"]},
                   var_scan={f"{T:.4g}": var[T] for T in CFG["grid"]},
                   tmin_candidates=results, tstar=(Ts if branch == "factor" else None),
                   chosen=dict(tmin=tmin, M=bestM, T=fin["T"], acc=fin["acc"], f=fin["f"],
                               round_trips=fin["round_trips"], sweeps_total=fin["sweeps_total"],
                               cost=fin["cost"], p0_tmin=fin["p0_tmin"]),
                   M_scan={str(m): dict(cost=r["cost"], min_acc=min(r["acc"]), round_trips=r["round_trips"]) for m, r in allM.items()}),
              open(f"ladders/{branch}_d{d:02d}.json", "w"), indent=1, default=float)
