#!/usr/bin/env python3
"""Per-d temperature-ladder optimization for parallel tempering.

Steps for one (branch, d):
  1. T_max  = smallest T on a geometric grid at which the integrated
              autocorrelation time of E for a single replica is <= 1 sweep
              (median over training instances).
  2. ladder shape: feedback optimization (Katzgraber, Trebst, Huse, Troyer,
              J. Stat. Mech. P03018, 2006).  Replicas are labelled by the end
              they last visited; f(T) = n_up / (n_up + n_down).  New
              temperatures put equal shares of sqrt(df) in each interval,
              i.e. density eta'(T) ~ sqrt(df/dT / dT).  Iterated with damping
              and growing run length; f is averaged over training instances.
  3. M      = number of temperatures minimising the cost per round trip,
              M * sweeps / round_trips, over a small set of M, each ladder
              feedback-optimized; ladders whose minimum adjacent swap
              acceptance falls below 0.15 are rejected.
  4. T_min  = chosen by the caller (thermodynamic criterion and training
              time-to-solution); this module optimizes the ladder for a given
              T_min and reports the ground-state occupation p0 at T_min.
"""
import json, math, subprocess, statistics as st
from concurrent.futures import ThreadPoolExecutor
import numpy as np

PT = "/home/claude/pt/pt"
POOL = ThreadPoolExecutor(2)


def model_args(inst):
    if inst["branch"] == "cnf":
        return ["--model", "cnf", inst["cnf"], "--move", "uniform"]
    return ["--model", "factor", "--N", str(inst["N"]), "--np", str(inst["np"]), "--nq", str(inst["nq"])]


def run(args, timeout=None):
    return subprocess.run([PT] + args, capture_output=True, text=True, timeout=timeout).stdout


def single(inst, T, sweeps, seed=1):
    o = run(model_args(inst) + ["--mode", "single", "--tmin", f"{float(T):.10g}", "--sweeps", str(sweeps), "--burn", "2000", "--seed", str(seed)])
    x = o.split()
    return dict(T=T, meanE=float(x[x.index("meanE") + 1]), tau=float(x[x.index("tau_int") + 1]),
                acc=float(x[x.index("flip_acc") + 1]))


def tune_run(inst, temps, sweeps, burn, seed=1):
    o = run(model_args(inst) + ["--mode", "tune", "--temps", ",".join(f"{t:.6g}" for t in temps),
                                "--sweeps", str(sweeps), "--burn", str(burn), "--seed", str(seed)])
    lines = o.strip().splitlines()
    hdr = lines[0].split()
    rt = int(hdr[hdr.index("round_trips") + 1])
    rows = [l.split() for l in lines[2:]]
    f = np.array([float(r[3]) for r in rows]); acc = np.array([float(r[2]) for r in rows[:-1]])
    p0 = np.array([float(r[6]) for r in rows]); meanE = np.array([float(r[4]) for r in rows])
    return dict(rt=rt, f=f, acc=acc, p0=p0, meanE=meanE, sweeps=sweeps)


def pick_tmax(insts, grid, sweeps):
    jobs = [(i, T) for T in grid for i in insts]
    res = list(POOL.map(lambda a: single(a[0], a[1], sweeps), jobs))
    by_T = {}
    for r in res:
        by_T.setdefault(r["T"], []).append(r["tau"])
    for T in grid:
        if st.median(by_T[T]) <= 1.0:
            return T, {T: st.median(v) for T, v in by_T.items()}
    return grid[-1], {T: st.median(v) for T, v in by_T.items()}


def feedback_update(T, f, damp=0.5):
    """New temperatures with equal sqrt(df) per interval (KTHT 2006)."""
    M = len(T)
    f = np.array(f, float)
    f[0], f[-1] = 1.0, 0.0
    # fill NaN (temperatures never reached by labelled replicas) by interpolation
    good = ~np.isnan(f)
    f = np.interp(np.arange(M), np.arange(M)[good], f[good])
    # enforce monotone decrease (pool adjacent violators)
    f = -_isotonic(-f)
    df = np.maximum(f[:-1] - f[1:], 1e-4)
    w = np.sqrt(df)
    cum = np.concatenate([[0], np.cumsum(w)]); cum /= cum[-1]
    targets = np.linspace(0, 1, M)
    newT = np.interp(targets, cum, T)
    out = (1 - damp) * np.array(T) + damp * newT
    out[0], out[-1] = T[0], T[-1]
    return np.sort(out)


def _isotonic(y):
    """Non-decreasing least-squares fit (pool adjacent violators)."""
    y = list(map(float, y)); blocks = [[v, 1] for v in y]
    i = 0
    while i < len(blocks) - 1:
        if blocks[i][0] > blocks[i + 1][0]:
            v = (blocks[i][0] * blocks[i][1] + blocks[i + 1][0] * blocks[i + 1][1]) / (blocks[i][1] + blocks[i + 1][1])
            blocks[i] = [v, blocks[i][1] + blocks[i + 1][1]]; del blocks[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return np.array([v for v, n in blocks for _ in range(n)])


def optimize_ladder(insts, tmin, tmax, M, sweeps0, iters=4, log=None):
    T = np.geomspace(tmin, tmax, M)
    hist = []
    sweeps = sweeps0
    for it in range(iters + 1):
        res = list(POOL.map(lambda i: tune_run(i, T, sweeps, max(1000, sweeps // 20), seed=11 + it), insts))
        f = np.nanmean(np.array([r["f"] for r in res]), axis=0)
        acc = np.mean(np.array([r["acc"] for r in res]), axis=0)
        rt = sum(r["rt"] for r in res)
        p0 = float(np.mean([r["p0"][0] for r in res]))
        rec = dict(iter=it, T=[float(t) for t in T], f=[float(x) for x in f], acc=[float(a) for a in acc],
                   round_trips=rt, sweeps_total=sweeps * len(insts), p0_tmin=p0,
                   cost=(M * sweeps * len(insts) / rt) if rt else math.inf)
        hist.append(rec)
        if log:
            log(f"    M={M} iter {it}: round trips {rt}, min acc {acc.min():.3f}, p0(Tmin) {p0:.3f}, "
                f"cost/rt {rec['cost']:.3g}")
        if it == iters:
            break
        T = feedback_update(T, f)
        sweeps *= 2 if it < 2 else 1
    return hist[-1], hist


def choose_M(insts, tmin, tmax, Ms, sweeps0, iters=4, log=None):
    best, allres = None, {}
    for M in Ms:
        final, hist = optimize_ladder(insts, tmin, tmax, M, sweeps0, iters, log)
        allres[M] = final
        ok = min(final["acc"]) >= 0.15 and final["round_trips"] > 0
        if ok and (best is None or final["cost"] < allres[best]["cost"]):
            best = M
    return best, allres
