#!/usr/bin/env python3
"""
analyze_rerun.py -- medians, quartiles and scaling fits for the solver re-run.

Reads manifest.tsv (every planned (instance, solver) task) and the per-task
result files written by sat_benchmark.  Every planned task is accounted for:

  solved   result file has a row with SAT and a passing check (F or C)
  timeout  result is TIMEOUT (the harness's per-solver limit was reached)
  killed   result file missing or header-only (Slurm ended the task first)
  error    any other result (ERROR, failed check, ...)

Timeouts and killed runs are right-censored: their runtime is at least the
timeout.  The median at a given d is reported only if more than half of the
runs at that d were solved; otherwise only the lower bound ">= timeout" is
known and the point is shown as a censored marker and excluded from the fit.
Quartiles follow the same rule.  Errors are listed and excluded.

Usage:
  python3 analyze_rerun.py [--manifest manifest.tsv] [--min-n 10] [--logs logs]

Outputs: rerun_summary.tsv, rerun_fits.tsv, median_runtime_vs_bits_rerun.{pdf,png}
"""
import argparse
import csv
import glob
import math
import os
import re
from collections import defaultdict

import numpy as np


def read_result(path):
    """Return (status, time_s) for one result file."""
    if not os.path.exists(path):
        return "killed", None
    with open(path) as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if len(rows) < 2:
        return "killed", None
    header, row = rows[0], rows[1]
    t_col = next(i for i, h in enumerate(header) if h.endswith("_time"))
    r_col = next(i for i, h in enumerate(header) if h.endswith("_result"))
    t = float(row[t_col])
    tag = row[r_col].split()
    if tag and tag[0] == "TIMEOUT":
        return "timeout", t
    if tag and tag[0] == "SAT" and (len(tag) == 1 or tag[1] in ("F", "C")):
        return "solved", t
    return "error", t


def censored_quantile(times_solved, n_total, q):
    """q-quantile when unsolved runs are known to exceed every solved time.

    Returns None if the quantile falls among the censored runs."""
    k = q * (n_total - 1)           # numpy 'linear' convention on the full sample
    hi = math.ceil(k)
    if hi >= len(times_solved):
        return None
    s = sorted(times_solved)
    lo = math.floor(k)
    return s[lo] + (k - lo) * (s[hi] - s[lo])


def versions_from_logs(logdir):
    seen = defaultdict(set)
    for f in glob.glob(os.path.join(logdir, "out_*.txt")):
        with open(f, errors="ignore") as fh:
            for line in fh:
                m = re.match(r"(kissat|cadical)_version: (.*)", line)
                if m:
                    seen[m[1]].add(m[2].strip())
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="manifest.tsv")
    ap.add_argument("--logs", default="logs")
    ap.add_argument("--min-n", type=int, default=10,
                    help="minimum number of runs at a d for it to enter the fit")
    args = ap.parse_args()

    with open(args.manifest) as f:
        tasks = list(csv.DictReader(f, delimiter="\t"))
    base = os.path.dirname(os.path.abspath(args.manifest))

    by = defaultdict(lambda: defaultdict(list))   # solver -> d -> [(status, t)]
    errors = []
    for t in tasks:
        status, time_s = read_result(os.path.join(base, t["result"]))
        by[t["solver"]][int(t["d"])].append((status, time_s, float(t["timeout"])))
        if status == "error":
            errors.append((t["solver"], t["cnf"]))

    summary, fits = [], []
    for solver in sorted(by):
        xs, ys = [], []
        for d in sorted(by[solver]):
            runs = by[solver][d]
            solved = [t for s, t, _ in runs if s == "solved"]
            n_err = sum(1 for s, _, _ in runs if s == "error")
            n = len(runs) - n_err
            n_cens = n - len(solved)
            med = censored_quantile(solved, n, 0.5) if n else None
            q1 = censored_quantile(solved, n, 0.25) if n else None
            q3 = censored_quantile(solved, n, 0.75) if n else None
            tmax = max(to for _, _, to in runs)
            summary.append(dict(solver=solver, d=d, n=n, solved=len(solved), censored=n_cens,
                                errors=n_err, q1=q1, median=med, q3=q3,
                                median_lower_bound=(None if med is not None else tmax)))
            if med is not None and n >= args.min_n:
                xs.append(d)
                ys.append(math.log10(med))
        if len(xs) >= 3:
            a, b = np.polyfit(xs, ys, 1)
            fits.append(dict(solver=solver, slope_log10_per_bit=a, intercept=b,
                             factor_per_bit=10 ** a, d_min=min(xs), d_max=max(xs), points=len(xs)))

    def fmt(v):
        return "" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v))

    with open("rerun_summary.tsv", "w") as f:
        keys = list(summary[0].keys())
        f.write("\t".join(keys) + "\n")
        for r in summary:
            f.write("\t".join(fmt(r[k]) for k in keys) + "\n")
    if fits:
        with open("rerun_fits.tsv", "w") as f:
            keys = list(fits[0].keys())
            f.write("\t".join(keys) + "\n")
            for r in fits:
                f.write("\t".join(fmt(r[k]) for k in keys) + "\n")

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"kissat": "#2a78d6", "cadical": "#eb6834", "glucose": "#1baa7c"}
    names = {"kissat": "Kissat", "cadical": "CaDiCaL", "glucose": "Glucose"}
    plt.rcParams.update({"font.size": 11})
    fig, ax = plt.subplots(figsize=(5.2, 3.6), constrained_layout=True)
    for solver in sorted(by):
        rows = [r for r in summary if r["solver"] == solver]
        ok = [r for r in rows if r["median"] is not None]
        cens = [r for r in rows if r["median"] is None]
        c = colors.get(solver, "#555555")
        fit = next((f for f in fits if f["solver"] == solver), None)
        lab = names.get(solver, solver)
        if fit:
            lab += f" (slope {fit['slope_log10_per_bit']:.3f} decades/bit)"
        lo = [r["median"] - (r["q1"] if r["q1"] is not None else r["median"]) for r in ok]
        hi = [((r["q3"] if r["q3"] is not None else r["median"]) - r["median"]) for r in ok]
        ax.errorbar([r["d"] for r in ok], [r["median"] for r in ok], yerr=[lo, hi],
                    fmt="o", ms=5, color=c, capsize=2.5, lw=1.4, label=lab)
        if cens:
            ax.plot([r["d"] for r in cens], [r["median_lower_bound"] for r in cens], "^",
                    color=c, ms=6, mfc="none")
        if fit:
            xx = np.array([fit["d_min"], fit["d_max"]], float)
            ax.plot(xx, 10 ** (fit["slope_log10_per_bit"] * xx + fit["intercept"]), "--",
                    color=c, lw=1)
    if any(r["median"] is None for r in summary):
        ax.plot([], [], "^", color="#555555", mfc="none", ms=6,
                label="median censored (more than half unsolved within the timeout)")
    ax.set_yscale("log")
    ax.set_xlabel("Bit-length $d$ of each factor")
    ax.set_ylabel("Median runtime (s)")
    ax.grid(True, color="#e6e5e0", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    fig.savefig("median_runtime_vs_bits_rerun.pdf")
    fig.savefig("median_runtime_vs_bits_rerun.png", dpi=200)

    # ---- report ----
    print(f"{'solver':8} {'d':>3} {'n':>3} {'solved':>6} {'cens':>4} {'err':>3} "
          f"{'Q1':>10} {'median':>10} {'Q3':>10}")
    for r in summary:
        med = fmt(r["median"]) or f">= {r['median_lower_bound']:.0f}"
        print(f"{r['solver']:8} {r['d']:>3} {r['n']:>3} {r['solved']:>6} {r['censored']:>4} "
              f"{r['errors']:>3} {fmt(r['q1']):>10} {med:>10} {fmt(r['q3']):>10}")
    for f in fits:
        print(f"fit {f['solver']}: slope {f['slope_log10_per_bit']:.3f} log10(s)/bit "
              f"(x{f['factor_per_bit']:.2f} per bit), d={f['d_min']}..{f['d_max']}, "
              f"{f['points']} points")
    if errors:
        print(f"{len(errors)} runs ended in ERROR or failed verification (excluded):")
        for s, c in errors[:20]:
            print("   ", s, c)
    vers = versions_from_logs(os.path.join(base, args.logs))
    for s, v in vers.items():
        print(f"{s} versions seen in logs: {sorted(v)}")
    if not vers:
        print("no solver versions found in logs (check the logs folder)")


if __name__ == "__main__":
    main()
