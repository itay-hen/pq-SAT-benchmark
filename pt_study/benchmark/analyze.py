#!/usr/bin/env python3
"""Summarize results.tsv: censored medians and quartiles per (method, d),
log-linear fits of the median, and a comparison figure."""
import csv, math, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TO = 300.0
R = list(csv.DictReader(open("results.tsv"), delimiter="\t"))
methods = ["kissat", "pt_factor", "pt_focused", "pt_uniform", "random_factor"]
label = {"kissat": "Kissat 4.0.4 (CDCL)",
         "pt_factor": "PT, factor bits (E = Hamming distance)",
         "pt_focused": "PT, residual CNF, focused moves",
         "pt_uniform": "PT, residual CNF, uniform moves",
         "random_factor": "random walk, factor bits (T = inf)"}


def q(vals, f):
    """Quantile of censored data (timeouts as +inf); None if it is censored."""
    v = sorted(vals)
    x = np.quantile(np.array([t if t < math.inf else 1e30 for t in v]), f, method="linear")
    return None if x >= TO else float(x)


summ = {}
with open("summary.tsv", "w") as out:
    out.write("method\td\tn\tsolved\tq1\tmedian\tq3\n")
    for m in methods:
        for d in range(5, 19):
            rs = [r for r in R if r["method"] == m and int(r["d"]) == d]
            if not rs:
                continue
            t = [float(r["time"]) if r["status"] == "SAT" else math.inf for r in rs]
            # resolution floor: times are reported to 0.1 ms
            t = [max(x, 1e-4) for x in t]
            s = (len(rs), sum(x < math.inf for x in t), q(t, .25), q(t, .5), q(t, .75))
            summ[(m, d)] = s
            f = lambda x: "censored" if x is None else f"{x:.4g}"
            out.write(f"{m}\t{d}\t{s[0]}\t{s[1]}\t{f(s[2])}\t{f(s[3])}\t{f(s[4])}\n")

fits = {}
for m in methods:
    pts = [(d, summ[(m, d)][3]) for d in range(5, 19) if (m, d) in summ and summ[(m, d)][3]]
    # fit over the range where medians are above the 1 ms timer/startup floor
    pts = [(d, t) for d, t in pts if t >= 1e-3]
    if len(pts) >= 3:
        x = np.array([d for d, _ in pts]); y = np.log10([t for _, t in pts])
        b, a = np.polyfit(x, y, 1)
        fits[m] = (b, a, int(x.min()), int(x.max()))
with open("fits.tsv", "w") as out:
    out.write("method\tslope_decades_per_bit\tintercept\td_min\td_max\n")
    for m, (b, a, lo, hi) in fits.items():
        out.write(f"{m}\t{b:.3f}\t{a:.3f}\t{lo}\t{hi}\n")

colors = {"kissat": "#1f77b4", "pt_factor": "#2ca02c", "pt_focused": "#d62728",
          "pt_uniform": "#ff7f0e", "random_factor": "#7f7f7f"}
marks = {"kissat": "o", "pt_factor": "s", "pt_focused": "^", "pt_uniform": "v", "random_factor": "x"}
fig, ax = plt.subplots(figsize=(5.6, 4.6))
for m in methods:
    ds = [d for d in range(5, 19) if (m, d) in summ and summ[(m, d)][3] is not None]
    if not ds:
        continue
    med = np.array([summ[(m, d)][3] for d in ds])
    lo = np.array([summ[(m, d)][2] or med[i] for i, d in enumerate(ds)])
    hi = np.array([summ[(m, d)][4] if summ[(m, d)][4] is not None else TO for d in ds])
    ax.errorbar(ds, med, yerr=[med - lo, hi - med], fmt=marks[m], color=colors[m], ms=5,
                capsize=2, lw=1, label=label[m])
    if m in fits:
        b, a, l, h = fits[m]
        xx = np.linspace(l, h, 20)
        ax.plot(xx, 10 ** (a + b * xx), color=colors[m], lw=0.8, ls="--")
    # mark the first d with a censored median
    cens = [d for d in range(5, 19) if (m, d) in summ and summ[(m, d)][3] is None]
    if cens:
        ax.plot([cens[0]], [TO], marker=marks[m], mfc="none", color=colors[m], ms=7, ls="none")
ax.axhline(TO, color="0.6", lw=0.6, ls=":")
ax.text(18.2, TO * 1.15, "300 s limit", ha="right", va="bottom", fontsize=8, color="0.4")
ax.set_yscale("log")
ax.set_xlabel("bits per factor, d")
ax.set_ylabel("median time to solution (s)")
ax.set_xticks(range(5, 19))
ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)
ax.text(5, 1.4e-4, "0.1 ms timer resolution", fontsize=7, color="0.4", va="bottom")
ax.grid(alpha=0.25, which="both", lw=0.4)
fig.tight_layout()
fig.savefig("pt_vs_kissat.pdf"); fig.savefig("pt_vs_kissat.png", dpi=200)
print(open("summary.tsv").read()); print(open("fits.tsv").read())
