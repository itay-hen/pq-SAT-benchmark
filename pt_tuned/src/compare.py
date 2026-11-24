#!/usr/bin/env python3
"""Tuned vs untuned PT vs Kissat on the test set; ladder summary table and figures."""
import csv, glob, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TO = 300.0
def med_q(ts):
    a = np.array(sorted(ts))
    q = lambda f: (lambda x: None if x >= TO else float(max(x, 1e-4)))(np.quantile(np.where(np.isinf(a), 1e30, a), f))
    return q(.25), q(.5), q(.75)

data = {}
for r in csv.DictReader(open("/home/claude/ptbench/results.tsv"), delimiter="\t"):
    t = float(r["time"]) if r["status"] == "SAT" else math.inf
    data.setdefault((r["method"], int(r["d"])), []).append(t)
for br, f in [("factor", "prod_factor.tsv"), ("cnf", "prod_cnf.tsv")]:
    for r in csv.DictReader(open(f), delimiter="\t"):
        t = float(r["time"]) if r["status"] == "SAT" else math.inf
        data.setdefault((f"tuned_{br}", int(r["d"])), []).append(t)

methods = ["kissat", "pt_factor", "tuned_factor", "pt_uniform", "tuned_cnf", "pt_focused"]
label = {"kissat": "Kissat 4.0.4",
         "pt_factor": "factor-bit PT, untuned ladder", "tuned_factor": "factor-bit PT, tuned ladder",
         "pt_uniform": "CNF PT, untuned ladder", "tuned_cnf": "CNF PT, tuned ladder",
         "pt_focused": "CNF PT, focused moves (untuned, no detailed balance)"}
color = {"kissat": "#1f77b4", "pt_factor": "#98df8a", "tuned_factor": "#2ca02c",
         "pt_uniform": "#ffbb78", "tuned_cnf": "#d62728", "pt_focused": "#9467bd"}
mk = {"kissat": "o", "pt_factor": "s", "tuned_factor": "s", "pt_uniform": "v", "tuned_cnf": "v", "pt_focused": "^"}

rows = []
for m in methods:
    for d in range(5, 19):
        if (m, d) not in data: continue
        ts = data[(m, d)]
        q1, md, q3 = med_q(ts)
        rows.append((m, d, len(ts), sum(t < math.inf for t in ts), q1, md, q3))
with open("compare_summary.tsv", "w") as f:
    f.write("method\td\tn\tsolved\tq1\tmedian\tq3\n")
    for r in rows:
        f.write("\t".join(str(x if x is not None else "censored") if not isinstance(x, float) else f"{x:.4g}" for x in r) + "\n")

fits = {}
for m, lo, hi in [("kissat", 10, 18), ("pt_factor", 12, 18), ("tuned_factor", 12, 18), ("pt_uniform", 5, 7), ("tuned_cnf", 5, 8)]:
    pts = [(d, md) for (mm, d, n, s, q1, md, q3) in rows if mm == m and lo <= d <= hi and md]
    if len(pts) >= 3:
        x = np.array([p[0] for p in pts]); y = np.log10([p[1] for p in pts]); b, a = np.polyfit(x, y, 1)
        fits[m] = (b, a, x.min(), x.max())
with open("compare_fits.tsv", "w") as f:
    f.write("method\tslope_decades_per_bit\tintercept\td_min\td_max\n")
    for m, (b, a, l, h) in fits.items(): f.write(f"{m}\t{b:.3f}\t{a:.3f}\t{l}\t{h}\n")

fig, ax = plt.subplots(figsize=(5.8, 4.8))
for m in methods:
    pts = [(d, q1, md, q3) for (mm, d, n, s, q1, md, q3) in rows if mm == m and md is not None]
    if not pts: continue
    d_, q1, md, q3 = map(np.array, zip(*[(p[0], p[1] or p[2], p[2], p[3] or TO) for p in pts]))
    ls = "--" if m.startswith("pt_") else "-"
    ax.errorbar(d_, md, yerr=[md - q1, q3 - md], fmt=mk[m], color=color[m], ms=5, capsize=2, lw=1,
                mfc="none" if m.startswith("pt_") else color[m], label=label[m])
    if m in fits:
        b, a, l, h = fits[m]; xx = np.linspace(l, h, 20); ax.plot(xx, 10 ** (a + b * xx), color=color[m], lw=0.8, ls=ls)
    cens = [d for (mm, d, n, s, q1, md, q3) in rows if mm == m and md is None]
    if cens: ax.plot([min(cens)], [TO], marker=mk[m], mfc="none", color=color[m], ms=8, ls="none")
ax.axhline(TO, color="0.6", lw=0.6, ls=":"); ax.text(18.2, TO * 1.15, "300 s limit", ha="right", fontsize=8, color="0.4")
ax.set_yscale("log"); ax.set_xlabel("bits per factor, d"); ax.set_ylabel("median time to solution (s)")
ax.set_xticks(range(5, 19)); ax.grid(alpha=0.25, which="both", lw=0.4)
ax.text(5, 1.4e-4, "0.1 ms timer resolution", fontsize=7, color="0.4", va="bottom")
ax.legend(fontsize=7.2, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
fig.tight_layout(); fig.savefig("tuned_vs_untuned.pdf"); fig.savefig("tuned_vs_untuned.png", dpi=200)

# ladder figure: T_min, T_max, T*, M per d for both branches
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2))
for br, c in [("factor", "#2ca02c"), ("cnf", "#d62728")]:
    J = [json.load(open(f)) for f in sorted(glob.glob(f"ladders/{br}_d*.json"))]
    d = [j["d"] for j in J]
    axs[0].plot(d, [j["chosen"]["tmin"] for j in J], "o-", color=c, ms=4, label=f"{br}: $T_{{min}}$")
    axs[0].plot(d, [j["tmax"] for j in J], "s--", color=c, ms=4, mfc="none", label=f"{br}: $T_{{max}}$")
    if br == "factor":
        axs[0].plot(d, [j["tstar"] for j in J], "k:", lw=1, label="factor: $T^*$ (exact $P_0=0.5$)")
    axs[1].plot(d, [j["chosen"]["M"] for j in J], "o-", color=c, ms=4, label=br)
axs[0].set_yscale("log"); axs[0].set_xlabel("d"); axs[0].set_ylabel("temperature")
axs[0].legend(fontsize=6.5, frameon=False); axs[1].set_xlabel("d"); axs[1].set_ylabel("number of temperatures M")
axs[1].legend(fontsize=7, frameon=False)
for a in axs: a.grid(alpha=0.25, lw=0.4); a.set_xticks(range(5, 19, 1)); a.tick_params(axis="x", labelsize=7)
fig.tight_layout(); fig.savefig("ladders.pdf"); fig.savefig("ladders.png", dpi=200)
print(open("compare_summary.tsv").read()); print(open("compare_fits.tsv").read())
