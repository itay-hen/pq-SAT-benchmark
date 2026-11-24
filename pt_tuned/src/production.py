#!/usr/bin/env python3
"""Run the test instances with the tuned per-d ladders.

usage: production.py factor|cnf
Test set = the instances of the first benchmark (d = 5..18), same seed (1),
same 300 s limit.  A branch stops after the first d with a censored median.
Results appended to prod_<branch>.tsv."""
import csv, json, math, os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
PT = "/home/claude/pt/pt"; TO = 300
branch = sys.argv[1]
out = f"prod_{branch}.tsv"
rows = []
seen = set()
for man, roots in [("/home/claude/ptbench/small_manifest.tsv", ["/home/claude/ptbench/s5", "/home/claude/ptbench/s69"]),
                   ("/home/claude/rerun/manifest.tsv", ["/home/claude/rerun"])]:
    for r in csv.DictReader(open(man), delimiter="\t"):
        key = (int(r["d"]), int(r["p"]), int(r["q"]))
        if key in seen: continue
        seen.add(key)
        path = r["cnf"].replace("CNFROOT/", "cnf/")
        for b in roots:
            if os.path.exists(os.path.join(b, path)): path = os.path.join(b, path); break
        rows.append(dict(d=key[0], p=key[1], q=key[2], N=int(r["N"]), cnf=path))
if not os.path.exists(out):
    open(out, "w").write("d\tp\tq\tstatus\ttime\twork\tround_trips\tmin_swap_acc\tM\tTmin\tTmax\n")
done = {(int(l.split("\t")[0]), int(l.split("\t")[1]), int(l.split("\t")[2])) for l in list(open(out))[1:]}

def one(r, lad):
    temps = ",".join(f"{t:.6g}" for t in lad["chosen"]["T"])
    if branch == "cnf":
        args = ["--model", "cnf", r["cnf"], "--move", "uniform"]
    else:
        args = ["--model", "factor", "--N", str(r["N"]), "--np", str(r["p"].bit_length()), "--nq", str(r["q"].bit_length())]
    o = subprocess.run([PT] + args + ["--mode", "solve", "--temps", temps, "--time", str(TO), "--seed", "1"],
                       capture_output=True, text=True).stdout.split()
    g = lambda k: o[o.index(k) + 1]
    return o[1], float(g("time")), int(g("work")), int(g("round_trips")), float(g("min_swap_acc"))

pool = ThreadPoolExecutor(2)
for d in sorted({r["d"] for r in rows}):
    f = f"ladders/{branch}_d{d:02d}.json"
    if not os.path.exists(f):
        continue
    lad = json.load(open(f))
    batch = [r for r in rows if r["d"] == d and (r["d"], r["p"], r["q"]) not in done]
    res = list(pool.map(lambda r: (r, one(r, lad)), batch))
    with open(out, "a") as fo:
        for r, (stt, t, w, rt, ma) in res:
            fo.write(f"{d}\t{r['p']}\t{r['q']}\t{stt}\t{t:.4f}\t{w}\t{rt}\t{ma:.3f}\t{lad['chosen']['M']}\t"
                     f"{lad['chosen']['tmin']:.4g}\t{lad['tmax']:.4g}\n")
    allr = [l.split("\t") for l in list(open(out))[1:] if int(l.split("\t")[0]) == d]
    nun = sum(x[3] != "SAT" for x in allr)
    print(f"{branch} d={d}: {len(allr)} runs, {nun} unsolved", flush=True)
    if nun >= len(allr) / 2:
        print(f"{branch}: stopped after d={d} (median censored)", flush=True); break
