#!/usr/bin/env python3
"""
make_benchmark_set.py -- build the CNF benchmark set and the Slurm command
list for the re-run of the solver-scaling study (Fig. 3 of the paper).

Design (fixes the problems of the March 2026 run):
  * symmetric instances: p and q are both d-bit primes, p != q;
  * the same number M of instances at every d in the requested range;
  * instances made by the fixed generator (forced-zero columns removed);
  * one Slurm task per (instance, solver), so the two solvers do not share
    a wall-time budget, each with its own result file;
  * an explicit per-solver timeout, recorded in the manifest.

Usage (from this folder, with the fixed generator in ../../../pq python):
  python3 make_benchmark_set.py --dmin 10 --dmax 28 --per-d 20 --seed 2026 \
      --timeout 144000 --generator "../../../pq python"

Outputs:
  cnf/d<d>/factor_<k>_p<p>_q<q>.cnf          (+ _cnf_solution.txt sidecar)
  manifest.tsv     one row per (instance, solver): d, p, q, N, vars, clauses,
                   cnf path, solver, timeout, result path
  commands.txt     one command per Slurm array task (same order as manifest)
"""
import argparse
import csv
import os
import random
import sys


def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:  # deterministic for n < 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def random_prime_bits(d, rng):
    lo, hi = 2 ** (d - 1), 2 ** d - 1
    while True:
        n = rng.randint(lo, hi) | (1 if d > 2 else 0)
        if lo <= n <= hi and is_prime(n):
            return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dmin", type=int, default=10)
    ap.add_argument("--dmax", type=int, default=28)
    ap.add_argument("--per-d", type=int, default=12, help="instances per bit-length")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--timeout", type=int, default=144000,
                    help="per-solver timeout in seconds (default 40 h)")
    ap.add_argument("--solvers", default="kissat,cadical")
    ap.add_argument("--generator", default="../../../pq python",
                    help="folder containing factorization_sat.py")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    sys.path.insert(0, os.path.abspath(args.generator))
    import factorization_sat as fsat  # the fixed generator

    rng = random.Random(args.seed)
    solvers = [s.strip() for s in args.solvers.split(",") if s.strip()]
    os.makedirs(os.path.join(args.outdir, "results"), exist_ok=True)
    rows = []

    for d in range(args.dmin, args.dmax + 1):
        ddir = os.path.join(args.outdir, "cnf", f"d{d:02d}")
        os.makedirs(ddir, exist_ok=True)
        seen = set()
        k = 0
        while k < args.per_d:
            p, q = random_prime_bits(d, rng), random_prime_bits(d, rng)
            if p == q:
                continue
            p, q = min(p, q), max(p, q)
            if p * q in seen:
                continue
            seen.add(p * q)

            vm, and_cl, xor_cl, pins, planted, pv, qv = fsat.generate_clauses(p, q)
            fsat.preprocess(vm, and_cl, xor_cl, pins, planted)
            cnf, var_map, num_vars = fsat.to_cnf(vm, and_cl, xor_cl)
            ok, _ = fsat.verify_cnf(cnf, var_map, planted, vm)
            if not ok:
                raise RuntimeError(f"planted solution fails CNF for {p}x{q}")

            stem = f"factor_{k:04d}_p{p}_q{q}"
            path = os.path.join(ddir, stem + ".cnf")
            comments = [f"Factorization SAT instance: N = {p*q} = {p} x {q}",
                        f"p, q: {d}-bit primes; seed={args.seed}",
                        "Generator: make_benchmark_set.py with factorization_sat.py (Nov 2026 fix)"]
            with open(path, "w") as f:
                fsat.write_dimacs(cnf, num_vars, file=f, comments=comments)
            # planted solution + factor-bit metadata, read by sat_benchmark for
            # the factorization check (result tag "F")
            fsat.write_cnf_solution(os.path.join(ddir, stem + "_cnf_solution.txt"),
                                    p, q, vm, var_map, num_vars, planted, pv, qv)

            for s in solvers:
                rel_cnf = os.path.relpath(path, args.outdir)
                res = os.path.join("results", s, f"d{d:02d}", stem + ".tsv")
                rows.append(dict(d=d, p=p, q=q, N=p * q, vars=num_vars, clauses=len(cnf),
                                 cnf=rel_cnf, solver=s, timeout=args.timeout, result=res))
            k += 1
        print(f"d={d}: {args.per_d} instances", flush=True)

    with open(os.path.join(args.outdir, "manifest.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(args.outdir, "commands.txt"), "w") as f:
        for r in rows:
            os.makedirs(os.path.join(args.outdir, os.path.dirname(r["result"])), exist_ok=True)
            f.write(f"./sat_benchmark --solver {r['solver']} --timeout {r['timeout']} "
                    f"--tsv '{r['result']}' '{r['cnf']}'\n")
    print(f"{len(rows)} tasks written to commands.txt "
          f"(use --array=0-{len(rows)-1}%<max concurrent> in the job file)")


if __name__ == "__main__":
    main()
