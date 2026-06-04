#!/usr/bin/env python3
"""
generate_instances.py — Batch generator of factorization-based planted-solution
SAT and Ising benchmark instances.

Given a digit count d and an instance count k, generates k instances,
each encoding N = p * q for two randomly chosen d-digit primes p and q.
By default, both DIMACS CNF and Ising TSV files are produced.

Output files:  factor_<serial>_p<p>_q<q>.cnf
               factor_<serial>_p<p>_q<q>_ising_E<E0>.tsv
               factor_<serial>_p<p>_q<q>_planted.txt

Usage:
  python generate_instances.py 3 10                  # 10 instances with 3-digit primes
  python generate_instances.py 5 20 -o ./instances   # output to directory
  python generate_instances.py 4 5 --no-ising        # CNF only, skip Ising files
  python generate_instances.py 3 10 --seed 42        # reproducible generation
  python generate_instances.py 6 1 --distinct-N      # ensure all N are distinct
"""

import argparse
import os
import sys
import random
from math import isqrt

# ---------------------------------------------------------------------------
# Prime generation utilities
# ---------------------------------------------------------------------------

def _is_prime_small(n):
    """Deterministic trial-division primality test (for n < ~10^7 or so)."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def _miller_rabin(n, a):
    """Single round of Miller-Rabin with witness a."""
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(r - 1):
        x = pow(x, 2, n)
        if x == n - 1:
            return True
    return False


def is_prime(n):
    """Primality test: deterministic for n < 3.3*10^24, probabilistic beyond."""
    if n < 2:
        return False
    # Small primes
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p:
            return True
        if n % p == 0:
            return False
    # Deterministic Miller-Rabin witnesses for n < 3,317,044,064,679,887,385,961,981
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if not _miller_rabin(n, a):
            return False
    return True


def random_d_digit_prime(d, rng):
    """Return a random d-digit prime (uniformly sampled by rejection)."""
    if d < 1:
        raise ValueError("d must be >= 1")
    lo = 10 ** (d - 1)
    hi = 10 ** d - 1
    if d == 1:
        lo = 2  # smallest 1-digit prime
    while True:
        n = rng.randint(lo, hi)
        # Make odd (except when d=1 we might want 2)
        if n > 2 and n % 2 == 0:
            n += 1
            if n > hi:
                continue
        if is_prime(n):
            return n


# ---------------------------------------------------------------------------
# Instance generation (wraps factorization_sat.py)
# ---------------------------------------------------------------------------

def generate_batch(d, k, output_dir=".", seed=None, also_ising=True,
                   distinct_n=False, verify=True):
    """Generate k factorization SAT instances with d-digit primes.

    Parameters
    ----------
    d : int
        Number of decimal digits for each prime factor.
    k : int
        Number of instances to generate.
    output_dir : str
        Directory for output files.
    seed : int or None
        RNG seed for reproducibility.
    also_ising : bool
        If True (default), also write Ising and planted-spin files.
    distinct_n : bool
        If True, ensure every instance has a distinct product N.
    verify : bool
        If True (default), verify each planted solution.

    Returns
    -------
    list of dict
        Metadata for each generated instance.
    """
    # Import the core module
    try:
        import factorization_sat as fsat
    except ImportError:
        # Try importing from the same directory as this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, script_dir)
        import factorization_sat as fsat

    os.makedirs(output_dir, exist_ok=True)
    rng = random.Random(seed)
    seen_N = set()
    results = []

    for idx in range(k):
        # Generate two d-digit primes
        attempts = 0
        while True:
            p = random_d_digit_prime(d, rng)
            q = random_d_digit_prime(d, rng)
            # Canonically order so p <= q (avoids trivial duplicates)
            if p > q:
                p, q = q, p
            N = p * q
            if distinct_n and N in seen_N:
                attempts += 1
                if attempts > 10000:
                    print(f"Warning: could not find distinct N after {attempts} "
                          f"attempts at index {idx}; reusing.", file=sys.stderr)
                    break
                continue
            break
        seen_N.add(N)

        # Stage 1: generate clauses
        vm, and_clauses, xor_clauses, pins, planted, p_vars, q_vars = \
            fsat.generate_clauses(p, q)

        # Stage 2: preprocess
        fsat.preprocess(vm, and_clauses, xor_clauses, pins, planted)

        # Stage 3: convert to CNF
        cnf_clauses, var_map, num_vars = fsat.to_cnf(vm, and_clauses, xor_clauses)

        # Verify planted solution
        if verify:
            ok, bad_idx = fsat.verify_cnf(cnf_clauses, var_map, planted, vm)
            if not ok:
                print(f"ERROR: planted solution fails clause {bad_idx} for "
                      f"p={p}, q={q}. Skipping.", file=sys.stderr)
                continue

        # Build filename
        basename = f"factor_{idx:04d}_p{p}_q{q}"
        cnf_path = os.path.join(output_dir, f"{basename}.cnf")

        comments = [
            f"Factorization SAT instance #{idx}: N = {N} = {p} x {q}",
            f"p = {p} ({p.bit_length()} bits), q = {q} ({q.bit_length()} bits)",
            f"N = {N} ({N.bit_length()} bits)",
            f"Generator: generate_instances.py  d={d}, k={k}, seed={seed}",
        ]

        # Write DIMACS CNF
        with open(cnf_path, "w") as f:
            fsat.write_dimacs(cnf_clauses, num_vars, file=f, comments=comments)

        info = {
            "index": idx,
            "p": p,
            "q": q,
            "N": N,
            "num_vars": num_vars,
            "num_clauses": len(cnf_clauses),
            "cnf_file": cnf_path,
        }

        # Optionally write Ising files
        if also_ising:
            ising_data = fsat.to_ising(vm, and_clauses, xor_clauses, planted)
            h, J, E0, spin_map, planted_spins, num_spins = ising_data

            if verify:
                ok_ising, energy = fsat.verify_ising(h, J, E0, planted_spins, num_spins)
                if not ok_ising:
                    print(f"ERROR: Ising verification failed for p={p}, q={q}, "
                          f"energy={energy}", file=sys.stderr)

            ising_path = os.path.join(output_dir, f"{basename}_ising_E{int(E0)}.tsv")
            with open(ising_path, "w") as f:
                fsat.write_ising(h, J, E0, num_spins, planted_spins, file=f)

            planted_path = os.path.join(output_dir, f"{basename}_planted.txt")
            with open(planted_path, "w") as f:
                fsat.write_planted(planted_spins, num_spins, file=f)

            info["ising_file"] = ising_path
            info["planted_file"] = planted_path
            info["num_spins"] = num_spins
            info["E0"] = int(E0)

        results.append(info)
        print(f"[{idx+1}/{k}] {basename}.cnf  "
              f"({num_vars} vars, {len(cnf_clauses)} clauses, "
              f"N={N})", file=sys.stderr)

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Batch-generate factorization SAT and Ising benchmark instances.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s 3 10                        # 10 instances, 3-digit primes (CNF + Ising)
  %(prog)s 5 20 -o ./benchmarks        # write to ./benchmarks/
  %(prog)s 4 5 --seed 42               # reproducible
  %(prog)s 3 10 --no-ising             # CNF only, skip Ising + planted files
  %(prog)s 6 1 --distinct-N            # all products N are distinct
        """)
    parser.add_argument("d", type=int, help="number of decimal digits per prime")
    parser.add_argument("k", type=int, help="number of instances to generate")
    parser.add_argument("-o", "--output-dir", default=".",
                        help="output directory (default: current directory)")
    parser.add_argument("--seed", type=int, default=None,
                        help="RNG seed for reproducibility")
    parser.add_argument("--no-ising", action="store_true",
                        help="skip Ising and planted-spin file output")
    parser.add_argument("--distinct-N", action="store_true",
                        help="ensure all instances have distinct products N")
    parser.add_argument("--no-verify", action="store_true",
                        help="skip planted-solution verification")
    parser.add_argument("--summary", action="store_true",
                        help="print summary table after generation")
    args = parser.parse_args()

    if args.d < 1:
        print("Error: d must be >= 1", file=sys.stderr)
        sys.exit(1)
    if args.k < 1:
        print("Error: k must be >= 1", file=sys.stderr)
        sys.exit(1)

    results = generate_batch(
        d=args.d,
        k=args.k,
        output_dir=args.output_dir,
        seed=args.seed,
        also_ising=not args.no_ising,
        distinct_n=args.distinct_N,
        verify=not args.no_verify,
    )

    if args.summary and results:
        print(f"\n{'='*72}")
        print(f"Generated {len(results)} instances  (d={args.d}, k={args.k}, "
              f"seed={args.seed})")
        print(f"{'='*72}")
        print(f"{'#':>5}  {'p':>12}  {'q':>12}  {'N':>24}  {'vars':>6}  {'clauses':>8}")
        print(f"{'-'*5}  {'-'*12}  {'-'*12}  {'-'*24}  {'-'*6}  {'-'*8}")
        for r in results:
            print(f"{r['index']:5d}  {r['p']:12d}  {r['q']:12d}  "
                  f"{r['N']:24d}  {r['num_vars']:6d}  {r['num_clauses']:8d}")
        print()


if __name__ == "__main__":
    main()
