#!/usr/bin/env python3
"""
verify_factorization_cnf.py — Verify a factorization SAT instance using p and q
extracted from the filename.

Given a .cnf file whose name contains p and q (e.g. factor_0005_p203767_q767317.cnf),
this script:
  1. Extracts p and q from the filename.
  2. Regenerates the SAT instance via factorization_sat.py (same pipeline that
     created the file) to recover the variable mapping.
  3. Constructs the planted assignment in the CNF variable space.
  4. Checks every clause in the *original file on disk* against that assignment.

This confirms end-to-end that knowing the factors is sufficient to satisfy the
exact CNF instance stored on disk.

Usage:
  python verify_factorization_cnf.py factor_0005_p203767_q767317.cnf
  python verify_factorization_cnf.py *.cnf          # batch verify
  python verify_factorization_cnf.py file.cnf -v    # verbose output
"""

import argparse
import os
import re
import sys


# ---------------------------------------------------------------------------
# Extract p, q from filename
# ---------------------------------------------------------------------------

def extract_pq_from_filename(filepath):
    """Parse p and q from a filename like factor_0005_p203767_q767317.cnf.

    Returns (p, q) as integers, or raises ValueError.
    """
    basename = os.path.basename(filepath)
    # Match _p<digits>_ and _q<digits> anywhere in the filename
    m_p = re.search(r'_p(\d+)', basename)
    m_q = re.search(r'_q(\d+)', basename)
    if not m_p or not m_q:
        raise ValueError(
            f"Cannot extract p and q from filename '{basename}'. "
            f"Expected pattern like factor_XXXX_p<number>_q<number>.cnf"
        )
    return int(m_p.group(1)), int(m_q.group(1))


# ---------------------------------------------------------------------------
# Read DIMACS CNF from file
# ---------------------------------------------------------------------------

def read_dimacs(filepath):
    """Read a DIMACS CNF file. Returns (num_vars, clauses).

    clauses is a list of lists of signed integers.
    """
    clauses = []
    num_vars = 0
    num_clauses_declared = 0
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("c"):
                continue
            if line.startswith("p "):
                parts = line.split()
                # p cnf <num_vars> <num_clauses>
                num_vars = int(parts[2])
                num_clauses_declared = int(parts[3])
                continue
            # Clause line: integers terminated by 0
            lits = list(map(int, line.split()))
            if lits and lits[-1] == 0:
                lits = lits[:-1]
            if lits:
                clauses.append(lits)
    return num_vars, num_clauses_declared, clauses


# ---------------------------------------------------------------------------
# Build planted assignment from factors
# ---------------------------------------------------------------------------

def build_planted_assignment(p, q):
    """Regenerate the factorization SAT instance for (p, q) and return the
    planted assignment as a dict {cnf_var: bool}.

    Uses factorization_sat.py's pipeline identically to how generate_instances.py
    creates the CNF, ensuring the variable mapping matches.
    """
    try:
        import factorization_sat as fsat
    except ImportError:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, script_dir)
        import factorization_sat as fsat

    # Reproduce the exact generation pipeline
    vm, and_clauses, xor_clauses, pins, planted, p_vars, q_vars = \
        fsat.generate_clauses(p, q)

    fsat.preprocess(vm, and_clauses, xor_clauses, pins, planted)

    cnf_clauses, var_map, num_vars = fsat.to_cnf(vm, and_clauses, xor_clauses)

    # Build assignment: cnf_variable_index -> True/False
    assignment = {}
    for old_var, new_var in var_map.items():
        canon = vm.find(old_var)
        var = abs(canon)
        neg = (canon < 0)
        if var in planted:
            val = planted[var] ^ neg
        elif var in vm.pinned:
            val = vm.pinned[var] ^ neg
        else:
            val = False
        assignment[new_var] = bool(val)

    return assignment, num_vars, len(cnf_clauses)


# ---------------------------------------------------------------------------
# Verify clauses against assignment
# ---------------------------------------------------------------------------

def verify_clauses(clauses, assignment):
    """Check every clause against the assignment.

    Returns (all_satisfied, num_satisfied, first_failing_index).
    """
    num_satisfied = 0
    for i, clause in enumerate(clauses):
        satisfied = False
        for lit in clause:
            var = abs(lit)
            if var not in assignment:
                # Variable not in assignment — treat as unassigned;
                # this shouldn't happen for a well-formed instance
                continue
            val = assignment[var]
            if lit < 0:
                val = not val
            if val:
                satisfied = True
                break
        if satisfied:
            num_satisfied += 1
        else:
            return False, num_satisfied, i
    return True, num_satisfied, -1


# ---------------------------------------------------------------------------
# Main verification logic
# ---------------------------------------------------------------------------

def verify_file(filepath, verbose=False):
    """Verify a single CNF file. Returns True if verification passes."""
    basename = os.path.basename(filepath)

    # Step 1: extract p, q
    try:
        p, q = extract_pq_from_filename(filepath)
    except ValueError as e:
        print(f"SKIP  {basename}: {e}", file=sys.stderr)
        return None

    N = p * q

    if verbose:
        print(f"File:    {basename}")
        print(f"Factors: p={p}, q={q}, N={N}")

    # Step 2: read the CNF from disk
    num_vars_file, num_clauses_declared, clauses = read_dimacs(filepath)

    if verbose:
        print(f"DIMACS:  {num_vars_file} vars, {num_clauses_declared} clauses "
              f"(read {len(clauses)} clauses)")

    # Step 3: regenerate the planted assignment
    assignment, num_vars_gen, num_clauses_gen = build_planted_assignment(p, q)

    if verbose:
        print(f"Regen:   {num_vars_gen} vars, {num_clauses_gen} clauses")

    # Sanity check: regenerated instance should match file dimensions
    dims_match = (num_vars_file == num_vars_gen and
                  len(clauses) == num_clauses_gen)
    if not dims_match:
        print(f"WARNING  {basename}: dimension mismatch — "
              f"file has {num_vars_file} vars / {len(clauses)} clauses, "
              f"regenerated has {num_vars_gen} vars / {num_clauses_gen} clauses. "
              f"Verifying against file clauses anyway.", file=sys.stderr)

    # Step 4: verify every clause in the file
    all_ok, num_sat, fail_idx = verify_clauses(clauses, assignment)

    if all_ok:
        print(f"PASS  {basename}  "
              f"(p={p}, q={q}, N={N}, "
              f"{len(clauses)} clauses all satisfied)")
        return True
    else:
        print(f"FAIL  {basename}  "
              f"(p={p}, q={q}, N={N}, "
              f"clause {fail_idx} unsatisfied, "
              f"{num_sat}/{len(clauses)} satisfied)")
        if verbose:
            print(f"  Failing clause: {clauses[fail_idx]}")
            for lit in clauses[fail_idx]:
                var = abs(lit)
                val = assignment.get(var, "MISSING")
                print(f"    lit {lit:+d}: var {var} = {val}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Verify factorization SAT instances using p, q from filename.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s factor_0005_p203767_q767317.cnf
  %(prog)s *.cnf
  %(prog)s instances/*.cnf -v
        """)
    parser.add_argument("files", nargs="+", help="CNF file(s) to verify")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="print detailed information")
    args = parser.parse_args()

    passed = 0
    failed = 0
    skipped = 0

    for filepath in args.files:
        result = verify_file(filepath, verbose=args.verbose)
        if result is True:
            passed += 1
        elif result is False:
            failed += 1
        else:
            skipped += 1
        if args.verbose and len(args.files) > 1:
            print()

    # Summary for batch runs
    if len(args.files) > 1:
        total = passed + failed + skipped
        print(f"\nSummary: {passed} passed, {failed} failed, "
              f"{skipped} skipped out of {total}")

    sys.exit(1 if failed > 0 else 0)


if __name__ == "__main__":
    main()
