#!/usr/bin/env python3
"""
make_dwave_instances.py -- build the Ising instances of the quantum
annealing study (Sec. 5.5 / Table 3 of paper_v17) with the released
generator, and check each one before it goes to the QPU.

For every prime pair below it writes, into ./dwave_instances/:
  <p>x<q>_ising_E<E0>.tsv   fields and couplings, "i<TAB>j<TAB>value"
                            (i == j is a field h_i, i < j a coupling J_ij);
                            spins are numbered 1..n and include the XOR
                            auxiliary spins
  <p>x<q>_planted.txt       planted value (+1/-1) of every spin, "i s_i"
and a summary table instances.tsv with the spin count as run in the
summer study, the spin count now, the number of couplings and E0.

Energy convention: with s_i = +-1, the stored Hamiltonian
    H(s) = sum_i h_i s_i + sum_{i<j} J_ij s_i s_j
has its minimum at the planted state, where H = -E0.  A read is a success
when its energy equals -E0 (and, as a second check, when it agrees with
the planted file spin by spin).  Every pair below has factors of different
bit lengths or p = q, so the ground state is unique.

With --gadget combined, each XOR relation is compiled together with its
half-adder AND partner as (1/2)(s_a + s_b - s_s - 2 s_t - 1)^2, which needs no
auxiliary spin.  Same problems, same unique planted solutions, fewer spins.
These are for the gadget comparison (item 2b of the instructions).

Usage:
    python3 make_dwave_instances.py
    python3 make_dwave_instances.py \
        --gadget combined --out dwave_instances_combined
"""
import argparse
import itertools
import os
import sys

# (p, q, spins as run in summer 2026), the 20 instances of Table 3
PAIRS = [
    (2, 2, 5), (5, 11, 10), (5, 19, 17), (7, 13, 18), (5, 5, 23),
    (7, 17, 23), (7, 41, 41), (5, 23, 43), (5, 43, 48), (7, 29, 49),
    (5, 13, 70), (7, 31, 76), (5, 17, 78), (7, 19, 92), (7, 23, 101),
    (11, 17, 113), (7, 43, 115), (7, 37, 128), (7, 47, 147), (13, 17, 152),
]


def energy(h, J, s):
    return sum(v * s[i] for i, v in h.items()) + sum(v * s[i] * s[j] for (i, j), v in J.items())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--generator",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."),
                    help="folder containing the released factorization_sat.py "
                         "(default: the repository root)")
    ap.add_argument("--out", default="dwave_instances")
    ap.add_argument("--gadget", choices=("separate", "combined"), default="separate",
                    help="separate = released benchmark (default); combined = "
                         "half-adder penalty without XOR auxiliaries")
    args = ap.parse_args()
    sys.path.insert(0, os.path.abspath(args.generator))
    import factorization_sat as fsat

    os.makedirs(args.out, exist_ok=True)
    rows = []
    for p, q, old_spins in PAIRS:
        vm, and_cl, xor_cl, pins, planted, pv, qv = fsat.generate_clauses(p, q)
        fsat.preprocess(vm, and_cl, xor_cl, pins, planted)
        h, J, E0, spin_map, planted_spins, n = fsat.to_ising(vm, and_cl, xor_cl, planted,
                                                              gadget=args.gadget)
        E0i = int(round(E0))
        stem = os.path.join(args.out, f"{p}x{q}")
        with open(f"{stem}_ising_E{E0i}.tsv", "w") as f:
            fsat.write_ising(h, J, E0, n, planted_spins, file=f)
        with open(f"{stem}_planted.txt", "w") as f:
            fsat.write_planted(planted_spins, n, file=f)

        # checks: planted state has energy -E0, and no single spin flip lowers it
        s = {i: planted_spins.get(i, 1) for i in range(1, n + 1)}
        e = energy(h, J, s)
        assert abs(e + E0) < 1e-9, (p, q, e, E0)
        for i in range(1, n + 1):
            s[i] = -s[i]
            assert energy(h, J, s) > e, (p, q, "single flip lowers energy")
            s[i] = -s[i]
        # exhaustive ground-state check for small instances
        if n <= 22:
            best = min(energy(h, J, dict(zip(range(1, n + 1), c)))
                       for c in itertools.product((1, -1), repeat=n))
            assert abs(best - e) < 1e-9, (p, q, "planted state is not the ground state")

        rows.append((f"{p}x{q}", p * q, old_spins, n, sum(1 for v in J.values() if v),
                     E0i))
        print(f"{p:>3} x {q:<3} N={p*q:<5} spins as run {old_spins:>4}  now {n:>4}  "
              f"couplings {rows[-1][4]:>4}  E0 {E0i}")

    with open(os.path.join(args.out, "instances.tsv"), "w") as f:
        f.write("instance\tN\tspins_as_run\tspins_now\tcouplings\tE0\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
