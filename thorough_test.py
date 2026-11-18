"""
Thorough test suite for factorization_sat.py

Tests:
1. Correctness: SAT solver finds valid factorization for many (p,q) pairs
2. Uniqueness: only valid factorizations satisfy the CNF
3. Edge cases: p=q, p=2, small primes, one factor much larger
4. Consistency: pre-reduction counts match analytical formula C(d)=d^2(d-1)^2/2
5. Preprocessing: reduced instance is equivalent to original
6. No-preprocess vs preprocess: both produce correct results
7. Planted solution actually satisfies every clause
8. Ising compilation: planted spins give zero penalty
9. Exhaustive check for small N: enumerate ALL satisfying assignments
"""
import sys
import traceback
from pysat.solvers import Glucose4
from factorization_sat import (
    generate_clauses, preprocess, to_cnf, to_ising, verify_ising,
    VariableManager
)

PASS = 0
FAIL = 0

def check(condition, msg):
    global PASS, FAIL
    if condition:
        PASS += 1
    else:
        FAIL += 1
        print(f"  FAIL: {msg}")

def recover_pq(p_vars, q_vars, vm, planted, var_map, assignment):
    """Extract p and q from a SAT assignment."""
    def get_val(pv):
        canon = vm.find(pv)
        var = abs(canon)
        neg = (canon < 0)
        if var in vm.pinned:
            return vm.pinned[var] ^ neg
        if var in var_map:
            new_v = var_map[var]
            if new_v in assignment:
                return assignment[new_v] ^ neg
        # Try finding through equivalences
        for old_v, new_v in var_map.items():
            c = vm.find(old_v)
            if abs(c) == var:
                val = assignment.get(new_v, False)
                if c < 0:
                    val = not val
                return val ^ neg
        return None

    p_bits = [get_val(pv) for pv in p_vars]
    q_bits = [get_val(qv) for qv in q_vars]

    if None in p_bits or None in q_bits:
        return None, None

    p_rec = sum((1 if b else 0) << i for i, b in enumerate(p_bits))
    q_rec = sum((1 if b else 0) << i for i, b in enumerate(q_bits))
    return p_rec, q_rec


def test_basic(p, q, preprocess_on=True, label=""):
    """Test that SAT solver finds valid factorization."""
    N = p * q
    tag = f"N={N}={p}x{q}" + (f" [{label}]" if label else "")

    try:
        vm, and_cl, xor_cl, pins, planted, pv, qv = generate_clauses(p, q)

        if preprocess_on:
            preprocess(vm, and_cl, xor_cl, pins, planted)
        else:
            for lit, val in pins:
                vm.pin(lit, val)

        cnf_clauses, var_map, num_vars = to_cnf(vm, and_cl, xor_cl)

        # Check planted solution satisfies all clauses
        inv_map = {v: k for k, v in var_map.items()}
        planted_assignment = {}
        for new_var in range(1, num_vars + 1):
            old_var = inv_map.get(new_var)
            if old_var is not None:
                canon = vm.find(old_var)
                var = abs(canon)
                neg = (canon < 0)
                if var in planted:
                    val = planted[var] ^ neg
                elif var in vm.pinned:
                    val = vm.pinned[var] ^ neg
                else:
                    val = False
                planted_assignment[new_var] = val

        for i, clause in enumerate(cnf_clauses):
            sat = False
            for lit in clause:
                v = abs(lit)
                val = planted_assignment.get(v, False)
                if lit < 0:
                    val = not val
                if val:
                    sat = True
                    break
            check(sat, f"{tag}: planted solution fails clause {i}: {clause}")
            if not sat:
                return False

        if num_vars == 0:
            check(True, f"{tag}: trivially solved")
            return True

        # Solve with SAT solver
        solver = Glucose4()
        for clause in cnf_clauses:
            solver.add_clause(clause)
        sat = solver.solve()
        check(sat, f"{tag}: instance is UNSAT")
        if not sat:
            solver.delete()
            return False

        model = solver.get_model()
        solver.delete()

        assignment = {abs(lit): (lit > 0) for lit in model}
        p_rec, q_rec = recover_pq(pv, qv, vm, planted, var_map, assignment)

        if p_rec is None or q_rec is None:
            check(False, f"{tag}: could not recover p,q from assignment")
            return False

        check(p_rec * q_rec == N,
              f"{tag}: product mismatch {p_rec}*{q_rec}={p_rec*q_rec} != {N}")
        return p_rec * q_rec == N

    except Exception as e:
        check(False, f"{tag}: exception: {e}")
        traceback.print_exc()
        return False


def test_exhaustive_small(p, q):
    """For small instances, enumerate ALL satisfying assignments and verify
    each one corresponds to a valid factorization of N."""
    N = p * q
    tag = f"exhaustive N={N}={p}x{q}"

    vm, and_cl, xor_cl, pins, planted, pv, qv = generate_clauses(p, q)
    preprocess(vm, and_cl, xor_cl, pins, planted)
    cnf_clauses, var_map, num_vars = to_cnf(vm, and_cl, xor_cl)

    if num_vars == 0:
        check(True, f"{tag}: trivially solved, no vars")
        return

    if num_vars > 20:
        print(f"  SKIP exhaustive for {tag}: {num_vars} vars too many")
        return

    # Enumerate all solutions using solver
    solver = Glucose4()
    for clause in cnf_clauses:
        solver.add_clause(clause)

    solutions = []
    while solver.solve():
        model = solver.get_model()
        assignment = {abs(lit): (lit > 0) for lit in model}
        solutions.append(assignment)
        # Block this solution
        blocking = [-lit for lit in model]
        solver.add_clause(blocking)

    solver.delete()

    check(len(solutions) > 0, f"{tag}: no solutions found")

    factorizations = set()
    for assignment in solutions:
        p_rec, q_rec = recover_pq(pv, qv, vm, planted, var_map, assignment)
        if p_rec is not None and q_rec is not None:
            if p_rec * q_rec == N:
                factorizations.add((min(p_rec, q_rec), max(p_rec, q_rec)))
            else:
                check(False, f"{tag}: solution gives {p_rec}*{q_rec}={p_rec*q_rec} != {N}")

    check((min(p,q), max(p,q)) in factorizations,
          f"{tag}: planted solution not among found factorizations")

    # Every factorization should be valid
    for (a, b) in factorizations:
        check(a * b == N, f"{tag}: invalid factorization {a}*{b}")

    print(f"  {tag}: {len(solutions)} total SAT assignments, "
          f"{len(factorizations)} distinct factorizations: {factorizations}")


def test_formula_match(d):
    """Verify pre-reduction counts match C(d) = d^2(d-1)^2/2."""
    from sympy import nextprime
    import random
    rng = random.Random(42 + d)
    lo = 1 << (d - 1)
    hi = (1 << d) - 1
    p = nextprime(rng.randint(lo, hi) - 1)
    if p > hi:
        p = nextprime(lo)
    q = nextprime(rng.randint(lo, hi) - 1)
    if q > hi or q == p:
        q = nextprime(p)

    n_p = p.bit_length()
    n_q = q.bit_length()

    vm, and_cl, xor_cl, pins, planted, pv, qv = generate_clauses(p, q)

    n_pp = n_p * n_q
    n_contr = len(xor_cl)
    n_and = len(and_cl)
    total_vars = len(vm.all_vars())

    # Check counts
    check(n_and == n_pp + n_contr,
          f"d={d}: AND clauses {n_and} != pp({n_pp}) + contr({n_contr})")
    check(len(xor_cl) == n_contr,
          f"d={d}: XOR clauses mismatch")
    check(total_vars == 2*max(n_p,n_q) + n_pp + 2*n_contr or
          total_vars == n_p + n_q + n_pp + 2*n_contr,
          f"d={d}: total vars {total_vars} unexpected")

    if n_p == n_q:
        C_pred = n_p**2 * (n_p - 1)**2 // 2
        check(n_contr == C_pred,
              f"d={d}: C={n_contr} != predicted {C_pred}")


def test_ising(p, q):
    """Test Ising compilation."""
    N = p * q
    tag = f"Ising N={N}={p}x{q}"

    vm, and_cl, xor_cl, pins, planted, pv, qv = generate_clauses(p, q)
    preprocess(vm, and_cl, xor_cl, pins, planted)

    if len(and_cl) == 0 and len(xor_cl) == 0:
        return  # trivially solved

    h, J, E0, spin_map, planted_spins, num_spins = to_ising(
        vm, and_cl, xor_cl, planted)

    ok, energy = verify_ising(h, J, E0, planted_spins, num_spins)
    check(ok, f"{tag}: Ising energy = {energy}, expected 0")

    # Also check that flipping any single spin increases energy
    for s in range(1, min(num_spins + 1, 50)):  # check first 50
        flipped = dict(planted_spins)
        flipped[s] = -flipped.get(s, 1)
        e_flip = E0
        for i, hi in h.items():
            e_flip += hi * flipped.get(i, 1)
        for (i, j), Jij in J.items():
            e_flip += Jij * flipped.get(i, 1) * flipped.get(j, 1)
        check(e_flip >= -1e-10,
              f"{tag}: flipping spin {s} gives energy {e_flip} < 0")


# ===========================================================================
# RUN ALL TESTS
# ===========================================================================

print("=" * 70)
print("TEST 1: Basic correctness (with preprocessing)")
print("=" * 70)
basic_cases = [
    (2, 3), (2, 5), (2, 7), (2, 11), (2, 13),
    (3, 5), (3, 7), (3, 11), (3, 13), (3, 17),
    (5, 7), (5, 11), (5, 13), (5, 17), (5, 19),
    (7, 11), (7, 13), (7, 17), (7, 19), (7, 23),
    (11, 13), (11, 17), (11, 19), (11, 23), (11, 29),
    (13, 17), (13, 19), (13, 23), (13, 29), (13, 31),
    (23, 29), (29, 31), (31, 37), (37, 41), (41, 43),
    (101, 103), (127, 131), (251, 257),
]
for p, q in basic_cases:
    test_basic(p, q, preprocess_on=True)

print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
print("TEST 2: Without preprocessing")
print("=" * 70)
for p, q in [(3, 5), (7, 11), (11, 13), (23, 29), (101, 103)]:
    test_basic(p, q, preprocess_on=False, label="no-preprocess")

print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
print("TEST 3: Edge case p = q (perfect square)")
print("=" * 70)
for p in [3, 5, 7, 11, 13]:
    test_basic(p, p, label="p=q")

print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
print("TEST 4: Edge case p = 2")
print("=" * 70)
for q in [3, 5, 7, 11, 13, 17, 101]:
    test_basic(2, q, label="p=2")

print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
print("TEST 5: Asymmetric sizes")
print("=" * 70)
for p, q in [(3, 127), (5, 251), (7, 1009), (11, 127), (3, 1021)]:
    test_basic(p, q, label="asymmetric")

print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
print("TEST 6: Formula C(d) = d^2(d-1)^2/2 verification")
print("=" * 70)
for d in range(2, 16):
    test_formula_match(d)

print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
print("TEST 7: Ising compilation")
print("=" * 70)
for p, q in [(5, 7), (7, 11), (11, 13), (13, 17), (23, 29), (101, 103)]:
    test_ising(p, q)

print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
print("TEST 8: Exhaustive (all solutions) for small instances")
print("=" * 70)
for p, q in [(3, 5), (3, 7), (5, 7), (7, 11), (7, 13), (5, 11)]:
    test_exhaustive_small(p, q)


def test_tail_forced_zero(p, q):
    """Columns k >= n_p + n_q can only hold zeros: any 1 there would make the
    weighted column sum exceed p'q' < 2^(n_p+n_q).  Preprocessing must find
    this, i.e. pin every such variable to False."""
    vm, and_cl, xor_cl, pins, planted, pv, qv = generate_clauses(p, q)
    width = p.bit_length() + q.bit_length()
    tail = [v for v, k in vm.column.items() if k >= width]
    preprocess(vm, and_cl, xor_cl, pins, planted)
    check(all(vm.value_of(v) is False for v in tail),
          f"tail not eliminated for {p}x{q}")
    # The residual system must be fully simplified: the only fixed literal
    # allowed in a residual relation is the 0-output of an AND, which leaves
    # the binary constraint NOT(x AND y) on two free inputs.  (That output
    # literal may be represented by any variable of the class pinned to 0,
    # including a tail variable.)
    ok = True
    for (o, a, b) in and_cl:
        if vm.value_of(a) is not None or vm.value_of(b) is not None:
            ok = False
        if vm.value_of(o) not in (None, False):
            ok = False
    for c in xor_cl:
        if any(vm.value_of(l) is not None for l in c):
            ok = False
    check(ok, f"residual relation contains a fixed literal for {p}x{q}")


print()
print("=" * 70)
print("TEST 9: High-order columns are eliminated by preprocessing")
print("=" * 70)
import random as _r
from sympy import nextprime as _np
_rng = _r.Random(2026)          # fixed seed: same instances on every run
def _prime_bits(nbits):
    while True:
        cand = _np(_rng.randint(2 ** (nbits - 1), 2 ** nbits - 1))
        if cand.bit_length() == nbits:
            return cand
for _a in range(2, 13):
    for _b in range(2, 13):         # symmetric and asymmetric widths
        test_tail_forced_zero(_prime_bits(_a), _prime_bits(_b))
print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print()
print("=" * 70)
print("TEST 10: Combined half-adder gadget (exhaustive spectrum, small cases)")
print("=" * 70)
import itertools as _it
def test_combined_gadget(p, q):
    """With gadget='combined' the planted state must be the only zero-energy
    state, every other state must have H >= 2, all coefficients are integers,
    and the spin count must not exceed that of the separate gadgets."""
    vm, a_cl, x_cl, pins, planted, pv, qv = generate_clauses(p, q)
    preprocess(vm, a_cl, x_cl, pins, planted)
    n_sep = to_ising(vm, a_cl, x_cl, planted)[5]
    h, J, E0, smap, ps, n = to_ising(vm, a_cl, x_cl, planted, gadget="combined")
    check(n <= n_sep, f"combined gadget has more spins for {p}x{q}")
    check(all(float(v).is_integer() for v in list(h.values()) + list(J.values()) + [E0]),
          f"non-integer coefficient for {p}x{q}")
    if n_sep > 20:
        return
    import numpy as _npy
    def _spectrum(h_, J_, E0_, n_):
        S = _npy.array(list(_it.product((-1, 1), repeat=n_)), dtype=_npy.int8).reshape(2 ** n_, n_)
        E = _npy.full(len(S), float(E0_))
        for i, v in h_.items():
            E += v * S[:, i - 1]
        for (i, j), v in J_.items():
            E += v * S[:, i - 1] * S[:, j - 1]
        return S, E
    S, E = _spectrum(h, J, E0, n)
    zero = _npy.abs(E) < 1e-9
    check(E.min() > -1e-9, f"negative energy for {p}x{q}")
    check(_npy.all(E[~zero] >= 2 - 1e-9), f"excited state below 2 for {p}x{q}")
    planted_row = _npy.array([ps[i] for i in range(1, n + 1)], dtype=_npy.int8)
    check(n == 0 or bool(_npy.any(_npy.all(S[zero] == planted_row, axis=1))),
          f"planted state not a ground state for {p}x{q}")
    # Same number of ground states as the separate-gadget model (= models of
    # the residual CNF; the swap degeneracy can sit in a variable that no
    # residual relation mentions, in which case both counts are 1).
    hs, Js, E0s, _, _, ns_ = to_ising(vm, a_cl, x_cl, planted)
    _, Es = _spectrum(hs, Js, E0s, ns_)
    check(int(zero.sum()) == int((_npy.abs(Es) < 1e-9).sum()),
          f"ground-state count differs between gadgets for {p}x{q}")
for _p, _q in [(2, 2), (2, 3), (3, 5), (5, 7), (5, 11), (7, 13), (5, 19), (7, 17),
               (5, 5), (11, 13), (13, 17), (7, 31), (151, 227)]:
    test_combined_gadget(_p, _q)
print(f"\nPassed so far: {PASS}, Failed: {FAIL}")

print("\n" + "=" * 70)
if FAIL == 0:
    print(f"ALL {PASS} TESTS PASSED")
else:
    print(f"FAILED: {FAIL} out of {PASS + FAIL} tests")
    sys.exit(1)
print("=" * 70)
