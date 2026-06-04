#!/usr/bin/env python3
"""
factorization_sat.py — Factorization-based planted-solution SAT/Ising benchmark generator.

Given two integers p and q, constructs SAT and Ising instances encoding
the binary factorization N = p * q.

Pipeline:
  1. Clause generation from binary multiplication (AND + XOR constraints)
  2. Boolean preprocessing (iterative logical reduction)
  3. Output as DIMACS CNF and Ising TSV (both by default)

Usage:
  python factorization_sat.py 7 11                     # both CNF + Ising to files
  python factorization_sat.py 7 11 -o inst             # inst.cnf + inst_ising_E<E0>.tsv + ...
  python factorization_sat.py 7 11 --cnf-only          # CNF to stdout
  python factorization_sat.py 7 11 --ising-only        # Ising TSV to stdout
  python factorization_sat.py 7 11 --no-preprocess     # skip preprocessing
  python factorization_sat.py 7 11 --stats             # print statistics only

Reference:
  Hen, "Factorization-based planted-solution Ising benchmarks" (2025/2026).
"""

import argparse
import os
import sys
from collections import defaultdict, deque


# =============================================================================
# Variable manager with union-find for equivalence classes
# =============================================================================

class VariableManager:
    """Manages Boolean variables with union-find equivalence tracking.

    Variables are positive integers.  Literals are signed integers:
    +v means variable v, -v means NOT(variable v).

    Two variables can be merged as equal (v == w) or anti-equal (v == NOT w).
    After merging, find(v) returns the canonical literal for v.
    """

    def __init__(self):
        self._parent = {}   # var -> (parent_var, flipped)
        self._rank = {}
        self._next_id = 1
        # Pinned values: canonical_var -> bool
        self.pinned = {}

    def new_var(self):
        v = self._next_id
        self._next_id += 1
        self._parent[v] = (v, False)
        self._rank[v] = 0
        return v

    def _find_root(self, v):
        """Returns (root_var, flipped) where flipped indicates whether v
        is negated relative to root."""
        flip = False
        cur = v
        while True:
            pv, pf = self._parent[cur]
            if pv == cur:
                break
            flip ^= pf
            cur = pv
        return cur, flip

    def find(self, lit):
        """Resolve literal to canonical form.  Returns a literal (signed int)."""
        neg = (lit < 0)
        v = abs(lit)
        if v not in self._parent:
            self._parent[v] = (v, False)
            self._rank[v] = 0
        root, flip = self._find_root(v)
        flip ^= neg
        return -root if flip else root

    def merge_equal(self, lit_a, lit_b):
        """Assert that lit_a == lit_b."""
        ca = self.find(lit_a)
        cb = self.find(lit_b)
        if ca == cb:
            return True  # already known
        if ca == -cb:
            return False  # contradiction
        va, fa = abs(ca), (ca < 0)
        vb, fb = abs(cb), (cb < 0)
        # Union by rank
        ra, rb = self._rank.get(va, 0), self._rank.get(vb, 0)
        need_flip = (fa != fb)  # if both same sign -> same; if different -> flip
        if ra < rb:
            va, vb = vb, va
            # flip stays the same
        self._parent[vb] = (va, need_flip)
        if ra == rb:
            self._rank[va] = ra + 1
        return True

    def merge_anti(self, lit_a, lit_b):
        """Assert that lit_a == NOT(lit_b)."""
        return self.merge_equal(lit_a, -lit_b)

    def pin(self, lit, value):
        """Pin a literal to a Boolean value.
        Returns True if new info, False if redundant, raises ValueError on contradiction."""
        canon = self.find(lit)
        var = abs(canon)
        val = value ^ (canon < 0)
        if var in self.pinned:
            if self.pinned[var] != val:
                raise ValueError(f"Contradiction: var {var} pinned to both True and False")
            return False  # no new info
        self.pinned[var] = val
        return True  # new pin

    def value_of(self, lit):
        """If the literal is pinned, return its Boolean value; else None."""
        canon = self.find(lit)
        var = abs(canon)
        if var not in self.pinned:
            return None
        val = self.pinned[var]
        if canon < 0:
            val = not val
        return val

    def all_vars(self):
        """Return all variable IDs that have been created."""
        return set(self._parent.keys())


# =============================================================================
# Stage 1: Clause generation from binary multiplication
# =============================================================================

def generate_clauses(p, q):
    """Generate AND/XOR/pin constraints from the multiplication N = p * q.

    Returns:
        vm: VariableManager with all variables
        and_clauses: list of (out_lit, in1_lit, in2_lit) meaning out = in1 AND in2
        xor_clauses: list of (out_lit, in1_lit, in2_lit) meaning out = in1 XOR in2
        pins: list of (lit, bool_value)
        planted: dict var -> bool_value for the planted solution
        p_vars: list of variable IDs for bits of p
        q_vars: list of variable IDs for bits of q
    """
    N = p * q
    bits_p = []
    tmp = p
    while tmp > 0:
        bits_p.append(tmp & 1)
        tmp >>= 1
    bits_q = []
    tmp = q
    while tmp > 0:
        bits_q.append(tmp & 1)
        tmp >>= 1
    bits_N = []
    tmp = N
    while tmp > 0:
        bits_N.append(tmp & 1)
        tmp >>= 1

    n_p = len(bits_p)
    n_q = len(bits_q)

    vm = VariableManager()
    and_clauses = []
    xor_clauses = []
    pins = []
    planted = {}  # var -> bool

    # Create variables for p and q bits
    p_vars = [vm.new_var() for _ in range(n_p)]
    q_vars = [vm.new_var() for _ in range(n_q)]
    for i, v in enumerate(p_vars):
        planted[v] = bool(bits_p[i])
    for j, v in enumerate(q_vars):
        planted[v] = bool(bits_q[j])

    # Partial products: a_{ij} = p_i AND q_j, placed in column i+j
    columns = defaultdict(deque)  # column_index -> deque of literals
    for i in range(n_p):
        for j in range(n_q):
            a_ij = vm.new_var()
            and_clauses.append((a_ij, p_vars[i], q_vars[j]))
            planted[a_ij] = bool(bits_p[i] and bits_q[j])
            columns[i + j].append(a_ij)

    # Column contraction: pairwise half-adder (XOR for sum, AND for carry)
    max_col = max(columns.keys()) if columns else 0
    # Process columns; carries can extend beyond max_col
    k = 0
    while True:
        if k not in columns or len(columns[k]) == 0:
            if k > max_col + 1:
                break
            k += 1
            continue
        col = columns[k]
        if len(col) <= 1:
            k += 1
            continue
        # Pop two entries and contract
        x = col.popleft()
        y = col.popleft()
        # sum = x XOR y (stays in column k)
        s = vm.new_var()
        xor_clauses.append((s, x, y))
        planted[s] = bool(planted.get(abs(x), False) ^ planted.get(abs(y), False))
        col.append(s)
        # carry = x AND y (goes to column k+1)
        c = vm.new_var()
        and_clauses.append((c, x, y))
        planted[c] = bool(planted.get(abs(x), False) and planted.get(abs(y), False))
        columns[k + 1].append(c)
        if k + 1 > max_col:
            max_col = k + 1
        # Don't advance k; there may be more entries in this column

    # Pinning: each column's single remaining entry must match N_k
    for k in sorted(columns.keys()):
        col = columns[k]
        if len(col) == 1:
            lit = col[0]
            if k < len(bits_N):
                pins.append((lit, bool(bits_N[k])))
            else:
                pins.append((lit, False))
        elif len(col) == 0:
            pass
        else:
            raise RuntimeError(f"Column {k} has {len(col)} entries after contraction")

    return vm, and_clauses, xor_clauses, pins, planted, p_vars, q_vars


# =============================================================================
# Stage 2: Boolean preprocessing
# =============================================================================

def preprocess(vm, and_clauses, xor_clauses, pins, planted, max_iters=1000):
    """Iterative Boolean reduction.

    Propagates pins, simplifies AND/XOR clauses, merges equivalent variables.
    Modifies lists in place and returns residual (and_clauses, xor_clauses).

    Raises ValueError if a contradiction is detected (indicating a bug in
    the deduction rules or an unsatisfiable instance).
    """
    # Apply initial pins
    for lit, val in pins:
        vm.pin(lit, val)

    def canon_lit(lit):
        return vm.find(lit)

    def lit_val(lit):
        return vm.value_of(lit)

    _preprocess_loop(vm, and_clauses, xor_clauses, canon_lit, lit_val, max_iters)

    return and_clauses, xor_clauses


def _preprocess_loop(vm, and_clauses, xor_clauses, canon_lit, lit_val, max_iters):
    for iteration in range(max_iters):
        progress = False

        # --- Simplify AND clauses ---
        new_and = []
        for (out, a, b) in and_clauses:
            out, a, b = canon_lit(out), canon_lit(a), canon_lit(b)
            vo, va, vb = lit_val(out), lit_val(a), lit_val(b)

            # Case: output is known
            if vo is not None:
                if vo:  # out = True => a = True AND b = True
                    if vm.pin(a, True): progress = True
                    if vm.pin(b, True): progress = True
                    continue
                else:  # out = False => at least one input is False
                    if va is True:  # a=T, out=F => b=F
                        if vm.pin(b, False): progress = True
                        continue
                    elif vb is True:  # b=T, out=F => a=F
                        if vm.pin(a, False): progress = True
                        continue
                    elif va is False or vb is False:
                        continue  # already satisfied
                    elif a == -b:  # a AND NOT(a) = False, already known
                        continue
                    # else: can't resolve further, but out=0 is still useful
                    # keep as reduced clause
                    new_and.append((out, a, b))
                    continue

            # Case: an input is known
            if va is not None or vb is not None:
                if va is False or vb is False:  # out = False
                    if vm.pin(out, False): progress = True
                    continue
                if va is True and vb is True:
                    if vm.pin(out, True): progress = True
                    continue
                if va is True:  # out = b
                    if vm.merge_equal(out, b): progress = True
                    continue
                if vb is True:  # out = a
                    if vm.merge_equal(out, a): progress = True
                    continue

            # Case: inputs are equal or anti-equal
            if a == b:  # out = a AND a = a
                if vm.merge_equal(out, a): progress = True
                continue
            if a == -b:  # out = a AND NOT(a) = False
                if vm.pin(out, False): progress = True
                continue

            # Case: out equals one of the inputs
            if out == a:  # a = a AND b => a implies b (a <= b), so a => b=T or a=F
                # a = a AND b means: if a=T then b=T; if a=F then OK.
                # This is equivalent to: NOT(a) OR b, i.e., a => b
                # Can't resolve without more info, keep it
                new_and.append((out, a, b))
                continue
            if out == b:
                new_and.append((out, a, b))
                continue

            new_and.append((out, a, b))

        and_clauses[:] = new_and

        # --- Simplify XOR clauses ---
        new_xor = []
        for (out, a, b) in xor_clauses:
            out, a, b = canon_lit(out), canon_lit(a), canon_lit(b)
            vo, va, vb = lit_val(out), lit_val(a), lit_val(b)

            # Count known values
            known = sum(x is not None for x in [vo, va, vb])

            if known >= 2:
                # Can determine the third
                if vo is not None and va is not None:
                    # b = out XOR a
                    if vm.pin(b, vo ^ va): progress = True
                    continue
                if vo is not None and vb is not None:
                    if vm.pin(a, vo ^ vb): progress = True
                    continue
                if va is not None and vb is not None:
                    if vm.pin(out, va ^ vb): progress = True
                    continue

            if known == 1:
                if vo is not None:
                    if vo:  # out=T => a XOR b = T => a = NOT b
                        if vm.merge_anti(a, b): progress = True
                    else:   # out=F => a = b
                        if vm.merge_equal(a, b): progress = True
                    continue
                if va is not None:
                    if va:  # a=T => out = NOT b
                        if vm.merge_anti(out, b): progress = True
                    else:   # a=F => out = b
                        if vm.merge_equal(out, b): progress = True
                    continue
                if vb is not None:
                    if vb:
                        if vm.merge_anti(out, a): progress = True
                    else:
                        if vm.merge_equal(out, a): progress = True
                    continue

            # Check for equal/anti-equal variables
            if a == b:  # out = a XOR a = 0
                if vm.pin(out, False): progress = True
                continue
            if a == -b:  # out = a XOR NOT(a) = 1
                if vm.pin(out, True): progress = True
                continue
            if out == a:  # a = a XOR b => b = 0
                if vm.pin(b, False): progress = True
                continue
            if out == -a:  # NOT(a) = a XOR b => b = 1
                if vm.pin(b, True): progress = True
                continue
            if out == b:
                if vm.pin(a, False): progress = True
                continue
            if out == -b:
                if vm.pin(a, True): progress = True
                continue

            new_xor.append((out, a, b))

        xor_clauses[:] = new_xor

        # --- Cross-clause inference: XOR pairs sharing variables ---
        # For each pair of XOR clauses sharing two literals, derive new info
        # (This is the symmetric-difference / Gaussian elimination step)
        if len(xor_clauses) < 500:  # only for manageable sizes
            xor_by_var = defaultdict(list)
            for idx, (out, a, b) in enumerate(xor_clauses):
                for lit in [out, a, b]:
                    xor_by_var[abs(lit)].append(idx)

            inferred = []
            seen_pairs = set()
            for var, indices in xor_by_var.items():
                for i in range(len(indices)):
                    for j in range(i + 1, len(indices)):
                        idx_i, idx_j = indices[i], indices[j]
                        pair = (min(idx_i, idx_j), max(idx_i, idx_j))
                        if pair in seen_pairs:
                            continue
                        seen_pairs.add(pair)
                        if idx_i >= len(xor_clauses) or idx_j >= len(xor_clauses):
                            continue
                        c1 = set(xor_clauses[idx_i])
                        c2 = set(xor_clauses[idx_j])
                        # XOR clauses: out = a XOR b means {out, a, b} with
                        # even parity.  Symmetric difference of two such sets
                        # (after accounting for signs) gives a new XOR relation.
                        # This is complex with signed literals; skip for now.

        if not progress:
            break

    return and_clauses, xor_clauses


# =============================================================================
# Stage 3a: Convert to DIMACS CNF
# =============================================================================

def to_cnf(vm, and_clauses, xor_clauses):
    """Convert residual AND and XOR clauses to CNF clauses.

    AND clause (out = a AND b) in CNF:
        (NOT a OR NOT b OR out)      — if both inputs T, output must be T
        (a OR NOT out)               — if output T, input a must be T
        (b OR NOT out)               — if output T, input b must be T

    XOR clause (out = a XOR b) in CNF:
        (NOT out OR NOT a OR NOT b)  — not all three true
        (NOT out OR a OR b)          — if out=T, at least one input T
        (out OR NOT a OR b)          — if out=F, inputs must match
        (out OR a OR NOT b)          — if out=F, inputs must match

    Returns:
        cnf_clauses: list of lists of signed integers (DIMACS-style)
        var_map: dict old_var -> new_contiguous_var (1-indexed)
        num_vars: int
    """
    # Collect all active (unpinned) variables
    active_lits = set()
    for (out, a, b) in and_clauses:
        for lit in [out, a, b]:
            active_lits.add(abs(vm.find(lit)))
    for (out, a, b) in xor_clauses:
        for lit in [out, a, b]:
            active_lits.add(abs(vm.find(lit)))

    # Remove pinned variables from active set
    active_vars = sorted(v for v in active_lits if v not in vm.pinned)
    var_map = {v: i + 1 for i, v in enumerate(active_vars)}
    num_vars = len(active_vars)

    def map_lit(lit):
        """Map a literal to the new variable space, substituting pinned values."""
        canon = vm.find(lit)
        var = abs(canon)
        neg = (canon < 0)
        if var in vm.pinned:
            # Return a trivially true/false indicator
            val = vm.pinned[var] ^ neg
            return ('const', val)
        if var not in var_map:
            # Variable was eliminated by equivalence; shouldn't happen after find()
            raise ValueError(f"Variable {var} not in active set")
        new_var = var_map[var]
        return ('lit', -new_var if neg else new_var)

    cnf_clauses = []

    def add_clause(lits_raw):
        """Add a clause, handling constant literals."""
        clause = []
        for lr in lits_raw:
            if lr[0] == 'const':
                if lr[1]:  # True constant -> clause is trivially satisfied
                    return
                else:  # False constant -> literal contributes nothing
                    continue
            clause.append(lr[1])
        if clause:
            cnf_clauses.append(clause)
        # Empty clause = unsatisfiable (shouldn't happen with correct input)

    # AND clauses: out = a AND b
    for (out, a, b) in and_clauses:
        o = map_lit(out)
        x = map_lit(a)
        y = map_lit(b)
        # Negate helper
        def neg(lr):
            if lr[0] == 'const':
                return ('const', not lr[1])
            return ('lit', -lr[1])
        # (NOT a OR NOT b OR out)
        add_clause([neg(x), neg(y), o])
        # (a OR NOT out)
        add_clause([x, neg(o)])
        # (b OR NOT out)
        add_clause([y, neg(o)])

    # XOR clauses: out = a XOR b
    for (out, a, b) in xor_clauses:
        o = map_lit(out)
        x = map_lit(a)
        y = map_lit(b)
        def neg(lr):
            if lr[0] == 'const':
                return ('const', not lr[1])
            return ('lit', -lr[1])
        # Four clauses for XOR:
        # (NOT out OR NOT a OR NOT b)
        add_clause([neg(o), neg(x), neg(y)])
        # (NOT out OR a OR b)
        add_clause([neg(o), x, y])
        # (out OR NOT a OR b)
        add_clause([o, neg(x), y])
        # (out OR a OR NOT b)
        add_clause([o, x, neg(y)])

    cnf_clauses = cleanup_cnf(cnf_clauses, num_vars)

    return cnf_clauses, var_map, num_vars


def cleanup_cnf(cnf_clauses, num_vars):
    """Post-process CNF: remove redundancies and propagate easy deductions.

    Iteratively applies:
      1. Deduplicate literals within each clause
      2. Remove tautological clauses (x OR NOT x)
      3. Remove duplicate clauses
      4. Unit propagation (single-literal clauses force variable values)
      5. Pure literal elimination (variables appearing in only one polarity)
      6. Subsumption (clause A subsumes clause B if A ⊂ B; remove B)
    """
    changed = True
    while changed:
        changed = False

        # 1. Deduplicate literals & remove tautologies
        new_clauses = []
        for clause in cnf_clauses:
            lits = set(clause)
            # Tautology check: x and -x in same clause
            if any(-l in lits for l in lits):
                changed = True
                continue
            deduped = sorted(lits, key=lambda x: (abs(x), x))
            if len(deduped) < len(clause):
                changed = True
            new_clauses.append(deduped)
        cnf_clauses = new_clauses

        # 2. Remove duplicate clauses
        seen = set()
        new_clauses = []
        for clause in cnf_clauses:
            key = tuple(clause)
            if key not in seen:
                seen.add(key)
                new_clauses.append(clause)
            else:
                changed = True
        cnf_clauses = new_clauses

        # 3. Subsumption (fast): only check unit and binary clauses as subsumers
        # Unit clause {L} subsumes any clause containing L.
        # Binary clause {L1, L2} subsumes any clause containing both L1 and L2.
        unit_set = set()
        binary_list = []
        for clause in cnf_clauses:
            if len(clause) == 1:
                unit_set.add(clause[0])
            elif len(clause) == 2:
                binary_list.append(frozenset(clause))

        if unit_set or binary_list:
            # Index: for each literal, which binary clauses contain it
            binary_set = set(binary_list)
            lit_to_binaries = {}
            for bc in binary_set:
                for lit in bc:
                    lit_to_binaries.setdefault(lit, []).append(bc)

            new_clauses = []
            for clause in cnf_clauses:
                cs = frozenset(clause)
                # Check unit subsumption
                if any(l in unit_set for l in clause) and len(clause) > 1:
                    changed = True
                    continue
                # Check binary subsumption
                if len(clause) > 2:
                    subsumed = False
                    for l in clause:
                        for bc in lit_to_binaries.get(l, []):
                            if bc.issubset(cs):
                                subsumed = True
                                break
                        if subsumed:
                            break
                    if subsumed:
                        changed = True
                        continue
                new_clauses.append(clause)
            cnf_clauses = new_clauses

    return cnf_clauses


def write_dimacs(cnf_clauses, num_vars, file=sys.stdout, comments=None):
    """Write CNF in DIMACS format."""
    if comments:
        for line in comments:
            file.write(f"c {line}\n")
    file.write(f"p cnf {num_vars} {len(cnf_clauses)}\n")
    for clause in cnf_clauses:
        file.write(" ".join(str(l) for l in clause) + " 0\n")


# =============================================================================
# Stage 3b: Compile to Ising Hamiltonian
# =============================================================================

def to_ising(vm, and_clauses, xor_clauses, planted):
    """Compile residual clauses into an Ising Hamiltonian H = E0 + sum h_i s_i + sum J_ij s_i s_j.

    Returns:
        h: dict var -> field
        J: dict (i,j) -> coupling (i < j)
        E0: float constant offset
        spin_map: dict old_var -> new_spin_index
        planted_spins: dict new_spin_index -> +1/-1
        num_spins: int
    """
    # Collect active variables
    active_lits = set()
    for (out, a, b) in and_clauses:
        for lit in [out, a, b]:
            active_lits.add(abs(vm.find(lit)))
    for (out, a, b) in xor_clauses:
        for lit in [out, a, b]:
            active_lits.add(abs(vm.find(lit)))

    active_vars = sorted(v for v in active_lits if v not in vm.pinned)
    spin_map = {v: i + 1 for i, v in enumerate(active_vars)}
    # XOR clauses need auxiliary spins
    aux_base = len(active_vars) + 1
    num_spins = len(active_vars)

    h = defaultdict(float)
    J = defaultdict(float)
    E0 = 0.0
    planted_spins = {}

    # Map planted solution
    for v in active_vars:
        canon = vm.find(v)
        var = abs(canon)
        neg = (canon < 0)
        if var in planted:
            val = planted[var] ^ neg
        else:
            val = False
        planted_spins[spin_map[v]] = 1 if val else -1

    def resolve_lit(lit):
        """Returns (spin_index, sign) or ('const', +1/-1) for pinned variables.
        sign = +1 if literal is positive, -1 if negated.
        """
        canon = vm.find(lit)
        var = abs(canon)
        neg = (canon < 0)
        if var in vm.pinned:
            val = vm.pinned[var] ^ neg
            return ('const', 1 if val else -1)  # Ising value
        si = spin_map[var]
        sign = -1 if neg else 1
        return ('spin', si, sign)

    def add_ising_term(coeff, *factors):
        """Add a term coeff * product(factors) to the Hamiltonian.
        Each factor is either a constant or (spin_index, sign)."""
        resolved = []
        c = coeff
        for f in factors:
            if f[0] == 'const':
                c *= f[1]
            else:
                resolved.append((f[1], f[2]))  # (spin_idx, sign)

        if len(resolved) == 0:
            E0_ref[0] += c
        elif len(resolved) == 1:
            si, sg = resolved[0]
            h[si] += c * sg
        elif len(resolved) == 2:
            si, sgi = resolved[0]
            sj, sgj = resolved[1]
            if si == sj:
                E0_ref[0] += c * sgi * sgj  # s_i^2 = 1
            else:
                i, j = (min(si, sj), max(si, sj))
                # Need to handle signs correctly
                J[(i, j)] += c * sgi * sgj
        else:
            raise ValueError("Ising gadgets should not produce terms with >2 spins")

    E0_ref = [0.0]

    # AND gadget: E = 3 - s1 - s2 + s1*s2 + 2*s3 - 2*s1*s3 - 2*s2*s3
    # where s1, s2 are inputs, s3 is output
    # Note: the gadget in the paper has a sign pattern:
    # E(s1,s2,s3) = 3 - s1 - s2 + s1*s2 - 2(-1 + s1 + s2)*s3
    #             = 3 - s1 - s2 + s1*s2 + 2*s3 - 2*s1*s3 - 2*s2*s3
    for (out, a, b) in and_clauses:
        s1 = resolve_lit(a)
        s2 = resolve_lit(b)
        s3 = resolve_lit(out)

        add_ising_term(3.0)        # constant
        add_ising_term(-1.0, s1)   # -s1
        add_ising_term(-1.0, s2)   # -s2
        add_ising_term(1.0, s1, s2)  # +s1*s2
        add_ising_term(2.0, s3)    # +2*s3
        add_ising_term(-2.0, s1, s3)  # -2*s1*s3
        add_ising_term(-2.0, s2, s3)  # -2*s2*s3

    # XOR gadget: needs auxiliary spin s_a
    # E(s1,s2,s3,sa) = 4 - s1 - s2 + s1*s2 + s3 - s1*s3 - s2*s3
    #                  + 2*sa - 2*s1*sa - 2*s2*sa + 2*s3*sa
    for (out, a, b) in xor_clauses:
        s1 = resolve_lit(a)
        s2 = resolve_lit(b)
        s3 = resolve_lit(out)

        # Create auxiliary spin
        num_spins += 1
        sa_idx = num_spins
        # Planted value: sa = +1 if s1 = s2 = +1, else sa = -1
        s1_planted = _planted_ising_val(a, vm, planted)
        s2_planted = _planted_ising_val(b, vm, planted)
        planted_spins[sa_idx] = 1 if (s1_planted == 1 and s2_planted == 1) else -1
        sa = ('spin', sa_idx, 1)

        add_ising_term(4.0)
        add_ising_term(-1.0, s1)
        add_ising_term(-1.0, s2)
        add_ising_term(1.0, s1, s2)
        add_ising_term(1.0, s3)
        add_ising_term(-1.0, s1, s3)
        add_ising_term(-1.0, s2, s3)
        add_ising_term(2.0, sa)
        add_ising_term(-2.0, s1, sa)
        add_ising_term(-2.0, s2, sa)
        add_ising_term(2.0, s3, sa)

    E0 = E0_ref[0]
    return dict(h), dict(J), E0, spin_map, planted_spins, num_spins


def _planted_ising_val(lit, vm, planted):
    """Get the planted Ising value (+1/-1) for a literal."""
    canon = vm.find(lit)
    var = abs(canon)
    neg = (canon < 0)
    if var in vm.pinned:
        val = vm.pinned[var] ^ neg
    elif var in planted:
        val = planted[var] ^ neg
    else:
        val = False
    return 1 if val else -1


def write_ising(h, J, E0, num_spins, planted_spins, file=sys.stdout, comments=None):
    """Write Ising instance as TSV triplets (i\\tj\\tvalue).

    Diagonal entries (i==j) encode local fields h_i.
    Off-diagonal entries (i<j) encode couplings J_ij.
    Pure data — no comment or header lines.
    The ground-state energy E0 is not stored in the file;
    it should be encoded in the filename by the caller.
    """
    for i in sorted(h.keys()):
        if h[i] != 0:
            file.write(f"{i}\t{i}\t{h[i]:g}\n")
    for (i, j) in sorted(J.keys()):
        if J[(i, j)] != 0:
            file.write(f"{i}\t{j}\t{J[(i,j)]:g}\n")


def write_planted(planted_spins, num_spins, file=sys.stdout):
    """Write planted spin assignment."""
    for i in range(1, num_spins + 1):
        s = planted_spins.get(i, 1)
        file.write(f"{i} {s}\n")


# =============================================================================
# Statistics
# =============================================================================

def print_stats(p, q, vm, and_clauses, xor_clauses, cnf_clauses, num_vars,
                ising_data=None, pre_and=None, pre_xor=None, pre_vars=None):
    """Print summary statistics."""
    N = p * q
    n_p = len(bin(p)) - 2
    n_q = len(bin(q)) - 2

    print(f"=== Factorization SAT instance: N = {N} = {p} x {q} ===")
    print(f"  p = {p} ({n_p} bits), q = {q} ({n_q} bits)")
    print(f"  N = {N} ({len(bin(N))-2} bits)")
    print()
    if pre_and is not None:
        print(f"  Pre-reduction:")
        print(f"    Boolean variables:  {pre_vars}")
        print(f"    AND clauses:        {pre_and}")
        print(f"    XOR clauses:        {pre_xor}")
        print(f"    Total constraints:  {pre_and + pre_xor}")
        print()
    print(f"  Post-reduction (residual):")
    print(f"    AND clauses:        {len(and_clauses)}")
    print(f"    XOR clauses:        {len(xor_clauses)}")
    n_pinned = len(vm.pinned)
    print(f"    Variables pinned:    {n_pinned}")
    print(f"    Free variables:     {num_vars}")
    print()
    print(f"  DIMACS CNF:")
    print(f"    Variables:          {num_vars}")
    print(f"    Clauses:            {len(cnf_clauses)}")
    if ising_data:
        h, J, E0, _, planted_spins, num_spins = ising_data
        n_fields = sum(1 for v in h.values() if v != 0)
        n_couplings = sum(1 for v in J.values() if v != 0)
        print()
        print(f"  Ising Hamiltonian:")
        print(f"    Spins:              {num_spins}")
        print(f"    Nonzero fields:     {n_fields}")
        print(f"    Nonzero couplings:  {n_couplings}")
        print(f"    Constant offset E0: {E0}")


# =============================================================================
# Verification
# =============================================================================

def verify_cnf(cnf_clauses, var_map, planted, vm):
    """Verify that the planted solution satisfies all CNF clauses."""
    # Build assignment: new_var -> True/False
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
        assignment[new_var] = val

    for i, clause in enumerate(cnf_clauses):
        satisfied = False
        for lit in clause:
            var = abs(lit)
            if var in assignment:
                val = assignment[var]
                if lit < 0:
                    val = not val
                if val:
                    satisfied = True
                    break
        if not satisfied:
            return False, i
    return True, -1


def verify_ising(h, J, E0, planted_spins, num_spins):
    """Verify planted spins achieve energy E0 (zero penalty)."""
    energy = E0
    for i, hi in h.items():
        energy += hi * planted_spins.get(i, 1)
    for (i, j), Jij in J.items():
        energy += Jij * planted_spins.get(i, 1) * planted_spins.get(j, 1)
    return abs(energy) < 1e-10, energy


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate SAT/Ising instances from integer factorization.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s 7 11                       # DIMACS CNF to stdout
  %(prog)s 7 11 -o inst.cnf           # CNF to file
  %(prog)s 7 11 --no-preprocess       # skip Boolean preprocessing
  %(prog)s 7 11 --ising -o inst.txt   # Ising format
  %(prog)s 7 11 --stats               # statistics only
  %(prog)s 7 11 --all -o inst         # all outputs: inst.cnf, inst_ising.txt, inst_planted.txt
        """)
    parser.add_argument("p", type=int, help="first factor")
    parser.add_argument("q", type=int, help="second factor")
    parser.add_argument("-o", "--output", help="output file prefix (default: factor_PxQ)")
    parser.add_argument("--no-preprocess", action="store_true",
                        help="skip Boolean preprocessing")
    parser.add_argument("--cnf-only", action="store_true",
                        help="output only DIMACS CNF (skip Ising)")
    parser.add_argument("--ising-only", action="store_true",
                        help="output only Ising TSV (skip CNF)")
    parser.add_argument("--stats", action="store_true",
                        help="print statistics only, no instance output")
    parser.add_argument("--verify", action="store_true",
                        help="verify planted solution satisfies the instance")
    args = parser.parse_args()

    p, q = args.p, args.q
    if p < 2 or q < 2:
        print("Error: p and q must be >= 2", file=sys.stderr)
        sys.exit(1)

    N = p * q

    # Stage 1: generate clauses
    vm, and_clauses, xor_clauses, pins, planted, p_vars, q_vars = generate_clauses(p, q)

    pre_and = len(and_clauses)
    pre_xor = len(xor_clauses)
    pre_vars = len(vm.all_vars())

    # Stage 2: preprocessing
    if not args.no_preprocess:
        preprocess(vm, and_clauses, xor_clauses, pins, planted)
    else:
        # Still apply pins
        for lit, val in pins:
            vm.pin(lit, val)

    # Stage 3: convert to CNF
    cnf_clauses, var_map, num_vars = to_cnf(vm, and_clauses, xor_clauses)

    comments = [
        f"Factorization SAT instance: N = {N} = {p} x {q}",
        f"p = {p} ({len(bin(p))-2} bits), q = {q} ({len(bin(q))-2} bits)",
        f"Preprocessing: {'off' if args.no_preprocess else 'on'}",
        f"Residual AND clauses: {len(and_clauses)}, XOR clauses: {len(xor_clauses)}",
    ]

    # Determine what to output
    do_cnf = not args.ising_only
    do_ising = not args.cnf_only

    # Ising compilation
    ising_data = None
    if do_ising:
        ising_data = to_ising(vm, and_clauses, xor_clauses, planted)

    # Verification (always performed — this is the planted-solution guarantee)
    ok_cnf, bad_idx = verify_cnf(cnf_clauses, var_map, planted, vm)
    if not ok_cnf:
        print(f"INTERNAL ERROR: planted solution fails clause {bad_idx}. "
              f"This is a bug.", file=sys.stderr)
        sys.exit(2)
    if ising_data:
        h, J, E0, _, planted_spins, num_spins = ising_data
        ok_ising, energy = verify_ising(h, J, E0, planted_spins, num_spins)
        if not ok_ising:
            print(f"INTERNAL ERROR: planted spins give energy {energy} != 0. "
                  f"This is a bug.", file=sys.stderr)
            sys.exit(2)
    if args.verify:
        print("Verification: PASSED (planted solution satisfies all clauses)",
              file=sys.stderr)
        if ising_data:
            print(f"Verification: PASSED (Ising energy = {energy})",
                  file=sys.stderr)

    # Stats
    if args.stats or (do_cnf and do_ising):
        print_stats(p, q, vm, and_clauses, xor_clauses, cnf_clauses, num_vars,
                    ising_data=ising_data, pre_and=pre_and, pre_xor=pre_xor,
                    pre_vars=pre_vars)
        if args.stats:
            return

    # Output
    prefix = args.output or f"factor_{p}x{q}"

    if do_cnf:
        if args.output or do_ising:
            # Writing to file
            cnf_path = f"{prefix}.cnf"
            with open(cnf_path, "w") as f:
                write_dimacs(cnf_clauses, num_vars, file=f, comments=comments)
            print(f"Wrote {cnf_path}", file=sys.stderr)
        else:
            # CNF-only with no -o: write to stdout
            write_dimacs(cnf_clauses, num_vars, comments=comments)

    if do_ising:
        h, J, E0, _, planted_spins, num_spins = ising_data
        if args.output or do_cnf:
            # Writing to file
            ising_path = f"{prefix}_ising_E{int(E0)}.tsv"
            with open(ising_path, "w") as f:
                write_ising(h, J, E0, num_spins, planted_spins, file=f)
            print(f"Wrote {ising_path}", file=sys.stderr)
            # Planted spins
            with open(f"{prefix}_planted.txt", "w") as f:
                write_planted(planted_spins, num_spins, file=f)
            print(f"Wrote {prefix}_planted.txt", file=sys.stderr)
        else:
            # Ising-only with no -o: write to stdout
            write_ising(h, J, E0, num_spins, planted_spins)

    if do_cnf and do_ising:
        # Also write planted CNF solution
        inv_map = {v: k for k, v in var_map.items()}
        with open(f"{prefix}_cnf_solution.txt", "w") as f:
            f.write(f"c Planted SAT solution for N = {N} = {p} x {q}\n")
            f.write(f"c Format: variable_index value(0/1)\n")
            f.write(f"c\n")
            f.write(f"c Metadata for factorization verification:\n")
            f.write(f"p {p}\n")
            f.write(f"q {q}\n")
            f.write(f"N {N}\n")
            f.write(f"p_bits {len(p_vars)}\n")
            for i, pv in enumerate(p_vars):
                canon = vm.find(pv)
                var = abs(canon)
                neg = (canon < 0)
                if var in vm.pinned:
                    val = vm.pinned[var] ^ neg
                    f.write(f"p_bit {i} pinned {1 if val else 0}\n")
                elif var in var_map:
                    cnf_var = var_map[var]
                    f.write(f"p_bit {i} var {cnf_var} {'neg' if neg else 'pos'}\n")
                else:
                    val = planted.get(var, False) ^ neg
                    f.write(f"p_bit {i} pinned {1 if val else 0}\n")
            f.write(f"q_bits {len(q_vars)}\n")
            for j, qv in enumerate(q_vars):
                canon = vm.find(qv)
                var = abs(canon)
                neg = (canon < 0)
                if var in vm.pinned:
                    val = vm.pinned[var] ^ neg
                    f.write(f"q_bit {j} pinned {1 if val else 0}\n")
                elif var in var_map:
                    cnf_var = var_map[var]
                    f.write(f"q_bit {j} var {cnf_var} {'neg' if neg else 'pos'}\n")
                else:
                    val = planted.get(var, False) ^ neg
                    f.write(f"q_bit {j} pinned {1 if val else 0}\n")
            f.write(f"c\n")
            f.write(f"c Planted assignment (all CNF variables):\n")
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
                    f.write(f"{new_var} {1 if val else 0}\n")
            print(f"Wrote {prefix}_cnf_solution.txt", file=sys.stderr)


if __name__ == "__main__":
    main()
