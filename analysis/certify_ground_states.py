#!/usr/bin/env python3
"""Requires python-sat (pip install python-sat).

Certify (1) uniqueness of the satisfying assignment / Ising ground state and
(2) the exact Ising gap Delta = 2, for released pq instances.

Ground states of H = sum of gadgets (all >= 0) are exactly the satisfying
assignments of the residual AND/XOR relations with every XOR auxiliary at its
planted value, so (1) reduces to counting CNF models.  For (2): every violated
XOR gadget can be brought to energy 2 by its auxiliary, every other non-zero
gadget value is >= 4, so H = 2 is attainable iff for some XOR relation g the
system "all relations except g hold, g is violated" is satisfiable.  Each
witness found is checked by evaluating the compiled Ising energy directly.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import factorization_sat as fs
from pysat.solvers import Cadical153


def build(p, q):
    vm, A, X, pins, planted, pv, qv = fs.generate_clauses(p, q)
    fs.preprocess(vm, A, X, pins, planted)
    cnf, var_map, nv = fs.to_cnf(vm, A, X)
    h, J, E0, spin_map, ps, ns = fs.to_ising(vm, A, X, planted)
    assert var_map == spin_map
    return vm, A, X, cnf, nv, h, J, E0, spin_map, ps, ns


def lit(vm, var_map, l):
    c = vm.find(l); v = abs(c); neg = c < 0
    if v in vm.pinned:
        return ('c', vm.pinned[v] ^ neg)
    return ('l', -var_map[v] if neg else var_map[v])


def rel_clauses(vm, vmap, kind, out, a, b, violate=False):
    """CNF for out = a op b (or its negation) with constants folded."""
    L = [lit(vm, vmap, x) for x in (a, b, out)]
    cls = []
    for va in (0, 1):
        for vb in (0, 1):
            vo = (va & vb) if kind == 'and' else (va ^ vb)
            for vout in (0, 1):
                ok = (vout == vo)
                if ok != violate:
                    continue
                # forbid the assignment (va, vb, vout): clause of negations
                clause, sat = [], False
                for (t, x), val in zip(L, (va, vb, vout)):
                    if t == 'c':
                        if x != val:
                            sat = True
                    else:
                        clause.append(-x if val else x)
                if sat:
                    continue
                cls.append(clause)
    return cls


def energy(h, J, E0, s):
    return E0 + sum(v * s[i] for i, v in h.items()) + sum(v * s[i] * s[j] for (i, j), v in J.items())


def certify(p, q):
    vm, A, X, cnf, nv, h, J, E0, vmap, ps, ns = build(p, q)
    # (1) uniqueness
    with Cadical153(bootstrap_with=cnf) as S:
        assert S.solve()
        m = S.get_model()
        planted_lits = [i if ps[i] == 1 else -i for i in range(1, nv + 1)]
        S.add_clause([-x for x in planted_lits])
        unique = not S.solve()
    # (2) gap = 2 witness
    base = []
    for (o, a, b) in A:
        base += rel_clauses(vm, vmap, 'and', o, a, b)
    witness = None
    for gi, (o, a, b) in enumerate(X):
        cls = list(base)
        for gj, (o2, a2, b2) in enumerate(X):
            cls += rel_clauses(vm, vmap, 'xor', o2, a2, b2, violate=(gj == gi))
        if any(len(c) == 0 for c in cls):
            continue
        with Cadical153(bootstrap_with=cls) as S:
            if S.solve():
                m = S.get_model()
                s = {abs(x): (1 if x > 0 else -1) for x in m if abs(x) <= nv}
                for i in range(1, nv + 1):
                    s.setdefault(i, 1)
                # auxiliaries: index nv+1+k for k-th XOR relation
                def val(l):
                    t, x = lit(vm, vmap, l)
                    return (1 if x else -1) if t == 'c' else (s[abs(x)] * (1 if x > 0 else -1))
                for k in range(len(X)):
                    s[nv + 1 + k] = -1
                for k, (o2, a2, b2) in enumerate(X):
                    best = None
                    for sa in (-1, 1):
                        s[nv + 1 + k] = sa
                        e = energy(h, J, E0, s)
                        if best is None or e < best[0]:
                            best = (e, sa)
                    s[nv + 1 + k] = best[1]
                witness = (gi, energy(h, J, E0, s))
                break
    return ns, unique, witness


if __name__ == "__main__":
    pairs = [tuple(map(int, a.split('x'))) for a in sys.argv[1:]]
    for p, q in pairs:
        ns, unique, w = certify(p, q)
        print(f"{p:>3}x{q:<3} spins {ns:4d}  unique ground state: {unique}  "
              f"H=2 witness: {'none' if w is None else 'XOR #%d, H = %g' % w}")
