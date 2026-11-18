#!/usr/bin/env python3
"""Requires python-sat (pip install python-sat).

Exact Ising gap (first excited energy) and ground-state count for released
pq instances, by weighted MaxSAT over the gadget tables.  Each gadget's energy
is encoded exactly as soft clauses forbidding each of its nonzero-energy
tuples, so the MaxSAT cost of an assignment equals its Ising energy H."""
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from certify_ground_states import build, lit, energy
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
from pysat.solvers import Cadical153

def Ex(s1,s2,s3,sa): return 4-s1-s2+s1*s2+s3-s1*s3-s2*s3+2*sa-2*s1*sa-2*s2*sa+2*s3*sa
def Ea(s1,s2,s3): return 3-s1-s2+s1*s2+2*s3-2*s1*s3-2*s2*s3

def run(p, q):
    vm, A, X, cnf, nv, h, J, E0, vmap, ps, ns = build(p, q)
    # count ground states (CNF models)
    gs = []
    with Cadical153(bootstrap_with=cnf) as S:
        while S.solve():
            m = [x for x in S.get_model() if abs(x) <= nv]
            gs.append(m); S.add_clause([-x for x in m])
            if len(gs) > 4: break
    W = WCNF()
    def soft_for(L, Efun, nlits):
        for vals in itertools.product([-1, 1], repeat=nlits):
            e = Efun(*vals)
            if e == 0: continue
            clause, sat = [], False
            for (t, x), v in zip(L, vals):
                if t == 'c':
                    if (1 if x else -1) != v: sat = True
                else:
                    clause.append(-x if v == 1 else x)  # forbid literal x having spin v
            if sat: continue
            if clause: W.append(clause, weight=e)
            else: W.append([nv + 10**6], weight=e); W.append([-(nv + 10**6)])  # constant cost
    for (o, a, b) in A:
        soft_for([lit(vm, vmap, a), lit(vm, vmap, b), lit(vm, vmap, o)], Ea, 3)
    for k, (o, a, b) in enumerate(X):
        soft_for([lit(vm, vmap, a), lit(vm, vmap, b), lit(vm, vmap, o), ('l', nv + 1 + k)], Ex, 4)
    # exclude every ground state on the physical spins
    for m in gs:
        W.append([-x for x in m])
    with RC2(W) as r:
        m = r.compute()
        cost = r.cost
    s = {abs(x): (1 if x > 0 else -1) for x in m if abs(x) <= ns}
    for i in range(1, ns + 1): s.setdefault(i, 1)
    e = energy(h, J, E0, s)
    return ns, len(gs), cost, e

if __name__ == "__main__":
    for a in sys.argv[1:]:
        p, q = map(int, a.split('x'))
        ns, ngs, cost, e = run(p, q)
        print(f"{p:>3}x{q:<3} spins {ns:4d}  ground states {ngs}  gap (MaxSAT) {cost}  H(witness) {e:g}")
