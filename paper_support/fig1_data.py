"""Data for Fig. 1 (effect of preprocessing on instance size).

Run with PQ_REPO=<path to pq-SAT-benchmark>; writes fig1.json in the current
directory.  The prime draws are reproducible: they use a seeded
random.Random (seed 2026) and sympy.isprime, not sympy.randprime, whose
internal RNG is not controlled by random.seed().

For each d = 3..19 we make five ordered draws (p, q) of d-bit primes,
independently and uniformly over the d-bit primes, so squares (p = q) and
repeated or swapped pairs can occur.  At d = 2 the only 2-bit primes are 2
and 3; the five draws there are fixed to (3,3), (3,2), (3,3), (3,2), (3,3).
"""
import sys, json, random
import os; sys.path.insert(0, os.environ.get('PQ_REPO', '.'))  # path to pq-SAT-benchmark
import factorization_sat as F
from sympy import isprime

SEED = 2026
rng = random.Random(SEED)


def draw_prime(d):
    """Uniform random d-bit prime, by rejection sampling with the seeded rng."""
    while True:
        x = rng.randrange(2**(d-1), 2**d)
        if isprime(x):
            return x


rows=[]
for d in range(2,20):
    for t in range(5):
        if d==2:
            p,q=(3,2) if t%2 else (3,3)
        else:
            p=draw_prime(d); q=draw_prime(d)
        vm,a,x,pins,pl,pv,qv=F.generate_clauses(p,q)
        pre_vars=len(vm.all_vars()); pre_cls=3*len(a)+4*len(x)
        F.preprocess(vm,a,x,pins,pl)
        cnf,vmap,nv=F.to_cnf(vm,a,x)
        rows.append(dict(d=d,p=p,q=q,pre_vars=pre_vars,pre_cls=pre_cls,post_vars=nv,post_cls=len(cnf),AND=len(a),XOR=len(x)))
    print(d, rows[-1], flush=True)
json.dump(rows,open('fig1.json','w'))
