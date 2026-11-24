import sys, json, random
import os; sys.path.insert(0, os.environ.get('PQ_REPO', '.'))  # path to pq-SAT-benchmark
import factorization_sat as F
from sympy import randprime
random.seed(2026)
rows=[]
for d in range(2,20):
    for t in range(5):
        if d==2: p,q=3,2 if t%2 else 3
        else:
            p=randprime(2**(d-1),2**d); q=randprime(2**(d-1),2**d)
        if d==2: p,q=(3,2) if t%2 else (3,3)
        vm,a,x,pins,pl,pv,qv=F.generate_clauses(p,q)
        pre_vars=len(vm.all_vars()); pre_cls=3*len(a)+4*len(x)
        F.preprocess(vm,a,x,pins,pl)
        cnf,vmap,nv=F.to_cnf(vm,a,x)
        rows.append(dict(d=d,p=p,q=q,pre_vars=pre_vars,pre_cls=pre_cls,post_vars=nv,post_cls=len(cnf),AND=len(a),XOR=len(x)))
    print(d, rows[-1], flush=True)
json.dump(rows,open('fig1.json','w'))
