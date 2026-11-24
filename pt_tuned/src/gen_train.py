"""Training instances: 4 per d, disjoint from the test set (d = 5 reuses test
pairs: only 10 distinct pairs of 5-bit primes exist).  CNF files for
d = 5..10 (residual CNF, same generator), factor list for d = 5..18."""
import sys, csv, random, os
sys.path.insert(0, '/home/claude/repo'); import factorization_sat as fs
from sympy import primerange
test = set()
for man in ['/home/claude/ptbench/small_manifest.tsv', '/home/claude/rerun/manifest.tsv']:
    for r in csv.DictReader(open(man), delimiter='\t'):
        test.add((int(r['p']), int(r['q'])))
rng = random.Random(4242)
os.makedirs('train_cnf', exist_ok=True)
rows = []
for d in range(5, 19):
    pr = list(primerange(2 ** (d - 1), 2 ** d))
    pairs = [(p, q) for i, p in enumerate(pr) for q in pr[i + 1:]]
    pool = [x for x in pairs if x not in test] if d > 5 else pairs
    chosen = []
    while len(chosen) < 4:
        x = rng.choice(pool) if len(pool) > 5000 else None
        if x is None:
            rng.shuffle(pool); chosen = pool[:4]; break
        if x not in chosen: chosen.append(x)
    for p, q in chosen:
        f = ''
        if d <= 10:
            vm, A, X, pins, planted, pv, qv = fs.generate_clauses(p, q); fs.preprocess(vm, A, X, pins, planted)
            cnf, vmap, nv = fs.to_cnf(vm, A, X)
            f = f'train_cnf/d{d:02d}_{p}x{q}.cnf'
            open(f, 'w').write(f"p cnf {nv} {len(cnf)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cnf))
        rows.append((d, p, q, p * q, f))
with open('train.tsv', 'w') as fo:
    fo.write('d\tp\tq\tN\tcnf\n')
    for r in rows: fo.write('\t'.join(map(str, r)) + '\n')
print(len(rows), 'training instances;', sum(1 for r in rows if r[4]), 'with CNF')
