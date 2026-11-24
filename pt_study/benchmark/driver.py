#!/usr/bin/env python3
"""Run PT-CNF (uniform, focused), PT-factor and Kissat on d = 5..18 instances.
Bit sizes in increasing order; a method is dropped after a d at which a
majority of its runs time out.  Results appended to results.tsv."""
import csv, os, subprocess, sys, time, glob, re
from concurrent.futures import ThreadPoolExecutor
TO = 300
KISSAT = "/home/claude/rerun/kissat"
rows = []
def add(man, root):
    seen = set()
    for r in csv.DictReader(open(man), delimiter='\t'):
        key = (r['d'], r['p'], r['q'])
        if key in seen: continue
        seen.add(key)
        path = r['cnf'].replace('CNFROOT/', 'cnf/')
        if not os.path.isabs(path):
            for base in root:
                if os.path.exists(os.path.join(base, path)): path = os.path.join(base, path); break
        rows.append(dict(d=int(r['d']), p=int(r['p']), q=int(r['q']), N=int(r['N']), cnf=path))
add('small_manifest.tsv', ['s5', 's69'])
add('/home/claude/rerun/manifest.tsv', ['/home/claude/rerun'])
rows = [r for r in rows if 5 <= r['d'] <= 18]
assert all(os.path.exists(r['cnf']) for r in rows), "missing cnf"
out = 'results.tsv'
done = set()
if os.path.exists(out):
    for l in open(out):
        x = l.split('\t'); done.add((x[0], x[1], x[2], x[3]))
else:
    open(out, 'w').write('method\td\tp\tq\tstatus\ttime\textra\n')

def run(method, r):
    if method == 'kissat':
        t0 = time.time()
        pr = subprocess.run(['timeout', str(TO), KISSAT, '-q', r['cnf']], capture_output=True, text=True)
        el = time.time() - t0
        st = 'SAT' if 's SATISFIABLE' in pr.stdout else ('TIMEOUT' if pr.returncode == 124 else 'ERR%d' % pr.returncode)
        return st, el, ''
    if method in ('pt_factor', 'random_factor'):
        np_ = r['p'].bit_length(); nq = r['q'].bit_length()
        extra = [] if method == 'pt_factor' else ['--replicas', '1', '--tmin', '1e12', '--tmax', '1e12']
        o = subprocess.run(['/home/claude/pt/pt_factor', '--N', str(r['N']), '--np', str(np_), '--nq', str(nq),
                            '--time', str(TO), '--seed', '1'] + extra, capture_output=True, text=True).stdout
    else:
        o = subprocess.run(['./run_pt.sh', r['cnf'], method.split('_')[-1], str(TO), '1'], capture_output=True, text=True).stdout
    m = re.search(r'result (\S+) time (\S+)', o)
    if not m: return 'ERR', 0.0, o.strip()[:80]
    extra = re.search(r'(evals \d+|bestE \d+)', o)
    return m.group(1), float(m.group(2)), extra.group(1) if extra else ''

methods = sys.argv[1:] or ['kissat', 'pt_factor', 'pt_uniform', 'pt_focused']
for method in methods:
    for d in sorted({r['d'] for r in rows}):
        batch = [r for r in rows if r['d'] == d and (method, str(d), str(r['p']), str(r['q'])) not in done]
        with ThreadPoolExecutor(2) as ex:
            res = list(ex.map(lambda r: (r, run(method, r)), batch))
        with open(out, 'a') as f:
            for r, (st, el, ex_) in res:
                f.write(f"{method}\t{d}\t{r['p']}\t{r['q']}\t{st}\t{el:.4f}\t{ex_}\n")
        allr = [l.split('\t') for l in open(out)][1:]
        mine = [x for x in allr if x[0] == method and int(x[1]) == d]
        nto = sum(x[4] != 'SAT' for x in mine)
        print(f"{method} d={d}: {len(mine)} runs, {nto} unsolved", flush=True)
        if nto > len(mine) / 2:
            print(f"{method}: dropped after d={d}", flush=True); break
