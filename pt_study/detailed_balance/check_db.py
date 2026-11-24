"""Compare per-temperature histograms from pt_cnf --mode sample with exact
Boltzmann distributions exp(-E/T)/Z obtained by enumerating all 2^n states."""
import sys, numpy as np, itertools
cnf, hist = sys.argv[1], sys.argv[2]
cl=[]; n=0
for line in open(cnf):
    if line.startswith('p'): n=int(line.split()[2]); continue
    xs=[int(t) for t in line.split()][:-1]
    if xs: cl.append(xs)
# energy of every state; state index bit (v-1) = x_v
idx=np.arange(2**n, dtype=np.int64)
X=((idx[:,None]>>np.arange(n))&1).astype(bool)
E=np.zeros(2**n,dtype=np.int64)
for c in cl:
    sat=np.zeros(2**n,dtype=bool)
    for l in c: sat|= X[:,l-1] if l>0 else ~X[:,-l-1]
    E+=~sat
worst=0
for line in open(hist):
    if not line.startswith('T'): continue
    parts=line.split(); T=float(parts[1]); h=np.array(parts[2:],dtype=float)
    p=h/h.sum(); w=np.exp(-(E-E.min())/T); q=w/w.sum()
    tv=0.5*np.abs(p-q).sum()
    levels=np.unique(E); pe=np.array([p[E==e].sum() for e in levels]); qe=np.array([q[E==e].sum() for e in levels])
    tvE=0.5*np.abs(pe-qe).sum()
    meanE_emp=(p*E).sum(); meanE_ex=(q*E).sum()
    worst=max(worst,tv)
    print(f"T={T:6.3f}  TV(states)={tv:.4f}  TV(energy)={tvE:.4f}  <E> sampled {meanE_emp:7.4f} exact {meanE_ex:7.4f}")
print(f"max TV over replicas: {worst:.4f}")
