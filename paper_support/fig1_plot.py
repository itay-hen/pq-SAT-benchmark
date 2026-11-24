import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=json.load(open('fig1.json'))
ds=sorted({r['d'] for r in R})
BLUE,ORANGE,INK,MUTED='#2a78d6','#eb6834','#0b0b0b','#8a8984'
plt.rcParams.update({'font.size':10,'axes.edgecolor':MUTED,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
fig,axs=plt.subplots(1,2,figsize=(9.2,3.6),constrained_layout=True)
for ax,key,lab,tag in [(axs[0],'cls','Number of CNF clauses','(a)'),(axs[1],'vars','Number of Boolean variables','(b)')]:
    pre=np.array([np.median([r['pre_'+key] for r in R if r['d']==d]) for d in ds])
    post=[[r['post_'+key] for r in R if r['d']==d] for d in ds]
    med=np.array([np.median(v) for v in post]); lo=np.array([min(v) for v in post]); hi=np.array([max(v) for v in post])
    D=np.array(ds,float)
    m=D>=4
    ax.plot(D[m],pre[m],'-s',color=BLUE,ms=5,lw=1.6,label='Before preprocessing')
    ax.errorbar(D[m],med[m],yerr=[med[m]-lo[m],hi[m]-med[m]],fmt='-o',color=ORANGE,ms=5,lw=1.6,capsize=2.5,label='After preprocessing (median, range)')
    xs=np.linspace(4,19,100)
    ax.plot(xs,pre[-1]*(xs/19)**4,'--',color=MUTED,lw=1,label=r'$\propto d^4$')
    ax.plot(xs,med[-1]*(xs/19)**3,':',color=MUTED,lw=1.3,label=r'$\propto d^3$')
    ax.set_yscale('log'); ax.set_xlabel('Bit-length $d$'); ax.set_ylabel(lab)
    ax.set_xticks(range(4,20,2)); ax.grid(True,which='major',color='#e6e5e0',lw=0.6); ax.set_axisbelow(True)
    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.text(0.03,0.95,tag,transform=ax.transAxes,fontweight='bold',va='top')
axs[1].legend(loc='lower right',frameon=False,fontsize=8.5)
fig.savefig('preprocessing_reduction.pdf'); fig.savefig('preprocessing_reduction.png',dpi=200)
# print table for text
for d in ds:
    v=[r for r in R if r['d']==d]
    print(d, int(np.median([r['pre_cls'] for r in v])), int(np.median([r['post_cls'] for r in v])), int(np.median([r['post_vars'] for r in v])), round(np.median([r['post_vars'] for r in v])/d**3,2), round(1-np.median([r['post_cls'] for r in v])/np.median([r['pre_cls'] for r in v]),3))
