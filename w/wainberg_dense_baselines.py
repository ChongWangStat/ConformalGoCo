from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import binom, norm
from math import exp, log, sqrt

ROOT=Path('/mnt/data/conformalgoco_analysis_2026-09-06')
W=ROOT/'wainberg'; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
# conservative to liberal
thr=np.array(sorted(list(range(600,1001,25))+[1050,1100,1200,1400,1600,2000,2500,3000], reverse=True),float)
assert len(thr)==25
# analysis units from network_communities
genes=pd.read_csv(W/'network_communities.csv')['gene'].astype(str).tolist(); gi={g:i for i,g in enumerate(genes)}
# truth sets
D=pd.read_csv(W/'go_truth_direct.csv.gz',usecols=['gene','go_id'])
T=pd.read_csv(W/'go_truth_true_path.csv.gz',usecols=['gene','go_id'])
D['key']=D['gene'].astype(str)+'\t'+D['go_id'].astype(str); T['key']=T['gene'].astype(str)+'\t'+T['go_id'].astype(str)
ds=set(D['key']); ts=set(T['key']); del D,T
# load only dense-grid-relevant scores
parts=[]
for ch in pd.read_csv(W/'wainberg_gene_go_scores.csv.gz',usecols=['gene','go_id','score'],chunksize=500000):
    q=ch[ch.score>=thr.min()].copy()
    if len(q): parts.append(q)
s=pd.concat(parts,ignore_index=True); del parts
s['key']=s.gene.astype(str)+'\t'+s.go_id.astype(str)
s['direct']=s.key.isin(ds).astype(np.int8); s['truepath']=s.key.isin(ts).astype(np.int8)
# matrices
n=len(genes); k=len(thr)
lossD=np.zeros((n,k),dtype=np.float64); lossT=np.zeros((n,k),dtype=np.float64)
yieldD=np.zeros((n,k),dtype=np.int32); yieldT=np.zeros((n,k),dtype=np.int32); calls=np.zeros((n,k),dtype=np.int32)
# aggregate per threshold from 17k rows
s['i']=s['gene'].map(gi)
print('dense rows outside eligible units', int(s['i'].isna().sum()))
s=s[s['i'].notna()].copy(); s['i']=s['i'].astype(int)
for j,c in enumerate(thr):
    q=s[s.score>=c]
    g=q.groupby('i',sort=False).agg(calls=('go_id','size'),direct=('direct','sum'),truepath=('truepath','sum'))
    idx=g.index.to_numpy(dtype=int); cc=g.calls.to_numpy(); dd=g.direct.to_numpy(); tt=g.truepath.to_numpy()
    calls[idx,j]=cc; yieldD[idx,j]=dd; yieldT[idx,j]=tt
    lossD[idx,j]=1-dd/cc; lossT[idx,j]=1-tt/cc
np.savez_compressed(OUT/'wainberg_dense_matrices.npz', thresholds=thr, genes=np.array(genes), loss_direct=lossD, loss_truepath=lossT,yield_direct=yieldD,yield_truepath=yieldT,calls=calls)

def kl(a,b):
    eps=np.finfo(float).eps; a=np.clip(a,eps,1-eps); b=np.clip(b,eps,1-eps)
    return a*np.log(a/b)+(1-a)*np.log((1-a)/(1-b))
def p_hoeff(x,a=.1):
    r=x.mean(); d=max(a-r,0); return min(1.,np.exp(-2*len(x)*d*d))
def p_hb(x,a=.1):
    r=x.mean(); n=len(x)
    if r>=a:return 1.
    return min(1., np.exp(-n*kl(r,a)), np.e*binom.cdf(int(np.ceil(n*r)),n,a))
def p_norm(x,a=.1):
    r=x.mean(); v=np.sum((x-r)**2)/(len(x)**2)
    return (0. if r<a else 1.) if v<=0 else norm.cdf((r-a)/np.sqrt(v))
def eb_ucb(x,delta):
    n=len(x)
    if n<2:return 1.
    sd=x.std(ddof=1); z=np.log(2/delta)
    return min(1.,x.mean()+sd*np.sqrt(2*z/n)+7*z/(3*(n-1)))
def p_eb(x,a=.1):
    if x.mean()>=a:return 1.
    lo,hi=1e-15,1-1e-12
    if eb_ucb(x,hi)>a:return 1.
    for _ in range(100):
        mid=(lo+hi)/2
        if eb_ucb(x,mid)<=a: hi=mid
        else:lo=mid
    return hi
def select(mat,mask,pfun,delta=.1):
    ps=[pfun(mat[mask,j]) for j in range(mat.shape[1])]
    sel=None
    for j,p in enumerate(ps):
        if p<=delta: sel=j
        else: break
    return sel,ps

def splits(style='default_rng'):
    for seed in range(100):
        if style=='default_rng': rng=np.random.default_rng(seed); idx=rng.permutation(n)
        elif style=='RandomState': rng=np.random.RandomState(seed); idx=rng.permutation(n)
        elif style=='seedplus': rng=np.random.default_rng(20260905+seed); idx=rng.permutation(n)
        mask=np.zeros(n,bool); mask[idx[:round(.7*n)]]=True
        yield seed,mask

def eval_style(style):
    rows=[]
    for lossname,mat,y in [('Direct',lossD,yieldD),('TruePath',lossT,yieldT)]:
      for seed,mask in splits(style):
        ev=~mask
        for name,fun in [('IID-Hoeffding',p_hoeff),('HB',p_hb),('EmpBern',p_eb),('IID-Normal',p_norm)]:
          j,ps=select(mat,mask,fun)
          if j is None: rows.append([style,lossname,seed,name,np.nan,np.nan,np.nan,0]); continue
          rows.append([style,lossname,seed,name,thr[j],mat[ev,j].mean(),y[ev,j].sum(),1])
    return pd.DataFrame(rows,columns=['split_style','loss','seed','method','threshold','heldout_fdp','yield','certified'])
allr=pd.concat([eval_style(s) for s in ['default_rng','RandomState','seedplus']],ignore_index=True)
allr.to_csv(OUT/'wainberg_dense_iid_baselines_splitstyle_probe.csv',index=False)
summary=allr.groupby(['split_style','loss','method'],as_index=False).agg(mean_fdp=('heldout_fdp','mean'),mean_yield=('yield','mean'),gt10=('heldout_fdp',lambda x: np.mean(x>.1)),median_thr=('threshold','median'),cert=('certified','mean'))
summary.to_csv(OUT/'wainberg_dense_iid_baselines_splitstyle_summary.csv',index=False)
print('score rows >=600',len(s))
print(summary.to_string(index=False))
