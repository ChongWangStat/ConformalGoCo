from __future__ import annotations
import numpy as np, pandas as pd, scipy.sparse as sp, math, hashlib
from scipy.stats import norm


def build_source_fraction_matrices(pred: pd.DataFrame, unit_ids, source_ids, thresholds):
    """Build a_ir(lambda)=average across selected labels of equal-share supporting sources."""
    ui={str(g):i for i,g in enumerate(unit_ids)}; si={str(g):i for i,g in enumerate(source_ids)}
    rows=[]
    for r in pred[['gene','go_id','score','support_sources']].itertuples(index=False,name=None):
        g,go,score,supp=r
        i=ui.get(str(g))
        if i is None: continue
        ss=[si[s] for s in str(supp).split(';') if s in si]
        if ss: rows.append((i,str(go),float(score),ss))
    out=[]
    n=len(unit_ids); M=len(source_ids)
    for t in thresholds:
        sel=[x for x in rows if x[2]>=t]
        nt=np.zeros(n,np.int32)
        for i,go,s,ss in sel: nt[i]+=1
        rr=[];cc=[];dd=[]
        for i,go,s,ss in sel:
            w=1.0/(nt[i]*len(ss)); rr.extend([i]*len(ss));cc.extend(ss);dd.extend([w]*len(ss))
        A=sp.coo_matrix((dd,(rr,cc)),shape=(n,M),dtype=np.float64).tocsr();A.sum_duplicates();out.append(A)
    return out


def assign_folds(unit_ids_subset, K=5, rule='order', seed=0):
    m=len(unit_ids_subset)
    if rule=='order': return np.arange(m)%K
    if rule=='seedperm':
        p=np.random.default_rng(20260905+seed+99173).permutation(m);f=np.empty(m,int);f[p]=np.arange(m)%K;return f
    if rule=='genehash':
        return np.array([int(hashlib.sha256(str(g).encode()).hexdigest()[:8],16)%K for g in unit_ids_subset],int)
    raise ValueError(rule)


def goco_a_variance(A_full, losses_full, cal_idx, unit_ids, K=5,kappa=8.0,fold_rule='order',seed=0):
    A=A_full[cal_idx,:].tocsr(); y=np.asarray(losses_full)[cal_idx].astype(float); m=len(y);M=A.shape[1]
    folds=assign_folds(np.asarray(unit_ids)[cal_idx],K,fold_rule,seed)
    totalD=np.asarray(A.sum(axis=0)).ravel(); totalN=np.asarray(A.T.dot(y)).ravel(); totalY=float(y.sum())
    D=[];N=[];SQ=[]; idxs=[]
    for f in range(K):
        ix=np.flatnonzero(folds==f);idxs.append(ix);Af=A[ix,:]
        D.append(np.asarray(Af.sum(axis=0)).ravel());N.append(np.asarray(Af.T.dot(y[ix])).ravel());SQ.append(np.asarray(Af.multiply(Af).sum(axis=0)).ravel())
    # fold-excluded fitted means
    mu=np.zeros(m,float)
    for f in range(K):
        ix=idxs[f]
        if len(ix)==0:continue
        ntr=m-len(ix); glob=(totalY-y[ix].sum())/ntr if ntr else y.mean()
        den=totalD-D[f]; num=totalN-N[f]
        pr=(num+kappa*glob)/(den+kappa)
        mu[ix]=np.asarray(A[ix,:].dot(pr)).ravel()
        rowmass=np.asarray(A[ix,:].sum(axis=1)).ravel();mu[ix[rowmass<=1e-15]]=0.0
    diag=float(np.sum((y-mu)**2))
    off=0.0
    for f in range(K):
        for g in range(K):
            if f==g:
                nex=m-len(idxs[f]); sy=totalY-y[idxs[f]].sum(); den=totalD-D[f]; num=totalN-N[f]
            else:
                nex=m-len(idxs[f])-len(idxs[g]); sy=totalY-y[idxs[f]].sum()-y[idxs[g]].sum(); den=totalD-D[f]-D[g]; num=totalN-N[f]-N[g]
            if nex<=0:continue
            glob=sy/nex
            pr=(num+kappa*glob)/(den+kappa);nu=pr*(1-pr)
            cross=D[f]*D[g]
            if f==g:cross=cross-SQ[f]
            off += float(np.dot(cross,nu))
    off=max(off,0.0)
    return (diag+off)/(m*m), diag/(m*m), off/(m*m), mu


def pvalue_from_var(losses_cal, alpha, var):
    r=float(np.mean(losses_cal))
    if var<=0:return 0.0 if r<alpha else 1.0
    return float(norm.cdf((r-alpha)/math.sqrt(var)))


def fixed_select(pvals,delta=.10):
    last=None
    for j,p in enumerate(pvals):
        if p<=delta:last=j
        else:break
    return last
