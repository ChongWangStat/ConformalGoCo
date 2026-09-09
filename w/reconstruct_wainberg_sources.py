import pandas as pd, numpy as np, json, math, gzip
from collections import defaultdict
from pathlib import Path
W=Path('/mnt/data/conformalgoco_analysis_2026-09-06/wainberg')
O=Path('/mnt/data/conformalgoco_analysis_2026-09-06/results')
O.mkdir(exist_ok=True)
# eligible units/background
units=pd.read_csv(W/'network_communities.csv',usecols=['gene'])['gene'].astype(str).tolist()
bg=set(units); N0=len(bg)
# non-syntenic modules and memberships
mods=pd.read_csv(W/'modules.csv')
non=set(mods.loc[mods['synteny'].eq('Non-syntenic'),'module_id'].dropna().astype(int).astype(str))
mm=pd.read_csv(W/'module_memberships.csv.gz',usecols=['module_id','gene'])
mm['module_id']=mm['module_id'].astype(str)
mm['gene']=mm['gene'].astype(str)
mm=mm[mm.module_id.isin(non) & mm.gene.isin(bg)].drop_duplicates()
mod_genes={mid:set(g.gene) for mid,g in mm.groupby('module_id')}
gene_mods=defaultdict(list)
for mid,gs in mod_genes.items():
    for g in gs: gene_mods[g].append(mid)
# direct truth restricted to background
truth=pd.read_csv(W/'go_truth_direct.csv.gz',usecols=['gene','go_id'])
truth['gene']=truth.gene.astype(str); truth['go_id']=truth.go_id.astype(str)
truth=truth[truth.gene.isin(bg)].drop_duplicates()
term_genes={go:set(g.gene) for go,g in truth.groupby('go_id')}
# only high score rows relevant to dense grid
scores=pd.read_csv(W/'wainberg_gene_go_scores.csv.gz',usecols=['gene','go_id','score'])
scores=scores[scores.gene.astype(str).isin(bg) & (scores.score>=600)].copy()
scores['gene']=scores.gene.astype(str); scores['go_id']=scores.go_id.astype(str)
print('high rows',len(scores),'genes',scores.gene.nunique(),'terms',scores.go_id.nunique())
rows=[]; no_match=0; no_mod=0
for z,(target,term,score) in enumerate(scores[['gene','go_id','score']].itertuples(index=False,name=None),1):
    tg=term_genes.get(term,set())
    # historical term eligibility count >=5 in background before target removal
    if len(tg) < 5: continue
    best=None
    N=N0-1
    K=len(tg - {target})
    if K<=0: continue
    for mid in gene_mods.get(target,[]):
        gs=mod_genes[mid]-{target}
        n=len(gs)
        if n<=0: continue
        supp=gs & tg
        k=len(supp)
        if k<=0: continue
        s=(k/n)/(K/N)
        # exact winner by score, then deterministic module id tie-break lower numeric
        key=(s, -int(mid) if mid.isdigit() else 0)
        if best is None or key>best[0]:
            best=(key,mid,supp,s,n,k,K,N)
    if best is None:
        no_mod+=1; continue
    _,mid,supp,s,n,k,K,N=best
    rel=abs(s-score)/max(abs(score),1e-300)
    if rel>1e-8:
        # find closest matching candidate in case score construction has special eligibility
        cand=[]
        for mid2 in gene_mods.get(target,[]):
            gs=mod_genes[mid2]-{target}; n2=len(gs)
            if not n2: continue
            supp2=gs&tg; k2=len(supp2)
            if not k2: continue
            s2=(k2/n2)/(K/N)
            cand.append((abs(s2-score),mid2,supp2,s2,n2,k2))
        if cand:
            cand.sort(key=lambda x:(x[0],int(x[1]) if x[1].isdigit() else 10**9))
            diff,mid2,supp2,s2,n2,k2=cand[0]
            rel2=abs(s2-score)/max(abs(score),1e-300)
            # If score is exactly matched by non-max module, historical code might impose module filters beyond synteny.
            if rel2<rel:
                mid,supp,s,n,k=mid2,supp2,s2,n2,k2; rel=rel2
    if rel>1e-6: no_match+=1
    rows.append((target,term,float(score),mid,float(s),rel,n,k,K,';'.join(sorted(supp))))
    if z%5000==0: print('processed',z,'no_match',no_match)
out=pd.DataFrame(rows,columns=['gene','go_id','score','winning_module','reconstructed_score','relative_error','module_n_ex_target','support_k','background_K_ex_target','support_sources'])
out.to_csv(O/'wainberg_highscore_source_attribution.csv.gz',index=False,compression='gzip')
print('wrote',len(out),'no_match >1e-6',no_match,'no_mod',no_mod)
print(out.relative_error.describe(percentiles=[.5,.9,.99,.999]).to_string())
print('exact <1e-10', (out.relative_error<1e-10).mean())
