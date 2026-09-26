"""Restore any missing explicit split tables and verify all six against the v4 bundle.

Checksums cover decompressed CSV bytes, so gzip timestamps do not matter.
No analysis is fitted and no reported outcome is recomputed here.
"""
from pathlib import Path
import gzip, hashlib, json, sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; DEST=ROOT/'splits'
manifest=json.loads((DEST/'manifest.json').read_text())
truth=None
for name,item in manifest.items():
    target=DEST/(name+'_splits.csv.gz')
    if not target.exists():
        if name=='STRING':
            nb=pd.read_csv(ROOT/'string'/'string_top50_neighborhoods.csv.gz')
            truth_genes=set(pd.read_csv(ROOT/'funmap'/'go_truth_direct.csv.gz').gene)
            units=sorted(set(nb.iloc[:,0]) & truth_genes)
        elif name=='FunMap':
            raise FileNotFoundError('The already uploaded FunMap split table is required.')
        else:
            import goco_rerun as gr
            if truth is None: truth=gr.load_truth()
            units=list(gr.Dataset(name,*truth).units)
        n=len(units);assert n==item['n'],(name,n,item['n'])
        nt,nc=round(.1*n),round(.7*n)
        folds=np.empty((100,n),dtype='U1')
        for s in range(100):
            perm=np.random.default_rng(item['seed0']+s).permutation(n)
            folds[s,perm[:nt]]='T';folds[s,perm[nt:nc]]='C';folds[s,perm[nc:]]='E'
        df=pd.DataFrame({'gene':np.tile(units,100),'split':np.repeat(np.arange(100),n),'fold':folds.ravel()})
        df.to_csv(target,index=False,compression={'method':'gzip','mtime':0})
    raw=gzip.decompress(target.read_bytes())
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=item['csv_sha256']: raise RuntimeError('Split table mismatch: '+name+' '+digest)
    print('PASS split assignments:',name,item['n'],'genes x 100 splits')
