"""Ranking-fold selection rule with K-fold cross-fitting inside T (rehearsal only).
Certified results are untouched; this only decides M vs N from T. Splits byte-identical to the main analysis."""
from __future__ import annotations
import os, sys, time, argparse, types
import numpy as np, pandas as pd
from scipy import sparse
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goco_rerun as gr, goco_learned as gl

ap = argparse.ArgumentParser()
ap.add_argument("--datasets", default="Wainberg,Sanger,DRIVE,HAP1"); ap.add_argument("--repeats", type=int, default=100)
ap.add_argument("--K", type=int, default=10); ap.add_argument("--out", default="wf_selection_kfold.csv")
args = ap.parse_args(); K = args.K
ALPHAS = [0.05, 0.10, 0.20]; MFRAC = gr.MFRAC_N

def fit_model_kfold(ds, train, kind):
    """Same as gl._fit_model but T calls are scored by the model fitted on the other K-1 hash folds of T."""
    rows = ds.rows; s = np.round(rows.score.to_numpy(float), 6); T = rows['T'].to_numpy(float); i_all = rows.i.to_numpy()
    srcs = rows.srcs.tolist(); terms = rows.go_id.to_numpy(); genes = rows.gene.to_numpy()
    train_set = np.zeros(ds.n, bool); train_set[train] = True
    fold = np.zeros(ds.n, np.int16) - 1
    fold[train] = np.minimum((gr.hashu(ds.units[train], salt='xfit') * K).astype(int), K - 1)
    keep = s >= ds.dense_floor; tr_rows = np.where(keep & train_set[i_all])[0]
    y = 1.0 - T; glob = float(y[tr_rows].mean()); calls_total = np.bincount(i_all[keep], minlength=ds.n).astype(float)
    g2t = gl._truth_by_gene(); cnt = Counter()
    for g in ds.units[train]:
        for t in g2t.get(g, ()): cnt[t] += 1
    ntr = float(len(train))
    def f6(idx):
        out = np.empty(len(idx))
        for q, r in enumerate(idx):
            c = cnt.get(terms[r], 0)
            if train_set[i_all[r]] and terms[r] in g2t.get(genes[r], ()): c -= 1
            out[q] = np.log((c + 0.5) / (ntr + 1.0))
        return out
    base = dict(f1=np.log(np.maximum(s, 1e-9)), f2=np.log1p(np.array([len(x) for x in srcs], float)), f3=np.log1p(calls_total[i_all]), f6=f6, x1=gl._extra_col(ds), terms=terms, srcs=srcs)
    folds = {None: gl._Fold(tr_rows, y, terms, srcs, glob, kind, None).fit(base)}
    for h in range(K):
        sub = tr_rows[fold[i_all[tr_rows]] != h]
        folds[h] = gl._Fold(sub, y, terms, srcs, glob, kind, None).fit(base)
    def predict(idx):
        out = np.empty(len(idx)); hidx = np.where(train_set[i_all[idx]], fold[i_all[idx]], -1)
        for h in [-1] + list(range(K)):
            sel = np.where(hidx == h)[0]
            if len(sel): out[sel] = folds[None if h == -1 else h].predict(idx[sel], base)
        return out
    return predict, s, i_all

dset, tset = gr.load_truth(); rows_out = []; t0 = time.time()
for name in args.datasets.split(","):
    ds = gr.Dataset(name, dset, tset); mg = ds.mats(ds.grid)
    rows = ds.rows; s_all = np.round(rows.score.to_numpy(float), 6); keep = s_all >= ds.dense_floor
    call_i = rows.i.to_numpy()[keep]; call_T = rows["T"].to_numpy(np.int32)[keep]; call_rows = np.where(keep)[0]; ncall = len(call_rows)
    hcall = gr.hashu(rows.gene.to_numpy()[keep], salt="half")
    for rep in range(args.repeats):
        perm = np.random.default_rng(gr.SEED0 + rep).permutation(ds.n); nt = round(0.10 * ds.n); train = perm[:nt]
        predict, _, _ = fit_model_kfold(ds, train, "learned")
        # monkeypatch so the M path's utility uses the K-fold model for this rehearsal
        gl._MODEL_CACHE.clear(); gl._fit_model = lambda ds_, tr_, kind_, _p=predict, _s=s_all, _i=rows.i.to_numpy(): (_p, _s, _i)
        phat = predict(call_rows)
        C, Y, L, edges = gr.topm_path(phat, hcall, call_i, call_T, ds.n, MFRAC)   # exact top-m sets (one-based ranks)
        for alpha in ALPHAS:
            selN = None
            for j in range(len(edges)):
                if float(L[train, j].mean()) > alpha: break
                selN = j
            yN = int(Y[train, selN].sum()) if selN is not None else 0
            bf = lambda j: ds.Aadd(mg["thr"][j + 1], mg["thr"][j])
            lastM = gr.run_path(ds, mg, train, train, train, alpha, 0.5, "TruePath", "goco", gr.QGRID_E, True, bf, {})
            yM = int(gr.metrics(ds, mg, lastM, train, train)["go_yield"]) if lastM is not None else 0
            rows_out.append(dict(dataset=name, repeat=rep, alpha=alpha, K=K, yield_M_T=yM, yield_N_T=yN, choose_N=int(yN > yM)))
        if rep % 20 == 0: print("  [%s] rep %d %.0fs" % (name, rep, time.time() - t0), flush=True)
    pd.DataFrame(rows_out).to_csv(args.out, index=False)
print("done %.0fs" % (time.time() - t0))
