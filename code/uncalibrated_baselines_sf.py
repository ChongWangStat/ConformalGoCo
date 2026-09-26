"""Release rules in current use, second predictor family (FunMap / STRING), same 100 splits.
Rules and conventions identical to wf/uncalibrated_baselines.py."""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from goco_second_family import FunMap, SEED0, hashu, P as DATA_DIR

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=100)
ap.add_argument("--out", required=True)
args = ap.parse_args()

t0 = time.time()
ds = FunMap()
s = ds.rows.score.to_numpy(float); gi = ds.rows.i.to_numpy(); T = ds.rows["T"].to_numpy(np.int32)
# refined global grid of run_second_family_v4.py (a function of frozen scores alone), ascending
COARSE = np.unique(np.round(np.concatenate([np.geomspace(1.0, 120.0, 30), [1.30103, 2.0]]), 6))
levels = np.unique(np.round(np.concatenate([COARSE, np.geomspace(1.0, 120.0, 300)]), 6))
K = len(levels)
# cumulative counts at each level, most liberal last: a call is released at level j iff score >= levels[j]
b = np.searchsorted(levels, np.round(s, 6), side="right") - 1
keep = b >= 0; b = b[keep]; gi = gi[keep]; T = T[keep]; s_keep = s[keep]
from scipy import sparse
def cum(v):
    m = sparse.coo_matrix((v, (gi, b)), shape=(ds.n, K)).tocsr().toarray()
    return np.cumsum(m[:, ::-1], axis=1)[:, ::-1]        # column j = calls with score >= levels[j]
C = cum(np.ones(len(b), np.int32)).astype(np.int32); Y = cum(T).astype(np.int32)
L = np.divide(C - Y, C, out=np.zeros(C.shape), where=C > 0)
print("units %d calls %d levels %d  %.0fs" % (ds.n, ds.ncall, K, time.time() - t0), flush=True)

# complete TruePath annotation set size per gene (CAFA recall denominator)
tp = pd.read_csv(os.path.join(str(DATA_DIR), "go_truth_true_path.csv.gz"), usecols=["gene", "go_id"]).drop_duplicates()
ntrue_full = tp.groupby("gene").size().reindex(pd.Index(ds.units)).fillna(0).to_numpy(float)
ntrue = Y[:, 0].astype(float)                            # supported calls inside the candidate space

# per-gene top-k
gi_all = ds.rows.i.to_numpy(); T_all = ds.rows["T"].to_numpy(np.int32)
h = hashu(ds.rows.gene.to_numpy() + "|" + ds.rows.go_id.to_numpy(), "topk")
order = np.lexsort((h, -s, gi_all)); rank = np.empty(len(s), np.int64)
starts = np.searchsorted(gi_all[order], np.arange(ds.n), side="left"); rank[order] = np.arange(len(s)) - starts[gi_all[order]]
topk = {}
for k in (1, 3):
    m = rank < k
    c = np.bincount(gi_all[m], minlength=ds.n).astype(float); y = np.bincount(gi_all[m], weights=T_all[m], minlength=ds.n)
    topk[k] = (c, y, np.where(c > 0, 1 - y / np.maximum(c, 1), 0.0))

rows = []
for rep in range(args.repeats):
    perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
    nc = round(0.70 * ds.n); cal, ev = perm[:nc], perm[nc:]
    Ccal, Ycal = C[cal], Y[cal]
    npred = (Ccal > 0).sum(axis=0)
    pr = np.where(npred > 0, np.where(Ccal > 0, Ycal / np.maximum(Ccal, 1), 0.0).sum(axis=0) / np.maximum(npred, 1), 0.0)
    for tag, nt in (("Fmax", ntrue), ("Fmax-full", ntrue_full)):
        ct = cal[nt[cal] > 0]
        rc = (Y[ct] / nt[ct][:, None]).mean(axis=0)
        F = np.where(pr + rc > 0, 2 * pr * rc / np.maximum(pr + rc, 1e-12), 0.0)
        j = int(np.argmax(F))
        rows.append(dict(repeat=rep, rule=tag, label=f"score>={levels[j]:g}", fdp=float(L[ev, j].mean()),
                         correct=float(Y[ev, j].sum()), total=float(C[ev, j].sum()), units=int((C[ev, j] > 0).sum()),
                         supp_units=int((Y[ev, j] > 0).sum()), fmax=float(F[j]), prec=float(pr[j]), rec=float(rc[j])))
    for k in (1, 3):
        c, y, l = topk[k]
        rows.append(dict(repeat=rep, rule=f"Top-{k}", label=f"top{k}", fdp=float(l[ev].mean()), correct=float(y[ev].sum()),
                         total=float(c[ev].sum()), units=int((c[ev] > 0).sum()), supp_units=int((y[ev] > 0).sum()),
                         fmax=np.nan, prec=np.nan, rec=np.nan))
    if rep % 25 == 0: print("  repeat %3d  %.0fs" % (rep, time.time() - t0), flush=True)
pd.DataFrame(rows).to_csv(args.out, index=False)
print("wrote", args.out, "%.0fs" % (time.time() - t0), flush=True)
