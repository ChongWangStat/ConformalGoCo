"""Does the call-resolution path (goco-rank) beat GoCo-L on Wainberg?

PREREGISTERED PREDICTION (written before running): it should LOSE.

The scope condition derived in DERIVATION_parallel_to_wainberg.md says the call-resolution path
dominates when (a) the granularity term is already exhausted by a dense global grid, and (b) the
base unsupported rate is low enough that spreading releases across genes is affordable.  Wainberg
fails both.  Its ties are exact, so no global grid reaches part (a) -- that is the entire reason
GoCo exists.  And its base unsupported rate is ~0.79, close to ProteomeHD truepath (0.71), where
goco-rank lost 1,174 supported calls and won 0/100 splits.  Per-gene FDP with abstention costing
zero rewards concentration, and at a high base error rate spreading is unaffordable.

If goco-rank WINS here, the derivation is wrong and the main paper's method choice needs revisiting.

Nothing frozen is modified.  goco_learned._fit_model supplies p_hat per call using the published
f1-f7 features and cross-fitting; only the PATH is new.  Splits are byte-identical to the frozen
runs (gr.SEED0), so GoCo-L is read from the frozen results rather than recomputed.
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
from scipy import sparse

R = os.environ.get("GOCO_CODE", r"C:/Users/chwang/Box/Research/26 09 ConformalGOCo/GoCo_BIB_revision_v2/analysis_rerun")
sys.path.insert(0, R)
import goco_rerun as gr          # frozen
import goco_learned as gl        # frozen

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=50)
ap.add_argument("--out", default=r"C:/g1/wainberg_rank.csv")
args = ap.parse_args()
ALPHAS = [0.05, 0.10, 0.20]
DELTAS = [0.50, 0.10]
MFRAC = np.unique(np.round(np.geomspace(0.0005, 1.0, 60), 6))   # identical to the FunMap/STRING runs

dset, tset = gr.load_truth()
ds = gr.Dataset("Wainberg", dset, tset)
rows = ds.rows
s_all = np.round(rows.score.to_numpy(float), 6)
keep = s_all >= ds.dense_floor
call_i = rows.i.to_numpy()[keep]
call_T = rows["T"].to_numpy(np.int32)[keep]
call_rows = np.where(keep)[0]
ncall = len(call_rows)
base_unsup = 1.0 - call_T.mean()
print("units %d | candidate calls %d | base unsupported rate %.3f"
      % (ds.n, ncall, base_unsup), flush=True)

hcall = gr.hashu(rows.gene.to_numpy()[keep], salt="half")


def rank_path(phat):
    """Per-unit cumulative counts along the top-m path ordered by predicted support."""
    key = np.lexsort((hcall, phat))
    rank = np.empty(ncall, np.int64); rank[key] = np.arange(ncall)
    edges = np.unique(np.maximum(1, (MFRAC * ncall).astype(np.int64)))
    b = np.searchsorted(edges, rank, side="left")
    ok = b < len(edges)
    shape = (ds.n, len(edges))
    cum = lambda v: np.cumsum(sparse.coo_matrix(
        (v, (call_i[ok], b[ok])), shape=shape).tocsr().toarray(), axis=1)
    C = cum(np.ones(int(ok.sum()), np.int32)).astype(np.int32)
    Y = cum(call_T[ok]).astype(np.int32)
    L = np.divide(C - Y, C, out=np.zeros(shape), where=C > 0)
    return C, Y, L, edges


out, t0 = [], time.time()
for rep in range(args.repeats):
    perm = np.random.default_rng(gr.SEED0 + rep).permutation(ds.n)
    nt, nc = round(0.10 * ds.n), round(0.70 * ds.n)
    train, cert, ev = perm[:nt], perm[nt:nc], perm[nc:]

    predict, s_m, i_m = gl._fit_model(ds, train, "learned")     # frozen model, published features
    phat = predict(call_rows)
    C, Y, L, edges = rank_path(phat)

    for alpha in ALPHAS:
        for delta in DELTAS:
            sel = None
            for j in range(len(edges)):
                p, _, _ = gr.clt_p(L[cert, j], alpha)
                if p > delta:
                    break
                sel = j
            if sel is None:
                out.append(dict(repeat=rep, alpha=alpha, delta=delta, arm="goco-rank",
                                label="abstain", fdp=0.0, correct=0, total=0, units=0))
                continue
            out.append(dict(repeat=rep, alpha=alpha, delta=delta, arm="goco-rank",
                            label="top%d" % edges[sel], fdp=float(L[ev, sel].mean()),
                            correct=int(Y[ev, sel].sum()), total=int(C[ev, sel].sum()),
                            units=int((C[ev, sel] > 0).sum())))
    if rep % 5 == 0:
        print("  repeat %2d  %.0fs" % (rep, time.time() - t0), flush=True)

pd.DataFrame(out).to_csv(args.out, index=False)
print("wrote %s in %.0fs" % (args.out, time.time() - t0), flush=True)
