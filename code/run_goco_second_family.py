"""GoCo on FunMap against global-threshold baselines, identical splits and certificate.

Arms (all certified by fixed-sequence LTT with the finite-population normal p-value, so any
difference in yield is attributable to the candidate path alone):
  global-coarse : geometric grid, global thresholds only        -- the GoDag analogue
  global-dense  : 10x denser grid, global thresholds only       -- free densification
  goco-score    : coarse grid + partial admission ordered by the added-call score (control)
  goco          : coarse grid + partial admission ordered by the learned per-call model
  goco-rank     : path defined directly by the learned per-call ordering (top-m calls)

goco-rank is the structure-specific step beyond GoCo as published.  GoCo uses the model only to
break ties inside a step, because Wainberg's enrichment scores tie in blocks of thousands and the
block is the natural unit.  FunMap's BH-adjusted scores carry no material ties, so the model can
order every candidate call and define the path itself: policies are the top-m calls by predicted
support, nested in m, and measurable with respect to the ranking fold, so the same fixed-sequence
certificate applies without modification.
"""
from __future__ import annotations
import sys, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, r"C:/g1")
from goco_second_family import (FunMap, clt_p, hashu, QGRID, TRIGGER_FRAC, MIN_RANK_UNITS, SEED0)
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=20)
ap.add_argument("--alphas", default="0.05,0.10,0.20")
ap.add_argument("--delta", type=float, default=0.10)
ap.add_argument("--out", default=r"C:/g1/goco_funmap_results.csv")
args = ap.parse_args()
ALPHAS = [float(x) for x in args.alphas.split(",")]

COARSE = np.unique(np.round(np.concatenate([np.geomspace(1.0, 120.0, 30), [1.30103, 2.0]]), 6))
DENSE = np.unique(np.round(np.concatenate([COARSE, np.geomspace(1.0, 120.0, 300)]), 6))

print("loading frozen predictor objects ...", flush=True)
ds = FunMap()
Mc, Md = ds.mats(COARSE), ds.mats(DENSE)
print("units %d | candidate calls %d | coarse %d | dense %d"
      % (ds.n, ds.ncall, len(COARSE), len(DENSE)), flush=True)

CALL_I = ds.rows.i.to_numpy()
CALL_S = ds.rows.score.to_numpy(float)
CALL_T = ds.rows["T"].to_numpy(np.int32)
MFRAC = np.unique(np.round(np.geomspace(0.0005, 1.0, 60), 6))   # prespecified admitted-call fractions


def rank_path(ds, phat, ev_all):
    """Per-unit counts along the path that admits calls in increasing predicted-unsupported order."""
    from scipy import sparse
    key = np.lexsort((ds.hcall, phat))          # most-likely-supported first, hash breaks ties
    rank = np.empty(ds.ncall, np.int64); rank[key] = np.arange(ds.ncall)
    edges = np.unique(np.maximum(1, (MFRAC * ds.ncall).astype(np.int64)))
    b = np.searchsorted(edges, rank, side="left")       # bin index of each call
    ok = b < len(edges)
    shape = (ds.n, len(edges))
    cum = lambda v: np.cumsum(sparse.coo_matrix(
        (v, (CALL_I[ok], b[ok])), shape=shape).tocsr().toarray(), axis=1)
    C = cum(np.ones(int(ok.sum()), np.int32)).astype(np.int32)
    Y = cum(CALL_T[ok]).astype(np.int32)
    L = np.divide(C - Y, C, out=np.zeros(shape), where=C > 0)
    return dict(C=C, Y=Y, L=L, grid=edges.astype(float))


def rank_arm(M, cert, alpha, delta):
    L, grid = M["L"], M["grid"]
    return fixed_sequence([{"loss": L[:, j], "j": j, "label": f"top{int(grid[j])}"}
                           for j in range(len(grid))], cert, alpha, delta)


def fixed_sequence(cands, cert, alpha, delta):
    last = None
    for c in cands:
        if clt_p(c["loss"][cert], alpha) > delta:
            break
        last = c
    return last


def global_arm(M, cert, alpha, delta):
    L, grid = M["L"], M["grid"]
    return fixed_sequence([{"loss": L[:, j], "j": j, "label": f"global_{grid[j]:g}"}
                           for j in range(len(grid) - 1, -1, -1)], cert, alpha, delta)


def fit_phat(ds, train):
    """Cross-fitted per-call probability of being unsupported, fitted on ranking-fold calls only."""
    X, intr, unsup = ds.call_features(train)
    y = unsup.astype(int)
    p = np.full(ds.ncall, np.nan)
    for half in (0, 1):
        m = intr & ((ds.hcall < 0.5) if half == 0 else (ds.hcall >= 0.5))
        if m.sum() < 50 or len(np.unique(y[m])) < 2:
            continue
        s = StandardScaler().fit(X[m])
        clf = LogisticRegression(C=0.5, max_iter=1000).fit(s.transform(X[m]), y[m])
        tgt = (ds.hcall >= 0.5) if half == 0 else (ds.hcall < 0.5)   # score the opposite half
        p[tgt] = clf.predict_proba(s.transform(X[tgt]))[:, 1]
    if np.isnan(p).all():
        return None
    return np.where(np.isnan(p), np.nanmean(p), p)


def goco_arm(ds, M, cert, train, alpha, delta, phat, mode="learned"):
    L, grid = M["L"], M["grid"]
    cands = []
    for j in range(len(grid) - 1, 0, -1):
        cands.append({"loss": L[:, j], "j": j, "label": f"global_{grid[j]:g}"})
        if float(L[train, j - 1].mean()) <= TRIGGER_FRAC * alpha:
            continue
        lo, hi = float(grid[j - 1]), float(grid[j])
        sel = (CALL_S >= lo) & (CALL_S < hi)
        if not sel.any():
            continue
        bi = CALL_I[sel]
        aff = np.unique(bi)
        if len(np.intersect1d(aff, train)) < MIN_RANK_UNITS:
            continue
        pv = phat[sel] if mode == "learned" else 1.0 - CALL_S[sel] / max(CALL_S[sel].max(), 1e-9)
        df = pd.DataFrame({"i": bi, "p": pv}).groupby("i").agg(
            dmean=("p", "mean"), ssum=("p", lambda v: float((1 - v).sum())))
        u = np.full(ds.n, np.inf)
        u[df.index.to_numpy()] = df.dmean.to_numpy() / np.maximum(df.ssum.to_numpy(), 1e-3)
        ranked = aff[np.lexsort((ds.hh[aff], u[aff]))]
        tset = set(train.tolist())
        tpos = [k for k, i in enumerate(ranked) if i in tset]
        for q in QGRID:
            kq = max(1, int(np.floor(q * len(tpos))))
            if kq > len(tpos):
                continue
            adm = ranked[: tpos[kq - 1] + 1]
            loss = L[:, j].copy()
            loss[adm] = L[adm, j - 1]
            cands.append({"loss": loss, "j": j, "jlo": j - 1, "adm": adm,
                          "label": f"adapt_{hi:g}_to_{lo:g}_q{q:g}"})
    return fixed_sequence(cands, cert, alpha, delta)


def report(res, M, ev):
    if res is None:
        return dict(label="abstain", fdp=0.0, correct=0, total=0, units=0)
    C, Y = M["C"][:, res["j"]].copy(), M["Y"][:, res["j"]].copy()
    if "adm" in res:
        C[res["adm"]] = M["C"][res["adm"], res["jlo"]]
        Y[res["adm"]] = M["Y"][res["adm"], res["jlo"]]
    return dict(label=res["label"], fdp=float(res["loss"][ev].mean()),
                correct=int(Y[ev].sum()), total=int(C[ev].sum()), units=int((C[ev] > 0).sum()))


rows, t0 = [], time.time()
for rep in range(args.repeats):
    perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
    nt, nc = round(0.10 * ds.n), round(0.70 * ds.n)
    train, cert, ev = perm[:nt], perm[nt:nc], perm[nc:]
    phat = fit_phat(ds, train)
    Mr = rank_path(ds, phat, ev)
    for alpha in ALPHAS:
        arms = {
            "global-coarse": (global_arm(Mc, cert, alpha, args.delta), Mc),
            "global-dense": (global_arm(Md, cert, alpha, args.delta), Md),
            "goco-score": (goco_arm(ds, Mc, cert, train, alpha, args.delta, phat, "score"), Mc),
            "goco": (goco_arm(ds, Mc, cert, train, alpha, args.delta, phat, "learned"), Mc),
            "goco-rank": (rank_arm(Mr, cert, alpha, args.delta), Mr),
        }
        for name, (res, M) in arms.items():
            rows.append(dict(repeat=rep, alpha=alpha, arm=name, **report(res, M, ev)))
    print("  repeat %2d  %.0fs" % (rep, time.time() - t0), flush=True)

pd.DataFrame(rows).to_csv(args.out, index=False)
print("wrote", args.out, "in %.0fs" % (time.time() - t0))
