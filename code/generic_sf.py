"""Generic Direct-loss calibrators for the neighbourhood-structure family (FunMap / STRING).

Same panel logic as the module-structure family: every rule selects ONE global threshold on the
coarse grid, walking from conservative to liberal and stopping at the first grid point whose
p-value (or plug-in predicted risk) exceeds the budget.  All rules see Direct labels only -- no GO
DAG, no neighbourhood evidence -- and every rule uses the same 60% certification fold and is
evaluated on the same 30% held-out fold as GoDag and GoCo, with the same split streams.

Rules
  Multilabel / IID-normal   finite-population normal p-value (the released multilabel rule)
  Hoeffding/LTT             Hoeffding p-value
  Hoeffding-Bentkus         min(KL, e*Binom) Bentkus bound
  Empirical Bernstein       Maurer-Pontil empirical-Bernstein UCB inverted in delta
  Platt / isotonic          delta-free plug-in: fit score -> P(supported) on the certification fold
                            and take the most liberal grid point whose predicted per-gene FDP <= alpha
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
from scipy.stats import norm, binom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from goco_second_family import FunMap, SEED0

from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=100)
ap.add_argument("--start", type=int, default=0)
ap.add_argument("--alphas", default="0.05,0.10,0.20")
ap.add_argument("--deltas", default="0.50,0.10")
ap.add_argument("--grid", default="coarse", choices=["coarse","dense"])
ap.add_argument("--out", required=True)
args = ap.parse_args()
ALPHAS = [float(x) for x in args.alphas.split(",")]
DELTAS = [float(x) for x in args.deltas.split(",")]
_COARSE = np.unique(np.round(np.concatenate([np.geomspace(1.0, 120.0, 30), [1.30103, 2.0]]), 6))
_DENSE = np.unique(np.round(np.concatenate([_COARSE, np.geomspace(1.0, 120.0, 300)]), 6))
COARSE = _COARSE if args.grid == "coarse" else _DENSE

t0 = time.time()
print("loading ...", flush=True)
ds = FunMap()
P = Path(os.environ.get("GOCO_DATA", Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) / "funmap"))
direct = pd.read_csv(P / "go_truth_direct.csv.gz", usecols=["gene", "go_id"]).drop_duplicates()
ds.rows["D"] = pd.MultiIndex.from_frame(ds.rows[["gene", "go_id"]]).isin(pd.MultiIndex.from_frame(direct))
print("units %d | calls %d | direct-supported %.3f | %.0fs"
      % (ds.n, ds.ncall, ds.rows["D"].mean(), time.time() - t0), flush=True)


def mats(grid, flag):
    """Per-gene x grid cumulative calls C, supported Y and loss L for the given label vector."""
    grid = np.asarray(grid, float)
    mb = np.searchsorted(grid, ds.rows.score.to_numpy(), side="right") - 1
    keep = mb >= 0
    r, m = ds.rows.i.to_numpy()[keep], mb[keep]
    t = ds.rows[flag].to_numpy(np.int32)[keep]
    shape = (ds.n, len(grid))
    cum = lambda v: np.cumsum(sparse.coo_matrix((v, (r, m)), shape=shape).tocsr().toarray()[:, ::-1],
                              axis=1)[:, ::-1]
    C = cum(np.ones(int(keep.sum()), np.int32)).astype(np.int32)
    Y = cum(t).astype(np.int32)
    L = np.divide(C - Y, C, out=np.zeros(shape), where=C > 0)
    return dict(C=C, Y=Y, L=L, grid=grid)


MD = mats(COARSE, "D")   # Direct loss: what the generic calibrators see
MT = mats(COARSE, "T")   # TruePath: what everything is evaluated against
NG = len(COARSE)


def kl(a, b):
    eps = np.finfo(float).eps
    a = np.clip(a, eps, 1 - eps); b = np.clip(b, eps, 1 - eps)
    return a * np.log(a / b) + (1 - a) * np.log((1 - a) / (1 - b))


def p_hoeff(x, a):
    r = float(x.mean()); d = max(a - r, 0.0)
    return float(min(1.0, np.exp(-2 * len(x) * d * d)))


def p_hb(x, a):
    r = float(x.mean()); n = len(x)
    if r >= a:
        return 1.0
    return float(min(1.0, np.exp(-n * kl(r, a)), np.e * binom.cdf(int(np.ceil(n * r)), n, a)))


def p_norm(x, a):
    """Finite-population one-sided normal p-value (the released multilabel rule)."""
    r = float(x.mean()); m = len(x)
    s = float(np.std(x))
    if s <= 0:
        return 0.0 if r < a else 1.0
    return float(norm.cdf((r - a) / (s / np.sqrt(m))))


def eb_ucb(x, delta):
    n = len(x)
    if n < 2:
        return 1.0
    sd = float(x.std(ddof=1)); z = np.log(2 / delta)
    return float(min(1.0, x.mean() + sd * np.sqrt(2 * z / n) + 7 * z / (3 * (n - 1))))


def p_eb(x, a):
    if x.mean() >= a:
        return 1.0
    lo, hi = 1e-15, 1 - 1e-12
    if eb_ucb(x, hi) > a:
        return 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if eb_ucb(x, mid) <= a:
            hi = mid
        else:
            lo = mid
    return hi


RULES = [("Multilabel", p_norm), ("Hoeffding/LTT", p_hoeff), ("Hoeffding-Bentkus", p_hb),
         ("IID-normal", p_norm), ("Empirical Bernstein", p_eb)]


def walk(pvals, delta):
    """Fixed sequence, conservative (high threshold) to liberal; last index before the first failure."""
    sel = None
    for j in range(NG - 1, -1, -1):
        if pvals[j] > delta:
            break
        sel = j
    return sel


def plugin_threshold(cert, alpha, phat_call):
    """Most liberal grid point whose plug-in predicted per-gene FDP on the certification fold <= alpha."""
    mb = np.searchsorted(COARSE, ds.rows.score.to_numpy(), side="right") - 1
    keep = mb >= 0
    r, m = ds.rows.i.to_numpy()[keep], mb[keep]
    shape = (ds.n, NG)
    cum = lambda v: np.cumsum(sparse.coo_matrix((v, (r, m)), shape=shape).tocsr().toarray()[:, ::-1],
                              axis=1)[:, ::-1]
    C = cum(np.ones(int(keep.sum()), np.float64))
    S = cum(phat_call[keep].astype(np.float64))
    pred = np.divide(C - S, C, out=np.zeros(shape), where=C > 0)
    sel = None
    for j in range(NG - 1, -1, -1):
        if float(pred[cert, j].mean()) > alpha:
            break
        sel = j
    return sel


def report(j, ev, pool):
    if j is None:
        return dict(threshold=np.nan, fdp=0.0, pool_risk=0.0, supp_genes=0, supp_terms=0, calls=0, genes=0)
    Lt, Yt, Ct = MT["L"], MT["Y"], MT["C"]
    return dict(threshold=float(COARSE[j]), fdp=float(Lt[ev, j].mean()), pool_risk=float(Lt[pool, j].mean()),
                supp_genes=int((Yt[ev, j] > 0).sum()), supp_terms=int(Yt[ev, j].sum()),
                calls=int(Ct[ev, j].sum()), genes=int((Ct[ev, j] > 0).sum()))


score_call = ds.rows.score.to_numpy(float)
D_call = ds.rows["D"].to_numpy(int)
rows = []
for rep in range(args.start, args.start + args.repeats):
    perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
    nt, nc = round(0.10 * ds.n), round(0.70 * ds.n)
    cert, ev, pool = perm[nt:nc], perm[nc:], perm[nt:]
    in_cert = np.isin(ds.rows.i.to_numpy(), cert)

    # delta-free plug-in maps, fitted on certification-fold calls
    x = np.log1p(score_call[in_cert]).reshape(-1, 1); y = D_call[in_cert]
    try:
        platt = LogisticRegression(max_iter=1000).fit(x, y).predict_proba(
            np.log1p(score_call).reshape(-1, 1))[:, 1]
    except Exception:
        platt = np.full(ds.ncall, float(y.mean()))
    try:
        iso = IsotonicRegression(out_of_bounds="clip").fit(np.log1p(score_call[in_cert]), y).predict(
            np.log1p(score_call))
    except Exception:
        iso = np.full(ds.ncall, float(y.mean()))

    for alpha in ALPHAS:
        pv = {name: np.array([fun(MD["L"][cert, j], alpha) for j in range(NG)]) for name, fun in RULES}
        jp = plugin_threshold(cert, alpha, platt)
        ji = plugin_threshold(cert, alpha, iso)
        for delta in DELTAS:
            for name, _ in RULES:
                rows.append(dict(repeat=rep, alpha=alpha, delta=delta, method=name,
                                 **report(walk(pv[name], delta), ev, pool)))
            rows.append(dict(repeat=rep, alpha=alpha, delta=delta, method="Platt calibration",
                             **report(jp, ev, pool)))
            rows.append(dict(repeat=rep, alpha=alpha, delta=delta, method="Isotonic calibration",
                             **report(ji, ev, pool)))
    if (rep + 1) % 10 == 0:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print("  repeat %3d  %.0fs" % (rep, time.time() - t0), flush=True)
pd.DataFrame(rows).to_csv(args.out, index=False)
print("wrote", args.out, "%.0fs" % (time.time() - t0), flush=True)
