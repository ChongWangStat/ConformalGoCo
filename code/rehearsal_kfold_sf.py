"""Ranking-fold selection rule for the neighbourhood family (FunMap / STRING) with K-fold cross-fitting inside T.
Rehearsal only: certified results are those of run_goco_second_family.py and are untouched. For every split the
T genes are hashed into K folds; the calls of fold h are scored (label-derived features AND classifier) by a fit
that uses only the T genes of the other K-1 folds, so no ranking-fold call sees its own label. Both paths are then
walked on T alone and stopped at the last candidate with T-mean loss <= alpha; the larger T yield wins.
Splits, grids, trigger rules and candidate construction are identical to run_goco_second_family.py (Eq. (6) statistic,
exact tie order, top-m releases exactly m calls; revision of 2026-09-25)."""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from goco_second_family import (FunMap, QGRID, TRIGGER_FRAC, MIN_RANK_UNITS, SEED0, COARSE, MFRAC_N, hashu, topm_path, admitted_prefix)
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=100); ap.add_argument("--start", type=int, default=0)
ap.add_argument("--K", type=int, default=10); ap.add_argument("--out", required=True)
args = ap.parse_args(); K = args.K
ALPHAS = [0.05, 0.10, 0.20]

print("loading ...", flush=True); t0 = time.time()
ds = FunMap(); Mc = ds.mats(COARSE)
CALL_I = ds.rows.i.to_numpy(); CALL_S = ds.rows.score.to_numpy(float); CALL_T = ds.rows["T"].to_numpy(np.int32)
gfold_all = np.minimum((hashu(ds.units, "xfit") * K).astype(int), K - 1)      # gene-level fold id, frozen
print("units %d | calls %d | %.0fs" % (ds.n, ds.ncall, time.time() - t0), flush=True)


def fit_phat_kfold(train):
    """p-hat for every call: T calls K-fold cross-fitted inside T (label-derived features AND classifier from the other
    folds); all other calls scored by the fit on the whole ranking fold, as pool genes are in the module family."""
    tr_set = np.zeros(ds.n, bool); tr_set[train] = True; inT = tr_set[CALL_I]
    X, intr, unsup = ds.call_features(train); y = unsup.astype(int)
    s = StandardScaler().fit(X[intr]); clf = LogisticRegression(C=0.5, max_iter=1000).fit(s.transform(X[intr]), y[intr])
    p = clf.predict_proba(s.transform(X))[:, 1]
    pT = np.full(ds.ncall, np.nan)
    for h in range(K):
        other = train[gfold_all[train] != h]; mine = train[gfold_all[train] == h]
        if len(other) == 0 or len(mine) == 0: continue
        Xh, intrh, unsuph = ds.call_features(other); yh = unsuph.astype(int)
        if intrh.sum() < 50 or len(np.unique(yh[intrh])) < 2: continue
        sh = StandardScaler().fit(Xh[intrh]); ch = LogisticRegression(C=0.5, max_iter=1000).fit(sh.transform(Xh[intrh]), yh[intrh])
        tgt = np.isin(CALL_I, mine); pT[tgt] = ch.predict_proba(sh.transform(Xh[tgt]))[:, 1]
    fill = np.nanmean(pT[inT]) if np.isfinite(pT[inT]).any() else 0.5
    pT[inT & ~np.isfinite(pT)] = fill
    p[inT] = pT[inT]
    return p, inT


def rank_path_T(phat, train):
    """GoCo-N path over ALL calls (top-m by p-hat, exact ties, same size grid as the main runner), read on T genes."""
    C, Y, L, edges = topm_path(phat, ds.hcall, CALL_I, CALL_T, ds.n, MFRAC_N)
    return [(float(L[train, j].mean()), int(Y[train, j].sum()), f"top{int(edges[j])}") for j in range(len(edges))]


def goco_M_T(train, alpha, Pw):
    """GoCo-M candidates evaluated on T only (same construction as cand_goco in the main runner)."""
    L, grid, C = Mc["L"], Mc["grid"], Mc["C"]; out = []
    for j in range(len(grid) - 1, 0, -1):
        out.append((float(L[train, j].mean()), int(Mc["Y"][train, j].sum()), f"global_{grid[j]:g}"))
        if float(L[train, j - 1].mean()) <= TRIGGER_FRAC * alpha: continue
        lo, hi = float(grid[j - 1]), float(grid[j])
        affected = C[:, j - 1] > C[:, j]
        if int(affected[train].sum()) < MIN_RANK_UNITS: continue
        nL, nH = C[:, j - 1].astype(float), C[:, j].astype(float)
        LhatL = np.divide(Pw[:, j - 1], nL, out=np.zeros(ds.n), where=nL > 0); LhatH = np.divide(Pw[:, j], nH, out=np.zeros(ds.n), where=nH > 0)
        tau = (nL - nH) - (Pw[:, j - 1] - Pw[:, j])
        u = np.where(affected, (LhatL - LhatH) / np.maximum(tau, 1e-3), np.inf)
        for q in QGRID:
            adm = admitted_prefix(u, ds.hh, affected, train, float(q))
            loss = np.where(adm, L[:, j - 1], L[:, j])[train]; Yv = np.where(adm, Mc["Y"][:, j - 1], Mc["Y"][:, j])[train]
            out.append((float(loss.mean()), int(Yv.sum()), f"adapt_{hi:g}_to_{lo:g}_q{q:g}"))
    return out


def rehearse(cands, alpha):
    last = None
    for c in cands:
        if c[0] > alpha: break
        last = c
    return (last[1], last[2]) if last is not None else (0, "abstain")


rows = []
for rep in range(args.start, args.start + args.repeats):
    perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
    nt = round(0.10 * ds.n); train = perm[:nt]
    phat, inT = fit_phat_kfold(train)
    Pw = ds.weighted_mats(COARSE, phat)
    candN = rank_path_T(phat, train)
    for alpha in ALPHAS:
        yN, labN = rehearse(candN, alpha)
        yM, labM = rehearse(goco_M_T(train, alpha, Pw), alpha)
        rows.append(dict(repeat=rep, alpha=alpha, K=K, yield_M_T=yM, yield_N_T=yN, choose_N=int(yN > yM), label_M=labM, label_N=labN))
    print("  repeat %3d  %.0fs" % (rep, time.time() - t0), flush=True)
    if (rep + 1) % 10 == 0: pd.DataFrame(rows).to_csv(args.out, index=False)
pd.DataFrame(rows).to_csv(args.out, index=False)
print("wrote", args.out, "%.0fs" % (time.time() - t0), flush=True)
