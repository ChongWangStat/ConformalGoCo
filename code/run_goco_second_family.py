"""Neighbourhood-structure family runner (FunMap / STRING through the same view): all arms, both delta, plus
(i) pool risk and held-out exceedance, (ii) matched-risk frontier of the refined global path,
(iii) the ranking-fold selection rule between GoCo-M and GoCo-N (two-half cross-fitting; the ten-fold version
used in the paper is rehearsal_kfold_sf.py).

Arms: global-coarse (GoDag on the 32-threshold grid), global-dense (GoDag on the refined grid), goco-M (grid path
with partial admission ordered by Eq. (6)), goco-N (nested top-m calls by cross-fitted p-hat).  Every arm certifies
on the 60% certification fold and is evaluated on the 30% evaluation fold; the pool of the certificate is their union.

Revision of 2026-09-25 (pre-submission review): p-hat is cross-fitted in features and coefficients (FunMap.fit_phat);
the GoCo-M statistic is exactly Eq. (6) (Delta_hat_i = L_hat(S_i^L) - L_hat(S_i^H), tau_hat_i = expected supported
added calls) with the exact (u, hash, index) tie order; top-m releases exactly m calls; the unused 'goco-score' arm
is dropped.  Splits: default_rng(20260910 + s).
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from goco_second_family import (FunMap, clt_p, QGRID, TRIGGER_FRAC, MIN_RANK_UNITS, SEED0, COARSE, DENSE, MFRAC_N,
                                topm_path, admitted_prefix)

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=100)
ap.add_argument("--start", type=int, default=0)
ap.add_argument("--alphas", default="0.05,0.10,0.20")
ap.add_argument("--deltas", default="0.50,0.10")
ap.add_argument("--leak", action="store_true", help="include term_background_count (original v3 feature set)")
ap.add_argument("--out", required=True)
args = ap.parse_args()
ALPHAS = [float(x) for x in args.alphas.split(",")]; DELTAS = [float(x) for x in args.deltas.split(",")]

print("loading ...", flush=True); t0 = time.time()
ds = FunMap(); Mc, Md = ds.mats(COARSE), ds.mats(DENSE)
CALL_I = ds.rows.i.to_numpy(); CALL_S = ds.rows.score.to_numpy(float); CALL_T = ds.rows["T"].to_numpy(np.int32)
print("units %d | calls %d | %.0fs" % (ds.n, ds.ncall, time.time() - t0), flush=True)


def cand_global(M):
    L, grid = M["L"], M["grid"]
    return [{"loss": L[:, j], "C": M["C"][:, j], "Y": M["Y"][:, j], "label": f"global_{grid[j]:g}"} for j in range(len(grid) - 1, -1, -1)]


def cand_rank(R):
    C, Y, L, edges = R
    return [{"loss": L[:, j], "C": C[:, j], "Y": Y[:, j], "label": f"top{int(edges[j])}"} for j in range(len(edges))]


def cand_goco(M, Pw, train, alpha):
    """GoCo-M candidates on the coarse grid: after each global threshold grid[j], partial admission of the block
    grid[j] -> grid[j-1] in the order of Eq. (6), u_i = Delta_hat_i / max(tau_hat_i, 1e-3), on the ranking-fold
    quantile grid QGRID.  Pw[:, k] is the per-unit sum of p-hat over the calls with score >= grid[k]."""
    L, grid, C = M["L"], M["grid"], M["C"]; cands = []
    for j in range(len(grid) - 1, 0, -1):
        cands.append({"loss": L[:, j], "C": C[:, j], "Y": M["Y"][:, j], "label": f"global_{grid[j]:g}"})
        if float(L[train, j - 1].mean()) <= TRIGGER_FRAC * alpha: continue      # ranking-fold band rule
        lo, hi = float(grid[j - 1]), float(grid[j])
        affected = C[:, j - 1] > C[:, j]                                        # B_j: units with a call in [lo, hi)
        if int(affected[train].sum()) < MIN_RANK_UNITS: continue
        nL, nH = C[:, j - 1].astype(float), C[:, j].astype(float)
        LhatL = np.divide(Pw[:, j - 1], nL, out=np.zeros(ds.n), where=nL > 0)
        LhatH = np.divide(Pw[:, j], nH, out=np.zeros(ds.n), where=nH > 0)
        dhat = LhatL - LhatH                                                    # Delta_hat_i = L_hat(S^L) - L_hat(S^H)
        tau = (nL - nH) - (Pw[:, j - 1] - Pw[:, j])                             # sum over added calls of (1 - p-hat)
        u = np.where(affected, dhat / np.maximum(tau, 1e-3), np.inf)
        for q in QGRID:
            adm = admitted_prefix(u, ds.hh, affected, train, float(q))
            loss = np.where(adm, L[:, j - 1], L[:, j]); Cq = np.where(adm, C[:, j - 1], C[:, j]); Yq = np.where(adm, M["Y"][:, j - 1], M["Y"][:, j])
            cands.append({"loss": loss, "C": Cq, "Y": Yq, "label": f"adapt_{hi:g}_to_{lo:g}_q{q:g}"})
    return cands


def certify(cands, cert, alpha, delta):
    last = None
    for c in cands:
        if clt_p(c["loss"][cert], alpha) > delta: break
        last = c
    return last


def rehearse(cands, train, alpha):
    """Ranking-fold walk: last candidate whose T-mean loss <= alpha; return its supported yield on T."""
    last = None
    for c in cands:
        if float(c["loss"][train].mean()) > alpha: break
        last = c
    return (int(last["Y"][train].sum()) if last is not None else 0), (last["label"] if last is not None else "abstain")


def frontier(Le, Ye, risk):
    """Best refined-global supported yield on the evaluation fold at realised risk <= risk."""
    ok = Le <= risk + 1e-12
    return int(Ye[ok].max()) if ok.any() else 0


def report(res, ev, pool, Le, Ye):
    if res is None: return dict(label="abstain", fdp=0.0, pool_risk=0.0, correct=0, total=0, units=0, supp_units=0, frontier=0)
    fdp = float(res["loss"][ev].mean())
    return dict(label=res["label"], fdp=fdp, pool_risk=float(res["loss"][pool].mean()),
                correct=int(res["Y"][ev].sum()), total=int(res["C"][ev].sum()), units=int((res["C"][ev] > 0).sum()),
                supp_units=int((res["Y"][ev] > 0).sum()), frontier=frontier(Le, Ye, fdp))


rows, sel_rows = [], []
for rep in range(args.start, args.start + args.repeats):
    perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
    nt, nc = round(0.10 * ds.n), round(0.70 * ds.n)
    train, cert, ev = perm[:nt], perm[nt:nc], perm[nc:]; pool = perm[nt:]
    phat, _ = ds.fit_phat(train, include_background_count=args.leak)
    Pw = ds.weighted_mats(COARSE, phat)
    cg, cd, cr = cand_global(Mc), cand_global(Md), cand_rank(topm_path(phat, ds.hcall, CALL_I, CALL_T, ds.n, MFRAC_N))
    Le = Md["L"][ev].mean(axis=0); Ye = Md["Y"][ev].sum(axis=0)
    for alpha in ALPHAS:
        cm = cand_goco(Mc, Pw, train, alpha)
        yM, labM = rehearse(cm, train, alpha); yN, labN = rehearse(cr, train, alpha)
        sel_rows.append(dict(repeat=rep, alpha=alpha, yield_M_T=yM, yield_N_T=yN, choose_N=int(yN > yM), label_M=labM, label_N=labN))
        for delta in DELTAS:
            for name, cands in [("global-coarse", cg), ("global-dense", cd), ("goco-M", cm), ("goco-N", cr)]:
                rows.append(dict(repeat=rep, alpha=alpha, delta=delta, arm=name, **report(certify(cands, cert, alpha, delta), ev, pool, Le, Ye)))
    print("  repeat %3d  %.0fs" % (rep, time.time() - t0), flush=True)
    if (rep + 1) % 10 == 0:
        pd.DataFrame(rows).to_csv(args.out, index=False); pd.DataFrame(sel_rows).to_csv(args.out.replace(".csv", "_selection.csv"), index=False)
pd.DataFrame(rows).to_csv(args.out, index=False); pd.DataFrame(sel_rows).to_csv(args.out.replace(".csv", "_selection.csv"), index=False)
print("wrote", args.out, "%.0fs" % (time.time() - t0), flush=True)
