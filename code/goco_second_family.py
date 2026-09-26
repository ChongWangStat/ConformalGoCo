"""GoCo for the neighbourhood-structure family (FunMap development data; STRING validation through the same view).

Port of the GoCo construction -- global grid, plus a ranking-fold-measurable refinement inside a
step, certified by fixed-sequence LTT with the finite-population normal p-value -- to a predictor
whose evidence object is a top-50 random-walk neighbourhood rather than a co-essential module.

Why the GoCo certificate and not the pair machinery: GoCo's guarantee is finite-population. The
only randomness is the split of a fixed panel; independence across genes is never assumed. Network
dependence among genes therefore does not enter, which is exactly the situation here.

Two structure-driven departures from GoCo as published:
  (1) Sources are RANKED by random-walk proximity, so attribution is proximity-weighted rather
      than uniform over module members (GoCo Eq. 4 gives every member 1/|R| of the evidence).
  (2) FunMap scores carry no material ties, so there is no tie-block to specialise a model to.
      One per-call model is fitted per split over all candidate calls, rather than one per block.

Revision of 2026-09-25 (pre-submission review), aligning the implementation with the module-family code
and with Supplementary S3.3 / S4.2:
  * bh_rank is the call's rank among the gene's calls WITHIN ITS ASPECT (the ontology the enrichment was run in).
  * fit_phat(): the classifier AND the two label-derived rate features are cross-fitted between two
    hash-defined halves of the ranking fold (each ranking-fold call is scored by a model whose rate tables and
    coefficients come from the other half only); every other call (certification and evaluation genes) is scored
    by the model fitted on the whole ranking fold, as in goco_learned.py.
  * weighted_mats(): per-unit cumulative sums of p-hat along a grid, so that the GoCo-M ordering statistic is
    exactly Eq. (6): Delta_hat_i = L_hat(S_i^L) - L_hat(S_i^H) with L_hat(S) the mean p-hat over S (0 if empty),
    tau_hat_i = sum over the added calls of (1 - p-hat).
  * topm_path(): the GoCo-N candidate family with one-based ranks, so that top-m releases exactly m calls.
  * Ties in every ordering are broken by the exact lexicographic key (statistic, frozen hash, index); no tolerance.
Seeds: split s of this family uses default_rng(20260910 + s); the module family uses 20260905 + s.
"""
from __future__ import annotations
import os
import numpy as np, pandas as pd, hashlib
from pathlib import Path
from scipy import sparse
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = Path(os.environ.get("GOCO_DATA", ROOT / "funmap"))
SEED0 = 20260910
FLOOR = 1.0
QGRID = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20,
                  0.225, 0.25, 0.275, 0.30, 0.35, 0.40, 0.50, 0.60, 0.75])
TRIGGER_FRAC = 0.8
MIN_RANK_UNITS = 5
MFRAC_N = np.unique(np.round(np.geomspace(0.0005, 1.0, 60), 6))   # GoCo-N candidate sizes as fractions of the call pool
COARSE = np.unique(np.round(np.concatenate([np.geomspace(1.0, 120.0, 30), [1.30103, 2.0]]), 6))
DENSE = np.unique(np.round(np.concatenate([COARSE, np.geomspace(1.0, 120.0, 300)]), 6))


def hashu(names, salt=""):
    return np.array([int(hashlib.sha256((salt + str(x)).encode()).hexdigest()[:12], 16) / 16**12
                     for x in names])


def clt_p(y, alpha):
    """One-sided finite-population normal p-value (GoCo Eq. 9)."""
    m = len(y)
    if m == 0:
        return 1.0
    r, s = float(np.mean(y)), float(np.std(y))
    return float(norm.cdf((r - alpha) / (max(s, 1e-6) / np.sqrt(m))))


def topm_path(phat, hcall, call_i, call_T, n, mfrac=MFRAC_N):
    """GoCo-N candidate family: nested top-m sets of calls in the total order (p-hat, frozen call hash, call index).

    Column j holds exactly m_j = edges[j] calls: ranks are one-based and a call of rank r belongs to the top-m set
    iff r <= m, so the label 'top{m}' names the number of released calls.  Returns the per-unit released count C,
    supported count Y, unit loss L (0 when the unit releases nothing) and the size vector edges.
    Byte-identical to goco_rerun.topm_path of the module-family code."""
    ncall = len(phat)
    key = np.lexsort((np.arange(ncall), hcall, phat)); rank = np.empty(ncall, np.int64); rank[key] = np.arange(1, ncall + 1)
    edges = np.unique(np.maximum(1, np.floor(np.asarray(mfrac) * ncall).astype(np.int64)))
    b = np.searchsorted(edges, rank, side="left"); ok = b < len(edges); shape = (n, len(edges))
    cum = lambda v: np.cumsum(sparse.coo_matrix((v, (call_i[ok], b[ok])), shape=shape).tocsr().toarray(), axis=1)
    C = cum(np.ones(int(ok.sum()), np.int32)).astype(np.int32); Y = cum(np.asarray(call_T, np.int32)[ok]).astype(np.int32)
    L = np.divide(C - Y, C, out=np.zeros(shape), where=C > 0)
    return C, Y, L, edges


def admitted_prefix(u, hh, affected, train, q):
    """Q_{j,q}: affected units whose key (u, hash, index) is at or below the key of the k-th affected ranking-fold unit,
    k = max(1, floor(q |B_j ∩ T|)), in the exact total order.  Same rule as goco_rerun.selmask; q >= 1 admits B_j."""
    a = train[affected[train]]
    if len(a) == 0: return np.zeros(len(u), bool)
    if q >= 1: return affected.copy()
    idx = np.arange(len(u))
    oo = np.lexsort((idx[a], hh[a], u[a])); k = max(1, int(np.floor(q * len(a)))); i = a[oo[k - 1]]; r = u[i]; h = hh[i]
    return affected & ((u < r) | ((u == r) & (hh < h)) | ((u == r) & (hh == h) & (idx <= i)))


class FunMap:
    """Frozen predictor output plus the objects GoCo needs: per-call sources and truth."""

    def __init__(self, floor=FLOOR):
        sc = pd.read_csv(P / "funmap_gene_go_scores.csv.gz",
                         usecols=["gene", "go_id", "aspect", "score", "neighbor_hits",
                                  "annotated_neighbors", "term_background_count"])
        sc = sc.dropna(subset=["gene", "go_id", "score"])
        sc = sc.sort_values(["gene", "score", "go_id"], ascending=[True, False, True])
        sc = sc.drop_duplicates(["gene", "go_id"], keep="first").reset_index(drop=True)

        tpath = pd.read_csv(P / "go_truth_true_path.csv.gz", usecols=["gene", "go_id"]).drop_duplicates()
        direct = pd.read_csv(P / "go_truth_direct.csv.gz", usecols=["gene", "go_id"]).drop_duplicates()
        edges = pd.read_csv(P / "funmap_edges.csv.gz")
        net = pd.unique(pd.concat([edges.gene_a, edges.gene_b], ignore_index=True))
        self.units = np.array(sorted(set(net) & set(direct.gene.unique())))
        self.n = len(self.units)
        self.uidx = {u: i for i, u in enumerate(self.units)}

        sc["T"] = pd.MultiIndex.from_frame(sc[["gene", "go_id"]]).isin(pd.MultiIndex.from_frame(tpath))
        sc = sc[sc.gene.isin(self.uidx) & (sc.score >= floor)].reset_index(drop=True)
        sc["i"] = sc.gene.map(self.uidx).to_numpy()
        # rank of the call among the gene's calls within its aspect (the ontology the BH adjustment was run in)
        sc["bh_rank"] = sc.groupby(["gene", "aspect"]).score.rank(ascending=False, method="first")
        sc["call"] = np.arange(len(sc))
        self.rows = sc
        self.ncall = len(sc)

        # ---- proximity-weighted source attribution, built once as a sparse calls x neighbours matrix
        nb = pd.read_csv(P / "author_top50_neighborhoods_lopo.csv.gz",
                         usecols=["focal_gene", "neighbor_gene", "author_rank"])
        nb = nb[nb.focal_gene.isin(self.uidx)]
        nbr_names = np.array(sorted(set(nb.neighbor_gene) | set(direct.gene)))
        nix = {g: k for k, g in enumerate(nbr_names)}
        self.nbr_names, self.nnbr = nbr_names, len(nbr_names)

        key = sc.gene.astype(str) + "\x00" + sc.go_id.astype(str)
        cmap = dict(zip(key, sc.call.to_numpy()))
        ann = direct.rename(columns={"gene": "neighbor_gene"})
        chunks = []
        for _, part in nb.groupby(np.arange(len(nb)) // 200000):
            s = part.merge(ann, on="neighbor_gene")
            k = s.focal_gene.astype(str) + "\x00" + s.go_id.astype(str)
            cid = k.map(cmap)
            s = s.loc[cid.notna()]
            chunks.append(pd.DataFrame({"call": cid[cid.notna()].to_numpy(np.int64),
                                        "nbr": s.neighbor_gene.map(nix).to_numpy(np.int64),
                                        "w": 1.0 / np.log2(1.0 + s.author_rank.to_numpy(float)),
                                        "rank": s.author_rank.to_numpy(float)}))
        S = pd.concat(chunks, ignore_index=True)
        self.A = sparse.csr_matrix((S.w.to_numpy(), (S.call.to_numpy(), S.nbr.to_numpy())),
                                   shape=(self.ncall, self.nnbr))
        self.wmass = np.asarray(self.A.sum(axis=1)).ravel()
        self.nsrc = np.asarray((self.A > 0).sum(axis=1)).ravel()
        self.minrank = (S.groupby("call")["rank"].min()
                        .reindex(np.arange(self.ncall)).fillna(50.0).to_numpy())
        self.hh = hashu(self.units)                                  # frozen unit hash (tie-breaker for genes)
        assert len(np.unique(self.hh)) == self.n, "unit hash collision: the index tie-break would be reached"
        self.hcall = hashu(self.rows.gene.to_numpy(), "half")        # frozen gene-level hash: cross-fitting halves and call ties
        self.half_unit = (hashu(self.units, "half") < 0.5)           # the same hash at unit level
        # frozen (label-free) part of the per-call design matrix, built once
        r = self.rows
        self._ncalls = r.groupby("i").size()
        self._frozen = np.column_stack([
            np.log1p(r.score.to_numpy(float)),
            np.log1p(r.neighbor_hits.to_numpy(float)),
            np.log1p(r.annotated_neighbors.to_numpy(float)),
            np.log1p(r.bh_rank.to_numpy(float)),
            np.log1p(self.nsrc.astype(float)),
            self.wmass,
            self.minrank,
            np.log1p(r.i.map(self._ncalls).to_numpy(float)),
        ])
        self._unsup = (~r["T"].to_numpy(bool)).astype(float)

    def set_truth(self, T):
        """Replace the per-call TruePath indicator (used by the measurability test to scramble non-ranking-fold labels)."""
        self.rows["T"] = np.asarray(T, bool)
        self._unsup = (~self.rows["T"].to_numpy(bool)).astype(float)

    def mats(self, grid):
        grid = np.asarray(grid, float)
        mb = np.searchsorted(grid, self.rows.score.to_numpy(), side="right") - 1
        keep = mb >= 0
        r, m = self.rows.i.to_numpy()[keep], mb[keep]
        t = self.rows["T"].to_numpy(np.int32)[keep]
        shape = (self.n, len(grid))
        cum = lambda v: np.cumsum(sparse.coo_matrix((v, (r, m)), shape=shape).tocsr().toarray()[:, ::-1],
                                  axis=1)[:, ::-1]
        C = cum(np.ones(int(keep.sum()), np.int32)).astype(np.int32)
        Y = cum(t).astype(np.int32)
        L = np.divide(C - Y, C, out=np.zeros(shape), where=C > 0)
        return dict(C=C, Y=Y, L=L, grid=grid)

    def weighted_mats(self, grid, w):
        """Per-unit sum of a per-call weight w over the calls with score >= grid[k], for every k (same layout as mats)."""
        grid = np.asarray(grid, float)
        mb = np.searchsorted(grid, self.rows.score.to_numpy(), side="right") - 1
        keep = mb >= 0
        r, m = self.rows.i.to_numpy()[keep], mb[keep]
        shape = (self.n, len(grid))
        return np.cumsum(sparse.coo_matrix((np.asarray(w, float)[keep], (r, m)), shape=shape).tocsr().toarray()[:, ::-1],
                         axis=1)[:, ::-1]

    # ---------------------------------------------------------------- features
    def rate_features(self, fit_units):
        """The two label-derived features, computed from the calls of fit_units ONLY and applied to every call:
        the term's shrunk unsupported rate and the proximity-weighted shrunk unsupported rate of the call's sources.
        Measurable with respect to the ranking fold whenever fit_units is a subset of it."""
        r = self.rows
        infit = np.isin(r.i.to_numpy(), fit_units)
        unsup = self._unsup
        glob = float(unsup[infit].mean()) if infit.any() else 0.5
        # term-level unsupported rate, shrunk toward the global rate by two pseudo-calls
        g = pd.DataFrame({"t": r.go_id.to_numpy()[infit], "u": unsup[infit]}).groupby("t").u.agg(["sum", "count"])
        trate_map = ((g["sum"] + 2 * glob) / (g["count"] + 2))
        trate = r.go_id.map(trate_map).fillna(glob).to_numpy()
        # proximity-weighted neighbour reliability: how often this neighbour's terms went unsupported on fit calls
        Afit = self.A[infit]
        num = np.asarray(Afit.T.dot(unsup[infit])).ravel()
        den = np.asarray(Afit.sum(axis=0)).ravel()
        srate = (num + 2 * glob) / (den + 2)
        wsrate = np.asarray(self.A.dot(srate)).ravel() / np.maximum(self.wmass, 1e-9)
        wsrate[self.wmass <= 0] = glob
        return np.column_stack([wsrate, trate]), infit

    def call_features(self, fit_units, include_background_count=False):
        """Per-call design matrix: frozen predictor output plus the rate features of fit_units.
        term_background_count is EXCLUDED by default: with the focal gene removed from the enrichment background it
        equals K_g - 1{focal gene Direct-annotated with g}, i.e. it encodes the focal gene's own label."""
        R, infit = self.rate_features(fit_units)
        cols = [self._frozen]
        if include_background_count:
            cols.append(np.log1p(self.rows.term_background_count.to_numpy(float))[:, None])
        cols.append(R)
        return np.column_stack(cols), infit, self._unsup

    def fit_phat(self, train, include_background_count=False, C=0.5):
        """Cross-fitted per-call probability of 'call unsupported' (TruePath).

        Whole-ranking-fold model -> every call of a certification or evaluation gene.  Half h of the ranking fold
        (frozen gene hash) -> scored by the model fitted on the other half, whose rate features are computed from
        that other half only, so no ranking-fold call sees its own label in features or coefficients.
        Same scheme as goco_learned._fit_model of the module family."""
        y = self._unsup.astype(int)
        tr_set = np.zeros(self.n, bool); tr_set[train] = True
        inT = tr_set[self.rows.i.to_numpy()]
        X, infit, _ = self.call_features(train, include_background_count)
        s = StandardScaler().fit(X[infit]); clf = LogisticRegression(C=C, max_iter=1000).fit(s.transform(X[infit]), y[infit])
        p = clf.predict_proba(s.transform(X))[:, 1]
        pT = np.full(self.ncall, np.nan)
        for h in (0, 1):
            mine = train[self.half_unit[train] == bool(h)]; other = train[self.half_unit[train] != bool(h)]
            if len(mine) == 0 or len(other) == 0: continue
            Xh, infith, _ = self.call_features(other, include_background_count)
            if infith.sum() < 50 or len(np.unique(y[infith])) < 2: continue
            sh = StandardScaler().fit(Xh[infith]); ch = LogisticRegression(C=C, max_iter=1000).fit(sh.transform(Xh[infith]), y[infith])
            tgt = np.isin(self.rows.i.to_numpy(), mine)
            pT[tgt] = ch.predict_proba(sh.transform(Xh[tgt]))[:, 1]
        fill = np.nanmean(pT[inT]) if np.isfinite(pT[inT]).any() else float(p[inT].mean())
        pT[inT & ~np.isfinite(pT)] = fill
        p[inT] = pT[inT]
        return p, inT
