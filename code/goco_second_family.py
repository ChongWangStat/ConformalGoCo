"""GoCo for the FunMap random-walk-neighbourhood predictor.

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
"""
from __future__ import annotations
import os
import numpy as np, pandas as pd, hashlib
from pathlib import Path
from scipy import sparse
from scipy.stats import norm

ROOT = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = Path(os.environ.get("GOCO_DATA", ROOT / "funmap"))
SEED0 = 20260910
FLOOR = 1.0
QGRID = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20,
                  0.225, 0.25, 0.275, 0.30, 0.35, 0.40, 0.50, 0.60, 0.75])
TRIGGER_FRAC = 0.8
MIN_RANK_UNITS = 5


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


class FunMap:
    """Frozen predictor output plus the objects GoCo needs: per-call sources and truth."""

    def __init__(self, floor=FLOOR):
        sc = pd.read_csv(P / "funmap_gene_go_scores.csv.gz",
                         usecols=["gene", "go_id", "score", "neighbor_hits",
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
        sc["bh_rank"] = sc.groupby("gene").score.rank(ascending=False, method="first")
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
        self.hh = hashu(self.units)
        self.hcall = hashu(self.rows.gene.to_numpy(), "half")

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

    # ---------------------------------------------------------------- features
    def call_features(self, train_units):
        """Per-call design matrix. Frozen predictor output plus ranking-fold-label rates only.

        Everything label-derived is computed from calls of ranking-fold genes, so the matrix is
        measurable with respect to the ranking fold and Corollary 1 applies unchanged.
        """
        r = self.rows
        intr = np.isin(r.i.to_numpy(), train_units)
        unsup = (~r["T"].to_numpy(bool)).astype(float)
        glob = float(unsup[intr].mean()) if intr.any() else 0.5

        # term-level unsupported rate, shrunk toward the global rate by two pseudo-calls
        g = pd.DataFrame({"t": r.go_id.to_numpy()[intr], "u": unsup[intr]}).groupby("t").u.agg(["sum", "count"])
        trate_map = ((g["sum"] + 2 * glob) / (g["count"] + 2))
        trate = r.go_id.map(trate_map).fillna(glob).to_numpy()

        # proximity-weighted neighbour reliability: how often this neighbour's terms went unsupported
        # on ranking-fold genes.  A^T restricted to ranking-fold calls, then shrunk.
        Atr = self.A[intr]
        num = np.asarray(Atr.T.dot(unsup[intr])).ravel()
        den = np.asarray(Atr.sum(axis=0)).ravel()
        srate = (num + 2 * glob) / (den + 2)
        wsrate = np.asarray(self.A.dot(srate)).ravel() / np.maximum(self.wmass, 1e-9)
        wsrate[self.wmass <= 0] = glob

        ncalls = r.groupby("i").size()
        return np.column_stack([
            np.log1p(r.score.to_numpy(float)),
            np.log1p(r.neighbor_hits.to_numpy(float)),
            np.log1p(r.annotated_neighbors.to_numpy(float)),
            np.log1p(r.term_background_count.to_numpy(float)),
            np.log1p(r.bh_rank.to_numpy(float)),
            np.log1p(self.nsrc.astype(float)),
            self.wmass,
            self.minrank,
            wsrate,
            trate,
            np.log1p(r.i.map(ncalls).to_numpy(float)),
        ]), intr, unsup
