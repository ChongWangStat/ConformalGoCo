"""Operational ranking-fold-measurability test for the FunMap ordering.

Validity of the fixed-sequence certificate requires the candidate path to be a function of the
ranking fold and frozen predictor output alone.  That is testable: scramble the TruePath labels of
every gene outside the ranking fold and confirm the per-call ordering statistic is bit-identical.

Note on what is legitimately used.  A call's source set is the annotated neighbours carrying the
term.  Those neighbour annotations are frozen predictor INPUT -- FunMap's own enrichment computes
the score from them -- and the author pipeline is leave-one-protein-out, so a gene's own
annotations never enter its own sources.  This is the property whose absence the ProteomeHD audit
documented, so it is worth testing rather than asserting.
"""
import sys
import numpy as np
sys.path.insert(0, r"C:/g1")
from goco_second_family import FunMap, SEED0
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


def fit_phat(ds, train):
    """Identical to the runner's: cross-fitted per-call model, ranking-fold calls only."""
    X, intr, unsup = ds.call_features(train)
    y = unsup.astype(int)
    p = np.full(ds.ncall, np.nan)
    for half in (0, 1):
        m = intr & ((ds.hcall < 0.5) if half == 0 else (ds.hcall >= 0.5))
        if m.sum() < 50 or len(np.unique(y[m])) < 2:
            continue
        s = StandardScaler().fit(X[m])
        clf = LogisticRegression(C=0.5, max_iter=1000).fit(s.transform(X[m]), y[m])
        tgt = (ds.hcall >= 0.5) if half == 0 else (ds.hcall < 0.5)
        p[tgt] = clf.predict_proba(s.transform(X[tgt]))[:, 1]
    return np.where(np.isnan(p), np.nanmean(p), p)

ds = FunMap()
rep = 0
perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
nt = round(0.10 * ds.n)
train, pool = perm[:nt], perm[nt:]

p_before = fit_phat(ds, train)

# scramble the labels of every non-ranking-fold gene's calls
rng = np.random.default_rng(12345)
i = ds.rows.i.to_numpy()
outside = ~np.isin(i, train)
T = ds.rows["T"].to_numpy().copy()
T[outside] = rng.permutation(T[outside])
ds.rows["T"] = T

p_after = fit_phat(ds, train)

pool_calls = np.isin(i, pool)
same = np.array_equal(p_before[pool_calls], p_after[pool_calls])
mx = float(np.max(np.abs(p_before[pool_calls] - p_after[pool_calls])))
print("pool calls compared      :", int(pool_calls.sum()))
print("ordering identical       :", same)
print("max |difference|         : %.3e" % mx)
print("PASS" if same else "FAIL")
sys.exit(0 if same else 1)
