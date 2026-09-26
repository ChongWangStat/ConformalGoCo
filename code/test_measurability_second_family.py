"""Operational ranking-fold-measurability test for the neighbourhood-family ordering (FunMap / STRING view).

Validity of the fixed-sequence certificate requires the candidate path to be a function of the ranking fold and
frozen predictor output alone.  That is testable: scramble the TruePath labels of every gene outside the ranking
fold and confirm that (i) the cross-fitted per-call statistic p-hat of every pool call, hence the GoCo-N order, and
(ii) the GoCo-M ordering statistic u_i of every affected pool gene in every block are bit-identical.

Note on what is legitimately used.  A call's source set is the annotated neighbours carrying the term.  Those
neighbour annotations are frozen predictor INPUT -- FunMap's own enrichment computes the score from them -- and the
author pipeline is leave-one-protein-out, so a gene's own annotations never enter its own sources.

Usage: GOCO_DATA=<funmap or string view dir> python test_measurability_second_family.py [split=0]
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from goco_second_family import FunMap, SEED0, COARSE, TRIGGER_FRAC, MIN_RANK_UNITS

rep = int(sys.argv[1]) if len(sys.argv) > 1 else 0
ds = FunMap()
perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
nt = round(0.10 * ds.n); train, pool = perm[:nt], perm[nt:]
i = ds.rows.i.to_numpy(); pool_calls = np.isin(i, pool)
Mc = ds.mats(COARSE); L, C = Mc["L"], Mc["C"]


def statistics():
    phat, _ = ds.fit_phat(train)
    Pw = ds.weighted_mats(COARSE, phat)
    us = {}
    for j in range(len(COARSE) - 1, 0, -1):
        if float(L[train, j - 1].mean()) <= TRIGGER_FRAC * 0.10: continue
        affected = C[:, j - 1] > C[:, j]
        if int(affected[train].sum()) < MIN_RANK_UNITS: continue
        nL, nH = C[:, j - 1].astype(float), C[:, j].astype(float)
        LhatL = np.divide(Pw[:, j - 1], nL, out=np.zeros(ds.n), where=nL > 0); LhatH = np.divide(Pw[:, j], nH, out=np.zeros(ds.n), where=nH > 0)
        tau = (nL - nH) - (Pw[:, j - 1] - Pw[:, j])
        u = np.where(affected, (LhatL - LhatH) / np.maximum(tau, 1e-3), np.inf)
        us[j] = u[pool[affected[pool]]].copy()
    return phat[pool_calls].copy(), us


p_before, u_before = statistics()
rng = np.random.default_rng(12345)
T = ds.rows["T"].to_numpy().copy(); T[pool_calls] = rng.permutation(T[pool_calls]); ds.set_truth(T)
p_after, u_after = statistics()

same_p = np.array_equal(p_before, p_after)
same_u = all(np.array_equal(u_before[j], u_after[j]) for j in u_before)
print("pool calls compared              :", int(pool_calls.sum()))
print("p-hat of pool calls identical    :", same_p, "(max |diff| = %.3e)" % float(np.max(np.abs(p_before - p_after))))
print("GoCo-M u_i of pool genes identical:", same_u, "over %d blocks" % len(u_before))
print("PASS" if (same_p and same_u) else "FAIL")
sys.exit(0 if (same_p and same_u) else 1)
