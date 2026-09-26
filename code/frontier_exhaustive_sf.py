"""Exhaustive matched-risk frontier of the global path (revision of 2026-09-25, review item 10).

For every split and every certified arm in a run_second_family_v4.py output, the best correct-call yield on the
evaluation fold that ANY global score threshold could have achieved at a realised held-out risk not exceeding the
risk the arm realised.  Candidates are every distinct score value of the frozen score file (not only the refined
grid) plus abstention, so the frontier is the exact oracle over all global cut-offs.  Adds the column
frontier_all next to the grid-based frontier and writes <out>.

Usage: GOCO_DATA=... python frontier_exhaustive_sf.py --run funmap_v4_100.csv --out funmap_v4_100.csv
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from goco_second_family import FunMap, SEED0

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True); ap.add_argument("--out", required=True)
args = ap.parse_args()

t0 = time.time(); ds = FunMap()
CALL_I = ds.rows.i.to_numpy(); CALL_S = ds.rows.score.to_numpy(float); CALL_T = ds.rows["T"].to_numpy(np.int64)
order = np.lexsort((np.arange(ds.ncall), -CALL_S))          # score descending; ties are immaterial at distinct-score cuts
S_o, I_o, T_o = CALL_S[order], CALL_I[order], CALL_T[order]
# last position of every distinct score value in the sorted order
last_of_value = np.flatnonzero(np.r_[S_o[1:] != S_o[:-1], True])
print("units %d | calls %d | distinct scores %d | %.0fs" % (ds.n, ds.ncall, len(last_of_value), time.time() - t0), flush=True)

run = pd.read_csv(args.run)
reps = sorted(run.repeat.unique())
curves = {}
for rep in reps:
    perm = np.random.default_rng(SEED0 + rep).permutation(ds.n)
    nt, nc = round(0.10 * ds.n), round(0.70 * ds.n); ev = perm[nc:]
    inev = np.zeros(ds.n, bool); inev[ev] = True; ne = float(len(ev))
    C = np.zeros(ds.n, np.int64); Y = np.zeros(ds.n, np.int64)
    m = inev[I_o]                                           # only calls of evaluation genes move the held-out statistics
    risk_sum = 0.0; ysum = 0; risks = [0.0]; yields = [0]
    pos = 0
    idx_ev = np.flatnonzero(m)                              # positions (in sorted order) of evaluation-gene calls
    ptr = 0
    for cut in last_of_value:                               # threshold = this distinct score: release every call at or above it
        while ptr < len(idx_ev) and idx_ev[ptr] <= cut:
            k = idx_ev[ptr]; i = I_o[k]
            old = 0.0 if C[i] == 0 else (C[i] - Y[i]) / C[i]
            C[i] += 1; Y[i] += T_o[k]
            new = (C[i] - Y[i]) / C[i]
            risk_sum += new - old; ysum += T_o[k]; ptr += 1
        risks.append(risk_sum / ne); yields.append(ysum)
    r = np.asarray(risks); y = np.asarray(yields, float)
    # frontier(risk) = max yield over cut-offs with realised risk <= risk: sort by risk, running max of yield
    o = np.argsort(r, kind="stable"); r_s = r[o]; y_max = np.maximum.accumulate(y[o])
    curves[rep] = (r_s, y_max)
    if rep % 10 == 0: print("  repeat %3d  %.0fs" % (rep, time.time() - t0), flush=True)

def frontier_all(rep, fdp):
    r_s, y_max = curves[rep]; k = np.searchsorted(r_s, fdp + 1e-12, side="right")
    return int(y_max[k - 1]) if k > 0 else 0

run["frontier_all"] = [frontier_all(int(a), float(b)) if lab != "abstain" else 0 for a, b, lab in zip(run.repeat, run.fdp, run.label)]
run.to_csv(args.out, index=False)
print("wrote", args.out, "| frontier_all >= frontier (grid) in %.4f of rows | %.0fs" % ((run.frontier_all >= run.frontier).mean(), time.time() - t0), flush=True)
