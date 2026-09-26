"""GoCo-N on the four dependency data sets: supported-gene counts (held-out genes with >=1 supported released call)
for the certified stopping points, plus ranking-fold diagnostics (AUC of score vs cross-fitted p-hat on T calls).
Splits byte-identical to the main analysis; the candidate family is gr.topm_path (top-m releases exactly m calls)."""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import rankdata
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goco_rerun as gr, goco_learned as gl

ap = argparse.ArgumentParser()
ap.add_argument("--datasets", default="Wainberg,Sanger,DRIVE,HAP1")
ap.add_argument("--repeats", type=int, default=100)
ap.add_argument("--out", default="wf_gocoN_detail.csv")
args = ap.parse_args()
ALPHAS = [0.05, 0.10, 0.20]; DELTAS = [0.50, 0.10]
MFRAC = gr.MFRAC_N

def auc(score, y):
    y = np.asarray(y, bool); n1 = y.sum(); n0 = len(y) - n1
    if n1 == 0 or n0 == 0: return np.nan
    r = rankdata(score)
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

dset, tset = gr.load_truth()
out, diag, t0 = [], [], time.time()
for name in args.datasets.split(","):
    ds = gr.Dataset(name, dset, tset)
    rows = ds.rows; s_all = np.round(rows.score.to_numpy(float), 6); keep = s_all >= ds.dense_floor
    call_i = rows.i.to_numpy()[keep]; call_T = rows["T"].to_numpy(np.int32)[keep]; call_rows = np.where(keep)[0]; ncall = len(call_rows)
    s_call = s_all[keep]
    hcall = gr.hashu(rows.gene.to_numpy()[keep], salt="half")
    for rep in range(args.repeats):
        perm = np.random.default_rng(gr.SEED0 + rep).permutation(ds.n)
        nt, nc = round(0.10 * ds.n), round(0.70 * ds.n)
        train, cert, ev = perm[:nt], perm[nt:nc], perm[nc:]; pool = perm[nt:]
        tr_set = np.zeros(ds.n, bool); tr_set[train] = True
        predict, _, _ = gl._fit_model(ds, train, "learned")
        phat = predict(call_rows)
        # ranking-fold diagnostics on T calls (p-hat cross-fitted inside T by gl._fit_model)
        mT = tr_set[call_i]; yT = call_T[mT] == 1
        diag.append(dict(dataset=name, repeat=rep, n_calls_T=int(mT.sum()), auc_score_T=auc(s_call[mT], yT), auc_phat_T=auc(-phat[mT], yT),
                         auc_score_E=auc(s_call[~tr_set[call_i]], call_T[~tr_set[call_i]] == 1), auc_phat_E=auc(-phat[~tr_set[call_i]], call_T[~tr_set[call_i]] == 1)))
        C, Y, L, edges = gr.topm_path(phat, hcall, call_i, call_T, ds.n, MFRAC)   # exact top-m sets (one-based ranks)
        for alpha in ALPHAS:
            for delta in DELTAS:
                sel = None
                for j in range(len(edges)):
                    p, _, _ = gr.clt_p(L[cert, j], alpha)
                    if p > delta: break
                    sel = j
                if sel is None:
                    out.append(dict(dataset=name, repeat=rep, alpha=alpha, delta=delta, arm="goco-N", label="abstain", fdp=0.0, pool_risk=0.0, correct=0, total=0, units=0, supp_units=0)); continue
                out.append(dict(dataset=name, repeat=rep, alpha=alpha, delta=delta, arm="goco-N", label="top%d" % edges[sel],
                                fdp=float(L[ev, sel].mean()), pool_risk=float(L[pool, sel].mean()),
                                correct=int(Y[ev, sel].sum()), total=int(C[ev, sel].sum()), units=int((C[ev, sel] > 0).sum()), supp_units=int((Y[ev, sel] > 0).sum())))
        if rep % 10 == 0:
            print("  [%s] repeat %3d  %.0fs" % (name, rep, time.time() - t0), flush=True)
            pd.DataFrame(out).to_csv(args.out, index=False); pd.DataFrame(diag).to_csv(args.out.replace(".csv", "_diag.csv"), index=False)
    pd.DataFrame(out).to_csv(args.out, index=False); pd.DataFrame(diag).to_csv(args.out.replace(".csv", "_diag.csv"), index=False)
print("done %.0fs" % (time.time() - t0), flush=True)
