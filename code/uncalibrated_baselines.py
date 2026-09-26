"""Release rules in current use, on the four dependency data sets, over the same 100 splits.

Three rules, none of which certifies anything:
  Fmax      the global score threshold maximising the protein-centric F-measure (CAFA convention) on the
            same 70% of the panel the global calibrators use, applied to the evaluation fold;
  Top-1     the single highest-scoring candidate term per gene (ties broken by the frozen gene--term hash);
  Top-3     the three highest-scoring candidate terms per gene.
Reported on the evaluation fold: mean gene-level TruePath FDP (abstention 0), the fraction of splits whose
held-out FDP exceeds alpha, genes with a call, released calls and reference-supported released calls.
Candidate pool, splits and truth are those of the main analysis.
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goco_rerun as gr

ap = argparse.ArgumentParser()
ap.add_argument("--datasets", default="Wainberg,Sanger,DRIVE,HAP1")
ap.add_argument("--repeats", type=int, default=100)
ap.add_argument("--out", default="wf_uncalibrated.csv")
args = ap.parse_args()
ALPHAS = [0.05, 0.10, 0.20]

dset, tset = gr.load_truth()
# |T_i|: number of TruePath-annotated terms of each gene in the frozen reference
tp_by_gene = pd.Series(list(tset)).apply(lambda x: x[0]).value_counts() if False else None
_tp = pd.DataFrame(list(tset), columns=["gene", "go_id"])
NT_ALL = _tp.groupby("gene").size()

rows, t0 = [], time.time()
for name in args.datasets.split(","):
    ds = gr.Dataset(name, dset, tset)
    m = ds.mats(ds.dense)                       # every distinct score level >= the reporting floor
    C, YT, LT, thr = m["C"], m["YT"], m["LT"], m["thr"]
    # Recall denominator. Primary: the supported terms inside the predictor's own candidate space (the terms it
    # scores at or above its reporting floor), i.e. YT at the most liberal level -- the convention favourable to the
    # baseline. Sensitivity: the gene's complete TruePath annotation set, which drives the Fmax threshold to the floor.
    ntrue = YT[:, -1].astype(float)
    ntrue_full = NT_ALL.reindex(ds.units).fillna(0).to_numpy(float)
    has_true = ntrue > 0; has_true_full = ntrue_full > 0
    # ---- per-gene top-k sets (k = 1, 3): rank a gene's candidate calls by score, ties by the frozen hash
    r = ds.rows
    s = np.round(r.score.to_numpy(float), 6); gi = r.i.to_numpy(); T = r["T"].to_numpy(np.int32)
    h = gr.hashu(r.gene.to_numpy() + "|" + r.go_id.to_numpy(), salt="topk")
    order = np.lexsort((h, -s, gi))             # within gene: score desc, then hash
    rank_in_gene = np.empty(len(r), np.int64)
    starts = np.searchsorted(gi[order], np.arange(ds.n), side="left")
    rank_in_gene[order] = np.arange(len(r)) - starts[gi[order]]
    topk = {}
    for k in (1, 3):
        sel = rank_in_gene < k
        c = np.bincount(gi[sel], minlength=ds.n).astype(float)
        y = np.bincount(gi[sel], weights=T[sel], minlength=ds.n)
        topk[k] = (c, y, np.where(c > 0, 1 - y / np.maximum(c, 1), 0.0))
    for rep in range(args.repeats):
        perm = np.random.default_rng(gr.SEED0 + rep).permutation(ds.n)
        nc = round(0.70 * ds.n)
        cal, ev = perm[:nc], perm[nc:]          # the 70% the global calibrators use, and the held-out fold
        # ---- Fmax on the calibration genes (protein-centric precision/recall over the score grid)
        Ccal, Ycal = C[cal], YT[cal]
        with np.errstate(invalid="ignore", divide="ignore"):
            prec = np.where(Ccal > 0, Ycal / np.maximum(Ccal, 1), np.nan)
        npred = (Ccal > 0).sum(axis=0)
        pr = np.where(npred > 0, np.nansum(np.where(Ccal > 0, prec, 0.0), axis=0) / np.maximum(npred, 1), 0.0)
        cal_t = cal[has_true[cal]]
        rc = (YT[cal_t] / ntrue[cal_t][:, None]).mean(axis=0)
        F = np.where(pr + rc > 0, 2 * pr * rc / np.maximum(pr + rc, 1e-12), 0.0)
        kstar = int(np.argmax(F))
        cal_tf = cal[has_true_full[cal]]
        rcf = (YT[cal_tf] / ntrue_full[cal_tf][:, None]).mean(axis=0)
        Ff = np.where(pr + rcf > 0, 2 * pr * rcf / np.maximum(pr + rcf, 1e-12), 0.0)
        kfull = int(np.argmax(Ff))
        for tag, kk in (("Fmax", kstar), ("Fmax-full", kfull)):
            rows.append(dict(dataset=name, repeat=rep, rule=tag, label=f"score>={thr[kk]:g}",
                             fdp=float(LT[ev, kk].mean()), correct=float(YT[ev, kk].sum()), total=float(C[ev, kk].sum()),
                             units=int((C[ev, kk] > 0).sum()), supp_units=int((YT[ev, kk] > 0).sum()),
                             fmax=float((F if tag == "Fmax" else Ff)[kk]), prec=float(pr[kk]),
                             rec=float((rc if tag == "Fmax" else rcf)[kk])))
        rows.append(dict(dataset=name, repeat=rep, rule="Floor", label=f"score>={thr[-1]:g}",
                         fdp=float(LT[ev, -1].mean()), correct=float(YT[ev, -1].sum()), total=float(C[ev, -1].sum()),
                         units=int((C[ev, -1] > 0).sum()), supp_units=int((YT[ev, -1] > 0).sum()),
                         fmax=np.nan, prec=np.nan, rec=np.nan))
        for k in (1, 3):
            c, y, l = topk[k]
            rows.append(dict(dataset=name, repeat=rep, rule=f"Top-{k}", label=f"top{k}",
                             fdp=float(l[ev].mean()), correct=float(y[ev].sum()), total=float(c[ev].sum()),
                             units=int((c[ev] > 0).sum()), supp_units=int((y[ev] > 0).sum()),
                             fmax=np.nan, prec=np.nan, rec=np.nan))
        if rep % 25 == 0: print("  [%s] repeat %3d  %.0fs" % (name, rep, time.time() - t0), flush=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
pd.DataFrame(rows).to_csv(args.out, index=False)
print("wrote", args.out, "%.0fs" % (time.time() - t0), flush=True)
