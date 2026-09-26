"""Build every number-bearing table of the paper (Tables 3-5, S1, S3, S8-S15 and key_numbers.json) from the split-level runs.

Inputs (all 100 splits, byte-identical split streams):
  repo/results/all_methods_harmonised_100splits.csv   Multilabel ('Boger'), GoDag, GoCo-M ('GoCo') on the four dependency data sets
  wf/wf_gocoN_detail.csv                              GoCo-N on the four dependency data sets (both delta, three alpha)
  WF_SEL (env; default wf/wf_selection_k10.csv)       ranking-fold selection rule, module family
  sf/funmap_v4_100.csv, sf/string_v4_100.csv          second family, five arms, both delta
  SF_SEL_FUNMAP / SF_SEL_STRING (env; default the embedded *_selection.csv)   ranking-fold selection rule, second family
Outputs: tables/*.tex and tables/key_numbers.json
"""
import json, os, sys
import numpy as np, pandas as pd

D = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(D)                                   # repository root
RES = os.path.join(ROOT, "results")                         # module family (goco_rerun.py -> make_tables.py)
WF = os.path.join(RES, "wf")                                # GoCo-N, ten-fold selection, Near-d, current-practice rules
SF = os.path.join(RES, "sf")                                # neighbourhood family
OUT = os.path.join(ROOT, "paper", "tables"); os.makedirs(OUT, exist_ok=True)
WF_SEL = os.environ.get("WF_SEL", f"{WF}/wf_selection_k10.csv")
def _sel(name):  # prefer the ten-fold rehearsal, matching the module family; fall back to the two-half one
    k10 = f"{SF}/{name}_v4_100_selection_k10.csv"
    return k10 if os.path.exists(k10) else f"{SF}/{name}_v4_100_selection.csv"
SF_SEL = {"FunMap": os.environ.get("SF_SEL_FUNMAP", _sel("funmap")),
          "STRING": os.environ.get("SF_SEL_STRING", _sel("string"))}
DS = ["Wainberg", "Sanger", "DRIVE", "HAP1"]
NE = {"Wainberg": 4738, "Sanger": 4465, "DRIVE": 1885, "HAP1": 4506}
ALPHAS = [0.05, 0.10, 0.20]
K = {}

# ------------------------------------------------------------------ load
har = pd.read_csv(f"{RES}/all_methods_harmonised_100splits.csv")
har = har.rename(columns={"seed": "repeat"})
wn = pd.read_csv(f"{WF}/wf_gocoN_detail.csv")
wsel = pd.read_csv(WF_SEL)
fm = pd.read_csv(f"{SF}/funmap_v4_100.csv"); fm["data"] = "FunMap"
st = pd.read_csv(f"{SF}/string_v4_100.csv"); st["data"] = "STRING"
sf = pd.concat([fm, st], ignore_index=True)
sfsel = pd.concat([pd.read_csv(SF_SEL[k]).assign(data=k) for k in SF_SEL], ignore_index=True)
NREP_SF = int(sf.groupby("data").repeat.nunique().min()); NREP_W = int(wn.groupby("dataset").repeat.nunique().min())
print("splits: second family", sf.groupby("data").repeat.nunique().to_dict(), "| GoCo-N module family", wn.groupby("dataset").repeat.nunique().to_dict(),
      "| selection rows", len(wsel), len(sfsel), "| WF_SEL =", os.path.basename(WF_SEL))
if os.environ.get("ALLOW_PARTIAL") != "1": assert NREP_W == 100, "GoCo-N run incomplete"

# ------------------------------------------------------------------ one long frame for the module family: arm x (dataset, repeat, alpha, delta)
def module_frame():
    rows = []
    for m, arm in [("Boger", "Multilabel"), ("GoDag", "GoDag"), ("GoCo", "GoCo-M")]:
        x = har[har.method == m]
        # pool_risk: TruePath risk over the certificate's own pool (calibration fold + evaluation fold), recorded by goco_rerun.metrics;
        # pool_risk_direct: the same on the Direct loss, the loss the Multilabel certificate bounds
        rows.append(pd.DataFrame(dict(dataset=x.dataset, repeat=x.repeat, alpha=x.alpha, delta=x.delta, arm=arm, fdp=x.unit_fdp, pool_risk=x.pool_risk,
                                      pool_risk_direct=x.pool_risk_direct,
                                      supp_genes=x.unit_yield, supp_terms=x.go_yield, calls=x.total_calls, genes_with_calls=x.genes_with_calls)))
    x = wn
    rows.append(pd.DataFrame(dict(dataset=x.dataset, repeat=x.repeat, alpha=x.alpha, delta=x.delta, arm="GoCo-N", fdp=x.fdp, pool_risk=x.pool_risk,
                                  pool_risk_direct=np.nan,
                                  supp_genes=x.supp_units, supp_terms=x.correct, calls=x.total, genes_with_calls=x.units)))
    f = pd.concat(rows, ignore_index=True)
    f["alpha"] = f.alpha.round(2); f["delta"] = f.delta.round(2)
    # ranking-fold-selected arm
    s = wsel[["dataset", "repeat", "alpha", "choose_N"]].copy(); s["alpha"] = s.alpha.round(2)
    M = f[f.arm == "GoCo-M"].merge(s, on=["dataset", "repeat", "alpha"]); N = f[f.arm == "GoCo-N"].merge(s, on=["dataset", "repeat", "alpha"])
    if os.environ.get("ALLOW_PARTIAL") != "1": assert len(M) == len(N) == len(f[f.arm == "GoCo-M"]), "selection rows missing"
    sel = pd.concat([M[M.choose_N == 0], N[N.choose_N == 1]], ignore_index=True).drop(columns="choose_N"); sel["arm"] = "GoCo-sel"
    return pd.concat([f, sel], ignore_index=True)
mf = module_frame()
mf.to_csv(f"{OUT}/module_family_long.csv", index=False)
ARM_LABEL = {"Multilabel": "Multilabel", "GoDag": "GoDag", "GoCo-M": r"\textbf{GoCo-M}", "GoCo-N": r"\textbf{GoCo-N}", "GoCo-sel": r"\textbf{GoCo}, $\mathcal T$-selected"}
ARM_LABEL_SWEEP = dict(ARM_LABEL); ARM_LABEL_SWEEP["Multilabel"] = "Multilabel (Direct loss)"

def cell(ds, alpha, delta, arm):
    return mf[(mf.dataset == ds) & (mf.alpha == alpha) & (mf.delta == delta) & (mf.arm == arm)].set_index("repeat").sort_index()
def ci(x):
    x = np.asarray(x, float); return x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))
def fmt_d(x, nd=1, plus=True):
    """Paired mean difference over splits. No interval: every quantity here is evaluated against
    the frozen reference on the held-out fold, so it is a realised difference, not an estimate."""
    m = float(np.mean(np.asarray(x, float))); f = ("{:+.%df}" if plus else "{:.%df}") % nd
    return f.format(m).replace("-", r"$-$")
def pct(a, b): return 100 * (a / b - 1)

# ------------------------------------------------------------------ Table 1: primary, alpha = 0.10, both delta
L = [r"\begin{table*}[!t]", r"\centering",
     r"\caption{Primary comparison at $\alpha=0.10$ on the module-structure family: means over the same 100 split streams at $\delta=0.50$ (the default of the released Multilabel procedure) and $\delta=0.10$ (the substantive guarantee of Proposition~\ref{prop:validity}). Wainberg is the development data set; Sanger, DRIVE and HAP1 are independent dependency resources on which the instantiation fixed on Wainberg is applied unchanged. Gene FDP: held-out gene-level TruePath FDP $\widehat R_{\mathcal E}$. Supp.\ genes: held-out genes with at least one TruePath-supported released term. Supp.\ terms: supported released gene--GO pairs. $\Delta$: paired mean difference in supported terms from GoDag over the same splits. Multilabel is calibrated on the Direct loss; GoDag and GoCo-M use the TruePath loss and share the 25-unit grid. The full panel of generic Direct-loss calibrators is in Table~\ref{tab:panelmodule}; exceedance and pool-level failure frequencies and released-call counts in Supplementary Table~S1.}",
     r"\label{tab:primary}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}", r"\resizebox{\linewidth}{!}{\begin{tabular}{ll rrrr rrrr}", r"\toprule",
     r" & & \multicolumn{4}{c}{$\delta=0.50$} & \multicolumn{4}{c}{$\delta=0.10$} \\", r"\cmidrule(lr){3-6}\cmidrule(lr){7-10}",
     r"Data set & Method & Gene FDP & Supp.\ genes & Supp.\ terms & $\Delta$ & Gene FDP & Supp.\ genes & Supp.\ terms & $\Delta$ \\", r"\midrule"]
for ds in DS:
    for k, arm in enumerate(["Multilabel", "GoDag", "GoCo-M"]):
        parts = []
        for delta in [0.50, 0.10]:
            a = cell(ds, 0.10, delta, arm); g = cell(ds, 0.10, delta, "GoDag")
            d = a.supp_terms - g.supp_terms.loc[a.index]
            parts.append(f"{a.fdp.mean():.4f} & {a.supp_genes.mean():.1f} & {a.supp_terms.mean():.1f} & " + ("--" if arm == "GoDag" else fmt_d(d)))
            K[f"primary_{ds}_{delta}_{arm}"] = dict(fdp=a.fdp.mean(), supp_genes=a.supp_genes.mean(), supp_terms=a.supp_terms.mean(), d_godag=d.mean(),
                                                     pct_godag=pct(a.supp_terms.mean(), g.supp_terms.mean()), p_godag=(d > 0).mean(), exceed=(a.fdp > 0.10).mean(),
                                                     pool_exceed=(a.pool_risk > 0.10).mean(), calls=a.calls.mean(), min_d_godag=d.min())
        first = (r"%s ($n_{\mathcal E}$=%s)" % (ds, format(NE[ds], ","))) if k == 0 else ""
        L.append(f"{first} & {ARM_LABEL[arm]} & " + " & ".join(parts) + r" \\")
    if ds != "HAP1": L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}}", r"\end{table*}"]
#  (the three-arm primary table is superseded by the module-structure panel; key numbers are still recorded below)
for ds in DS:                      # record every arm, including those now shown only in the supplement
    for arm in ["Multilabel", "GoDag", "GoCo-M", "GoCo-N", "GoCo-sel"]:
        for delta in [0.50, 0.10]:
            a = cell(ds, 0.10, delta, arm); g = cell(ds, 0.10, delta, "GoDag")
            d = a.supp_terms - g.supp_terms.loc[a.index]
            K[f"primary_{ds}_{delta}_{arm}"] = dict(fdp=a.fdp.mean(), supp_genes=a.supp_genes.mean(), supp_terms=a.supp_terms.mean(), d_godag=d.mean(),
                                                    pct_godag=pct(a.supp_terms.mean(), g.supp_terms.mean()), p_godag=(d > 0).mean(), exceed=(a.fdp > 0.10).mean(),
                                                    pool_exceed=(a.pool_risk > 0.10).mean(), pool_exceed_direct=(a.pool_risk_direct > 0.10).mean() if arm == "Multilabel" else np.nan,
                                                    calls=a.calls.mean(), min_d_godag=d.min())

# ------------------------------------------------------------------ alpha-sweep key numbers (the table is built after the second family is loaded)
for ds in DS:
    for arm in ["Multilabel", "GoDag", "GoCo-M", "GoCo-N", "GoCo-sel"]:
        for delta in [0.50, 0.10]:
            for alpha in ALPHAS:
                a = cell(ds, alpha, delta, arm)
                K[f"sweep_{ds}_{alpha}_{delta}_{arm}"] = dict(fdp=a.fdp.mean(), supp_terms=a.supp_terms.mean(), exceed=(a.fdp > alpha).mean())

# ------------------------------------------------------------------ Table S1: detail
def _ml_direct(delta):
    return "/".join(f"{(cell(ds, 0.10, delta, 'Multilabel').pool_risk_direct > 0.10).mean():.2f}" for ds in DS)
ML_DIRECT_NOTE = (f"{_ml_direct(0.50)} (Wainberg/Sanger/DRIVE/HAP1) at $\\delta=0.50$ and {_ml_direct(0.10)} at $\\delta=0.10$; "
                  "the Direct loss is stricter than TruePath, so a policy certified on it is conservative on the TruePath loss every table reports.")
K["ml_direct_pool_exceed"] = {f"{ds}_{delta}": float((cell(ds, 0.10, delta, 'Multilabel').pool_risk_direct > 0.10).mean()) for ds in DS for delta in [0.50, 0.10]}
L = [r"\begin{table*}[htbp]", r"\centering",
     r"\caption{Split-level statistics for Multilabel, GoDag, GoCo-M, GoCo-N and the $\mathcal T$-selected GoCo at $\alpha=0.10$ and $\delta\in\{0.50,0.10\}$, 100 splits: held-out exceedance $\Pr(\widehat R_{\mathcal E}>\alpha)$ and pool-level failure $\Pr(R>\alpha)$ frequencies, released calls, correct calls and the genes receiving them, paired differences from GoDag and the fraction of splits ahead of GoDag.}",
     r"\label{tab:detail}", r"\scriptsize", r"\setlength{\tabcolsep}{3.2pt}", r"\resizebox{\linewidth}{!}{\begin{tabular}{lll rrr rrr r r}", r"\toprule",
     r"$\delta$ & Data set & Method & Gene FDP & $\Pr(\widehat R_{\mathcal E}>\alpha)$ & $\Pr(R>\alpha)$ & Genes w/ correct call & Calls & Correct calls & $\Delta$ correct calls vs GoDag & $\Pr(>\text{GoDag})$ \\", r"\midrule"]
for di, delta in enumerate([0.50, 0.10]):
    for j, ds in enumerate(DS):
        for k, arm in enumerate(["Multilabel", "GoDag", "GoCo-M", "GoCo-N", "GoCo-sel"]):
            a = cell(ds, 0.10, delta, arm); g = cell(ds, 0.10, delta, "GoDag"); d = a.supp_terms - g.supp_terms.loc[a.index]
            dd = "--" if arm == "GoDag" else fmt_d(d); pw = "--" if arm in ("GoDag", "Multilabel") else f"{(d > 0).mean():.2f}"
            L.append(f"{f'{delta:.2f}' if (j == 0 and k == 0) else ''} & {ds if k == 0 else ''} & {ARM_LABEL_SWEEP[arm]} & {a.fdp.mean():.4f} & {(a.fdp > 0.10).mean():.2f} & {(a.pool_risk > 0.10).mean():.2f} & {a.supp_genes.mean():.1f} & {a.calls.mean():.0f} & {a.supp_terms.mean():.1f} & {dd} & {pw} \\\\")
    if di == 0: L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}}",
      r"\par\vspace{3pt}\begin{minipage}{\linewidth}\footnotesize Means over the same 100 split streams at $\alpha=0.10$. Gene FDP: held-out gene-level TruePath FDP $\widehat R_{\mathcal E}$. $\Pr(\widehat R_{\mathcal E}>\alpha)$: fraction of splits whose held-out FDP exceeded $\alpha$. $\Pr(R>\alpha)$: fraction of splits in which the TruePath risk of the certified policy over the pool of its own certificate---its calibration fold (70\% for Multilabel and GoDag, the 60\% certification fold for GoCo) together with the evaluation fold---exceeded $\alpha$; for GoDag and GoCo this is the event Proposition~1 bounds. Multilabel certifies the Direct loss instead, and on that loss its pool failure frequencies were " + ML_DIRECT_NOTE + r" Correct calls: released gene--GO pairs the frozen reference confirms; genes w/ correct call: held-out genes receiving at least one of them; Calls: released calls of any kind. Multilabel is the released multilabel calibration of Boger et al.\ on the Direct loss, calibrated on 70\% of the panel; GoDag and GoCo-M use the TruePath loss and share the 25-unit grid; GoCo-N walks the top-$m$ calls of the ranking-fold model; GoCo, $\mathcal T$-selected is whichever of the two the ranking-fold rule chose in that split (Table~S10). $\Delta$ correct calls: paired mean difference in correct calls from GoDag; $\Pr(>\text{GoDag})$: fraction of splits in which the method released more correct calls than GoDag.\end{minipage}",
      r"\end{table*}"]
open(f"{OUT}/TableS1_detail.tex", "w").write("\n".join(L) + "\n")

# ------------------------------------------------------------------ Table S8: GoCo-N and the selected GoCo on the module family, all alpha, both delta
L = [r"\begin{table*}[htbp]", r"\centering",
     r"\caption{Both instantiations on the four dependency data sets, 100 splits byte-identical to the main analysis. GoCo-N uses the per-call model of Section~2.3 unchanged---the same features and cross-fitting as GoCo-M---so only the candidate path differs. Correct calls---reference-supported released gene--GO pairs---on the evaluation fold for GoDag, GoCo-M, GoCo-N and the $\mathcal T$-selected GoCo; the paired difference GoCo-N minus GoCo-M and the fraction of splits in which GoCo-N was ahead; and the held-out gene-level FDP of GoCo-N. The prespecified structural mapping---shared evidence to GoCo-M---is what the main text reports, and it is not optimal on every predictor: Wainberg and DRIVE favour GoCo-M, Sanger and HAP1 GoCo-N, and the ranking fold sees which before certification (Table~S10).}",
     r"\label{tab:modulefamily}", r"\scriptsize", r"\begin{tabular}{l r r r r r r r r r}", r"\toprule",
     r"Data set & $\alpha$ & $\delta$ & GoDag & GoCo-M & GoCo-N & GoCo, $\mathcal T$-sel. & N $-$ M & Splits N ahead & FDP, GoCo-N \\", r"\midrule"]
for ds in DS:
    for alpha in ALPHAS:
        for delta in [0.50, 0.10]:
            g, m, n, s = (cell(ds, alpha, delta, a) for a in ["GoDag", "GoCo-M", "GoCo-N", "GoCo-sel"]); d = n.supp_terms - m.supp_terms
            L.append(f"{ds if (alpha == 0.05 and delta == 0.50) else ''} & {alpha:.2f} & {delta:.2f} & {g.supp_terms.mean():.1f} & {m.supp_terms.mean():.1f} & {n.supp_terms.mean():.1f} & {s.supp_terms.mean():.1f} & {fmt_d(d)} & {(d > 0).mean():.2f} & {n.fdp.mean():.4f} \\\\")
            K[f"module_{ds}_{alpha}_{delta}"] = dict(godag=g.supp_terms.mean(), M=m.supp_terms.mean(), N=n.supp_terms.mean(), sel=s.supp_terms.mean(), dNM=d.mean(), pN=(d > 0).mean(),
                                                    fdpN=n.fdp.mean(), exceedN=(n.fdp > alpha).mean(), pct_N_godag=pct(n.supp_terms.mean(), g.supp_terms.mean()),
                                                    pct_M_godag=pct(m.supp_terms.mean(), g.supp_terms.mean()), pct_sel_godag=pct(s.supp_terms.mean(), g.supp_terms.mean()),
                                                    min_sel_minus_godag=(s.supp_terms - g.supp_terms.loc[s.index]).min(), p_sel_ge_godag=((s.supp_terms - g.supp_terms.loc[s.index]) >= 0).mean())
    if ds != "HAP1": L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
open(f"{OUT}/TableS8_wainberg_control.tex", "w").write("\n".join(L) + "\n")

# ------------------------------------------------------------------ Table 5 (delta=0.10) and S13 (delta=0.50): second family
ARMS_SF = [("global-coarse", "Global grid"), ("global-dense", "Global grid, refined"), ("goco-M", r"\textbf{GoCo-M}"), ("goco-N", r"\textbf{GoCo-N}"), ("goco-sel", r"\textbf{GoCo}, $\mathcal T$-selected")]
def sfsel_frame():
    s = sfsel[["data", "repeat", "alpha", "choose_N"]].copy(); s["alpha"] = s.alpha.round(2)
    x = sf.copy(); x["alpha"] = x.alpha.round(2)
    M = x[x.arm == "goco-M"].merge(s, on=["data", "repeat", "alpha"]); N = x[x.arm == "goco-N"].merge(s, on=["data", "repeat", "alpha"])
    if os.environ.get("ALLOW_PARTIAL") != "1": assert len(M) == len(N) == len(x[x.arm == "goco-M"])
    sel = pd.concat([M[M.choose_N == 0], N[N.choose_N == 1]]).drop(columns="choose_N"); sel["arm"] = "goco-sel"
    return pd.concat([x, sel], ignore_index=True)
sfx = sfsel_frame()
sfx.to_csv(f"{OUT}/second_family_long.csv", index=False)
def scell(data, alpha, delta, arm):
    return sfx[(sfx.data == data) & (sfx.alpha == alpha) & (sfx.delta.round(2) == delta) & (sfx.arm == arm)].set_index("repeat").sort_index()
def second_family_table(delta, fname, label, caption, star):
    L = [r"\begin{table*}[!t]" if star else r"\begin{table*}[htbp]", r"\centering", r"\caption{" + caption + "}", r"\label{" + label + "}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
         r"\resizebox{\linewidth}{!}{\begin{tabular}{l l l rrr rrr r r r}", r"\toprule",
         r"Data & $\alpha$ & Candidate path & FDP & $\Pr(\widehat R_{\mathcal E}>\alpha)$ & $\Pr(R>\alpha)$ & Correct calls & Calls & Genes & vs.\ grid (\%) & $\Delta$ vs.\ refined & $\Pr(>\text{refined})$ \\", r"\midrule"]
    for data in ["FunMap", "STRING"]:
        for ai, alpha in enumerate(ALPHAS):
            base = scell(data, alpha, delta, "global-coarse"); ref = scell(data, alpha, delta, "global-dense")
            for k, (arm, name) in enumerate(ARMS_SF):
                a = scell(data, alpha, delta, arm).loc[ref.index]
                vs = pct(a.correct.mean(), base.correct.mean()); d = a.correct - ref.correct
                dstr = "--" if arm == "global-dense" else fmt_d(d, 0); pwin = "--" if arm in ("global-dense", "global-coarse") else f"{(d > 0).mean():.2f}"
                L.append(f"{data if (ai == 0 and k == 0) else ''} & {alpha:.2f} & {name} & {a.fdp.mean():.4f} & {(a.fdp > alpha).mean():.2f} & {(a.pool_risk > alpha).mean():.2f} & {a.correct.mean():.0f} & {a.total.mean():.0f} & {a.units.mean():.0f} & {vs:+.1f} & {dstr} & {pwin} \\\\")
                K[f"{data}_{alpha}_{delta}_{arm}"] = dict(fdp=a.fdp.mean(), exceed=(a.fdp > alpha).mean(), pool=(a.pool_risk > alpha).mean(), supported=a.correct.mean(), calls=a.total.mean(),
                                                         genes=a.units.mean(), vs_grid=vs, d_ref=d.mean(), p_ref=(d > 0).mean())
            if not (data == "STRING" and alpha == 0.20): L.append(r"\addlinespace[2pt]" if alpha != 0.20 else r"\midrule")
    L += [r"\bottomrule", r"\end{tabular}}", r"\end{table*}"]
    open(f"{OUT}/{fname}", "w").write("\n".join(L) + "\n")
cap_sf = (r"Neighbourhood-structure family in full: FunMap, the development data set, and STRING v12.0 (experiments channel), an independent network frozen eighteen months earlier that carries no Gene Ontology evidence. Means over the same 100 split streams. Every arm is certified by the same fixed-sequence rule on the same 60\% certification fold, so differences reflect the candidate path alone. FDP: held-out mean gene-level TruePath FDP; $\Pr(\widehat R_{\mathcal E}>\alpha)$: held-out exceedance; $\Pr(R>\alpha)$: pool-level failure (the event Proposition~1 bounds). Correct calls: released gene--GO pairs the frozen reference confirms; Calls: released calls of any kind; Genes: held-out genes receiving at least one call. vs.\ grid: change in correct calls relative to the global grid; $\Delta$ vs.\ refined: paired difference from the refined global grid, the fair baseline for a tie-free score since refinement is free, and the fraction of splits in which the arm exceeded it. GoCo-N is the instantiation prespecified for this structure; GoCo-M and the $\mathcal T$-selected arm are shown for comparison. ")
second_family_table(0.10, "TableS13_second_family_delta010.tex", "tab:secondfamily10", cap_sf + r"Operating point $\delta=0.10$.", False)
second_family_table(0.50, "TableS13b_second_family_delta050.tex", "tab:secondfamily50", cap_sf + r"Operating point $\delta=0.50$, the empirical check $\widehat R_{\mathcal C}\le\alpha$.", False)

# ================================================================== generic Direct-loss calibrator panel
# Module family: the frozen panel shipped with the reproducibility package (alpha = 0.10, both delta).
# Neighbourhood family: generic_sf.py, same rules, same split streams, same certification and evaluation folds.
GEN_ROWS = [("Hoeffding/LTT", "Hoeffding/LTT"), ("Hoeffding-Bentkus", "Hoeffding--Bentkus"), ("IID-normal", "IID-normal"),
            ("Empirical Bernstein", "Empirical Bernstein"), ("Binary-incidence McDiarmid", "McDiarmid (binary incidence)"),
            ("Janson dependency-graph", "Janson dependency graph"), ("Network-HAC (b=1)", r"Network-HAC ($b=1$)"),
            ("Platt calibration", "Platt calibration"), ("Isotonic calibration", "Isotonic calibration")]
genm = pd.read_csv(f"{RES}/first_draft_primary_method_summary_long.csv")
genm["delta"] = genm.delta.round(2)
def gcell_mod(ds, delta, method):
    x = genm[(genm.dataset == ds) & (genm.delta == delta) & (genm.method == method)]
    assert len(x) == 1, (ds, delta, method, len(x))
    return float(x.unit_fdp.iloc[0]), float(x.unit_yield.iloc[0]), float(x.go_yield.iloc[0])

GEN_SF_GRID = os.environ.get("GEN_SF_GRID", "dense")
GEN_SF_ROWS = [("Hoeffding/LTT", "Hoeffding/LTT"), ("Hoeffding-Bentkus", "Hoeffding--Bentkus"), ("IID-normal", "IID-normal"),
               ("Empirical Bernstein", "Empirical Bernstein"), ("Platt calibration", "Platt calibration"),
               ("Isotonic calibration", "Isotonic calibration")]
gensf = None
_gf = {"FunMap": f"{SF}/funmap_generic_{GEN_SF_GRID}.csv", "STRING": f"{SF}/string_generic_{GEN_SF_GRID}.csv"}
if all(os.path.exists(v) for v in _gf.values()):
    gensf = pd.concat([pd.read_csv(v).assign(data=k) for k, v in _gf.items()], ignore_index=True)
    gensf["alpha"] = gensf.alpha.round(2); gensf["delta"] = gensf.delta.round(2)
def gcell_sf(data, alpha, delta, method):
    x = gensf[(gensf.data == data) & (gensf.alpha == alpha) & (gensf.delta == delta) & (gensf.method == method)]
    return float(x.fdp.mean()), float(x.supp_terms.mean()), float(x.genes.mean())
_gfc = {"FunMap": f"{SF}/funmap_generic_coarse.csv", "STRING": f"{SF}/string_generic_coarse.csv"}
gensfc = None
if all(os.path.exists(v) for v in _gfc.values()):
    gensfc = pd.concat([pd.read_csv(v).assign(data=k) for k, v in _gfc.items()], ignore_index=True)
    gensfc["alpha"] = gensfc.alpha.round(2); gensfc["delta"] = gensfc.delta.round(2)
def gcell_sf_coarse(data, alpha, delta, method):
    """Same rules on the shared coarse grid, for rows that must be grid-matched with GoDag."""
    x = gensfc[(gensfc.data == data) & (gensfc.alpha == alpha) & (gensfc.delta == delta) & (gensfc.method == method)]
    return float(x.fdp.mean()), float(x.supp_terms.mean()), float(x.genes.mean())

HEAD_D = r"\emph{Direct loss, one global threshold (ontology-blind)}"
HEAD_T = r"\emph{TruePath loss, one global threshold (GO hierarchy)}"

# ------------------------------------------------------------------ Table 3 / S11: module-structure panel
def panel_module(delta, fname, label, caption, star):
    L = [r"\begin{table*}[!t]" if star else r"\begin{table*}[htbp]", r"\centering", r"\caption{" + caption + "}", r"\label{" + label + "}",
         r"\scriptsize", r"\setlength{\tabcolsep}{3pt}", r"\resizebox{\linewidth}{!}{\begin{tabular}{l rrr rrr rrr rrr}", r"\toprule",
         r"Method & \multicolumn{3}{c}{Wainberg ($n_{\mathcal E}$=4,738)} & \multicolumn{3}{c}{Sanger ($n_{\mathcal E}$=4,465)} & \multicolumn{3}{c}{DRIVE ($n_{\mathcal E}$=1,885)} & \multicolumn{3}{c}{HAP1 ($n_{\mathcal E}$=4,506)} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}\cmidrule(lr){11-13}",
         r" & Gene FDP & Genes with a & Correct & Gene FDP & Genes with a & Correct & Gene FDP & Genes with a & Correct & Gene FDP & Genes with a & Correct \\",
         r" & & correct call & calls & & correct call & calls & & correct call & calls & & correct call & calls \\", r"\midrule",
         r"\multicolumn{13}{l}{" + HEAD_D + r"} \\"]
    a = lambda ds, arm: cell(ds, 0.10, delta, arm)
    L.append("Multilabel & " + " & ".join(f"{a(ds,'Multilabel').fdp.mean():.4f} & {a(ds,'Multilabel').supp_genes.mean():.1f} & {a(ds,'Multilabel').supp_terms.mean():.1f}" for ds in DS) + r" \\")
    for src, lab in GEN_ROWS:
        cells = []
        for ds in DS:
            f_, g_, t_ = gcell_mod(ds, delta, src); cells.append(f"{f_:.4f} & {g_:.1f} & {t_:.1f}")
            K[f"gen_{ds}_{delta}_{src}"] = dict(fdp=f_, supp_genes=g_, supp_terms=t_)
        L.append(f"{lab} & " + " & ".join(cells) + r" \\")
    L.append(r"\addlinespace[2pt]\multicolumn{13}{l}{" + HEAD_T + r"} \\")
    L.append("GoDag & " + " & ".join(f"{a(ds,'GoDag').fdp.mean():.4f} & {a(ds,'GoDag').supp_genes.mean():.1f} & {a(ds,'GoDag').supp_terms.mean():.1f}" for ds in DS) + r" \\")
    L.append(r"\addlinespace[2pt]\multicolumn{13}{l}{\emph{TruePath loss, relational partial admission (GO hierarchy and module-shared evidence)}} \\")
    L.append(r"\textbf{GoCo-M} & " + " & ".join(f"\\textbf{{{a(ds,'GoCo-M').fdp.mean():.4f}}} & \\textbf{{{a(ds,'GoCo-M').supp_genes.mean():.1f}}} & \\textbf{{{a(ds,'GoCo-M').supp_terms.mean():.1f}}}" for ds in DS) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}}",
          r"\par\vspace{3pt}\begin{minipage}{\linewidth}\footnotesize Generic baselines calibrate on Direct labels only and have no GO-DAG/TruePath input; GoDag calibrates with TruePath; GoCo-M calibrates with TruePath plus the predictor's own module membership. All rows are evaluated against the same held-out TruePath reference on the same evaluation fold, over the same 100 split streams. Every global-threshold method calibrates over the shared 25-unit grid on the first 70\% of each permutation; GoCo-M certifies on the 60\% certification fold after fitting its ordering on the 10\% ranking fold. Platt and isotonic are $\delta$-free plug-in calibrators without a risk guarantee, so their rows are the same at both $\delta$. At $\delta=0.50$ every normal-type rule coincides with Multilabel (Remark~1 of the main text); at $\delta=0.10$ they separate, only the IID-normal CLT still agreeing with it.\end{minipage}",
          r"\end{table*}"]
    open(f"{OUT}/{fname}", "w").write("\n".join(L) + "\n")
panel_module(0.10, "Table3_panel_module.tex", "tab:panelmodule",
             r"Module-structure family, all calibrators, at $\alpha=0.10$ and $\delta=0.10$ (the substantive guarantee of Proposition~\ref{prop:validity}). Wainberg is the development data set; Sanger, DRIVE and HAP1 are independent dependency resources. Gene FDP: mean held-out gene-level TruePath false discovery proportion. Correct calls: released gene--GO pairs that the frozen TruePath reference confirms---the quantity called reference-supported elsewhere in the paper; genes with a correct call: held-out genes receiving at least one of them. Paired differences from GoDag, exceedance and pool-level failure frequencies are in Supplementary Table~S1; the $\delta=0.50$ operating point in Supplementary Table~S11.", True)
panel_module(0.50, "TableS11_panel_module_delta050.tex", "tab:panelmodule50",
             r"Module-structure family, all calibrators, at $\alpha=0.10$ and $\delta=0.50$ (the default of the released Multilabel procedure, equivalent to the empirical check $\widehat R_{\mathcal C}\le\alpha$); columns as in the module-structure panel of the main text.", False)

# ------------------------------------------------------------------ Table 4 / S12: neighbourhood-structure panel
def panel_nbr(delta, fname, label, caption, star):
    L = [r"\begin{table*}[!t]" if star else r"\begin{table*}[htbp]", r"\centering", r"\caption{" + caption + "}", r"\label{" + label + "}",
         r"\scriptsize", r"\setlength{\tabcolsep}{5pt}", r"\begin{tabular}{l rrr rrr}", r"\toprule",
         r"Method & \multicolumn{3}{c}{FunMap ($n_{\mathcal E}$=3,095)} & \multicolumn{3}{c}{STRING ($n_{\mathcal E}$=2,417)} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
         r" & Gene FDP & Correct & Genes with & Gene FDP & Correct & Genes with \\",
         r" & & calls & a call & & calls & a call \\", r"\midrule",
         r"\multicolumn{7}{l}{\emph{Direct loss, one global threshold on the refined grid (ontology-blind)}} \\"]
    for src, lab in [("Multilabel", "Multilabel")] + GEN_SF_ROWS:
        cells = []
        for data in ["FunMap", "STRING"]:
            f_, t_, g_ = gcell_sf(data, 0.10, delta, src); cells.append(f"{f_:.4f} & {t_:.0f} & {g_:.0f}")
            K[f"gensf_{data}_{delta}_{src}"] = dict(fdp=f_, supported=t_, genes=g_)
        L.append(f"{lab} & " + " & ".join(cells) + r" \\")
    L.append(r"\addlinespace[2pt]\multicolumn{7}{l}{" + HEAD_T + r"} \\")
    for arm, lab in [("global-coarse", "GoDag"), ("global-dense", "GoDag, refined grid")]:
        L.append(f"{lab} & " + " & ".join(
            (lambda x: f"{x.fdp.mean():.4f} & {x.correct.mean():.0f} & {x.units.mean():.0f}")(scell(data, 0.10, delta, arm))
            for data in ["FunMap", "STRING"]) + r" \\")
    L.append(r"\addlinespace[2pt]\multicolumn{7}{l}{\emph{TruePath loss, relational ordering of calls (GO hierarchy and neighbourhood evidence)}} \\")
    L.append(r"\textbf{GoCo-N} & " + " & ".join(
        (lambda x: f"\\textbf{{{x.fdp.mean():.4f}}} & \\textbf{{{x.correct.mean():.0f}}} & \\textbf{{{x.units.mean():.0f}}}")(scell(data, 0.10, delta, "goco-N"))
        for data in ["FunMap", "STRING"]) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\par\vspace{3pt}\begin{minipage}{\linewidth}\footnotesize Generic baselines calibrate on Direct labels only and have no GO-DAG input; GoDag calibrates with TruePath; GoCo-N adds the predictor's own neighbourhood evidence, ordering individual calls by the ranking-fold model, which is the instantiation this evidence structure calls for. GoCo-M, the instantiation for evidence shared between genes, is run on these data too and is in Supplementary Tables~S13 and~S14. Every Direct-loss rule selects one global threshold on the refined score grid---the fair global baseline for a tie-free score, since refinement costs nothing---and GoDag is shown on both the shared coarse grid and the refined grid. All rules certify on the same 60\% fold and are evaluated against the same held-out TruePath reference, over the same 100 split streams. Correct calls: released gene--GO pairs the frozen reference confirms; genes with a call: held-out genes receiving at least one released call, correct or not. The three network-variance rules of the module-structure panel (McDiarmid, Janson, Network-HAC) are defined through a gene-level dependency graph built from module membership and have no counterpart here; the module-structure panel shows that they buy no accuracy and cost yield. Platt and isotonic are $\delta$-free plug-in calibrators without a risk guarantee.\end{minipage}",
          r"\end{table*}"]
    open(f"{OUT}/{fname}", "w").write("\n".join(L) + "\n")
if gensf is not None:
    panel_nbr(0.10, "Table5_panel_neighbourhood.tex", "tab:panelnbr",
              r"Neighbourhood-structure family, all applicable calibrators, at $\alpha=0.10$ and $\delta=0.10$. FunMap is the development data set; STRING v12.0 (experiments channel), frozen eighteen months earlier and carrying no Gene Ontology evidence, is the independent validation. Gene FDP: mean held-out gene-level TruePath false discovery proportion. Correct calls: released gene--GO pairs the frozen reference confirms. Exceedance, pool-level failure and the matched-risk frontier are in Supplementary Tables~S13 and~S9; the $\delta=0.50$ operating point in Supplementary Table~S12.", True)
    panel_nbr(0.50, "TableS12_panel_neighbourhood_delta050.tex", "tab:panelnbr50",
              r"Neighbourhood-structure family, all applicable calibrators, at $\alpha=0.10$ and $\delta=0.50$; columns as in the neighbourhood-structure panel of the main text.", False)
else:
    print("WARNING: generic second-family panel missing; Table 4 not written")

# ------------------------------------------------------------------ Table 5 / S3: error-target sweep, every predictor
SWEEP_SF = [("Multilabel", "Multilabel (Direct loss, shared grid)"), ("GoDag", "GoDag, shared grid"), ("GoDag-refined", "GoDag, refined grid"), ("GoCo-N", r"\textbf{GoCo-N}")]
def sweep(delta, fname, label, caption, star):
    L = [r"\begin{table*}[!t]" if star else r"\begin{table*}[htbp]", r"\centering", r"\caption{" + caption + "}", r"\label{" + label + "}", r"\scriptsize",
         r"\begin{tabular}{ll rr rr rr}", r"\toprule",
         r"Data set & Method & \multicolumn{2}{c}{$\alpha=0.05$} & \multicolumn{2}{c}{$\alpha=0.10$} & \multicolumn{2}{c}{$\alpha=0.20$} \\",
         r" & & FDP & Correct & FDP & Correct & FDP & Correct \\",
         r" & & & calls & & calls & & calls \\ \midrule",
         r"\multicolumn{8}{l}{\emph{Module structure: GoCo-M (development: Wainberg)}} \\"]
    for ds in DS:
        for k, arm in enumerate(["Multilabel", "GoDag", "GoCo-M"]):
            parts = [f"{cell(ds, al, delta, arm).fdp.mean():.4f} & {cell(ds, al, delta, arm).supp_terms.mean():.1f}" for al in ALPHAS]
            L.append(f"{ds if k == 0 else ''} & {ARM_LABEL_SWEEP[arm]} & " + " & ".join(parts) + r" \\")
        L.append(r"\addlinespace[2pt]" if ds != "HAP1" else r"\midrule")
    L.append(r"\multicolumn{8}{l}{\emph{Neighbourhood structure: GoCo-N (development: FunMap)}} \\")
    for di, data in enumerate(["FunMap", "STRING"]):
        for k, (src, lab) in enumerate(SWEEP_SF):
            parts = []
            for al in ALPHAS:
                if src == "Multilabel":
                    f_, t_, _ = gcell_sf_coarse(data, al, delta, "Multilabel")
                else:
                    arm = {"GoDag": "global-coarse", "GoDag-refined": "global-dense", "GoCo-N": "goco-N"}[src]
                    x = scell(data, al, delta, arm); f_, t_ = x.fdp.mean(), x.correct.mean()
                parts.append(f"{f_:.4f} & {t_:.1f}")
            L.append(f"{data if k == 0 else ''} & {lab} & " + " & ".join(parts) + r" \\")
        if di == 0: L.append(r"\addlinespace[2pt]")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\par\vspace{3pt}\begin{minipage}{\linewidth}\footnotesize FDP: mean held-out gene-level TruePath false discovery proportion. Correct calls: released gene--GO pairs the frozen reference confirms. Each structure is shown with the instantiation prespecified for it. Within each block every rule uses the same grid, so the rows are comparable: GoCo-M and GoDag share the 25-unit grid in the module block, and Multilabel and GoDag the shared coarse grid in the neighbourhood block, where the refined grid---free for a tie-free score---is shown as a separate GoDag row.\end{minipage}",
          r"\end{table*}"]
    open(f"{OUT}/{fname}", "w").write("\n".join(L) + "\n")
if gensf is not None:
    sweep(0.50, "Table4_alpha_sweep_delta050.tex", "tab:sweep50",
          r"Error-target sweep at $\delta=0.50$: mean held-out gene-level FDP and supported GO-term yield at $\alpha\in\{0.05,0.10,0.20\}$ for every predictor, 100 splits.", True)
    sweep(0.10, "TableS3_alpha_sweep_delta010.tex", "tab:sweep10",
          r"Error-target sweep at $\delta=0.10$; layout as in the error-target sweep of the main text, 100 splits.", False)

# ------------------------------------------------------------------ Table S9: matched-risk frontier
L = [r"\begin{table*}[htbp]", r"\centering",
     r"\caption{Neighbourhood-structure family: gain over the exhaustive matched-risk frontier of the global path, 100 splits. For each certified arm the frontier value is the correct-call yield of the best global threshold---every distinct score value of the frozen score file, and abstention---whose realised held-out risk does not exceed the risk the arm realised: the exact oracle over all global cut-offs at that risk, so the difference isolates path quality from merely spending more of the budget. Paired mean differences and the fraction of splits in which the arm exceeded its frontier.}",
     r"\label{tab:frontier}", r"\scriptsize", r"\begin{tabular}{l l l rr rr}", r"\toprule",
     r"Data & $\alpha$ & $\delta$ & \multicolumn{2}{c}{GoCo-M} & \multicolumn{2}{c}{GoCo-N} \\", r" & & & Correct calls vs.\ frontier & Splits ahead & Correct calls vs.\ frontier & Splits ahead \\", r"\midrule"]
for data in ["FunMap", "STRING"]:
    for alpha in ALPHAS:
        for delta in [0.10, 0.50]:
            cells = []
            for arm in ["goco-M", "goco-N"]:
                a = scell(data, alpha, delta, arm)
                fr = a.frontier_all if "frontier_all" in a.columns else a.frontier      # exhaustive distinct-score oracle (frontier_exhaustive_sf.py)
                d = a.correct - fr; dg = a.correct - a.frontier
                cells += [fmt_d(d, 0), f"{(d > 0).mean():.2f}"]
                K[f"frontier_{data}_{alpha}_{delta}_{arm}"] = dict(d=d.mean(), p=(d > 0).mean(), d_grid=dg.mean(), p_grid=(dg > 0).mean(), exhaustive="frontier_all" in a.columns)
            L.append(f"{data if (alpha == 0.05 and delta == 0.10) else ''} & {alpha:.2f} & {delta:.2f} & " + " & ".join(cells) + r" \\")
    if data == "FunMap": L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
open(f"{OUT}/TableS9b_oracle_frontier.tex", "w").write("\n".join(L) + "\n")

# ------------------------------------------------------------------ Table S10: selection rule, every predictor, every alpha
RULE_DESC = os.environ.get("RULE_DESC", "ten-fold")
L = [r"\begin{table*}[htbp]", r"\centering",
     r"\caption{The ranking-fold selection rule of Section~2.3 on every predictor in this paper, 100 splits, $\delta=0.10$. In each split both instantiations are walked on the ranking fold $\mathcal T$ alone, with the per-call model cross-fitted inside $\mathcal T$, each stopping at the last candidate whose $\mathcal T$-mean loss is at most $\alpha$; the rule chooses the instantiation with the larger $\mathcal T$ supported yield. Columns: mean $\mathcal T$ yield of each instantiation; fraction of splits in which the rule chose GoCo-N; mean held-out supported yield of GoCo-M, GoCo-N and the instantiation the rule chose, when certified on $\mathcal C$; the fraction of splits in which GoCo-N was ahead of GoCo-M on the evaluation fold; and the fraction of splits in which the rule's choice was the held-out winner (ties counted as agreement). The rule is a function of $\mathcal T$ and frozen inputs only, so choosing by it leaves Proposition~1 unchanged (Corollary~1).}",
     r"\label{tab:selection}", r"\scriptsize", r"\begin{tabular}{l l l rr r rrr r r}", r"\toprule",
     r"Predictor & Evidence & $\alpha$ & \multicolumn{2}{c}{$\mathcal T$ yield} & Rule chooses & \multicolumn{3}{c}{Held-out supported} & $\Pr(\text{N}>\text{M})$ & Rule agrees \\",
     r" & & & GoCo-M & GoCo-N & GoCo-N & GoCo-M & GoCo-N & selected & held-out & with winner \\", r"\midrule"]
def held(data, alpha):
    if data in ("FunMap", "STRING"):
        m, n, s = (scell(data, alpha, 0.10, a).correct for a in ["goco-M", "goco-N", "goco-sel"])
        r = sfsel[(sfsel.data == data) & (sfsel.alpha.round(2) == alpha)].set_index("repeat").sort_index()
    else:
        m, n, s = (cell(data, alpha, 0.10, a).supp_terms for a in ["GoCo-M", "GoCo-N", "GoCo-sel"])
        r = wsel[(wsel.dataset == data) & (wsel.alpha.round(2) == alpha)].set_index("repeat").sort_index()
    idx = m.index.intersection(n.index).intersection(r.index); m, n, s, r = m.loc[idx], n.loc[idx], s.loc[idx], r.loc[idx]
    winner_N = (n > m); agree = np.where(n == m, True, (r.choose_N == 1) == winner_N).mean()
    return r.yield_M_T.mean(), r.yield_N_T.mean(), r.choose_N.mean(), m.mean(), n.mean(), s.mean(), winner_N.mean(), agree, np.maximum(m, n).mean()
for data, ev in [("Wainberg", "module-shared"), ("Sanger", "module-shared"), ("DRIVE", "module-shared"), ("HAP1", "module-shared"), ("FunMap", "neighbourhood"), ("STRING", "neighbourhood")]:
    for ai, alpha in enumerate(ALPHAS):
        yM, yN, cN, hM, hN, hS, pN, agree, orc = held(data, alpha)
        L.append(f"{data if ai == 0 else ''} & {ev if ai == 0 else ''} & {alpha:.2f} & {yM:.1f} & {yN:.1f} & {cN:.2f} & {hM:.1f} & {hN:.1f} & {hS:.1f} & {pN:.2f} & {agree:.2f} \\\\")
        K[f"select_{data}_{alpha}"] = dict(yM=yM, yN=yN, chooseN=cN, heldM=hM, heldN=hN, heldSel=hS, pN=pN, agree=agree, oracle=orc)
    if data != "STRING": L.append(r"\midrule" if data == "HAP1" else r"\addlinespace[2pt]")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
open(f"{OUT}/TableS10_selection_rule.tex", "w").write("\n".join(L) + "\n")

# ------------------------------------------------------------------ Table S8 (structures): what a cutoff misses, and the rule's verdict where tested
def verdict(*names):
    out = []
    for n in names:
        c = K[f"select_{n}_0.1"]["chooseN"]
        out.append(f"{n}: {'GoCo-N' if c > 0.5 else 'GoCo-M'} ({max(c, 1 - c)*100:.0f}\\%)")
    return "; ".join(out)
L = [r"\begin{table*}[htbp]", r"\centering",
     r"\caption{Source-object structures, what a cutoff misses on each, and the ranking-fold verdict where tested. The framework asks a predictor for a score and the evidence behind each call. Whether that evidence is a group shared between genes, which forces exact score ties, or is assembled per gene, which does not, decides what a global cutoff cannot reach: a tied block has no interior, so a cutoff takes or drops it whole, whereas a tie-free score lets a refined grid recover the granularity term of Lemma~1 for free. Both instantiations are available for every structure; which yields more is decided on the ranking fold (Section~2.3), and the last column gives the rule's majority verdict at $\alpha=0.10$ with the fraction of splits in which it was chosen (Table~S10). Rows without a verdict are classified by construction, which says what a cutoff misses and not which instantiation wins.}",
     r"\label{tab:structures}", r"\tiny", r"\setlength{\tabcolsep}{3pt}", r"\renewcommand{\arraystretch}{1.15}", r"\begin{tabular}{@{}p{0.145\textwidth} p{0.135\textwidth} p{0.05\textwidth} p{0.085\textwidth} p{0.105\textwidth} p{0.135\textwidth} p{0.225\textwidth}@{}}", r"\toprule",
     r"Predictor family & Source set $\mathcal R_{ig}$ & Evidence & Score ties & Structure & A cutoff misses & Ranking-fold verdict \\", r"\midrule",
     r"Co-essential modules (Wainberg; Sanger and DRIVE via the same catalogue) & annotated members of the gene's module & shared & exact, up to 3,361 pairs & module-shared & the interior of every tied block & " + verdict("Wainberg", "Sanger", "DRIVE") + r" \\",
     r"Protein complexes (CORUM, BioPlex) & annotated members of the complex & shared & exact & module-shared & the interior of every tied block & structural \\",
     r"Co-regulation clusters (ProteomeHD) & annotated members of the co-regulated group & shared & exact & module-shared & the interior of every tied block & structural \\",
     r"Community / partition detection (HAP1 Global-L2 clusters) & annotated members of the community & shared & exact, up to 14,245 pairs & module-shared & the interior of every tied block & " + verdict("HAP1") + r" \\",
     r"Random-walk neighbourhoods (FunMap) & annotated neighbours, proximity-ranked & gene-specific & none material & neighbourhood & little: a refined grid recovers granularity & " + verdict("FunMap") + r" \\",
     r"Interaction neighbourhoods (STRING) & annotated interactors, ranked & gene-specific & none material & neighbourhood & little: a refined grid recovers granularity & " + verdict("STRING") + r" \\",
     r"Network propagation / diffusion & diffusion-weighted neighbours & gene-specific & none material & neighbourhood & little & structural \\",
     r"Multi-network integration (GeneMANIA) & (network, neighbour) pairs with weights & gene-specific & none material & neighbourhood & little & structural \\",
     r"Sequence-homology hit lists & homologues with alignment scores & gene-specific & none material & neighbourhood & little & structural \\",
     r"Embedding $k$-nearest neighbours & nearest neighbours with distances & gene-specific & none material & neighbourhood & little & structural \\",
     r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
open(f"{OUT}/TableS9_structures.tex", "w").write("\n".join(L) + "\n")

# ------------------------------------------------------------------ ranking-fold AUC diagnostics (score vs cross-fitted model, T calls)
dg = pd.read_csv(f"{WF}/wf_gocoN_detail_diag.csv")
K["auc"] = {d: dict(score_T=x.auc_score_T.mean(), phat_T=x.auc_phat_T.mean(), score_E=x.auc_score_E.mean(), phat_E=x.auc_phat_E.mean(), p_phat_better_T=(x.auc_phat_T > x.auc_score_T).mean())
            for d, x in dg.groupby("dataset")}
json.dump(K, open(f"{OUT}/key_numbers.json", "w"), indent=1, default=float)
print("tables written to", OUT)

# ================================================================== Table 3: release rules in current use
UNC = {}
for f, n in [(f"{WF}/wf_uncalibrated.csv", None), (f"{SF}/funmap_uncalibrated.csv", "FunMap"), (f"{SF}/string_uncalibrated.csv", "STRING")]:
    if not os.path.exists(f): continue
    x = pd.read_csv(f)
    if n is not None: x["dataset"] = n
    for (dd, rr), y in x.groupby(["dataset", "rule"]):
        UNC[(dd, rr)] = y
        K[f"unc_{dd}_{rr}"] = dict(fdp=y.fdp.mean(), ex10=(y.fdp > 0.10).mean(), genes=y.units.mean(), calls=y.total.mean(), supp=y.correct.mean())
PRED = DS + ["FunMap", "STRING"]
# fixed score cutoffs on the second family, read from the raw score file (below the analysis floor as well)
FIX = {}
for f, n in [(f"{SF}/funmap_fixed_cutoff.csv", "FunMap"), (f"{SF}/string_fixed_cutoff.csv", "STRING")]:
    if not os.path.exists(f): continue
    x = pd.read_csv(f)
    for c, y in x.groupby("cutoff"):
        FIX[(n, round(float(c), 5))] = y
        K[f"fixed_{n}_{c:g}"] = dict(fdp=y.fdp.mean(), ex10=(y.fdp > 0.10).mean(), genes=y.units.mean(), calls=y.total.mean(), supp=y.correct.mean())
RULES = [("Fmax-full", r"Validation $F_{\max}$ threshold"),
         ("Top-1", "Top-1 GO term per gene"), ("Top-3", "Top-3 GO terms per gene")]
# A fixed cutoff is 0.5 when the predictor emits a probability. None of these six scores is on a 0--1 scale,
# so the fixed row uses each predictor's own reporting threshold -- the value at or above which it publishes.
FIXED_THR = {"Wainberg": "600", "Sanger": "100", "DRIVE": "100", "HAP1": "100", "FunMap": "1", "STRING": "1"}
FIXED_KEY = 1.0
have = all((p, r) in UNC or r == "Fixed" for p in PRED for r, _ in RULES)
if have and ("FunMap", FIXED_KEY) in FIX and ("Wainberg", "Floor") in UNC:
    def fixed_rows(p):
        return FIX[(p, FIXED_KEY)] if p in ("FunMap", "STRING") else UNC[(p, "Floor")]
    L = [r"\begin{table*}[htbp]", r"\centering",
         r"\caption{Release rules in current use, at their own operating points, on the same 100 split streams and the same held-out fold as every other table; means over splits. \emph{Validation $F_{\max}$}: the score threshold maximising the $F_1$ measure, in the protein-centric form CAFA uses---precision averaged over the genes with at least one released call, recall over the genes with at least one reference annotation, $F=2pr/(p+r)$---computed on the same 70\% of the panel the global calibrators use and applied to the held-out fold. \emph{Top-$k$}: the $k$ highest-scoring candidate terms of each gene. Validation $F$ is small in absolute terms here because a gene's complete TruePath reference contains many terms the predictor never scores, so recall stays low and $F_1$ is maximised by taking essentially everything the predictor offers: the maximiser sits at or just above each predictor's own reporting threshold. It is nonetheless the quantity the rule optimises, and Gene FDP is the held-out gene-level TruePath FDP it realises; the two are unrelated: neither rule in the upper part of each block estimates, bounds or exposes an error rate, and the FDP they realise runs from " + "%.2f to %.2f" % (min(min(UNC[(q, r)].fdp.mean() for r, _ in RULES) for q in PRED), max(max(UNC[(q, r)].fdp.mean() for r, _ in RULES) for q in PRED)) + r" across these six predictors. The last row of each block is GoCo certified at $\alpha=0.10$, $\delta=0.10$, with the instantiation prespecified for that predictor's evidence structure.}",
         r"\label{tab:current}", r"\scriptsize", r"\setlength{\tabcolsep}{4pt}",
         r"\resizebox{\linewidth}{!}{\begin{tabular}{l l r r r r r r r}", r"\toprule",
         r"Predictor & Release rule & Threshold & Validation $F$ & Gene FDP & $\Pr(\widehat R_{\mathcal E}>0.10)$ & Genes with a call & Calls & Correct calls \\", r"\midrule"]
    for p in PRED:
        for k, (r, lab) in enumerate(RULES):
            y = UNC[(p, r)]
            if r.startswith("Fmax"):
                thr = f"{y.label.str.replace('score>=', '', regex=False).astype(float).mean():.3g}"; vf = f"{y.fmax.mean():.3f}"
            else:
                thr = "--"; vf = "--"
            L.append(f"{p if k == 0 else ''} & {lab} & {thr} & {vf} & {y.fdp.mean():.3f} & {(y.fdp > 0.10).mean():.2f} & {y.units.mean():.0f} & {y.total.mean():.0f} & {y.correct.mean():.0f} \\\\")
            K[f"unc_{p}_{r}"] = dict(fdp=y.fdp.mean(), ex10=(y.fdp > 0.10).mean(), genes=y.units.mean(), calls=y.total.mean(), supp=y.correct.mean())
        if p in DS:
            a = cell(p, 0.10, 0.10, "GoCo-M"); fdp, gg, cc, ss = a.fdp.mean(), a.genes_with_calls.mean(), a.calls.mean(), a.supp_terms.mean()
        else:
            a = scell(p, 0.10, 0.10, "goco-N"); fdp, gg, cc, ss = a.fdp.mean(), a.units.mean(), a.total.mean(), a.correct.mean()
        L.append(r"& \textbf{GoCo, certified at $\alpha=0.10$} & nested policy & -- & " + f"\\textbf{{{fdp:.3f}}} & \\textbf{{{(a.fdp > 0.10).mean():.2f}}} & {gg:.0f} & {cc:.0f} & \\textbf{{{ss:.0f}}} \\\\")
        K[f"unc_{p}_GoCo"] = dict(fdp=fdp, ex10=(a.fdp > 0.10).mean(), genes=gg, calls=cc, supp=ss)
        if p != "STRING": L.append(r"\midrule")
    L += [r"\bottomrule", r"\end{tabular}}", r"\end{table*}"]
    open(f"{OUT}/TableS15_current_practice.tex", "w").write("\n".join(L) + "\n")
    print("Table S15 (current practice) written")

# ================================================================== Table S14: Near-d loss sensitivity
ND = f"{WF}/wf_neard.csv"
if os.path.exists(ND):
    nd_all = pd.read_csv(ND)
    nd = nd_all[nd_all.arm.isin(["GoDag", "GoCo-M", "GoCo-N"])]
    nrep = int(nd.groupby(["dataset", "loss", "arm"]).repeat.nunique().min())
    LOSSN = {"truepath": r"$d=0$ (TruePath)", "near1": r"$d=1$", "near2": r"$d=2$"}
    L = [r"\begin{table*}[htbp]", r"\centering",
         r"\caption{Hierarchy-aware loss sensitivity: the whole procedure rerun under the Near-$d$ family of Note~S3.2 at $\alpha=0.10$, $\delta=0.10$, " + str(nrep) + r" splits. A released call is acceptable when some term within undirected distance $d$ of it in the frozen GO DAG is TruePath-annotated for the gene; $d=0$ is the loss used throughout the paper. Splits, frozen scores, candidate paths, model class and certificate are unchanged---only the definition of an acceptable call moves, and the per-call model is refitted to it. Acceptable calls: released calls the loss accepts; genes w/ acceptable call: held-out genes receiving at least one. $\Delta$: paired mean difference in accepted calls from GoDag. For comparability across $d$ the GoCo-M path here is the uniform construction of Box~1 step~4a on all four data sets; at $d=0$ it reproduces the main-text Wainberg policy to within " + ("%.1f" % abs(nd_all[(nd_all.dataset == "Wainberg") & (nd_all.loss == "truepath") & (nd_all.arm == "GoCo-M")].go_yield.mean() - nd_all[(nd_all.dataset == "Wainberg") & (nd_all.loss == "truepath") & (nd_all.arm == "GoCo-M (main-text path)")].go_yield.mean())) + r" supported terms in 100 splits. Validity is untouched: the held-out FDP stays at or below target in every cell, and which instantiation yields more is the same under all three losses, so it is not an artefact of the strict reference.}",
         r"\label{tab:neard}", r"\scriptsize", r"\begin{tabular}{l l l r r r r r}", r"\toprule",
         r"Data set & Loss & Method & Gene FDP & $\Pr(R>\alpha)$ & Genes w/ acceptable call & Acceptable calls & $\Delta$ vs.\ GoDag \\", r"\midrule"]
    for ds_ in [x for x in DS if (nd.dataset == x).any()]:
        for li, loss in enumerate(["truepath", "near1", "near2"]):
            g = nd[(nd.dataset == ds_) & (nd.loss == loss) & (nd.arm == "GoDag")].set_index("repeat").sort_index()
            for k, arm in enumerate(["GoDag", "GoCo-M", "GoCo-N"]):
                a = nd[(nd.dataset == ds_) & (nd.loss == loss) & (nd.arm == arm)].set_index("repeat").sort_index()
                d_ = a.go_yield - g.go_yield.loc[a.index]
                nm = arm if arm == "GoDag" else r"\textbf{" + arm + "}"
                L.append(f"{ds_ if (li == 0 and k == 0) else ''} & {LOSSN[loss] if k == 0 else ''} & {nm} & {a.unit_fdp.mean():.4f} & {(a.pool_risk > 0.10).mean():.2f} & {a.unit_yield.mean():.0f} & {a.go_yield.mean():.0f} & " + ("--" if arm == "GoDag" else fmt_d(d_, 0)) + r" \\")
                K[f"neard_{ds_}_{loss}_{arm}"] = dict(fdp=a.unit_fdp.mean(), pool=(a.pool_risk > 0.10).mean(), genes=a.unit_yield.mean(),
                                                      acc=a.go_yield.mean(), d_godag=d_.mean(), pct_godag=pct(a.go_yield.mean(), g.go_yield.mean()), p_godag=(d_ > 0).mean())
            if loss != "near2": L.append(r"\addlinespace[2pt]")
        if ds_ != [x for x in DS if (nd.dataset == x).any()][-1]: L.append(r"\midrule")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    open(f"{OUT}/TableS14_neard.tex", "w").write("\n".join(L) + "\n")
    print("Table S14 written")
json.dump(K, open(f"{OUT}/key_numbers.json", "w"), indent=1, default=float)
