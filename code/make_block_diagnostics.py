import os
"""Supplementary diagnostics for the revised GoCo manuscript.
(1) Table S6: per data set x grid block used at alpha=0.10: affected pool/ranking genes, mean/sd Delta, fraction of affected
    genes with full source mass, out-of-sample Spearman correlations of Delta_hat, u (Eq. 6) and the added-call score with Delta.
(2) tie_counts.csv: exact-score tie statistics quoted in Results 3.1.
(3) Table S7: illustrative genes relaxed by GoCo (Eq. 6) but not by GoDag in one Wainberg split.
Run from the rerun directory (needs goco_rerun.py and the frozen inputs)."""
import numpy as np, pandas as pd, os, sys
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goco_rerun as gr
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep  # output location only
OUT = REV + 'tables/'; AN = os.path.join(HERE,'..','results') + os.sep
dset, tset = gr.load_truth()
ALPHA = 0.10
BLOCKS = {'Wainberg': [(800, 775)], 'Sanger': [(450, 425), (425, 400), (400, 375)], 'DRIVE': [(275, 250), (250, 225)], 'HAP1': [(475, 450)]}
rows = []; ties = []
for name in ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']:
    ds = gr.Dataset(name, dset, tset)
    mg = ds.mats(ds.grid); L = mg['LT']; C = mg['C']; YT = mg['YT']
    sc = np.round(ds.rows.score.to_numpy(float), 6)
    floor = 600 if name == 'Wainberg' else 200
    vals, cnts = np.unique(sc[sc >= floor], return_counts=True)
    ties.append(dict(dataset=name, score_floor=floor, rows_at_or_above_floor=int(cnts.sum()), distinct_scores=len(vals), largest_tie_rows=int(cnts.max()), largest_tie_score=float(vals[np.argmax(cnts)]),
                     rows_at_789_65=int(cnts[np.isclose(vals, 789.65)].sum()) if name == 'Wainberg' else np.nan))
    for hi, lo in BLOCKS[name]:
        j = int(np.where(np.isclose(ds.grid, hi))[0][0]); assert np.isclose(ds.grid[j + 1], lo)
        Aadd, meanscore = ds.Aadd(lo, hi)
        delta = (L[:, j + 1] - L[:, j]).astype(float); cadd = (C[:, j + 1] - C[:, j]).astype(float)
        affected = cadd > 0
        mass = np.asarray(Aadd.sum(axis=1)).ravel()
        gr.utility.extra = {'fsup': (YT[:, j + 1] - YT[:, j]).astype(float), 'lo': float(lo), 'hi': float(hi)}
        sp_d, sp_u, sp_s, sp_l, nT, pos = [], [], [], [], [], []
        for seed in range(100):
            perm = np.random.default_rng(gr.SEED0 + seed).permutation(ds.n); nt = round(.1 * ds.n); nc = round(.7 * ds.n)
            train = perm[:nt]; pool = perm[nt:]
            fit = train[affected[train]]; nT.append(len(fit))
            if len(fit) < 5: continue
            u, hh = gr.utility('knap', ds, Aadd, meanscore, delta, cadd, affected, train, ds.hh)
            ul, _ = gr.utility('learned', ds, Aadd, meanscore, delta, cadd, affected, train, ds.hh)
            dhat = gr.source_pred(Aadd, delta, fit, kappa=1.0)
            ap = pool[affected[pool]]
            sp_d.append(spearmanr(dhat[ap], delta[ap]).correlation); r = spearmanr(u[ap], delta[ap]).correlation; sp_u.append(r); pos.append(r > 0)
            sp_l.append(spearmanr(ul[ap], delta[ap]).correlation)
            ms = meanscore[ap]; ok = np.isfinite(ms); sp_s.append(spearmanr(-ms[ok], delta[ap][ok]).correlation)
        rows.append(dict(dataset=name, block=f'{hi:g}->{lo:g}', affected_genes=int(affected.sum()), affected_ranking_genes_mean=float(np.mean(nT)), mean_delta=float(delta[affected].mean()), sd_delta=float(delta[affected].std()),
                         frac_affected_full_source_mass=float(np.mean(np.isclose(mass[affected], 1.0, atol=1e-9))), frac_affected_zero_source_mass=float(np.mean(mass[affected] <= 1e-15)),
                         spearman_dhat_delta=float(np.mean(sp_d)), spearman_u_delta=float(np.mean(sp_u)), spearman_ulearned_delta=float(np.mean(sp_l)), spearman_negscore_delta=float(np.mean(sp_s)), frac_splits_u_positive=float(np.mean(pos))))
        print(rows[-1], flush=True)
    if name == 'Wainberg':
        # ---- Table S7: illustrative genes from split 0 at the GoCo-selected admission fraction ----
        seed = 0; perm = np.random.default_rng(gr.SEED0 + seed).permutation(ds.n); nt = round(.1 * ds.n); nc = round(.7 * ds.n)
        train, cert, ev = perm[:nt], perm[nt:nc], perm[nc:]
        j = int(np.where(np.isclose(ds.grid, 800))[0][0]); Aadd, meanscore = ds.Aadd(775, 800)
        delta = (L[:, j + 1] - L[:, j]).astype(float); cadd = (C[:, j + 1] - C[:, j]).astype(float); affected = cadd > 0
        gr.utility.extra = {'fsup': (YT[:, j + 1] - YT[:, j]).astype(float), 'lo': 775.0, 'hi': 800.0}
        u, hh = gr.utility('learned', ds, Aadd, meanscore, delta, cadd, affected, train, ds.hh)
        # replay the fixed-sequence selection for GoCo on this split to get the selected q
        sel_q = None
        for q in gr.QGRID_E:
            sel = gr.selmask(u, hh, affected, train, float(q)); y = np.where(sel, L[:, j + 1], L[:, j]).astype(float)
            p, r, sd = gr.clt_p(y[cert], ALPHA)
            if p > 0.5: break
            sel_q = float(q); sel_final = sel.copy()
        evsel = np.where(sel_final & np.isin(np.arange(ds.n), ev))[0]
        order = evsel[np.argsort(u[evsel])][:6]
        blk = ds.block_rows(775, 800)
        names = pd.read_csv(gr.W + 'wainberg_gene_go_scores.csv.gz', usecols=['gene', 'go_id', 'name'], dtype=str).drop_duplicates(['gene', 'go_id']).set_index(['gene', 'go_id'])['name']
        attr = pd.read_csv(gr.W + 'wainberg_highscore_source_attribution.csv.gz', usecols=['gene', 'go_id', 'winning_module', 'support_sources'], dtype={'gene': str, 'go_id': str}).drop_duplicates(['gene', 'go_id']).set_index(['gene', 'go_id'])
        lines = ['\\begin{table*}[!htbp]', '\\centering', f'\\caption{{Illustrative Wainberg genes released by GoCo and not by GoDag in split 0 (admission fraction $q={sel_q:g}$): the six evaluation-fold genes with the smallest ranking statistic $u_i$ among those admitted. Added terms are the gene--GO calls with scores in the step $800\\to775$; $\\checkmark$ marks TruePath support in the frozen reference. ``Module members supporting the term'' are the sources $\\mathcal R_{{ig}}$ that generated the call.}}', '\\label{tab:casestudy}', '\\scriptsize',
                 '\\resizebox{\\linewidth}{!}{\\begin{tabular}{l r r l p{0.42\\textwidth} l}', '\\toprule', 'Gene & $u_i$ & $c_i$ & Module & Added GO term & Module members supporting the term \\\\ \\midrule']
        case_rows = []
        for i in order:
            g = ds.units[i]; sub = blk[blk.i == i]
            first = True
            for _, r in sub.iterrows():
                nm = names.get((r.gene, r.go_id), r.go_id)
                mod = attr.loc[(r.gene, r.go_id), 'winning_module'] if (r.gene, r.go_id) in attr.index else ''
                src = attr.loc[(r.gene, r.go_id), 'support_sources'] if (r.gene, r.go_id) in attr.index else ''
                tick = ' $\\checkmark$' if r['T'] else ''
                lines.append(f"{g if first else ''} & {f'{u[i]:.3f}' if first else ''} & {int(cadd[i]) if first else ''} & {mod if first else ''} & {r.go_id} {nm}{tick} & \\texttt{{{str(src).replace(';', ', ')}}} \\\\")
                case_rows.append(dict(gene=g, u=float(u[i]), added_calls=int(cadd[i]), module=mod, go_id=r.go_id, term=nm, supported=int(r['T']), sources=src))
                first = False
            lines.append('\\addlinespace')
        lines += ['\\bottomrule', '\\end{tabular}}', '\\end{table*}']
        open(OUT + 'TableS7_case_study.tex', 'w', encoding='utf-8').write('\n'.join(lines)); pd.DataFrame(case_rows).to_csv(AN + 'case_study_split0.csv', index=False)
        print('case study written; q =', sel_q)
df = pd.DataFrame(rows); df.to_csv(AN + 'block_diagnostics_alpha010.csv', index=False)
pd.DataFrame(ties).to_csv(AN + 'tie_counts.csv', index=False)
L = ['\\begin{table*}[!htbp]', '\\centering', '\\caption{Block diagnostics at $\\alpha=0.10$ for the grid steps in which GoCo policies were certified. Affected genes: genes whose released set changes in the step; ranking genes: affected genes in the ranking fold (mean over 100 splits); $\\Delta$: realised loss increment over affected genes; full mass: fraction of affected genes with $\\sum_r a_{ir}=1$; Spearman columns: out-of-sample rank correlation over affected pool genes between the realised $\\Delta_i$ and, respectively, the source-smoothed $\\widehat\\Delta_i$, the ranking statistic $u_i$ of Equation~(6), the learned per-call statistic $u_i^{\\mathrm{L}}$ and the negative mean added-call score, averaged over splits; last column: fraction of splits with a positive correlation for $u_i$.}', '\\label{tab:blockdiag}', '\\scriptsize',
     '\\resizebox{\\linewidth}{!}{\\begin{tabular}{l l r r r r r r r r r r}', '\\toprule', 'Data set & Step & Affected & Ranking & mean $\\Delta$ & sd $\\Delta$ & full mass & $\\rho(\\widehat\\Delta,\\Delta)$ & $\\rho(u,\\Delta)$ & $\\rho(u^{\\mathrm L},\\Delta)$ & $\\rho(-\\text{score},\\Delta)$ & $\\Pr(\\rho_u>0)$ \\\\ \\midrule']
for _, r in df.iterrows():
    L.append(f"{r.dataset} & ${r.block.replace('->', '\\to ')}$ & {r.affected_genes} & {r.affected_ranking_genes_mean:.0f} & {r.mean_delta:.3f} & {r.sd_delta:.3f} & {r.frac_affected_full_source_mass:.2f} & {r.spearman_dhat_delta:.3f} & {r.spearman_u_delta:.3f} & {r.spearman_ulearned_delta:.3f} & {r.spearman_negscore_delta:.3f} & {r.frac_splits_u_positive:.2f} \\\\")
L += ['\\bottomrule', '\\end{tabular}}', '\\end{table*}']
open(OUT + 'TableS6_block_diagnostics.tex', 'w', encoding='utf-8').write('\n'.join(L))
print(df.round(3).to_string()); print(pd.DataFrame(ties).to_string())
