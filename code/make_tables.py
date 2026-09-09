import os
"""Build LaTeX/CSV tables for the revised GoCo manuscript (uniform grid-path definition).
Reads out/*_full, out/Wainberg_grid, out/*_knap, out/*_knap2 split-level CSVs.
Writes ../../GoCo_BIB_revision_v2/tables_v2/*.tex|csv and analysis_rerun/all_methods_harmonised_100splits.csv
"""
import pandas as pd, numpy as np, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep  # output location only
RESULTS = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results')) + os.sep
OUTDIR = REV + "tables/"; os.makedirs(OUTDIR, exist_ok=True)
DS = ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']
N = {'Wainberg': 15794, 'Sanger': 14885, 'DRIVE': 6283, 'HAP1': 15020}
NT = {k: round(.1 * v) for k, v in N.items()}; NC70 = {k: round(.7 * v) for k, v in N.items()}
NE = {k: N[k] - NC70[k] for k in N}; M60 = {k: NC70[k] - NT[k] for k in N}

frames = []
SPLITDIR = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results'))  # frozen split-level outputs shipped with the repo
for f in glob.glob(os.path.join(SPLITDIR, '*_100splits.csv')):
    if 'repro' in f: continue
    x = pd.read_csv(f); x['src'] = os.path.basename(f); frames.append(x)
d = pd.concat(frames, ignore_index=True)
# the oracle control was redefined (ordering by the realised Delta_i, Lemma 1b); keep only the 'fix1' rerun rows for it
d = d[~(d.method.str.contains('oracle') & ~d.src.str.contains('fix1'))]
# learned orderings: the audited implementation is the 'learn2' run; drop the superseded 'learn' rows (and knapEB / ens, not reported)
d = d[~(d.src.str.contains('_learn_'))]
for k, v in {' (paper path)': '', ' (tie-level path)': '', ' (global, all distinct scores)': '', ' (global grid)': '', ' (Direct, global grid)': '',
             ' (global grid, certification fold only)': '-cert60', ' (grid path, all blocks)': ''}.items():
    d['method'] = d.method.str.replace(k, v, regex=False)
# Harmonise: uniform grid-path definition. Wainberg: GoCoGrid-* are the GoCo rows; the original one-block rows become 'GoCo1block-*'.
w = d.dataset == 'Wainberg'
d.loc[w & d.method.str.startswith('GoCo-') & ~d.method.str.startswith('GoCo-dense'), 'method'] = d.loc[w & d.method.str.startswith('GoCo-') & ~d.method.str.startswith('GoCo-dense'), 'method'].str.replace('GoCo-', 'GoCo1block-', regex=False)
d['method'] = d.method.str.replace('GoCoGrid-', 'GoCo-', regex=False)
d = d.drop_duplicates(['dataset', 'seed', 'alpha', 'delta', 'method'])
d['exceed'] = (d.unit_fdp > d.alpha).astype(float)
d['supp_frac'] = 1 - d.term_fdp
# pool-level risk of the selected policy (the certified estimand): weighted mean of certification-fold mean and held-out mean
def pool(row):
    ds = row.dataset
    if row.method.startswith('GoCo') or row.method == 'GoDag-cert60':
        mC = M60[ds]
    else:
        mC = NC70[ds]
    return (mC * row.cal_mean + NE[ds] * row.unit_fdp) / (mC + NE[ds])
d['pool_risk'] = d.apply(pool, axis=1)
d['pool_exceed'] = (d.pool_risk > d.alpha).astype(float)
# average the three random salts split-wise
r = d[d.method.isin(['GoCo-random1', 'GoCo-random2', 'GoCo-random3'])].groupby(['dataset', 'seed', 'alpha', 'delta'], as_index=False).mean(numeric_only=True)
r['method'] = 'GoCo-random'; r['policy'] = ''
d = pd.concat([d, r], ignore_index=True)
d.to_csv(RESULTS + 'all_methods_harmonised_100splits.csv', index=False)

LABEL = {'Boger': 'Boger et al.\\ (Direct loss)', 'GoDag': 'GoDag', 'GoDag-dense': 'GoDag-dense',
         'GoCo-random': 'GoCo-random', 'GoCo-source': 'source-smoothed $\\widehat\\Delta_i/\\sqrt{c_i}$', 'GoCo-knap': 'source-smoothed knapsack',
         'GoCo-score': 'GoCo-score (added-call score)', 'GoCo-exp0': 'GoCo ($\\widehat\\Delta_i$, $e=0$)', 'GoCo-exp1': 'GoCo ($\\widehat\\Delta_i/c_i$, $e=1$)',
         'GoCo-oracle': 'GoCo-oracle (realised $\\Delta_i$; uses labels)', 'GoCo-dense-knap': 'GoCo-dense, source-smoothed', 'GoCo-dense-learned': 'GoCo-dense (every distinct score)', 'GoCo1block-source': 'GoCo, single-block construction of the first draft',
         'GoCo-knapk0.5': 'GoCo, $\\kappa=0.5$', 'GoCo-knapk2': 'GoCo, $\\kappa=2$', 'GoCo-knapk4': 'GoCo, $\\kappa=4$',
         'GoCo-knapEB': 'GoCo, $\\kappa$ by method of moments', 'GoCo-learned': '\\textbf{GoCo} (Eq.~\\ref{eq:utility})', 'GoCo-learnedgb': 'GoCo, gradient-boosted variant', 'GoCo-ens': 'GoCo, rank-average variant'}

def write_tex(L, fname):
    """Join, make wide table* fit the text width, and make supplementary files self-contained."""
    txt = '\n'.join(L)
    if 'table*' in txt:
        txt = txt.replace('\\begin{tabular}', '\\resizebox{\\linewidth}{!}{\\begin{tabular}').replace('\\end{tabular}', '\\end{tabular}}')
    # threeparttable measures the unscaled tabular; replace it by a plain minipage for the notes
    txt = txt.replace('\\begin{threeparttable}\n', '').replace('\\end{threeparttable}\n', '')
    txt = txt.replace('\\begin{tablenotes}[flushleft]\\footnotesize\n\\item ', '\\par\\vspace{3pt}\\begin{minipage}{\\linewidth}\\footnotesize ').replace('\\end{tablenotes}', '\\end{minipage}')
    if os.path.basename(fname).startswith('TableS'):
        txt = txt.replace('Eq.~\\ref{eq:utility}', 'Eq.~(6)').replace('Proposition~\\ref{prop:validity}', 'Proposition~1').replace('[!t]', '[!htbp]')
    open(fname, 'w', encoding='utf-8').write(txt)

def paired(sub, m, ref, col='go_yield'):
    a = sub[sub.method == m].set_index('seed')[col]; b = sub[sub.method == ref].set_index('seed')[col]
    ix = a.index.intersection(b.index); dif = (a.loc[ix] - b.loc[ix])
    if len(dif) == 0: return np.nan, np.nan, np.nan
    return dif.mean(), 1.96 * dif.std(ddof=1) / np.sqrt(len(dif)), (dif > 0).mean()

def main_table(delta, fname, caption, label):
    methods = ['Boger', 'GoDag', 'GoDag-dense', 'GoCo-random', 'GoCo-learned']
    rows = []
    L = ['\\begin{table*}[!t]', '\\centering', f'\\caption{{{caption}}}', f'\\label{{{label}}}', '\\scriptsize', '\\setlength{\\tabcolsep}{3.2pt}', '\\begin{threeparttable}',
         '\\begin{tabular}{ll rrr rrrr r}', '\\toprule',
         'Data set & Method & Gene FDP & $\\Pr(\\widehat R_{\\mathcal E}>\\alpha)$ & $\\Pr(R>\\alpha)$ & Supp.\\ genes & Calls & Supp.\\ terms & Supp.\\ frac. & $\\Delta$terms vs GoDag [95\\% CI] \\\\', '\\midrule']
    for ds in DS:
        sub = d[(d.dataset == ds) & (d.alpha == 0.10) & (d.delta == delta)]
        first = True
        for m in methods:
            s = sub[sub.method == m]
            if len(s) == 0: continue
            md, ci, p = paired(sub, m, 'GoDag')
            dstr = '--' if m == 'GoDag' else f'{md:+.1f} [{md-ci:+.1f}, {md+ci:+.1f}]'
            lab = LABEL[m]
            L.append(f'{(ds + " ($n_{\\mathcal E}$=" + format(NE[ds], ",") + ")") if first else ""} & {lab} & {s.unit_fdp.mean():.4f} & {s.exceed.mean():.2f} & {s.pool_exceed.mean():.2f} & {s.unit_yield.mean():.1f} & {s.total_calls.mean():.0f} & {s.go_yield.mean():.1f} & {s.supp_frac.mean():.3f} & {dstr} \\\\')
            first = False
            rows.append(dict(dataset=ds, delta=delta, method=m, unit_fdp=s.unit_fdp.mean(), heldout_exceed=s.exceed.mean(), pool_exceed=s.pool_exceed.mean(), genes_with_calls=s.genes_with_calls.mean(), supported_genes=s.unit_yield.mean(), calls=s.total_calls.mean(), supported_terms=s.go_yield.mean(), supported_fraction=s.supp_frac.mean(), dterms_vs_godag=md, ci95=ci, p_gt=p))
        L.append('\\midrule')
    L[-1] = '\\bottomrule'
    L += ['\\end{tabular}', '\\begin{tablenotes}[flushleft]\\footnotesize',
          '\\item Means over the same 100 split streams at $\\alpha=0.10$. Gene FDP: held-out gene-level TruePath FDP $\\widehat R_{\\mathcal E}$. $\\Pr(\\widehat R_{\\mathcal E}>\\alpha)$: fraction of splits whose held-out FDP exceeded $\\alpha$. $\\Pr(R>\\alpha)$: fraction of splits in which the pool risk of the certified policy, the quantity bounded by Proposition~\\ref{prop:validity}, exceeded $\\alpha$. Supp.\\ genes / terms: held-out genes with $\\ge1$ TruePath-supported released term / supported released gene--GO pairs. Supp.\\ frac.: supported fraction of all released calls (the number of held-out genes with a nonempty set is in the reproducibility tables). Boger et al.\\ is the released implementation on the Direct loss, calibrated on 70\\% of the panel; all other rows use the TruePath loss. GoDag and GoCo share the 25-unit grid; GoDag-dense uses every distinct score value as a threshold; GoCo-random uses GoCo\'s blocks and admission grid with a hash ordering (mean of three hash seeds). The last column is the paired mean difference in supported terms from GoDag with a descriptive 95\\% interval over splits.',
          '\\end{tablenotes}', '\\end{threeparttable}', '\\end{table*}']
    write_tex(L, OUTDIR + fname + '.tex'); pd.DataFrame(rows).to_csv(OUTDIR + fname + '.csv', index=False)

main_table(0.50, 'Table1_primary_alpha010_delta050', 'Primary comparison at $\\alpha=0.10$ and $\\delta=0.50$ (the operating point of the released Boger et al.\\ implementation).', 'tab:primary50')
main_table(0.10, 'Table2_primary_alpha010_delta010', 'Primary comparison at $\\alpha=0.10$ and $\\delta=0.10$ (the substantive guarantee of Proposition~\\ref{prop:validity}).', 'tab:primary10')

def ablation_table(delta, fname, caption, label):
    methods = ['GoCo-random', 'GoCo-score', 'GoCo-exp0', 'GoCo-source', 'GoCo-exp1', 'GoCo-knap', 'GoCo-learnedgb', 'GoCo-learned', 'GoCo-oracle']
    L = ['\\begin{table*}[!t]', '\\centering', f'\\caption{{{caption}}}', f'\\label{{{label}}}', '\\scriptsize', '\\setlength{\\tabcolsep}{3.5pt}', '\\begin{threeparttable}',
         '\\begin{tabular}{l ' + 'rr ' * 4 + '}', '\\toprule', 'Ordering of affected genes & ' + ' & '.join(f'\\multicolumn{{2}}{{c}}{{{ds}}}' for ds in DS) + ' \\\\',
         ' & ' + ' & '.join(['Terms & vs random [P]'] * 4) + ' \\\\ \\midrule']
    rows = []
    for m in methods:
        cells = []
        for ds in DS:
            sub = d[(d.dataset == ds) & (d.alpha == 0.10) & (d.delta == delta)]
            s = sub[sub.method == m]
            if len(s) == 0: cells.append('-- & --'); continue
            md, ci, p = paired(sub, m, 'GoCo-random')
            cells.append(f'{s.go_yield.mean():.1f} & ' + ('--' if m == 'GoCo-random' else f'{md:+.1f} [{p:.2f}]'))
            rows.append(dict(dataset=ds, delta=delta, method=m, terms=s.go_yield.mean(), unit_fdp=s.unit_fdp.mean(), diff_vs_random=md, ci95=ci, p_gt_random=p))
        L.append(f'{LABEL[m]} & ' + ' & '.join(cells) + ' \\\\')
    L += ['\\bottomrule', '\\end{tabular}', '\\begin{tablenotes}[flushleft]\\footnotesize',
          '\\item Same blocks, admission grid, certification fold and statistic; only the ordering of affected genes differs. Terms: mean held-out supported GO-term yield. vs random: paired mean difference from GoCo-random and, in brackets, the fraction of splits in which the ordering beat GoCo-random. $\\widehat\\Delta_i$ is the source-smoothed predicted loss increment; $e$ is the exponent of the added-call count in $\\widehat\\Delta_i/c_i^{\\,e}$. GoCo estimates the increment and the number of supported added calls with a per-call model fitted to ranking-fold labels using frozen predictor features (Methods); the source-smoothed rows replace that model by shrinkage over module members alone. GoCo-oracle orders by the realised $\\Delta_i$ of every gene, using the labels of all folds; it is a reference for the ordering effect, not a method. Mean held-out gene FDP was at most $\\alpha$ for every row.',
          '\\end{tablenotes}', '\\end{threeparttable}', '\\end{table*}']
    write_tex(L, OUTDIR + fname + '.tex'); pd.DataFrame(rows).to_csv(OUTDIR + fname + '.csv', index=False)

ablation_table(0.50, 'Table3_ordering_ablation_delta050', 'Ordering ablation at $\\alpha=0.10$, $\\delta=0.50$: the value of the ordering beyond partial admission itself.', 'tab:ablation50')
ablation_table(0.10, 'TableS2_ordering_ablation_delta010', 'Ordering ablation at $\\alpha=0.10$, $\\delta=0.10$.', 'tab:ablation10')

def sweep_table(delta, fname, caption, label):
    methods = ['GoDag', 'GoDag-dense', 'GoCo-random', 'GoCo-learned', 'GoCo-dense-learned']
    L = ['\\begin{table*}[!t]', '\\centering', f'\\caption{{{caption}}}', f'\\label{{{label}}}', '\\scriptsize',
         '\\begin{tabular}{ll rr rr rr}', '\\toprule', 'Data set & Method & \\multicolumn{2}{c}{$\\alpha=0.05$} & \\multicolumn{2}{c}{$\\alpha=0.10$} & \\multicolumn{2}{c}{$\\alpha=0.20$} \\\\',
         ' & & FDP & Terms & FDP & Terms & FDP & Terms \\\\ \\midrule']
    rows = []
    for ds in DS:
        first = True
        for m in methods:
            cells = []
            for a in [0.05, 0.10, 0.20]:
                s = d[(d.dataset == ds) & (d.alpha == a) & (d.delta == delta) & (d.method == m)]
                if len(s) == 0: cells.append('-- & --'); continue
                cells.append(f'{s.unit_fdp.mean():.4f} & {s.go_yield.mean():.1f}')
                rows.append(dict(dataset=ds, delta=delta, alpha=a, method=m, unit_fdp=s.unit_fdp.mean(), exceed=s.exceed.mean(), pool_exceed=s.pool_exceed.mean(), terms=s.go_yield.mean(), genes=s.unit_yield.mean(), calls=s.total_calls.mean()))
            L.append(f'{ds if first else ""} & {LABEL[m]} & ' + ' & '.join(cells) + ' \\\\'); first = False
        L.append('\\midrule')
    L[-1] = '\\bottomrule'
    L += ['\\end{tabular}', '\\end{table*}']
    write_tex(L, OUTDIR + fname + '.tex'); pd.DataFrame(rows).to_csv(OUTDIR + fname + '.csv', index=False)

sweep_table(0.50, 'Table4_alpha_sweep_delta050', 'Error-target sweep at $\\delta=0.50$: mean held-out unit FDP and supported GO-term yield at $\\alpha\\in\\{0.05,0.10,0.20\\}$. GoCo-dense inserts partial admission between every pair of adjacent distinct score values.', 'tab:sweep50')
sweep_table(0.10, 'TableS3_alpha_sweep_delta010', 'Error-target sweep at $\\delta=0.10$.', 'tab:sweep10')

# ---------------- Supplementary: single-block construction of the first draft (Wainberg) and cert-60 GoDag ----------------
L = ['\\begin{table*}[!t]', '\\centering', '\\caption{Wainberg: the single-block construction used in the first draft (partial admission only of the tied step $800\\to789.65$, path ending there) versus the uniform grid-path construction of this paper, at three targets, and GoDag calibrated on the certification fold only. The first-draft construction cannot move beyond its single block, which is why it coincides with GoDag at $\\alpha=0.05$ and loses to GoDag at $\\alpha=0.20$.}', '\\label{tab:oneblock}', '\\scriptsize',
     '\\begin{tabular}{l l l rrr}', '\\toprule', '$\\alpha$ & $\\delta$ & Method & Gene FDP & Supp.\\ genes & Supp.\\ terms \\\\ \\midrule']
for alpha, delta in [(0.05, 0.50), (0.10, 0.50), (0.20, 0.50), (0.10, 0.10)]:
    for m, lab in [('GoDag', 'GoDag (70\\% calibration)'), ('GoDag-cert60', 'GoDag (60\\% certification fold only)'), ('GoCo1block-source', 'GoCo, single block, $\\widehat\\Delta_i/\\sqrt{c_i}$ (first draft)'), ('GoCo-source', 'GoCo, grid path, $\\widehat\\Delta_i/\\sqrt{c_i}$'), ('GoCo-knap', 'GoCo, grid path, Eq.~\\ref{eq:utility}')]:
        s = d[(d.dataset == 'Wainberg') & (d.alpha == alpha) & (d.delta == delta) & (d.method == m)]
        if len(s) == 0: continue
        L.append(f'{alpha:.2f} & {delta:.2f} & {lab} & {s.unit_fdp.mean():.4f} & {s.unit_yield.mean():.1f} & {s.go_yield.mean():.1f} \\\\')
    L.append('\\midrule')
L[-1] = '\\bottomrule'; L += ['\\end{tabular}', '\\end{table*}']
write_tex(L, OUTDIR + 'TableS4_wainberg_oneblock.tex')

# ---------------- Supplementary: kappa sensitivity ----------------
L = ['\\begin{table*}[!t]', '\\centering', '\\caption{Sensitivity of GoCo (Eq.~\\ref{eq:utility}) to the shrinkage weight $\\kappa$ of Equation~(5) at $\\alpha=0.10$: mean held-out gene FDP and supported GO-term yield over 100 splits. $\\kappa=1$ is the prespecified value used throughout.}', '\\label{tab:kappa}', '\\scriptsize',
     '\\begin{tabular}{l l rr rr rr rr}', '\\toprule', '$\\delta$ & $\\kappa$ & \\multicolumn{2}{c}{Wainberg} & \\multicolumn{2}{c}{Sanger} & \\multicolumn{2}{c}{DRIVE} & \\multicolumn{2}{c}{HAP1} \\\\', ' & & FDP & Terms & FDP & Terms & FDP & Terms & FDP & Terms \\\\ \\midrule']
rows = []
for delta in [0.50, 0.10]:
    for m, klab in [('GoCo-knapk0.5', '0.5'), ('GoCo-knap', '1'), ('GoCo-knapk2', '2'), ('GoCo-knapk4', '4')]:
        cells = []
        for ds in DS:
            s = d[(d.dataset == ds) & (d.alpha == 0.10) & (d.delta == delta) & (d.method == m)]
            cells.append(f'{s.unit_fdp.mean():.4f} & {s.go_yield.mean():.1f}' if len(s) else '-- & --')
            if len(s): rows.append(dict(delta=delta, kappa=klab, dataset=ds, unit_fdp=s.unit_fdp.mean(), terms=s.go_yield.mean()))
        L.append(f'{delta:.2f} & {klab} & ' + ' & '.join(cells) + ' \\\\')
    L.append('\\midrule')
L[-1] = '\\bottomrule'; L += ['\\end{tabular}', '\\end{table*}']
write_tex(L, OUTDIR + 'TableS8_kappa_sensitivity.tex'); pd.DataFrame(rows).to_csv(OUTDIR + 'TableS8_kappa_sensitivity.csv', index=False)

# ---------------- summary CSV of everything ----------------
s = d.groupby(['dataset', 'alpha', 'delta', 'method']).agg(unit_fdp=('unit_fdp', 'mean'), heldout_exceed=('exceed', 'mean'), pool_exceed=('pool_exceed', 'mean'), mean_pool_risk=('pool_risk', 'mean'),
                                                          genes_with_calls=('genes_with_calls', 'mean'), supported_genes=('unit_yield', 'mean'), calls=('total_calls', 'mean'), supported_terms=('go_yield', 'mean'), supported_fraction=('supp_frac', 'mean'), sd_terms=('go_yield', 'std')).round(4)
s.to_csv(RESULTS + 'summary_by_method_harmonised.csv')
print(s.loc[(slice(None), 0.1, 0.5), :][['unit_fdp', 'heldout_exceed', 'pool_exceed', 'supported_terms', 'calls']].to_string())
print('tables written to', OUTDIR)
