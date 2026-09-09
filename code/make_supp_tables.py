"""Supplementary Table S5 (risk decomposition). Table S1 is produced by make_tables.py."""
import pandas as pd, os
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep
RESULTS = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results')) + os.sep
DS = ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']
def fit(L):
    return '\n'.join(L).replace('\\begin{tabular}', '\\resizebox{\\linewidth}{!}{\\begin{tabular}').replace('\\end{tabular}', '\\end{tabular}}')

d = pd.read_csv(RESULTS + 'all_methods_harmonised_100splits.csv')
L = ['\\begin{table*}[htbp]', '\\centering',
     '\\caption{Calibration mean ($\\widehat R_{\\mathcal C}$ over the 60\\% certification fold for GoCo; $\\widehat R_{\\mathcal T\\cup\\mathcal C}$ over the 70\\% calibration sample for GoDag), pool risk $R(\\widehat\\pi)$ and held-out mean $\\widehat R_{\\mathcal E}(\\widehat\\pi)$ of the certified policy, averaged over 100 splits at $\\alpha=0.10$, with the fractions of splits in which $R(\\widehat\\pi)>\\alpha$ (the event bounded by Proposition 1) and $\\widehat R_{\\mathcal E}(\\widehat\\pi)>\\alpha$; and the yield of the certified policy on the evaluation fold: the number of genes with a nonempty released set (of $n_{\\mathcal E}$ = 4,738 / 4,465 / 1,885 / 4,506), the number of released calls, the number of reference-supported gene--GO calls and the supported terms per gene carrying at least one. The pooled precision of the released call list is not the certified quantity and is not reported as an error rate (Note S3.7).}',
     '\\label{tab:riskdecomp}', '\\scriptsize', '\\begin{tabular}{l l l rrr rr rrrr}', '\\toprule',
     '$\\delta$ & Data set & Method & Calib.\\ mean & $R$ & $\\widehat R_{\\mathcal E}$ & $\\Pr(R>\\alpha)$ & $\\Pr(\\widehat R_{\\mathcal E}>\\alpha)$ & Genes w/ calls & Calls & Supp.\\ terms & Terms per supp.\\ gene \\\\ \\midrule']
NE = {'Wainberg': 4738, 'Sanger': 4465, 'DRIVE': 1885, 'HAP1': 4506}
for delta in [0.5, 0.1]:
    for ds in DS:
        for m, lab2 in [('GoDag', 'GoDag'), ('GoCo' if 'GoCo' in set(d.method) else 'GoCo-learned', 'GoCo')]:
            s = d[(d.dataset == ds) & (d.alpha == 0.1) & (d.delta == delta) & (d.method == m)]
            tpg = (s.go_yield / s.unit_yield).mean()
            L.append(f'{delta:.2f} & {ds} & {lab2} & {s.cal_mean.mean():.4f} & {s.pool_risk.mean():.4f} & {s.unit_fdp.mean():.4f} & {s.pool_exceed.mean():.2f} & {s.exceed.mean():.2f} & {s.genes_with_calls.mean():.0f} & {s.total_calls.mean():.0f} & {s.go_yield.mean():.1f} & {tpg:.2f} \\\\')
    L.append('\\midrule')
L[-1] = '\\bottomrule'; L += ['\\end{tabular}', '\\end{table*}']
open(REV + 'tables/TableS5_risk_decomposition.tex', 'w', encoding='utf-8').write(fit(L))
print('S5 written')
