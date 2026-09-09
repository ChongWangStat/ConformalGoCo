"""Supplementary Tables S1 (generic-calibrator panel of the first draft) and S5 (risk decomposition)."""
import pandas as pd, os
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep
RESULTS = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results')) + os.sep
BIB = RESULTS + 'first_draft_primary_method_summary_long.csv'
b = pd.read_csv(BIB)
order = ['Boger et al.', 'Hoeffding/LTT', 'Hoeffding-Bentkus', 'IID-normal', 'Empirical Bernstein', 'Binary-incidence McDiarmid', 'Janson dependency-graph', 'Network-HAC (b=1)', 'Platt calibration', 'Isotonic calibration', 'GoDag']
lab = {'Boger et al.': 'Boger et al.\\ (CLT, Direct)', 'Hoeffding/LTT': 'Hoeffding LTT', 'Hoeffding-Bentkus': 'Hoeffding--Bentkus', 'IID-normal': 'IID-normal CLT', 'Empirical Bernstein': 'Empirical Bernstein',
       'Binary-incidence McDiarmid': 'McDiarmid (binary source incidence)', 'Janson dependency-graph': 'Janson dependency graph', 'Network-HAC (b=1)': 'Network-HAC ($b=1$)',
       'Platt calibration': 'Platt calibration ($\\delta$-free)', 'Isotonic calibration': 'Isotonic calibration ($\\delta$-free)', 'GoDag': 'GoDag (TruePath)'}
DS = ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']
L = ['\\begin{table*}[!t]', '\\centering',
     '\\caption{Generic Direct-loss calibrators of the first draft at $\\alpha=0.10$ (same 100 splits, same held-out TruePath truth). Every generic method calibrates on Direct labels only, on the shared 25-unit grid, over the first 70\\% of each permutation; none exceeds GoDag. At $\\delta=0.50$ all normal-type members coincide with Boger et al.\\ (Remark 1).}',
     '\\label{tab:generic}', '\\scriptsize', '\\begin{tabular}{l rr rr rr rr}', '\\toprule',
     'Method & \\multicolumn{2}{c}{Wainberg} & \\multicolumn{2}{c}{Sanger} & \\multicolumn{2}{c}{DRIVE} & \\multicolumn{2}{c}{HAP1} \\\\',
     ' & FDP & Terms & FDP & Terms & FDP & Terms & FDP & Terms \\\\ \\midrule']
PUBLISHED_GODAG = {'Wainberg': (0.0832, 301.5), 'Sanger': (0.0790, 293.5), 'DRIVE': (0.0806, 193.1), 'HAP1': (0.0721, 130.5)}  # Hoeffding p-value, TruePath loss, delta=0.10, 70% calibration (earlier packages)
for delta in [0.5, 0.1]:
    L.append(f'\\multicolumn{{9}}{{l}}{{\\emph{{$\\delta={delta:.2f}$}}}} \\\\')
    for m in order:
        cells = []
        for ds in DS:
            r = b[(b.dataset == ds) & (b.delta == delta) & (b.method == m)]
            cells.append(f'{r.unit_fdp.iloc[0]:.4f} & {r.go_yield.iloc[0]:.1f}' if len(r) else '-- & --')
        L.append(f'{lab[m]} & ' + ' & '.join(cells) + ' \\\\')
    if delta == 0.1:
        L.append('GoDag as published (Hoeffding $p$-value, TruePath)\\tnote{a} & ' + ' & '.join(f'{PUBLISHED_GODAG[ds][0]:.4f} & {PUBLISHED_GODAG[ds][1]:.1f}' for ds in DS) + ' \\\\')
    L.append('\\midrule')
L[-1] = '\\bottomrule'; L += ['\\end{tabular}', '\\par\\vspace{3pt}\\begin{minipage}{\\linewidth}\\footnotesize $^{\\mathrm a}$The original ConformalGoDag certificate (Hoeffding $p$-value $\\exp\\{-2m(\\alpha-\\widehat R_k)_+^2\\}$ on the TruePath loss, 70\\% calibration, $\\delta=0.10$), taken from our earlier analysis packages on the same 100 splits; replacing the Hoeffding bound by the normal statistic accounts for the difference from the GoDag (TruePath) row above, and the GoCo path accounts for the further gain reported in the main text.\\end{minipage}', '\\end{table*}']
for k, line in enumerate(L):
    L[k] = line.replace('\\tnote{a}', '$^{\\mathrm a}$')
def fit(L):
    return '\n'.join(L).replace('\\begin{tabular}', '\\resizebox{\\linewidth}{!}{\\begin{tabular}').replace('\\end{tabular}', '\\end{tabular}}')
open(REV + 'tables/TableS1_generic_panel.tex', 'w', encoding='utf-8').write(fit(L))

d = pd.read_csv(RESULTS + 'all_methods_harmonised_100splits.csv')
L = ['\\begin{table*}[!t]', '\\centering',
     '\\caption{Calibration mean ($\\widehat R_{\\mathcal C}$ over the 60\\% certification fold for GoCo; $\\widehat R_{\\mathcal T\\cup\\mathcal C}$ over the 70\\% calibration sample for GoDag), pool risk $R(\\widehat\\pi)$ and held-out mean $\\widehat R_{\\mathcal E}(\\widehat\\pi)$ of the certified policy, averaged over 100 splits at $\\alpha=0.10$, with the fractions of splits in which $R(\\widehat\\pi)>\\alpha$ (the event bounded by Proposition 1) and $\\widehat R_{\\mathcal E}(\\widehat\\pi)>\\alpha$; the number of held-out genes with a nonempty released set (of $n_{\\mathcal E}$ = 4,738 / 4,465 / 1,885 / 4,506) and the mean FDP among those genes, $R^{+}=\\widehat R_{\\mathcal E}\\,n_{\\mathcal E}/(\\text{genes with calls})$.}',
     '\\label{tab:riskdecomp}', '\\scriptsize', '\\begin{tabular}{l l l rrr rr rr}', '\\toprule',
     '$\\delta$ & Data set & Method & Calib.\\ mean & $R$ & $\\widehat R_{\\mathcal E}$ & $\\Pr(R>\\alpha)$ & $\\Pr(\\widehat R_{\\mathcal E}>\\alpha)$ & Genes w/ calls & $R^{+}$ \\\\ \\midrule']
NE = {'Wainberg': 4738, 'Sanger': 4465, 'DRIVE': 1885, 'HAP1': 4506}
for delta in [0.5, 0.1]:
    for ds in DS:
        for m, lab2 in [('GoDag', 'GoDag'), ('GoCo-learned', 'GoCo')]:
            s = d[(d.dataset == ds) & (d.alpha == 0.1) & (d.delta == delta) & (d.method == m)]
            rplus = (s.unit_fdp * NE[ds] / s.genes_with_calls).mean()
            L.append(f'{delta:.2f} & {ds} & {lab2} & {s.cal_mean.mean():.4f} & {s.pool_risk.mean():.4f} & {s.unit_fdp.mean():.4f} & {s.pool_exceed.mean():.2f} & {s.exceed.mean():.2f} & {s.genes_with_calls.mean():.0f} & {rplus:.2f} \\\\')
    L.append('\\midrule')
L[-1] = '\\bottomrule'; L += ['\\end{tabular}', '\\end{table*}']
open(REV + 'tables/TableS5_risk_decomposition.tex', 'w', encoding='utf-8').write(fit(L))
print('S1, S5 written')
