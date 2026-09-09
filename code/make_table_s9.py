"""Supplementary Table S6 (file TableS9_estimands.tex): the two estimands on the same data (conformal selection vs GoCo)."""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep
RESULTS = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results')) + os.sep
B = chr(92)
NL = B + B

cs = pd.read_csv(RESULTS + 'conformal_selection_alpha0.1.csv')
g = pd.read_csv(RESULTS + 'all_methods_harmonised_100splits.csv')
GOCO = 'GoCo' if 'GoCo' in set(g.method) else 'GoCo-learned'
DS = ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']

cap = ("Two estimands on the same data, same 100 splits, same frozen scores, $" + B + "alpha=0.10$ and $" + B + "delta=0.50$. Conformal selection "
       "computes a conformal $p$-value for each gene--GO pair against the unsupported calibration pairs of the split and "
       "applies Benjamini--Hochberg, targeting the " + B + "emph{pooled} false discovery rate over released calls; GoCo targets "
       "the average over genes of the per-gene FDP. Conformal selection does what it promises---its realised pooled FDR is "
       "at or below the target in every data set---but reaching that target against an incomplete reference requires "
       "releasing almost nothing. Counts are means per evaluation fold; the realised pooled FDR is the mean over the 100 splits of each split's pooled false discovery proportion, not the ratio of the two count columns.")

L = [B + 'begin{table*}[!htbp]', B + 'centering', B + 'caption{' + cap + '}', B + 'label{tab:estimands}', B + 'scriptsize',
     B + 'begin{tabular}{l l rrrr}', B + 'toprule',
     'Data set & Target & Calls released & Supported calls & Realised pooled FDR & Genes annotated ' + NL + ' ' + B + 'midrule']
rows = []
for ds in DS:
    c = cs[cs.dataset == ds]
    G = g[(g.dataset == ds) & (g.alpha == 0.1) & (g.delta == 0.5) & (g.method == GOCO)]
    L.append(f'{ds} & pooled FDR over calls & {c.selected_calls.mean():.1f} & {c.supported_calls.mean():.1f} & '
             f'{c.pooled_fdr.mean():.3f} & {c.genes_with_calls.mean():.1f} ' + NL)
    L.append(f' & gene-level FDP (GoCo) & {G.total_calls.mean():.0f} & {G.go_yield.mean():.1f} & '
             f'{G.term_fdp.mean():.3f} & {G.genes_with_calls.mean():.0f} ' + NL)
    L.append(B + 'addlinespace')
    rows.append(dict(dataset=ds, cs_calls=c.selected_calls.mean(), cs_supported=c.supported_calls.mean(),
                     cs_pooled_fdr=c.pooled_fdr.mean(), goco_calls=G.total_calls.mean(),
                     goco_supported=G.go_yield.mean(), goco_pooled_fdp=G.term_fdp.mean()))
L = L[:-1] + [B + 'bottomrule', B + 'end{tabular}', B + 'end{table*}']
open(REV + 'tables/TableS9_estimands.tex', 'w', encoding='utf-8').write('\n'.join(L) + '\n')
pd.DataFrame(rows).to_csv(REV + 'tables/TableS9_estimands.csv', index=False)
print('estimand table written')
