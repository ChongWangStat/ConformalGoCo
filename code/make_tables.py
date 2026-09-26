"""Build the LaTeX/CSV tables of the GoCo manuscript from results/goco_results_100splits.csv and the frozen
generic-calibrator summary results/first_draft_primary_method_summary_long.csv."""
import pandas as pd, numpy as np, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep
RESULTS = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results')) + os.sep
OUTDIR = REV + "tables/"; os.makedirs(OUTDIR, exist_ok=True)
DS = ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']
N = {'Wainberg': 15794, 'Sanger': 14885, 'DRIVE': 6283, 'HAP1': 15020}
NT = {k: round(.1 * v) for k, v in N.items()}; NC70 = {k: round(.7 * v) for k, v in N.items()}
NE = {k: N[k] - NC70[k] for k in N}; M60 = {k: NC70[k] - NT[k] for k in N}

frames = []
SPLITDIR = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results'))  # frozen split-level outputs shipped with the repo
for f in glob.glob(os.path.join(SPLITDIR, '*_100splits.csv')):
    if 'repro' in f or 'prereview' in f or 'harmonised' in f: continue   # never re-read this script's own output
    x = pd.read_csv(f); x['src'] = os.path.basename(f); frames.append(x)
d = pd.concat(frames, ignore_index=True)
for k, v in {' (paper path)': '', ' (tie-level path)': '', ' (global, all distinct scores)': '', ' (global grid)': '', ' (Direct, global grid)': '',
             ' (global grid, certification fold only)': '-cert60', ' (grid path, all blocks)': ''}.items():
    d['method'] = d.method.str.replace(k, v, regex=False)
d = d.drop_duplicates(['dataset', 'seed', 'alpha', 'delta', 'method'])
d['exceed'] = (d.unit_fdp > d.alpha).astype(float)
d['supp_frac'] = 1 - d.term_fdp
# Pool-level risk of the certified policy.  goco_rerun.metrics() records it directly as the mean loss over the
# certificate's own pool (its calibration fold, 70% for the global methods and 60% for GoCo, plus the evaluation fold)
# on BOTH losses: pool_risk (TruePath, the loss every table reports) and pool_risk_direct (Direct, the loss the
# Boger et al. certificate itself bounds).  The pre-review script reconstructed a Direct calibration mean with a
# TruePath evaluation mean for Boger et al.; that hybrid is gone.
for c in ['pool_risk', 'pool_risk_direct', 'm_cal']:
    assert c in d.columns, f'{c} missing: regenerate the split-level files with the revised goco_rerun.py'
d['pool_exceed'] = (d.pool_risk > d.alpha).astype(float)               # TruePath pool risk above alpha
d['pool_exceed_direct'] = (d.pool_risk_direct > d.alpha).astype(float) # Direct pool risk above alpha
d['pool_exceed_certified'] = np.where(d.method.str.startswith('Boger'), d.pool_exceed_direct, d.pool_exceed)  # the loss each certificate bounds
d.to_csv(RESULTS + 'all_methods_harmonised_100splits.csv', index=False)

LABEL = {'Boger': 'Boger et al.\\ (Direct loss)', 'GoDag': 'GoDag', 'GoCo': '\\textbf{GoCo}', 'GoCo-learned': '\\textbf{GoCo}'}

def write_tex(L, fname):
    """Join, make wide table* fit the text width, and make supplementary files self-contained."""
    txt = '\n'.join(L)
    if 'table*' in txt:
        txt = txt.replace('\\begin{tabular}', '\\resizebox{\\linewidth}{!}{\\begin{tabular}').replace('\\end{tabular}', '\\end{tabular}}')
    # threeparttable measures the unscaled tabular; replace it by a plain minipage for the notes
    txt = txt.replace('\\begin{threeparttable}\n', '').replace('\\end{threeparttable}\n', '')
    txt = txt.replace('\\begin{tablenotes}[flushleft]\\footnotesize\n\\item ', '\\par\\vspace{3pt}\\begin{minipage}{\\linewidth}\\footnotesize ').replace('\\end{tablenotes}', '\\end{minipage}')
    if os.path.basename(fname).startswith('TableS'):
        txt = txt.replace('Eq.~\\ref{eq:utility}', 'Eq.~(6)').replace('Proposition~\\ref{prop:validity}', 'Proposition~1').replace('[!t]', '[htbp]')
    open(fname, 'w', encoding='utf-8').write(txt)


def paired(sub, m, ref, col='go_yield'):
    a = sub[sub.method == m].set_index('seed')[col]; b = sub[sub.method == ref].set_index('seed')[col]
    ix = a.index.intersection(b.index); dif = (a.loc[ix] - b.loc[ix])
    if len(dif) == 0: return np.nan, np.nan, np.nan
    return dif.mean(), 1.96 * dif.std(ddof=1) / np.sqrt(len(dif)), (dif > 0).mean()

GOCO = 'GoCo' if 'GoCo' in set(d.method) else 'GoCo-learned'
GENERIC = pd.read_csv(RESULTS + 'first_draft_primary_method_summary_long.csv')
GEN_ORDER = [('Hoeffding/LTT', 'Hoeffding/LTT'), ('Hoeffding-Bentkus', 'Hoeffding--Bentkus'), ('IID-normal', 'IID-normal'), ('Empirical Bernstein', 'Empirical Bernstein'),
             ('Binary-incidence McDiarmid', 'McDiarmid (binary incidence)'), ('Janson dependency-graph', 'Janson dependency graph'),
             ('Network-HAC (b=1)', 'Network-HAC ($b=1$)'), ('Platt calibration', 'Platt calibration'), ('Isotonic calibration', 'Isotonic calibration')]

def main_table(delta, fname, caption, label):
    """Tables 1-2: one row per calibrator, two columns per data set (held-out gene-level FDP; supported GO terms = power)."""
    rows = []
    L = ['\\begin{table*}[!t]', '\\centering', f'\\caption{{{caption}}}', f'\\label{{{label}}}', '\\scriptsize', '\\setlength{\\tabcolsep}{4pt}', '\\begin{threeparttable}',
         '\\begin{tabular}{l rr rr rr rr}', '\\toprule',
         'Method & \\multicolumn{2}{c}{Wainberg ($n_{\\mathcal E}$=4,738)} & \\multicolumn{2}{c}{Sanger ($n_{\\mathcal E}$=4,465)} & \\multicolumn{2}{c}{DRIVE ($n_{\\mathcal E}$=1,885)} & \\multicolumn{2}{c}{HAP1 ($n_{\\mathcal E}$=4,506)} \\\\',
         '\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}',
         ' & Unit FDP & Supp.\\ genes & Unit FDP & Supp.\\ genes & Unit FDP & Supp.\\ genes & Unit FDP & Supp.\\ genes \\\\ \\midrule']
    def harmonised(m, ds):
        s = d[(d.dataset == ds) & (d.alpha == 0.10) & (d.delta == delta) & (d.method == m)]
        return s.unit_fdp.mean(), s.unit_yield.mean()
    def generic(m, ds):
        r = GENERIC[(GENERIC.dataset == ds) & (GENERIC.delta == delta) & (GENERIC.method == m)]
        return (r.unit_fdp.iloc[0], r.unit_yield.iloc[0]) if len(r) else (np.nan, np.nan)
    def add(lab, key, getter, group):
        cells = []
        for ds in DS:
            f, y = getter(key, ds); cells.append(f'{f:.4f} & {y:.1f}')
            rows.append(dict(delta=delta, group=group, method=key, dataset=ds, unit_fdp=f, supported_genes=y))
        L.append(f'{lab} & ' + ' & '.join(cells) + ' \\\\')
    L.append('\\multicolumn{9}{l}{\\emph{Direct loss, one global threshold (ontology-blind)}} \\\\')
    add('Boger et al.', 'Boger', harmonised, 'direct')
    for key, lab in GEN_ORDER: add(lab, key, generic, 'direct')
    L.append('\\addlinespace[2pt]\\multicolumn{9}{l}{\\emph{TruePath loss, one global threshold (GO hierarchy)}} \\\\')
    add('GoDag', 'GoDag', harmonised, 'hierarchy')
    L.append('\\addlinespace[2pt]\\multicolumn{9}{l}{\\emph{TruePath loss, relational partial admission (GO hierarchy and co-essentiality)}} \\\\')
    add('\\textbf{GoCo}', GOCO, harmonised, 'relational')
    L += ['\\bottomrule', '\\end{tabular}', '\\begin{tablenotes}[flushleft]\\footnotesize',
          '\\item Generic baselines calibrate on Direct labels only and have no GO-DAG/TruePath input; GoDag calibrates with TruePath; GoCo calibrates with TruePath plus predictor-native co-essentiality source information. All rows are evaluated on the same held-out TruePath reference. Every global-threshold method calibrates on the first 70\\% of each permutation over the shared 25-unit grid; GoCo certifies on the 60\\% certification fold after fitting its ordering on the 10\\% ranking fold. Platt and isotonic are $\\delta$-free practical calibrators without a risk guarantee, so their rows are the same in both tables. At $\\delta=0.50$ every normal-type rule coincides with Boger et al.\\ (Remark~\\ref{rem:delta}); at $\\delta=0.10$ they separate, only the IID-normal CLT still agreeing with it. Exceedance frequencies, released calls, supported GO-term yields and paired differences from GoDag for Boger et al., GoDag and GoCo are in Supplementary Table~S1.',
          '\\end{tablenotes}', '\\end{threeparttable}', '\\end{table*}']
    write_tex(L, OUTDIR + fname + '.tex'); pd.DataFrame(rows).to_csv(OUTDIR + fname + '.csv', index=False)

main_table(0.50, 'Table1_primary_alpha010_delta050', 'Gene-level TruePath performance at $\\alpha=0.10$ and $\\delta=0.50$ (the operating point of the released Boger et al.\\ implementation). Unit FDP is the mean held-out gene-level false discovery proportion, the formal risk target. Supported genes are held-out genes with at least one TruePath-supported released GO term. Values are means over the same 100 split streams.', 'tab:primary50')
main_table(0.10, 'Table2_primary_alpha010_delta010', 'Gene-level TruePath performance at $\\alpha=0.10$ and $\\delta=0.10$ (the substantive guarantee of Proposition~\\ref{prop:validity}); columns as in Table~\\ref{tab:primary50}.', 'tab:primary10')

def detail_table(fname, caption, label):
    """Table S1: the full split-level statistics for the three main methods at both delta."""
    methods = ['Boger', 'GoDag', GOCO]
    rows = []
    L = ['\\begin{table*}[!t]', '\\centering', f'\\caption{{{caption}}}', f'\\label{{{label}}}', '\\scriptsize', '\\setlength{\\tabcolsep}{3.2pt}', '\\begin{threeparttable}',
         '\\begin{tabular}{lll rrr rrrr r}', '\\toprule',
         '$\\delta$ & Data set & Method & Gene FDP & $\\Pr(\\widehat R_{\\mathcal E}>\\alpha)$ & $\\Pr(R>\\alpha)$ & Supp.\\ genes & Calls & Supp.\\ terms & Supp.\\ frac. & $\\Delta$terms vs GoDag [95\\% CI] \\\\', '\\midrule']
    for delta in [0.50, 0.10]:
        for ds in DS:
            sub = d[(d.dataset == ds) & (d.alpha == 0.10) & (d.delta == delta)]
            first = True
            for m in methods:
                s = sub[sub.method == m]
                if len(s) == 0: continue
                md, ci, p = paired(sub, m, 'GoDag')
                dstr = '--' if m == 'GoDag' else f'{md:+.1f} [{md-ci:+.1f}, {md+ci:+.1f}]'
                dl = f'{delta:.2f}' if (first and ds == DS[0]) else ''
                L.append(f'{dl} & {ds if first else ""} & {LABEL[m]} & {s.unit_fdp.mean():.4f} & {s.exceed.mean():.2f} & {s.pool_exceed.mean():.2f} & {s.unit_yield.mean():.1f} & {s.total_calls.mean():.0f} & {s.go_yield.mean():.1f} & {s.supp_frac.mean():.3f} & {dstr} \\\\')
                first = False
                rows.append(dict(dataset=ds, delta=delta, method=m, unit_fdp=s.unit_fdp.mean(), heldout_exceed=s.exceed.mean(), pool_exceed=s.pool_exceed.mean(), pool_exceed_direct=s.pool_exceed_direct.mean(), pool_exceed_certified=s.pool_exceed_certified.mean(), genes_with_calls=s.genes_with_calls.mean(), supported_genes=s.unit_yield.mean(), calls=s.total_calls.mean(), supported_terms=s.go_yield.mean(), supported_fraction=s.supp_frac.mean(), dterms_vs_godag=md, ci95=ci, p_gt=p))
        L.append('\\midrule')
    L[-1] = '\\bottomrule'
    L += ['\\end{tabular}', '\\begin{tablenotes}[flushleft]\\footnotesize',
          '\\item Means over the same 100 split streams at $\\alpha=0.10$. Gene FDP: held-out gene-level TruePath FDP $\\widehat R_{\\mathcal E}$. $\\Pr(\\widehat R_{\\mathcal E}>\\alpha)$: fraction of splits whose held-out FDP exceeded $\\alpha$. $\\Pr(R>\\alpha)$: fraction of splits in which the TruePath risk of the certified policy over the pool of its own certificate (calibration fold plus evaluation fold) exceeded $\\alpha$; for Boger et al., whose certificate bounds the Direct loss, the corresponding Direct-loss frequency is in the summary CSV (pool\\_exceed\\_direct). Supp.\\ genes / terms: held-out genes with $\\ge1$ TruePath-supported released term / supported released gene--GO pairs. Supp.\\ frac.: supported fraction of all released calls. Boger et al.\\ is the released implementation on the Direct loss, calibrated on 70\\% of the panel; GoDag and GoCo use the TruePath loss and share the 25-unit grid. The last column is the paired mean difference in supported terms from GoDag with a descriptive 95\\% interval over splits.',
          '\\end{tablenotes}', '\\end{threeparttable}', '\\end{table*}']
    write_tex(L, OUTDIR + fname + '.tex'); pd.DataFrame(rows).to_csv(OUTDIR + fname + '.csv', index=False)

detail_table('TableS1_detail', 'Split-level statistics for Boger et al., GoDag and GoCo at $\\alpha=0.10$ and $\\delta\\in\\{0.50,0.10\\}$: exceedance frequencies, released calls, supported genes and terms, supported fraction and paired differences from GoDag.', 'tab:detail')

def sweep_table(delta, fname, caption, label):
    methods = ['Boger', 'GoDag', 'GoCo']
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

sweep_table(0.50, 'Table4_alpha_sweep_delta050', 'Error-target sweep at $\\delta=0.50$: mean held-out unit FDP and supported GO-term yield at $\\alpha\\in\\{0.05,0.10,0.20\\}$. ', 'tab:sweep50')
sweep_table(0.10, 'TableS3_alpha_sweep_delta010', 'Error-target sweep at $\\delta=0.10$.', 'tab:sweep10')

# ---------------- summary CSV of everything ----------------
s = d.groupby(['dataset', 'alpha', 'delta', 'method']).agg(unit_fdp=('unit_fdp', 'mean'), heldout_exceed=('exceed', 'mean'), pool_exceed=('pool_exceed', 'mean'), mean_pool_risk=('pool_risk', 'mean'), pool_exceed_direct=('pool_exceed_direct', 'mean'), mean_pool_risk_direct=('pool_risk_direct', 'mean'), pool_exceed_certified=('pool_exceed_certified', 'mean'), m_cal=('m_cal', 'first'),
                                                          genes_with_calls=('genes_with_calls', 'mean'), supported_genes=('unit_yield', 'mean'), calls=('total_calls', 'mean'), supported_terms=('go_yield', 'mean'), supported_fraction=('supp_frac', 'mean'), sd_terms=('go_yield', 'std')).round(4)
s.to_csv(RESULTS + 'summary_by_method_harmonised.csv')
print(s.loc[(slice(None), 0.1, 0.5), :][['unit_fdp', 'heldout_exceed', 'pool_exceed', 'supported_terms', 'calls']].to_string())
print('tables written to', OUTDIR)
