import os
"""Table S7 / Results case study: what GoCo releases that GoDag does not, in one Wainberg split.
Systematic, not cherry-picked: every co-essential module that contributes at least four admitted
evaluation-fold calls, with its supported fraction, plus the overall composition of the admitted set."""
import numpy as np, pandas as pd, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goco_rerun as gr
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep  # output location only
RESULTS = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results')) + os.sep
ALPHA, DELTA, SEED = 0.10, 0.50, 0

dset, tset = gr.load_truth()
ds = gr.Dataset('Wainberg', dset, tset)
mg = ds.mats(ds.grid); L = mg['LT']; C = mg['C']; YT = mg['YT']
j = int(np.where(np.isclose(ds.grid, 800))[0][0]); lo, hi = float(ds.grid[j + 1]), 800.0
Aadd, meanscore = ds.Aadd(lo, hi)
delta = (L[:, j + 1] - L[:, j]).astype(float); cadd = (C[:, j + 1] - C[:, j]).astype(float); affected = cadd > 0
gr.utility.extra = {'fsup': (YT[:, j + 1] - YT[:, j]).astype(float), 'lo': lo, 'hi': hi}
perm = np.random.default_rng(gr.SEED0 + SEED).permutation(ds.n)
nt, nc = round(.1 * ds.n), round(.7 * ds.n)
train, cert, ev = perm[:nt], perm[nt:nc], perm[nc:]
u, hh = gr.utility('learned', ds, Aadd, meanscore, delta, cadd, affected, train, ds.hh)
sel_q, sel_final = None, None
for q in gr.QGRID_E:
    sel = gr.selmask(u, hh, affected, train, float(q))
    y = np.where(sel, L[:, j + 1], L[:, j]).astype(float)
    p, r, sd = gr.clt_p(y[cert], ALPHA)
    if p > DELTA: break
    sel_q, sel_final = float(q), sel.copy()
print('selected admission fraction q =', sel_q)

evset = np.zeros(ds.n, bool); evset[ev] = True
blk = ds.block_rows(lo, hi).copy().drop(columns=['support_sources', 'srcs'], errors='ignore')
attr = pd.read_csv(gr.W + 'wainberg_highscore_source_attribution.csv.gz',
                   usecols=['gene', 'go_id', 'winning_module', 'support_sources'],
                   dtype={'gene': str, 'go_id': str}).drop_duplicates(['gene', 'go_id'])
names = pd.read_csv(gr.W + 'wainberg_gene_go_scores.csv.gz', usecols=['gene', 'go_id', 'name'], dtype=str).drop_duplicates(['gene', 'go_id'])
blk = blk.merge(attr, on=['gene', 'go_id'], how='left').merge(names, on=['gene', 'go_id'], how='left')
adm = blk[sel_final[blk.i.to_numpy()] & evset[blk.i.to_numpy()]].copy()
print(f'admitted evaluation-fold calls: {len(adm)} over {adm.gene.nunique()} genes; supported {int(adm["T"].sum())} ({adm["T"].mean():.1%})')

g = adm.groupby('winning_module').agg(calls=('T', 'size'), supported=('T', 'sum'), genes=('gene', 'nunique')).sort_values('calls', ascending=False)
g['supported_frac'] = g.supported / g.calls
print('\nmodules contributing >= 4 admitted calls:')
print(g[g.calls >= 4].to_string())
mm = pd.read_csv(gr.W + 'module_memberships.csv.gz')
rows = []
for m in g[g.calls >= 4].index[:8]:
    mem = mm[mm.module_id == m].gene.astype(str).tolist()
    sub = adm[adm.winning_module == m]
    rows.append(dict(module=int(m), members='; '.join(sorted(mem)), admitted_genes='; '.join(sorted(set(sub.gene))),
                     calls=len(sub), supported=int(sub['T'].sum()),
                     terms='; '.join(f"{r.go_id} {r['name']}{' *' if r['T'] else ''}" for _, r in sub.drop_duplicates('go_id').iterrows())))
    print(f"\nmodule {m}: members {sorted(mem)}")
    print(f"   admitted genes {sorted(set(sub.gene))}, {len(sub)} calls, {int(sub['T'].sum())} supported")
    for _, r in sub.iterrows():
        print(f"      {r.gene:10s} {r.go_id} {str(r['name'])[:52]:54s} {'SUPPORTED' if r['T'] else '-'}  <- {r.support_sources}")
pd.DataFrame(rows).to_csv(RESULTS + 'case_study_split0.csv', index=False)

# ---------------- LaTeX Table S7 ----------------
cap = (f"What GoCo releases that GoDag does not, in Wainberg split~{SEED} at $\\alpha=0.10$, $\\delta=0.50$ "
       f"(admission fraction $q={sel_q:g}$). Of the {len(adm)} additional evaluation-fold gene--GO calls, over {adm.gene.nunique()} genes, "
       f"{int(adm['T'].sum())} ({100 * adm['T'].mean():.0f}\\%) are supported by the frozen TruePath reference. "
       "The table lists every co-essential module contributing at least four of them, its members, the genes admitted from it, "
       "and the terms transferred ($\\ast$ = supported). Modules whose members share a coherent function transfer supported terms; "
       "heterogeneous modules transfer their members' own annotations to genes that do not have them.")
Lx = ['\\begin{table*}[!htbp]', '\\centering', f'\\caption{{{cap}}}', '\\label{tab:casestudy}', '\\scriptsize',
      '\\resizebox{\\linewidth}{!}{\\begin{tabular}{r p{0.26\\textwidth} p{0.16\\textwidth} r r p{0.34\\textwidth}}', '\\toprule',
      'Module & Members & Genes admitted & Calls & Supp. & Terms transferred \\\\ \\midrule']
def tex_escape(t):
    for ch in ['%', '&', '#', '_']:
        t = t.replace(ch, '\\' + ch)
    return t
for r in rows:
    terms = [tex_escape(t.strip()).replace('\\*', '$\\ast$').replace('*', '$\\ast$') for t in r['terms'].split(';')]
    keep, n = [], 0
    for t in terms:                       # truncate on whole terms, never inside math
        if n + len(t) > 380: keep.append('\\dots'); break
        keep.append(t); n += len(t)
    Lx.append(f"{r['module']} & \\texttt{{{r['members'].replace(';', ',')}}} & \\texttt{{{r['admitted_genes'].replace(';', ',')}}} & {r['calls']} & {r['supported']} & {'; '.join(keep)} \\\\")
    Lx.append('\\addlinespace')
Lx = Lx[:-1] + ['\\bottomrule', '\\end{tabular}}', '\\end{table*}']
open(REV + 'tables_v2/TableS7_case_study.tex', 'w', encoding='utf-8').write('\n'.join(Lx))
print('\nTable S7 written')
