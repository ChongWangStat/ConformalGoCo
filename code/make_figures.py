"""Figures of the GoCo manuscript. Reads results/all_methods_harmonised_100splits.csv (written by make_tables.py) and the
frozen Wainberg inputs; writes paper/figures/*.pdf and *.png."""
import pandas as pd, numpy as np, os, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.environ.get('GOCO_PAPER_DIR', os.path.join(HERE, '..', 'paper')) + os.sep
RESULTS = os.environ.get('GOCO_SPLITS', os.path.join(HERE, '..', 'results')) + os.sep
OUT = REV + 'figures/'; os.makedirs(OUT, exist_ok=True)
d = pd.read_csv(RESULTS + 'all_methods_harmonised_100splits.csv')
# v4: long frame with GoCo-N and the ranking-fold-selected GoCo (written by make_v4_tables.py); falls back to the harmonised file
LONG = os.environ.get('GOCO_LONG', '')
if LONG and os.path.exists(LONG):
    lf = pd.read_csv(LONG).rename(columns={'repeat': 'seed', 'arm': 'method', 'fdp': 'unit_fdp', 'supp_terms': 'go_yield'})
    lf['method'] = lf.method.replace({'Multilabel': 'Boger'})
    d = lf
    METHODS = ['Boger', 'GoDag', 'GoCo-M']
else:
    METHODS = ['Boger', 'GoDag', 'GoCo']
DS = ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
C = {'Boger': '#bdbdbd', 'GoDag': '#1f5fa8', 'GoCo': '#7b1fa2', 'GoCo-M': '#7b1fa2', 'GoCo-N': '#e67e22', 'GoCo-sel': '#c0392b'}
NAME = {'Boger': 'Multilabel (Direct)', 'GoDag': 'GoDag', 'GoCo': 'GoCo-M', 'GoCo-M': 'GoCo-M', 'GoCo-N': 'GoCo-N', 'GoCo-sel': 'GoCo, $\\mathcal{T}$-selected'}
MARK = {'Boger': 'o', 'GoDag': 'D', 'GoCo': '*', 'GoCo-M': '*', 'GoCo-N': '^', 'GoCo-sel': 's'}
GOCOM = 'GoCo-M' if 'GoCo-M' in METHODS else 'GoCo'

def save(fig, name):
    fig.savefig(OUT + name + '.pdf', bbox_inches='tight'); fig.savefig(OUT + name + '.png', dpi=300, bbox_inches='tight'); plt.close(fig)

# ------------------------------------------------------------------ Figure 3: primary comparison
def fig_primary(delta, name):
    methods = METHODS
    fig, axes = plt.subplots(2, 4, figsize=(11, 6.0), gridspec_kw={'hspace': 0.95, 'wspace': 0.35})
    for j, ds in enumerate(DS):
        sub = d[(d.dataset == ds) & (d.alpha == 0.1) & (d.delta == delta)]
        axes[0, j].set_xlim(0.065, 0.107); axes[0, j].xaxis.set_major_locator(plt.MaxNLocator(4))
        for row, (col, xlabel) in enumerate([('unit_fdp', 'held-out gene-level FDP'), ('go_yield', 'correct calls')]):
            ax = axes[row, j]
            for k, m in enumerate(methods):
                s = sub[sub.method == m][col]
                if len(s) == 0: continue
                q1, q3 = s.quantile([.25, .75])
                ax.plot([q1, q3], [k, k], color=C[m], lw=1.2, alpha=0.8)
                ax.plot(s.mean(), k, marker=MARK[m], ms=8 if MARK[m] == '*' else 5, color=C[m], mec='k', mew=0.4, zorder=3)
            if row == 0:
                ax.axvline(0.10, ls='--', color='k', lw=0.8); ax.set_title(ds, fontsize=9, fontweight='bold')
            ax.set_yticks(range(len(methods))); ax.set_yticklabels([NAME[m] for m in methods] if j == 0 else []); ax.invert_yaxis()
            ax.set_xlabel(xlabel, fontsize=7.5); ax.grid(axis='x', lw=0.3, alpha=0.5)
    axes[0, 0].text(-0.02, 1.28, 'A  Controlled risk (target 0.10)', transform=axes[0, 0].transAxes, fontsize=9, fontweight='bold')
    axes[1, 0].text(-0.02, 1.12, 'B  Correct calls released', transform=axes[1, 0].transAxes, fontsize=9, fontweight='bold')
    fig.suptitle(f'$\\alpha=0.10$, $\\delta={delta:.2f}$: mean over 100 splits (marker) and interquartile range (bar)', y=1.02, fontsize=9)
    save(fig, name)
fig_primary(0.50, 'Figure3_primary_delta050'); fig_primary(0.10, 'FigureS1_primary_delta010')

# ------------------------------------------------------------------ Figure 4: alpha sweep
def fig_sweep(delta, name):
    methods = METHODS
    fig, axes = plt.subplots(2, 4, figsize=(11, 5.4), gridspec_kw={'hspace': 0.85, 'wspace': 0.35})
    for j, ds in enumerate(DS):
        sub = d[(d.dataset == ds) & (d.delta == delta)]
        ref = sub[sub.method == 'Boger'].groupby('alpha').go_yield.mean()
        for m in methods:
            g = sub[sub.method == m].groupby('alpha').agg(fdp=('unit_fdp', 'mean'), terms=('go_yield', 'mean'))
            if len(g) == 0: continue
            ls = '--' if m == 'GoCo-sel' else '-'
            axes[0, j].plot(g.index, g.fdp, marker=MARK[m], ms=4 if MARK[m] != '*' else 6, color=C[m], label=NAME[m], ls=ls)
            axes[1, j].plot(g.index, 100 * (g.terms / ref.loc[g.index] - 1), marker=MARK[m], ms=4 if MARK[m] != '*' else 6, color=C[m], label=NAME[m], ls=ls)
        axes[0, j].plot([0.04, 0.21], [0.04, 0.21], ls='--', color='k', lw=0.8); axes[0, j].set_title(ds, fontsize=9, fontweight='bold')
        axes[1, j].axhline(0, color='k', lw=0.8)
        for r in range(2):
            axes[r, j].set_xticks([0.05, 0.10, 0.20]); axes[r, j].set_xlabel(r'target $\alpha$', fontsize=7.5); axes[r, j].grid(lw=0.3, alpha=0.5)
    axes[0, 0].set_ylabel('mean held-out gene-level FDP'); axes[1, 0].set_ylabel('correct calls, % change vs Multilabel')
    axes[0, 0].legend(fontsize=6.5, frameon=False, loc='upper left')
    axes[0, 0].text(-0.02, 1.25, 'A  Held-out risk tracks the target', transform=axes[0, 0].transAxes, fontsize=9, fontweight='bold')
    axes[1, 0].text(-0.02, 1.1, 'B  Yield relative to Multilabel at each target', transform=axes[1, 0].transAxes, fontsize=9, fontweight='bold')
    save(fig, name)
fig_sweep(0.50, 'Figure4_alpha_sweep_delta050')

# ------------------------------------------------------------------ Figure 2: ties and the risk staircase
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import goco_rerun as gr
    _dset, _tset = gr.load_truth()
    ds = gr.Dataset('Wainberg', _dset, _tset)
    sc = np.round(ds.rows.score.to_numpy(float), 6)
    md = ds.mats(ds.dense)
    have_w = True
except Exception as e:
    raise RuntimeError('Wainberg inputs required for Figure 2') from e
fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), gridspec_kw={'wspace': 0.35, 'width_ratios': [1, 1]})
if have_w:
    ax = axes[0]
    vals, cnts = np.unique(sc[(sc >= 600) & (sc <= 1000)], return_counts=True)
    ax.vlines(vals, 0, cnts, color='#555', lw=0.8)
    tie = 789.65; ntie = int(cnts[np.argmin(np.abs(vals - tie))])
    ax.vlines([tie], 0, ntie, color=C['GoCo'], lw=2.2)
    ax.annotate(f'{ntie:,} gene–GO pairs\nshare the score {tie:g}', xy=(tie, ntie), xytext=(830, 500), fontsize=7, arrowprops=dict(arrowstyle='->', lw=0.6)); ax.text(602, 3200, f'largest tie: {cnts.max():,} pairs at {vals[np.argmax(cnts)]:g}', fontsize=6.5, color='#333')
    for t in ds.grid[(ds.grid >= 600) & (ds.grid <= 1000)]:
        ax.axvline(t, color='#1f5fa8', lw=0.4, alpha=0.5)
    ax.set_xlabel('Wainberg enrichment score (600–1000)'); ax.set_ylabel('gene–GO pairs at this exact score'); ax.set_yscale('log'); ax.set_ylim(0.7, 6000)
    ax.set_title('A  Enrichment scores are heavily tied', fontsize=9, fontweight='bold', loc='left')
    ax = axes[1]
    thr = md['thr']; risk = md['LT'].mean(axis=0); calls = md['C'].sum(axis=0)
    keep = (thr >= 700) & (thr <= 900)
    ax.step(thr[keep], risk[keep], where='post', color='#1f5fa8', lw=1.2, label='global threshold path (all distinct scores)')
    ax.axhline(0.10, ls='--', color='k', lw=0.8)
    # risk of partial admission inside the 800->775 block interpolates between the two ends (Lemma 1a); GoCo's certified policy marked by the mean over splits
    ax.plot([800, 775], [risk[np.argmin(np.abs(thr - 800))], risk[np.argmin(np.abs(thr - 775))]], color='#8c8c8c', lw=1.2, ls=':', label='partial admission across the block (interpolated risk)')
    ax.scatter([800], [0.0923], color='#1f5fa8', zorder=3, s=25, label='GoDag stops at 800')
    star = float(d[(d.alpha == 0.1) & (d.delta == 0.5) & (d.dataset == 'Wainberg') & (d.method == GOCOM)].pool_risk.mean())
    ax.scatter([789.65], [star], color=C['GoCo'], marker='*', s=90, zorder=4, label='GoCo-M certified policy (mean pool risk)')
    ax.set_xlim(905, 695); ax.set_xlabel('score threshold (decreasing → more liberal)'); ax.set_ylabel('full-panel TruePath risk')
    ax.set_title('B  The 800→775 step overshoots α', fontsize=9, fontweight='bold', loc='left'); ax.legend(fontsize=6.3, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2)
save(fig, 'Figure2_ties')

from pathlib import Path
from make_workflow_figure import build
build(Path(OUT))
print('Five manuscript figures written to', OUT)
