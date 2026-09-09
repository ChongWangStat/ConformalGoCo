"""Comparator: conformal selection (Jin & Candes 2023) with Benjamini-Hochberg over gene-GO pairs.

This targets a DIFFERENT estimand from GoCo: the pooled false discovery rate among released
gene-GO calls, rather than the average over genes of the per-gene FDP. It is the alternative a
reader will ask about, so we report what it does on the same splits and the same frozen scores.

For each candidate pair (i,g) with score s_ig, a conformal p-value is computed against the
unsupported (null) calibration pairs of the same split,
    p_ig = (1 + #{j in cal : y_j = 0, s_j >= s_ig}) / (1 + #{j in cal : y_j = 0}),
which is valid for exchangeable pairs; BH at level alpha is then applied to the evaluation pairs.
Reported: the number of selected pairs, the realised pooled FDR among them, and the number of
reference-supported selected pairs, alongside GoCo at the same alpha.
Usage: python conformal_selection.py [alpha] [nsplits]
"""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goco_rerun as gr

ALPHA = float(sys.argv[1]) if len(sys.argv) > 1 else 0.10
NSPL = int(sys.argv[2]) if len(sys.argv) > 2 else 100
OUT = os.environ.get('GOCO_SPLITS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))


def bh(p, alpha):
    """Benjamini-Hochberg: return boolean mask of rejections."""
    n = len(p)
    if n == 0:
        return np.zeros(0, bool)
    o = np.argsort(p)
    thr = alpha * np.arange(1, n + 1) / n
    passed = p[o] <= thr
    k = np.nonzero(passed)[0]
    m = np.zeros(n, bool)
    if len(k):
        m[o[:k[-1] + 1]] = True
    return m


def run(name, dset, tset):
    ds = gr.Dataset(name, dset, tset)
    s = np.round(ds.rows.score.to_numpy(float), 6)
    y = ds.rows['T'].to_numpy(int)          # 1 = TruePath-supported
    gene = ds.rows.i.to_numpy()
    rows = []
    for seed in range(NSPL):
        perm = np.random.default_rng(gr.SEED0 + seed).permutation(ds.n)
        nc = round(.70 * ds.n)
        cal_g = np.zeros(ds.n, bool); cal_g[perm[:nc]] = True
        ev_g = ~cal_g
        cal = cal_g[gene]; ev = ev_g[gene]
        s_null = np.sort(s[cal & (y == 0)])          # calibration scores of unsupported pairs
        n0 = len(s_null)
        s_ev = s[ev]; y_ev = y[ev]
        # p = (1 + #{null cal scores >= s}) / (1 + n0)
        ge = n0 - np.searchsorted(s_null, s_ev, side='left')
        p = (1.0 + ge) / (1.0 + n0)
        sel = bh(p, ALPHA)
        k = int(sel.sum()); sup = int(y_ev[sel].sum())
        rows.append(dict(dataset=name, seed=seed, alpha=ALPHA, method='Conformal selection (BH over pairs)',
                         selected_calls=k, supported_calls=sup,
                         pooled_fdr=float((k - sup) / k) if k else 0.0,
                         genes_with_calls=int(len(np.unique(gene[ev][sel]))) if k else 0))
        if (seed + 1) % 25 == 0:
            print(f'[{name}] seed {seed+1}/{NSPL}', flush=True)
    d = pd.DataFrame(rows)
    print(f'[{name}] alpha={ALPHA}: selected {d.selected_calls.mean():.1f} calls, '
          f'{d.supported_calls.mean():.1f} supported, realised pooled FDR {d.pooled_fdr.mean():.3f}, '
          f'genes annotated {d.genes_with_calls.mean():.1f}', flush=True)
    return d


if __name__ == '__main__':
    dset, tset = gr.load_truth()
    out = pd.concat([run(n, dset, tset) for n in ['Wainberg', 'Sanger', 'DRIVE', 'HAP1']], ignore_index=True)
    os.makedirs(OUT, exist_ok=True)
    out.to_csv(os.path.join(OUT, f'conformal_selection_alpha{ALPHA:g}.csv'), index=False)
    print('written')
