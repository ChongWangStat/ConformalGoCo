"""Audit of the exact tie rule (revision of 2026-09-25).

Runs the full module-family analysis (Boger, GoDag, GoDag on the 60% fold, GoCo) through goco_rerun.run_dataset
while shadowing every selmask() call with the pre-revision rule (|u - r| <= 1e-14 ties) and counts the block
candidates on which the two rules admit different sets.  It also compares every certified GoCo policy with the
frozen pre-revision split-level file.  Outputs: out/<dataset>_run_100splits.csv (identical to `goco_rerun.py
--tag run`) and out/tie_audit_<dataset>.txt.

Usage: python code/tie_audit.py Wainberg,Sanger,DRIVE,HAP1 [seeds=100]
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goco_rerun as gr


def selmask_old(u, hh, affected, train, q):
    a = train[affected[train]]
    if len(a) == 0: return np.zeros(len(u), bool)
    if q >= 1: return affected.copy()
    oo = np.lexsort((hh[a], u[a])); k = max(1, int(np.floor(q * len(a)))); i = a[oo[k - 1]]; r = u[i]; h = hh[i]
    return affected & ((u < r) | ((np.abs(u - r) <= 1e-14) & (hh <= h)))


COUNT = dict(calls=0, diff=0, ties=0)
_new = gr.selmask


def shadow(u, hh, affected, train, q):
    a = _new(u, hh, affected, train, q); b = selmask_old(u, hh, affected, train, q)
    COUNT['calls'] += 1; COUNT['diff'] += int(not np.array_equal(a, b))
    ua = u[train[affected[train]]]; COUNT['ties'] += len(ua) - len(np.unique(ua))
    return a


gr.selmask = shadow

if __name__ == '__main__':
    names = sys.argv[1].split(','); seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    dset, tset = gr.load_truth()
    frozen = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'goco_results_100splits_prereview.csv')
    old = pd.read_csv(frozen) if os.path.exists(frozen) else None
    for name in names:
        for k in COUNT: COUNT[k] = 0
        ds = gr.Dataset(name, dset, tset)
        df = gr.run_dataset(ds, [0.05, 0.10, 0.20], [0.10, 0.50], range(seeds), ['Boger', 'GoDag', 'GoDag-cert60', 'GoCoGrid:goco'], 'run')
        lines = [f'{name}: selmask evaluations {COUNT["calls"]}, admitted sets differing between the exact and the 1e-14 rule {COUNT["diff"]}, '
                 f'exact u-ties among affected ranking-fold units {COUNT["ties"]}']
        if old is not None:
            o = old[(old.dataset == name) & (old.seed < seeds)]
            for m in ['Boger', 'GoDag', 'GoCo']:
                mm = m if m != 'Boger' else 'Boger'
                a = df[df.method.str.startswith(mm)]
                a = a[~a.method.str.contains('certification fold')]
                x = a.merge(o[o.method == m][['seed', 'alpha', 'delta', 'policy', 'unit_fdp', 'go_yield']], on=['seed', 'alpha', 'delta'], suffixes=('', '_old'))
                lines.append(f'  {m}: rows {len(x)}, identical policy {(x.policy == x.policy_old).mean():.4f}, max |dFDP| {(x.unit_fdp - x.unit_fdp_old).abs().max():.2e}, max |dYield| {(x.go_yield - x.go_yield_old).abs().max():.0f}')
        txt = '\n'.join(lines); print(txt, flush=True)
        open(gr.OUT + f'tie_audit_{name}.txt', 'w').write(txt + '\n')
