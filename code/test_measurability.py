"""Operational T-measurability test for the GoCo ordering.
For one split, scramble the TruePath/Direct labels of every non-ranking-fold gene (and the truth dictionary of pool
genes) and check that the ordering statistic u_i of affected POOL genes is unchanged.
Usage: python test_measurability.py [dataset] [seed]"""
import sys, numpy as np, goco_rerun as gr, goco_learned as gl
name = sys.argv[1] if len(sys.argv) > 1 else 'Wainberg'; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
dset, tset = gr.load_truth(); ds = gr.Dataset(name, dset, tset)
mg = ds.mats(ds.grid); L = mg['LT']; C = mg['C']; YT = mg['YT']
perm = np.random.default_rng(gr.SEED0 + seed).permutation(ds.n); nt = round(.1 * ds.n); train = perm[:nt]; pool = perm[nt:]
intrain = np.zeros(ds.n, bool); intrain[train] = True
blocks = [(int(np.where(np.isclose(ds.grid, 800))[0][0]),)] if name == 'Wainberg' else [(j,) for j in range(len(ds.grid) - 1) if float(L[train, j + 1].mean()) > 0.08 and (C[train, j + 1] > C[train, j]).sum() >= 5]
orders = ['goco']
def compute():
    out = {}
    for (j,) in blocks:
        lo, hi = float(ds.grid[j + 1]), float(ds.grid[j]); Aadd, meanscore = ds.Aadd(lo, hi)
        delta = (L[:, j + 1] - L[:, j]).astype(float); cadd = (C[:, j + 1] - C[:, j]).astype(float); affected = cadd > 0
        gr.utility.extra = {'fsup': (YT[:, j + 1] - YT[:, j]).astype(float), 'lo': lo, 'hi': hi}
        ap = pool[affected[pool]]
        for o in orders:
            u, hh = gr.utility(o, ds, Aadd, meanscore, delta, cadd, affected, train, ds.hh); out[(j, o)] = (u[ap].copy(), hh[ap].copy())
    return out
gl._MODEL_CACHE.clear(); before = compute()
rows = ds.rows; m = ~intrain[rows.i.to_numpy()]; rng = np.random.default_rng(12345)
for col in ['T', 'D']:
    v = rows[col].to_numpy().copy(); v[m] = rng.permutation(v[m]); rows[col] = v
mg2 = ds.mats(ds.grid); L[:] = mg2['LT']; C[:] = mg2['C']; YT[:] = mg2['YT']   # losses of pool genes now differ
gl._MODEL_CACHE.clear(); gl._TRUTH.clear(); g2t = gl._truth_by_gene()
poolgenes = set(ds.units[pool]); keys = [g for g in g2t if g in poolgenes]; vals = [g2t[g] for g in keys]; rng.shuffle(vals)
for g, v in zip(keys, vals): g2t[g] = v
after = compute()
ok = True
for k in before:
    same = np.array_equal(before[k][0], after[k][0]) and np.array_equal(before[k][1], after[k][1])
    print(f'block {k[0]} {k[1]:10s} identical under scrambling: {same} (max |diff| = {np.max(np.abs(before[k][0] - after[k][0])):.2e})'); ok &= same
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
