"""
GoCo re-analysis: reproduction gate + missing controls.

Reproduces the BIB-package GoDag / GoCo / Boger rows (alpha=.10, delta in {.10,.50}) from raw inputs and adds:
  GoDag-dense            global path over every distinct score value
  GoCo-random{k}         same blocks / q-grid, ordering by salted hash only (k salts)
  GoCo-score             ordering by mean score of the added calls (ties -> hash)
  GoCo-oracle            ordering by the TRUE incremental loss (uses labels; upper bound only)
  GoCo-exp0 / GoCo-exp1  utility exponent sensitivity (u = d_hat / c^e, e in {0, 1}; paper uses e = 1/2)
  GoCo-dense             GoCo along the dense (tie-level) path with the external block rules
for alpha in {.05,.10,.20} x delta in {.10,.50}, all four datasets, 100 common splits, with cost metrics.

Split rule (all methods): perm = default_rng(20260905+seed).permutation(n); T = perm[:round(.1n)];
C = perm[round(.1n):round(.7n)]; E = perm[round(.7n):]; global methods calibrate on perm[:round(.7n)].
"""
import numpy as np, pandas as pd, scipy.sparse as sp, hashlib, math, time, sys, os, argparse
from scipy.stats import norm

SEED0 = 20260905
HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(HERE, '..', 'w') + os.sep
X = os.path.join(HERE, '..', 'ext') + os.sep
OUT = os.path.join(HERE, 'out') + os.sep
os.makedirs(OUT, exist_ok=True)
QGRID_W = np.array([.01, .02, .03, .04, .05, .075, .10, .125, .15, .175, .20, .225, .25, .275, .30, .35, .40, .50, .60, .75, 1.])
QGRID_E = QGRID_W[:-1]
GRID_W = np.array(sorted(list(range(600, 1001, 25)) + [1050, 1100, 1200, 1400, 1600, 2000, 2500, 3000], reverse=True), float)
GRID_E = np.array([3000, 2500, 2000, 1600, 1400, 1200, 1100, 1050, 1000] + list(range(975, 99, -25)), float)


def hashu(units, salt=''):
    return np.array([int(hashlib.sha256((salt + str(g)).encode()).hexdigest()[:12], 16) / float(16 ** 12) for g in units])


def load_truth():
    td = pd.read_csv(W + 'go_truth_direct.csv.gz', usecols=['gene', 'go_id'], dtype=str)
    tt = pd.read_csv(W + 'go_truth_true_path.csv.gz', usecols=['gene', 'go_id'], dtype=str)
    return set(zip(td.gene, td.go_id)), set(zip(tt.gene, tt.go_id))


class Dataset:
    def __init__(self, name, direct_set, tp_set):
        self.name = name
        units_w = pd.read_csv(W + 'network_communities.csv')['gene'].astype(str).to_numpy()
        si = {g: i for i, g in enumerate(units_w)}  # source universe = Wainberg units (as in the original code)
        if name == 'Wainberg':
            self.units = units_w
            rows = pd.read_csv(W + 'wainberg_highscore_source_attribution.csv.gz', usecols=['gene', 'go_id', 'score', 'support_sources'], dtype={'gene': str, 'go_id': str})
            self.grid = GRID_W
            self.dense_floor = 600.0
            # reproduction check of the attribution file against the raw score file (rows >= 600 within units)
            parts = []
            for ch in pd.read_csv(W + 'wainberg_gene_go_scores.csv.gz', usecols=['gene', 'go_id', 'score'], chunksize=1000000, dtype={'gene': str, 'go_id': str}):
                q = ch[(ch.score >= 600) & ch.gene.isin(si)]
                if len(q): parts.append(q)
            raw = pd.concat(parts, ignore_index=True)
            a = set(zip(rows.gene, rows.go_id, rows.score.round(6)))
            b = set(zip(raw.gene, raw.go_id, raw.score.round(6)))
            print(f'[{name}] attribution rows {len(rows)} raw rows>=600 within units {len(raw)} identical sets: {a == b}', flush=True)
        else:
            z = np.load(X + f'{name.lower()}_external_direct_truepath_matrices_extended100.npz', allow_pickle=True)
            self.units = z['genes'].astype(str)
            self.grid = z['thresholds']; assert np.allclose(np.sort(self.grid), np.sort(GRID_E)), 'frozen external grid differs from GRID_E'.astype(float)
            self.ref = {k: z[k] for k in ['loss_direct', 'loss_truepath', 'yield_direct', 'yield_truepath', 'calls']}
            rows = pd.read_csv(X + f'{name.lower()}_external_go_scores_extended100.csv.gz', usecols=['gene', 'go_id', 'score', 'support_sources'], dtype={'gene': str, 'go_id': str})
            self.dense_floor = 100.0
        self.n = len(self.units)
        ui = {g: i for i, g in enumerate(self.units)}
        rows = rows[rows.gene.isin(ui)].copy()
        rows['i'] = rows.gene.map(ui).astype(int)
        key = list(zip(rows.gene, rows.go_id))
        rows['D'] = np.array([k in direct_set for k in key], dtype=np.int8)
        rows['T'] = np.array([k in tp_set for k in key], dtype=np.int8)
        rows['srcs'] = [[si[s] for s in str(x).split(';') if s in si] for x in rows.support_sources]
        self.rows = rows.reset_index(drop=True)
        self.M = len(units_w)
        self.dense = np.array(sorted(set(np.round(self.rows.score.to_numpy(float), 6)), reverse=True))
        self.dense = self.dense[self.dense >= self.dense_floor]
        self.hh = hashu(self.units)
        print(f'[{name}] n={self.n} rows={len(self.rows)} dense levels={len(self.dense)} grid={len(self.grid)}', flush=True)

    # ---- matrices on an arbitrary descending threshold vector ----
    def mats(self, thr):
        thr = np.asarray(thr, float); K = len(thr)
        s = np.round(self.rows.score.to_numpy(float), 6)
        k0 = np.searchsorted(-thr, -s, side='left')  # first k with thr[k] <= s
        keep = k0 < K
        i = self.rows.i.to_numpy()[keep]; k = k0[keep]
        def acc(w):
            m = np.zeros((self.n, K), np.float64); np.add.at(m, (i, k), w); return np.cumsum(m, axis=1)
        C = acc(np.ones(keep.sum()))
        YT = acc(self.rows['T'].to_numpy(float)[keep]); YD = acc(self.rows['D'].to_numpy(float)[keep])
        LT = np.where(C > 0, 1 - YT / np.maximum(C, 1), 0.0); LD = np.where(C > 0, 1 - YD / np.maximum(C, 1), 0.0)
        f = np.float32
        return dict(thr=thr, C=C.astype(f), YT=YT.astype(f), YD=YD.astype(f), LT=LT.astype(f), LD=LD.astype(f))

    def block_rows(self, lo, hi):
        s = np.round(self.rows.score.to_numpy(float), 6)
        return self.rows[(s >= lo) & (s < hi)]

    def Aadd(self, lo, hi):
        q = self.block_rows(lo, hi)
        counts = np.bincount(q.i.to_numpy(), minlength=self.n)
        rr, cc, dd = [], [], []
        for i, ss in zip(q.i.to_numpy(), q.srcs):
            if not ss: continue
            w = 1.0 / (counts[i] * len(ss)); rr.extend([i] * len(ss)); cc.extend(ss); dd.extend([w] * len(ss))
        A = sp.coo_matrix((dd, (rr, cc)), shape=(self.n, self.M), dtype=np.float64).tocsr(); A.sum_duplicates()
        meanscore = np.full(self.n, np.nan)
        if len(q):
            g = q.groupby('i').score.mean(); meanscore[g.index.to_numpy()] = g.to_numpy()
        return A, meanscore


def clt_p(y, alpha):
    y = np.asarray(y, float); m = len(y); r = float(y.mean()); sd = max(float(y.std(ddof=0)), 1e-6)
    return float(norm.cdf((r - alpha) / (sd / math.sqrt(m)))), r, sd


def source_pred(Aadd, delta, fit, kappa=1.0):
    yy = delta[fit]; glob = float(yy.mean()); Af = Aadd[fit, :]
    den = np.asarray(Af.sum(axis=0)).ravel(); num = np.asarray(Af.T.dot(yy)).ravel()
    pr = (num + kappa * glob) / (den + kappa)
    pred = np.asarray(Aadd.dot(pr)).ravel(); mass = np.asarray(Aadd.sum(axis=1)).ravel(); pred[mass <= 1e-15] = glob
    return pred


def selmask(u, hh, affected, train, q):
    a = train[affected[train]]
    if len(a) == 0: return np.zeros(len(u), bool)
    if q >= 1: return affected.copy()
    oo = np.lexsort((hh[a], u[a])); k = max(1, int(np.floor(q * len(a)))); i = a[oo[k - 1]]; r = u[i]; h = hh[i]
    return affected & ((u < r) | ((np.abs(u - r) <= 1e-14) & (hh <= h)))


def utility(order, ds, Aadd, meanscore, delta_true, cadd, affected, train, hh_default):
    """Return (u, hh) defining the admission order for a block. Lower u admitted first; hh breaks ties."""
    kap = 1.0
    if order.startswith('knapk'):
        kap = float(order[5:]); order = 'knap'
    hh = hh_default
    if order in ('knapEB', 'learned', 'learnedgb', 'ens'):
        from goco_learned import learned_utility
        return learned_utility(order, ds, Aadd, meanscore, delta_true, cadd, affected, train, hh_default, utility.extra)
    if order.startswith('random'):
        return np.zeros(ds.n), hashu(ds.units, salt=order)
    if order == 'score':
        return -np.nan_to_num(meanscore, nan=-1e9), hh
    fit = train[affected[train]]
    if order == 'oracle':
        # label-using reference ordering by the realised increment of every gene (not T-measurable; Lemma 1b optimum for risk)
        return delta_true.copy(), hh
    pred = source_pred(Aadd, delta_true, fit, kappa=kap) if len(fit) else np.full(ds.n, float(delta_true[affected].mean()) if affected.any() else 0.0)
    if order == 'source': return pred / np.sqrt(np.maximum(cadd, 1.0)), hh
    if order == 'exp0': return pred, hh
    if order == 'exp1': return pred / np.maximum(cadd, 1.0), hh
    if order == 'knap':
        # knapsack-style: predicted loss increment per predicted supported added call.
        # f_i = unsupported fraction of the added calls (nuisance fold), source-smoothed like Delta.
        fsup = utility.extra['fsup']  # supported count among added calls, all genes
        f_unsup = np.where(cadd > 0, 1.0 - fsup / np.maximum(cadd, 1.0), 0.0)
        fhat = source_pred(Aadd, f_unsup, fit, kappa=kap) if len(fit) else np.full(ds.n, float(f_unsup[affected].mean()) if affected.any() else 0.0)
        exp_supported = np.maximum(cadd, 1.0) * np.clip(1.0 - fhat, 0.0, 1.0)
        return pred / np.maximum(exp_supported, 1e-3), hh
    raise ValueError(order)


def run_path(ds, m, cert, ev, train, alpha, delta, loss, adaptive_order, qgrid, band_rule, block_fn, adapt_cache, adapt_blocks=None):
    """Fixed-sequence LTT along thresholds m['thr'] with optional adaptive insertions.
    block_fn(j) -> (Aadd, meanscore) for block thr[j] -> thr[j+1]; adapt_cache memoizes per block.
    adapt_blocks: optional set of block indices j at which adaptive candidates may be inserted (None = any).
    Returns the selected candidate tuple."""
    L = m['LT'] if loss == 'TruePath' else m['LD']
    thr = m['thr']; K = len(thr)
    # vectorised global p-values on the certification fold
    Lc = L[cert, :].astype(np.float64); means = Lc.mean(axis=0); sds = np.maximum(Lc.std(axis=0), 1e-6)
    pg = norm.cdf((means - alpha) / (sds / math.sqrt(len(cert))))
    last = None
    for j in range(K):
        if pg[j] > delta: break
        last = ('global', j, None, float(pg[j]), float(means[j]), float(sds[j]))
        if adaptive_order is None or j == K - 1: continue
        if adapt_blocks is not None and j not in adapt_blocks: continue
        if band_rule and float(L[train, j + 1].mean()) <= 0.8 * alpha: continue
        affected = m['C'][:, j + 1] > m['C'][:, j]
        if int(affected[train].sum()) < 5: continue
        if j not in adapt_cache: adapt_cache[j] = block_fn(j)
        Aadd, meanscore = adapt_cache[j]
        delta_true = (L[:, j + 1] - L[:, j]).astype(np.float64); cadd = (m['C'][:, j + 1] - m['C'][:, j]).astype(np.float64)
        utility.extra = {'fsup': (m['YT'][:, j + 1] - m['YT'][:, j]).astype(np.float64), 'lo': float(thr[j + 1]), 'hi': float(thr[j])}
        u, hh = utility(adaptive_order, ds, Aadd, meanscore, delta_true, cadd, affected, train, ds.hh)
        stop = False
        for q in qgrid:
            if q >= 1: continue  # q=1 is the next global candidate
            sel = selmask(u, hh, affected, train, float(q))
            y = np.where(sel, L[:, j + 1], L[:, j]).astype(np.float64)
            p, r, sd = clt_p(y[cert], alpha)
            if p > delta: stop = True; break
            last = ('adapt', j, sel, p, r, sd, float(q))
        if stop: break
    return last


def metrics(ds, m, last, ev):
    if last is None: return dict(policy='none', p_selected=float('nan'), cal_mean=0.0, cal_sd=0.0, unit_fdp=0.0, unit_yield=0, go_yield=0.0, total_calls=0.0, genes_with_calls=0, term_fdp=0.0, unit_fdp_direct=0.0, go_yield_direct=0.0, n_relaxed_eval=0)
    kind, j = last[0], last[1]
    if kind == 'global':
        cols = {k: m[k][:, j].astype(np.float64) for k in ['LT', 'LD', 'YT', 'YD', 'C']}; pol = f'global_{m["thr"][j]:g}'
    else:
        sel = last[2]
        cols = {k: np.where(sel, m[k][:, j + 1], m[k][:, j]).astype(np.float64) for k in ['LT', 'LD', 'YT', 'YD', 'C']}
        pol = f'adapt_{m["thr"][j]:g}_to_{m["thr"][j + 1]:g}_q{last[6]:g}'
    calls = float(cols['C'][ev].sum()); yt = float(cols['YT'][ev].sum()); yd = float(cols['YD'][ev].sum())
    return dict(policy=pol, p_selected=last[3], cal_mean=last[4], cal_sd=last[5],
                unit_fdp=float(cols['LT'][ev].mean()), unit_yield=int((cols['YT'][ev] > 0).sum()), go_yield=yt,
                total_calls=calls, genes_with_calls=int((cols['C'][ev] > 0).sum()), term_fdp=float((calls - yt) / max(calls, 1)),
                unit_fdp_direct=float(cols['LD'][ev].mean()), go_yield_direct=yd,
                n_relaxed_eval=int(sel[ev].sum()) if kind == 'adapt' else 0)


def run_dataset(ds, alphas, deltas, seeds, methods, tag):
    t0 = time.time()
    mg = ds.mats(ds.grid)
    if ds.name != 'Wainberg':  # reproduction gate for external matrices
        for a, b in [('LT', 'loss_truepath'), ('LD', 'loss_direct'), ('C', 'calls'), ('YT', 'yield_truepath')]:
            d = np.abs(mg[a] - ds.ref[b].astype(float)).max()
            print(f'[{ds.name}] matrix check {a}: max abs diff vs frozen npz = {d:.2e}', flush=True)
    md = ds.mats(ds.dense)
    # Wainberg original GoCo block: global grid down to 800 then partial admission of 800 -> 789.65
    if ds.name == 'Wainberg':
        j800 = int(np.where(ds.grid == 800)[0][0])
        thr_w = np.concatenate([ds.grid[:j800 + 1], [789.65]])
        mw = ds.mats(thr_w)
    rows = []
    cache_w, cache_g, cache_d = {}, {}, {}
    for seed in seeds:
        perm = np.random.default_rng(SEED0 + seed).permutation(ds.n)
        nt = round(.10 * ds.n); nc = round(.70 * ds.n)
        train, cert, ev, cal70 = perm[:nt], perm[nt:nc], perm[nc:], perm[:nc]
        for alpha in alphas:
            for delta in deltas:
                def add(method, last, m):
                    r = dict(dataset=ds.name, seed=seed, alpha=alpha, delta=delta, method=method); r.update(metrics(ds, m, last, ev)); rows.append(r)
                if 'Boger' in methods:
                    add('Boger (Direct, global grid)', run_path(ds, mg, cal70, ev, train, alpha, delta, 'Direct', None, None, False, None, {}), mg)
                if 'GoDag' in methods:
                    add('GoDag (global grid)', run_path(ds, mg, cal70, ev, train, alpha, delta, 'TruePath', None, None, False, None, {}), mg)
                if 'GoDag-dense' in methods:
                    add('GoDag-dense (global, all distinct scores)', run_path(ds, md, cal70, ev, train, alpha, delta, 'TruePath', None, None, False, None, {}), md)
                if 'GoDag-cert60' in methods:
                    add('GoDag (global grid, certification fold only)', run_path(ds, mg, cert, ev, train, alpha, delta, 'TruePath', None, None, False, None, {}), mg)
                orders = [o for o in methods if o.startswith('GoCo:')]
                if orders:
                    if ds.name == 'Wainberg':
                        cache = cache_w
                        bf = lambda j: ds.Aadd(789.65, 800.0)
                        for o in orders:
                            order = o.split(':')[1]
                            # original Wainberg construction: adaptive block only after the last grid threshold (800); q grid includes 1
                            last = run_path(ds, mw, cert, ev, train, alpha, delta, 'TruePath', order, QGRID_W, False, bf, cache, adapt_blocks={len(thr_w) - 2})
                            add(f'GoCo-{order} (paper path)', last, mw)
                    else:
                        cache = cache_g
                        bf = lambda j: ds.Aadd(mg['thr'][j + 1], mg['thr'][j])
                        for o in orders:
                            order = o.split(':')[1]
                            last = run_path(ds, mg, cert, ev, train, alpha, delta, 'TruePath', order, QGRID_E, True, bf, cache)
                            add(f'GoCo-{order} (paper path)', last, mg)
                gorders = [o for o in methods if o.startswith('GoCoGrid:')]
                if gorders:
                    # uniform definition: the same 25-unit grid as GoDag, partial admission inserted in every block that
                    # satisfies the nuisance band rule (mean loss at the liberal end > 0.8*alpha, >=5 affected nuisance genes)
                    cache = cache_g
                    bf = lambda j: ds.Aadd(mg['thr'][j + 1], mg['thr'][j])
                    for o in gorders:
                        order = o.split(':')[1]
                        last = run_path(ds, mg, cert, ev, train, alpha, delta, 'TruePath', order, QGRID_E, True, bf, cache)
                        add(f'GoCoGrid-{order} (grid path, all blocks)', last, mg)
                dorders = [o for o in methods if o.startswith('GoCoDense:')]
                if dorders:
                    cache = cache_d
                    bf = lambda j: ds.Aadd(md['thr'][j + 1], md['thr'][j])
                    for o in dorders:
                        order = o.split(':')[1]
                        last = run_path(ds, md, cert, ev, train, alpha, delta, 'TruePath', order, QGRID_E, True, bf, cache)
                        add(f'GoCo-dense-{order} (tie-level path)', last, md)
        if (seed + 1) % 10 == 0: print(f'[{ds.name}] seed {seed + 1}/{len(seeds)} elapsed {time.time() - t0:.0f}s', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT + f'{ds.name}_{tag}_100splits.csv', index=False)
    return df


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--datasets', default='Wainberg,Sanger,DRIVE,HAP1')
    ap.add_argument('--alphas', default='0.05,0.10,0.20')
    ap.add_argument('--deltas', default='0.10,0.50')
    ap.add_argument('--seeds', type=int, default=100)
    ap.add_argument('--methods', default='Boger,GoDag,GoDag-dense,GoDag-cert60,GoCo:source,GoCo:random1,GoCo:random2,GoCo:random3,GoCo:score,GoCo:oracle,GoCo:exp0,GoCo:exp1,GoCoDense:source,GoCoDense:random1')
    ap.add_argument('--tag', default='full')
    a = ap.parse_args()
    dset, tset = load_truth()
    alphas = [float(x) for x in a.alphas.split(',')]; deltas = [float(x) for x in a.deltas.split(',')]
    methods = a.methods.split(',')
    for name in a.datasets.split(','):
        ds = Dataset(name, dset, tset)
        df = run_dataset(ds, alphas, deltas, range(a.seeds), methods, a.tag)
        df['exceed'] = (df.unit_fdp > df.alpha).astype(float)
        s = df.groupby(['alpha', 'delta', 'method']).agg(unit_fdp=('unit_fdp', 'mean'), exceed=('exceed', 'mean'), unit_yield=('unit_yield', 'mean'), go_yield=('go_yield', 'mean'), calls=('total_calls', 'mean'), term_fdp=('term_fdp', 'mean'), genes_with_calls=('genes_with_calls', 'mean'))
        pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
        print(s.to_string(), flush=True)
