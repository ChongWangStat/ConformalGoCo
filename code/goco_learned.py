"""Learned, T-measurable ordering for GoCo partial admission (GoCo-L).

Every quantity used to rank a pool gene is a function of (i) frozen, label-free predictor output and (ii) labels of
the ranking fold T only. Labels of certification or evaluation genes never enter, so Proposition 1 / Corollary 1 apply.
The definitive operational check is test_measurability.py (scrambling all non-ranking-fold labels leaves the pool
ordering bit-identical).

Per-call features (call = gene i, term g, score >= lowest grid threshold):
  f1  log score s_ig                                   (frozen)
  f2  log(1 + number of supporting sources of the call) (frozen)
  f3  log(1 + number of the gene's calls above the floor) (frozen)
  f4  shrunk unsupported rate of term g among ranking-fold calls          (T labels; cross-fitted)
  f5  mean shrunk unsupported rate of the call's sources, ranking-fold calls (T labels; cross-fitted)
  f6  log frequency of g among TruePath annotations of ranking-fold genes (T labels; leave-one-gene-out for training rows)
  f7  module size excluding the target (Wainberg) / target-module correlation (Sanger, DRIVE) / cluster z (HAP1) (frozen, label-free)
NOT used: background_K_ex_target (equals the term's Direct count minus the target's own label -> would encode C/E labels),
n_support_sources (duplicate of f2), calls above the block's upper threshold and calls added by the block (constant on
training rows, hence uninformative).

Model: logistic regression (C = 0.5, standardized features) or histogram gradient boosting (max_iter 150, depth 3,
learning rate 0.06, min_samples_leaf 20, l2 = 1.0) for the indicator 'call unsupported' (TruePath). The classifier and
the rate features are cross-fitted between two hash-defined halves of the ranking fold: a ranking-fold call is scored by
the model fitted on the other half; pool calls are scored by the model fitted on the whole ranking fold.
Delta_hat_i = L_hat(S_i^L) - L_hat(S_i^H) with L_hat(S) = mean of p_hat over S (0 for empty S); expected supported added
calls = sum over D_i of (1 - p_hat); u_i^L = Delta_hat_i / max(expected supported, 1e-3), exactly as in Eq. (6).
Hyperparameters and the model class were fixed on 10 Wainberg development splits (seeds 0-9) before the external runs.
"""
import numpy as np, pandas as pd, hashlib
from collections import Counter
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
import goco_rerun as gr

_TRUTH = {}
_EXTRA = {}
_MODEL_CACHE = {}
RATE_KAPPA = 2.0


def _truth_by_gene():
    if 'g2t' not in _TRUTH:
        tt = pd.read_csv(gr.W + 'go_truth_true_path.csv.gz', usecols=['gene', 'go_id'], dtype=str)
        g2t = {}
        for g, t in zip(tt.gene, tt.go_id): g2t.setdefault(g, set()).add(t)
        _TRUTH['g2t'] = g2t
    return _TRUTH['g2t']


def _extra_col(ds):
    """Frozen, label-free per-call extra feature f7 aligned with ds.rows."""
    if ds.name in _EXTRA: return _EXTRA[ds.name]
    if ds.name == 'Wainberg':
        a = pd.read_csv(gr.W + 'wainberg_highscore_source_attribution.csv.gz', usecols=['gene', 'go_id', 'module_n_ex_target'], dtype={'gene': str, 'go_id': str}).rename(columns={'module_n_ex_target': 'x1'})
    else:
        f = gr.X + f'{ds.name.lower()}_external_go_scores_extended100.csv.gz'
        cols = pd.read_csv(f, nrows=1).columns.tolist()
        c1 = 'external_target_module_corr' if 'external_target_module_corr' in cols else 'cluster_mean_z'
        a = pd.read_csv(f, usecols=['gene', 'go_id', c1], dtype={'gene': str, 'go_id': str}).rename(columns={c1: 'x1'})
    a = a.drop_duplicates(['gene', 'go_id'])
    m = ds.rows[['gene', 'go_id']].merge(a, on=['gene', 'go_id'], how='left')
    assert len(m) == len(ds.rows)
    x1 = m.x1.to_numpy(float); x1 = np.where(np.isfinite(x1), x1, np.nanmedian(x1))
    _EXTRA[ds.name] = x1
    return x1


def _shrunk_rates(keys, y, kappa, glob):
    s = {}; n = {}
    for k, v in zip(keys, y):
        s[k] = s.get(k, 0.0) + v; n[k] = n.get(k, 0) + 1
    return lambda kl: glob if not kl else float(np.mean([(s.get(k, 0.0) + kappa * glob) / (n.get(k, 0) + kappa) for k in kl]))


class _Fold:
    """Rate tables and classifier fitted on a subset of ranking-fold rows."""
    def __init__(self, rows_fit, y, terms, srcs, glob, kind, X_fn):
        self.term_rate = _shrunk_rates(terms[rows_fit], y[rows_fit], RATE_KAPPA, glob)
        sk = [r for rr in rows_fit for r in srcs[rr]]; sy = [y[rr] for rr in rows_fit for _ in srcs[rr]]
        self.src_rate = _shrunk_rates(sk, sy, RATE_KAPPA, glob)
        self.kind = kind; self.X_fn = X_fn; self.clf = None; self.scaler = None; self.rows_fit = rows_fit; self.y = y

    def features(self, idx, base):
        f4 = np.array([self.term_rate([t]) for t in base['terms'][idx]])
        f5 = np.array([self.src_rate(base['srcs'][r]) for r in idx])
        return np.column_stack([base['f1'][idx], base['f2'][idx], base['f3'][idx], f4, f5, base['f6'](idx), base['x1'][idx]])

    def fit(self, base):
        X = self.features(self.rows_fit, base); yy = self.y[self.rows_fit]
        if self.kind == 'learnedgb':
            self.clf = HistGradientBoostingClassifier(max_iter=150, max_depth=3, learning_rate=0.06, min_samples_leaf=20, l2_regularization=1.0, random_state=0).fit(X, yy)
        else:
            self.scaler = StandardScaler().fit(X)
            self.clf = LogisticRegression(C=0.5, max_iter=1000).fit(self.scaler.transform(X), yy)
        return self

    def predict(self, idx, base):
        X = self.features(idx, base)
        if self.scaler is not None: X = self.scaler.transform(X)
        return self.clf.predict_proba(X)[:, 1]


def _fit_model(ds, train, kind):
    key = (ds.name, kind, hashlib.sha256(np.sort(train).astype(np.int64).tobytes()).hexdigest()[:16])
    if key in _MODEL_CACHE: return _MODEL_CACHE[key]
    rows = ds.rows; s = np.round(rows.score.to_numpy(float), 6); T = rows['T'].to_numpy(float); i_all = rows.i.to_numpy()
    srcs = rows.srcs.tolist(); terms = rows.go_id.to_numpy(); genes = rows.gene.to_numpy()
    train_set = np.zeros(ds.n, bool); train_set[train] = True
    half = np.zeros(ds.n, np.int8) - 1
    half[train] = (gr.hashu(ds.units[train], salt='xfit') < 0.5).astype(np.int8)
    keep = s >= ds.dense_floor
    tr_rows = np.where(keep & train_set[i_all])[0]
    y = 1.0 - T
    glob = float(y[tr_rows].mean())
    calls_total = np.bincount(i_all[keep], minlength=ds.n).astype(float)
    # term frequency among ranking-fold genes' TruePath truth; leave-one-gene-out for ranking-fold rows
    g2t = _truth_by_gene(); cnt = Counter()
    for g in ds.units[train]:
        for t in g2t.get(g, ()): cnt[t] += 1
    ntr = float(len(train))
    def f6(idx):
        out = np.empty(len(idx))
        for q, r in enumerate(idx):
            c = cnt.get(terms[r], 0)
            if train_set[i_all[r]] and terms[r] in g2t.get(genes[r], ()): c -= 1
            out[q] = np.log((c + 0.5) / (ntr + 1.0))
        return out
    base = dict(f1=np.log(np.maximum(s, 1e-9)), f2=np.log1p(np.array([len(x) for x in srcs], float)), f3=np.log1p(calls_total[i_all]), f6=f6, x1=_extra_col(ds), terms=terms, srcs=srcs)
    folds = {}
    folds[None] = _Fold(tr_rows, y, terms, srcs, glob, kind, None).fit(base)                       # whole ranking fold -> pool genes
    for h in (0, 1):
        sub = tr_rows[half[i_all[tr_rows]] != h]                                                # fitted without half h -> scores half h
        folds[h] = _Fold(sub, y, terms, srcs, glob, kind, None).fit(base)

    def predict(idx):
        out = np.empty(len(idx)); hidx = np.where(train_set[i_all[idx]], half[i_all[idx]], -1)
        for h in (-1, 0, 1):
            sel = np.where(hidx == h)[0]
            if len(sel): out[sel] = folds[None if h == -1 else h].predict(idx[sel], base)
        return out
    _MODEL_CACHE[key] = (predict, s, i_all)
    if len(_MODEL_CACHE) > 8:
        for k in list(_MODEL_CACHE)[:-4]: del _MODEL_CACHE[k]
    return _MODEL_CACHE[key]


def learned_utility(order, ds, Aadd, meanscore, delta_true, cadd, affected, train, hh, extra):
    lo, hi = extra['lo'], extra['hi']
    kind = 'learnedgb' if order == 'learnedgb' else 'learned'
    predict, s, i_all = _fit_model(ds, train, kind)
    aff_rows = np.where(affected[i_all] & (s >= lo))[0]          # calls of affected genes in S^L
    i = i_all[aff_rows]
    p_unsup = predict(aff_rows)
    in_H = s[aff_rows] >= hi
    sumL = np.bincount(i, weights=p_unsup, minlength=ds.n); nL = np.bincount(i, minlength=ds.n).astype(float)
    sumH = np.bincount(i[in_H], weights=p_unsup[in_H], minlength=ds.n); nH = np.bincount(i[in_H], minlength=ds.n).astype(float)
    LhatL = np.where(nL > 0, sumL / np.maximum(nL, 1), 0.0); LhatH = np.where(nH > 0, sumH / np.maximum(nH, 1), 0.0)
    dhat = LhatL - LhatH
    exp_sup = np.bincount(i[~in_H], weights=1.0 - p_unsup[~in_H], minlength=ds.n)
    u = np.where(affected, dhat / np.maximum(exp_sup, 1e-3), 0.0)
    if order in ('learned', 'learnedgb'): return u, hh
    raise ValueError(order)
