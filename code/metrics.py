"""metrics.py — contact-prediction metrics for the HySTraG-IP benchmark.

Per-protein scoring (precision@L/k, AUPRC) with helpers for aggregation
(bootstrap CI) and a paired significance test. Import and use from the notebook
or any evaluation script.
"""
import numpy as np

try:
    from scipy.stats import wilcoxon
except Exception:
    wilcoxon = None
try:
    from sklearn.metrics import average_precision_score
except Exception:
    average_precision_score = None


def sep_band(L, lo, hi=None):
    """Boolean (L,L) mask selecting residue pairs in a sequence-separation band."""
    s = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :])
    return (s >= lo) if hi is None else ((s >= lo) & (s <= hi))


def _score_mask(L, pair_valid, min_sep, extra=None):
    m = pair_valid & sep_band(L, min_sep) & np.triu(np.ones((L, L), bool), 1)
    if extra is not None:
        m = m & extra
    return m


def precision_at_Lk(prob, contact, pair_valid, L_eff, k, min_sep=24, extra=None):
    """Precision among the top floor(L_eff/k) ranked pairs in a separation band.
    prob, contact, pair_valid: (L,L). k in {1,2,5} -> top L, L/2, L/5."""
    L = prob.shape[0]
    m = _score_mask(L, pair_valid, min_sep, extra)
    ii, jj = np.where(m)
    if ii.size == 0:
        return np.nan
    top = np.argsort(-prob[ii, jj])[:max(1, int(round(L_eff / k)))]
    return float(contact[ii[top], jj[top]].mean())


def auprc(prob, contact, pair_valid, min_sep=24, extra=None):
    """Area under the precision-recall curve over scorable pairs in a band."""
    if average_precision_score is None:
        return np.nan
    L = prob.shape[0]
    m = _score_mask(L, pair_valid, min_sep, extra)
    y = contact[m].astype(int)
    p = prob[m]
    return float(average_precision_score(y, p)) if y.sum() > 0 else np.nan


def summarize(arr, n_boot=2000, seed=0):
    """mean, sd, and bootstrap 95% CI over a per-protein metric array (NaNs dropped)."""
    arr = np.asarray(arr, float)
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return dict(mean=np.nan, sd=np.nan, ci=(np.nan, np.nan), n=0)
    rng = np.random.default_rng(seed)
    boot = [rng.choice(arr, arr.size, replace=True).mean() for _ in range(n_boot)]
    return dict(mean=float(arr.mean()),
                sd=float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
                ci=(float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))),
                n=int(arr.size))


def paired_test(model_arr, base_arr):
    """Paired Wilcoxon signed-rank of model vs. baseline across proteins."""
    a = np.asarray(model_arr, float)
    b = np.asarray(base_arr, float)
    m = ~np.isnan(a) & ~np.isnan(b)
    if wilcoxon is None or m.sum() < 2:
        md = float(np.median((a - b)[m])) if m.any() else np.nan
        return dict(p_value=np.nan, median_delta=md)
    _stat, p = wilcoxon(a[m], b[m])
    return dict(p_value=float(p), median_delta=float(np.median(a[m] - b[m])))


def fmt(s):
    """Pretty 'mean +/- sd' from a summarize() dict."""
    if s["n"] > 1:
        return f"{s['mean']:.3f} +/- {s['sd']:.3f}"
    return f"{s['mean']:.3f}" if s["n"] else "nan"
