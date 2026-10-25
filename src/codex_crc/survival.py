# Cox, permutation p, cut-point correction
import numpy as np
from scipy.stats import chi2, norm


def _risk_sets(time, event):
    order = np.argsort(time, kind="stable")
    return order, time[order], event[order]


def cox1(x, time, event, iters=50):
    """one-covariate Cox PH, Breslow ties"""
    x, time, event = map(np.asarray, (x, time, event))
    x = x.astype(float)
    order, t, e = _risk_sets(time, event)
    xs = x[order]

    def loglik(b):
        eta = b * xs
        ll, g, h = 0.0, 0.0, 0.0
        for i in np.flatnonzero(e == 1):
            r = t >= t[i]
            w = np.exp(eta[r] - eta[r].max())
            s0, s1, s2 = w.sum(), (w * xs[r]).sum(), (w * xs[r] ** 2).sum()
            ll += eta[i] - (np.log(s0) + eta[r].max())
            g += xs[i] - s1 / s0
            h += s2 / s0 - (s1 / s0) ** 2
        return ll, g, h

    b = 0.0
    ll0 = loglik(0.0)[0]
    for _ in range(iters):
        ll, g, h = loglik(b)
        if h <= 0:
            break
        step = g / h
        b += step
        if abs(step) < 1e-10:
            break
    ll, g, h = loglik(b)
    se = 1 / np.sqrt(h) if h > 0 else np.inf
    return {"beta": b, "hr": np.exp(b), "se": se, "p_wald": 2 * norm.sf(abs(b / se)),
            "p_lr": chi2.sf(2 * (ll - ll0), 1)}


def score_z(X, time, event):
    """Cox score statistics at beta = 0, all columns of X"""
    X = np.asarray(X, float)
    ev = np.flatnonzero(np.asarray(event) == 1)
    R = (time[None, :] >= time[ev][:, None]).astype(float)  # risk set of each event
    s0 = R.sum(1)[:, None]
    m1, m2 = R @ X / s0, R @ X ** 2 / s0
    u = (X[ev] - m1).sum(0)
    info = (m2 - m1 ** 2).sum(0)
    return u / np.sqrt(info)


def westfall_young(X, time, event, n_perm=10000, seed=0):
    """max-T adjusted p, rows of X permuted"""
    X, time, event = np.asarray(X, float), np.asarray(time, float), np.asarray(event)
    z = np.abs(score_z(X, time, event))
    rng = np.random.default_rng(seed)
    maxz = np.array([np.abs(score_z(X[rng.permutation(len(X))], time, event)).max() for _ in range(n_perm)])
    return {"z": z, "p_raw": 2 * norm.sf(z), "p_adj": np.array([(1 + (maxz >= zj).sum()) / (n_perm + 1) for zj in z])}


def logrank_z(group, time, event):
    """standardized log-rank z, two groups"""
    order, t, e = _risk_sets(time, event)
    g = np.asarray(group)[order]
    o_e, v = 0.0, 0.0
    for tt in np.unique(t[e == 1]):
        r = t >= tt
        d = ((t == tt) & (e == 1))
        n, n1, dd, d1 = r.sum(), (r & g).sum(), d.sum(), (d & g).sum()
        o_e += d1 - dd * n1 / n
        if n > 1:
            v += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    return o_e / np.sqrt(v) if v > 0 else 0.0


def maxstat(x, time, event, lo=0.1, hi=0.9, n_perm=10000, seed=0):
    """best cut point + permutation p of max |z|"""
    x = np.asarray(x, float)
    cuts = np.unique(x[(x >= np.quantile(x, lo)) & (x < np.quantile(x, hi))])

    def best(tm, ev):
        zs = [abs(logrank_z(x > c, tm, ev)) for c in cuts]
        k = int(np.argmax(zs))
        return zs[k], cuts[k]
    z, cut = best(time, event)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for b in range(n_perm):
        p = rng.permutation(len(time))
        null[b] = best(time[p], event[p])[0]
    return {"cut": cut, "z": z, "p_naive": 2 * norm.sf(z), "p_adjusted": (1 + (null >= z).sum()) / (n_perm + 1)}


def icc_oneway(values, groups):
    """one-way ICC(1), spots within patients"""
    values, groups = np.asarray(values, float), np.asarray(groups)
    labs = np.unique(groups)
    k = np.array([(groups == g).sum() for g in labs])
    means = np.array([values[groups == g].mean() for g in labs])
    n, a = len(values), len(labs)
    msb = (k * (means - values.mean()) ** 2).sum() / (a - 1)
    msw = sum(((values[groups == g] - m) ** 2).sum() for g, m in zip(labs, means)) / (n - a)
    k0 = (n - (k ** 2).sum() / n) / (a - 1)
    return (msb - msw) / (msb + (k0 - 1) * msw)
