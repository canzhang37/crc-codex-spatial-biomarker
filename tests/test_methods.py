import numpy as np
from sklearn.metrics import roc_auc_score

from codex_crc.classify import loo_probs
from codex_crc.contacts import contact_counts, delaunay_edges
from codex_crc.features import otsu
from codex_crc.neighborhoods import windows
from codex_crc.survival import cox1, icc_oneway, logrank_z, maxstat, score_z, westfall_young


def sim(n=300, beta=0.7, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    t = rng.exponential(1 / np.exp(beta * x))
    c = rng.exponential(2.0, n)
    return x, np.minimum(t, c), (t <= c).astype(int)


def test_cox_recovers_effect():
    x, t, e = sim()
    f = cox1(x, t, e)
    assert abs(f["beta"] - 0.7) < 0.2 and f["p_lr"] < 1e-6


def test_score_equals_logrank_for_two_groups():
    x, t, e = sim(n=80, seed=3)
    g = x > 0
    assert np.isclose(score_z(g[:, None].astype(float), t, e)[0], logrank_z(g, t, e))


def test_westfall_young_is_conservative():
    x, t, e = sim(n=40, beta=0.0, seed=4)
    X = np.random.default_rng(5).normal(size=(40, 8))
    r = westfall_young(X, t, e, n_perm=300)
    assert (r["p_adj"] >= r["p_raw"] - 0.05).all() and r["p_adj"].min() > 0.01


def test_maxstat_corrects_the_search():
    x, t, e = sim(n=30, beta=0.0, seed=6)
    r = maxstat(x, t, e, n_perm=300)
    assert r["p_adjusted"] >= r["p_naive"]


def test_windows_rows_are_frequencies():
    rng = np.random.default_rng(0)
    xy, types = rng.uniform(size=(200, 2)), rng.integers(0, 4, 200)
    spots = np.repeat([0, 1], 100)
    W = windows(xy, types, spots, k=10)
    assert np.allclose(W.sum(1), 1) and W.shape == (200, 4)


def test_contacts_symmetric():
    xy = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [0.5, 0.5]])
    e = delaunay_edges(xy)
    m = contact_counts(np.array([0, 0, 1, 1, 2]), e, 3)
    assert np.allclose(m, m.T) and len(e) == 8


def test_icc_and_otsu():
    g = np.repeat(np.arange(10), 4)
    v = g * 10 + np.random.default_rng(1).normal(scale=0.1, size=40)
    assert icc_oneway(v, g) > 0.99
    bimodal = np.r_[np.random.default_rng(2).normal(0, 1, 500), np.random.default_rng(3).normal(8, 1, 500)]
    assert 2 < otsu(bimodal) < 6


def test_loo_classifier_on_separable_toy():
    rng = np.random.default_rng(7)
    y = np.repeat([0, 1], 15)
    X = rng.normal(size=(30, 3)) + y[:, None] * 2.5
    assert roc_auc_score(y, loo_probs(X, y)) > 0.9
