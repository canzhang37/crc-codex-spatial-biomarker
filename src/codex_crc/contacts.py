# Delaunay cell-type contacts vs label shuffles
import numpy as np
from scipy.spatial import Delaunay


def delaunay_edges(xy):
    """unique undirected Delaunay edges"""
    tri = Delaunay(xy).simplices
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]])
    e.sort(1)
    return np.unique(e, axis=0)


def contact_counts(types, edges, n_types):
    m = np.zeros((n_types, n_types))
    np.add.at(m, (types[edges[:, 0]], types[edges[:, 1]]), 1)
    return m + m.T - np.diag(np.diag(m))  # symmetric, same-type edges once


def contact_log_ratio(types, edges, n_types, n_shuffle=200, rng=None):
    """log enrichment vs shuffles, per type pair"""
    rng = rng or np.random.default_rng(0)
    obs = contact_counts(types, edges, n_types)
    exp = np.mean([contact_counts(rng.permutation(types), edges, n_types) for _ in range(n_shuffle)], axis=0)
    return np.log((obs + 0.5) / (exp + 0.5))
