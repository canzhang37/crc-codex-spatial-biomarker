# patient and spot-level features
import numpy as np
import pandas as pd


def otsu(v, bins=256):
    """Otsu threshold, 1-d array"""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    hist, edges = np.histogram(v, bins=bins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = ((hist * mids).sum() - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    between = w0 * w1 * (m0 - m1) ** 2
    return mids[int(np.argmax(between[:-1]))]


def positive(intensity, mask, cofactor=5.0):
    """Otsu call on arcsinh(x / cofactor)"""
    z = np.arcsinh(np.asarray(intensity, float) / cofactor)
    return z > otsu(z[mask])


def cn_frequency(df, flag, patient, cn, unit=None):
    """flagged share per CN, per patient or spot"""
    key = [unit] if unit else [patient]
    g = df.assign(_f=np.asarray(flag, float)).groupby(key + [cn])._f
    return (g.sum() / g.size()).unstack(cn).fillna(0.0)


def composition_clr(counts, pseudo=0.5):
    """centered log-ratio of cell-type counts"""
    lc = np.log(np.asarray(counts, float) + pseudo)
    return pd.DataFrame(lc - lc.mean(1, keepdims=True), index=counts.index, columns=counts.columns)


def spot_bootstrap(num_by_spot, den_by_spot, spot_patient, rng):
    """spot resample within patient"""
    out = {}
    for p, spots in spot_patient.groupby(spot_patient).groups.items():
        pick = rng.choice(np.asarray(spots), len(spots), replace=True)
        n, d = num_by_spot.loc[pick].sum(), den_by_spot.loc[pick].sum()
        out[p] = n / d if d > 0 else 0.0
    return pd.Series(out)
