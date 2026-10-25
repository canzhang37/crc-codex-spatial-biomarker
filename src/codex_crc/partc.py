# Part C: composition vs location features
# follicle CN, B cells, unclear clusters dropped
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score

from codex_crc.classify import boot_auc, loo_probs
from codex_crc.contacts import contact_log_ratio, delaunay_edges
from codex_crc.features import composition_clr
from codex_crc.prep import CN, CNS, load

DROP = ["dirt", "undefined", "tumor cells / immune cells", "immune cells / vasculature", "B cells"]
N_PERM, N_BOOT = 1000, 2000


def keep_cells(df):
    return df[~df["ClusterName"].isin(DROP) & df[CN].isin(CNS)].copy()


def set1(d):
    counts = d.groupby(["patients", "ClusterName"]).size().unstack(fill_value=0)
    return composition_clr(counts)


def set2_extra(d):
    g = d.groupby(["patients", CN, "ClusterName"]).size()
    share = (g / g.groupby(level=[0, 1]).transform("sum")).unstack([1, 2], fill_value=0.0)
    share.columns = [f"{t}|CN-{int(c)}" for c, t in share.columns]
    return share


def set3_extra(d, n_shuffle=200):
    types = pd.Categorical(d["ClusterName"])
    codes, names = types.codes, list(types.categories)
    k = len(names)
    iu = np.triu_indices(k)
    rows = {}
    for s, idx in d.groupby("spot").indices.items():
        if len(idx) < 10:
            continue
        xy = d[["X:X", "Y:Y"]].to_numpy(float)[idx]
        lr = contact_log_ratio(codes[idx], delaunay_edges(xy), k, n_shuffle, np.random.default_rng(len(idx)))
        rows[s] = lr[iu]
    spot = pd.DataFrame(rows, index=[f"{names[i]}~{names[j]}" for i, j in zip(*iu)]).T
    pat = d.groupby("spot")["patients"].first()
    return spot.groupby(pat.loc[spot.index].values).mean()


def evaluate(name, X, y):
    p = loo_probs(X, y)
    auc = roc_auc_score(y, p)
    rng = np.random.default_rng(2)
    perms = [rng.permutation(y) for _ in range(N_PERM)]
    null = Parallel(n_jobs=-1)(delayed(lambda yp: roc_auc_score(yp, loo_probs(X, yp)))(yp) for yp in perms)
    return p, {"set": name, "n_features": X.shape[1], "auc": auc, "p_perm": (1 + (np.array(null) >= auc).sum()) / (N_PERM + 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", type=Path)
    ap.add_argument("out", type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    d = keep_cells(load(a.raw / "CRC_clusters_neighborhoods_markers.csv"))
    y_ser = d.groupby("patients")["groups"].first()
    s1 = set1(d)
    sets = {"1 composition": s1, "2 composition + CN-specific": s1.join(set2_extra(d)),
            "3 composition + contacts": s1.join(set3_extra(d))}
    y = (y_ser.loc[s1.index].values == 2).astype(int)  # DII = 1
    rows, probs = [], {}
    for name, X in sets.items():
        X = X.loc[s1.index].fillna(0.0)
        probs[name], r = evaluate(name, X.values, y)
        rows.append(r)
        print(r, flush=True)
    # paired patient bootstrap
    rng = np.random.default_rng(3)
    idx = []
    while len(idx) < N_BOOT:
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(y):
            idx.append(i)
    names = list(sets)
    B = np.array([[roc_auc_score(y[i], probs[n][i]) for n in names] for i in idx])
    for k, r in enumerate(rows):
        r["auc_lo"], r["auc_hi"] = np.percentile(B[:, k], [2.5, 97.5])
        if k:
            dlt = B[:, k] - B[:, 0]
            r["diff_vs_set1"] = r["auc"] - rows[0]["auc"]
            r["diff_lo"], r["diff_hi"] = np.percentile(dlt, [2.5, 97.5])
            r["adds_by_rule"] = bool(r["diff_vs_set1"] >= 0.05 and r["diff_lo"] > 0)
    res = pd.DataFrame(rows)
    res.to_csv(a.out / "composition_vs_spatial.csv", index=False)
    pd.DataFrame(probs, index=s1.index).assign(dii=y).to_csv(a.out / "partc_loo_probs.csv")
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
