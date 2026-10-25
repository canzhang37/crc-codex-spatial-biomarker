# cellular neighborhoods, 10-NN + k-means
# ARI vs the authors' CNs
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

from codex_crc.data import COLS, load_cells, spot_key


def windows(xy, types, spots, k=10):
    """type frequencies of the k nearest cells, same spot"""
    n_types = types.max() + 1
    W = np.zeros((len(types), n_types))
    for s in np.unique(spots):
        idx = np.flatnonzero(spots == s)
        kk = min(k, len(idx))
        nn = NearestNeighbors(n_neighbors=kk).fit(xy[idx])
        _, nb = nn.kneighbors(xy[idx])
        for j in range(kk):
            np.add.at(W, (idx, types[idx[nb[:, j]]]), 1.0)
        W[idx] /= kk
    return W


def match(a, b):
    """1-to-1 label matching by shared cells"""
    ct = pd.crosstab(a, b)
    r, c = linear_sum_assignment(-ct.values)
    return pd.DataFrame({"ours": ct.index[r], "authors": ct.columns[c], "shared_cells": ct.values[r, c]})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--seeds", type=int, default=20)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    df = load_cells(a.csv)
    types = pd.Categorical(df[COLS["cluster"]]).codes
    xy = df[[COLS["x"], COLS["y"]]].to_numpy(float)
    W = windows(xy, types, spot_key(df).to_numpy())
    authors = df[COLS["cn"]].to_numpy()

    rows, best, labs = [], (-1, None), []
    for s in range(a.seeds):
        lab = MiniBatchKMeans(n_clusters=10, random_state=s, n_init=3).fit_predict(W)
        ari = adjusted_rand_score(authors, lab)
        rows.append({"seed": s, "ari": ari})
        labs.append(lab.astype(np.int8))
        if ari > best[0]:
            best = (ari, lab)
        print(f"seed {s}: ARI {ari:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(a.out / "ari_by_seed.csv", index=False)
    match(best[1], authors).to_csv(a.out / "best_seed_matching.csv", index=False)
    np.save(a.out / "best_labels.npy", best[1])
    np.save(a.out / "labels_all_seeds.npy", np.vstack(labs))  # post hoc seed check


if __name__ == "__main__":
    main()
