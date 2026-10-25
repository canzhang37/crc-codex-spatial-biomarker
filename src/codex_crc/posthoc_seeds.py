# post hoc: CN-9 from each k-means rebuild
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from codex_crc.neighborhoods import match
from codex_crc.prep import load, logf, survival_table
from codex_crc.survival import cox1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", type=Path)
    ap.add_argument("results", type=Path)
    a = ap.parse_args()
    df = load(a.raw / "CRC_clusters_neighborhoods_markers.csv")
    labs = np.load(a.results / "neighborhoods" / "labels_all_seeds.npy")
    ari = pd.read_csv(a.results / "neighborhoods" / "ari_by_seed.csv")
    surv = survival_table(a.raw / "supplementary" / "mmc2.xlsx")
    dii = surv[surv.group == 2]
    gran = df["neighborhood10"].to_numpy() == 9
    rows = []
    for s, lab in enumerate(labs):
        m = match(lab, df["neighborhood10"].to_numpy())
        k = int(m.loc[m.authors == 9, "ours"].iloc[0])
        in9 = lab == k
        share = df.loc[in9].groupby("patients")["CD4+PD-1+"].mean().reindex(dii.index).fillna(0.0)
        f = cox1(logf(share.values), dii.os.values, dii.death.values)
        rows.append({"seed": s, "ari": ari.ari[s], "cells_in_cn9": int(in9.sum()),
                     "overlap_with_published_cn9": (in9 & gran).sum() / gran.sum(), "hr": f["hr"], "p_lr": f["p_lr"]})
    res = pd.DataFrame(rows)
    res.to_csv(a.results / "posthoc_seed_survival.csv", index=False)
    print(res.round(4).to_string(index=False))
    print(f"\nseeds with p < 0.05: {(res.p_lr < 0.05).sum()} of {len(res)}; p range {res.p_lr.min():.4f} to {res.p_lr.max():.4f}")

    # excl. deaths < 3 months after surgery
    share = df[df["neighborhood10"] == 9].groupby("patients")["CD4+PD-1+"].mean().reindex(dii.index).fillna(0.0)
    late = dii.os >= 3
    f = cox1(logf(share[late].values), dii.os[late].values, dii.death[late].values)
    early = pd.DataFrame([{"subset": "OS >= 3 months", "n": int(late.sum()), "deaths": int(dii.death[late].sum()), **f}])
    early.to_csv(a.results / "posthoc_early_deaths.csv", index=False)
    print(early.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
