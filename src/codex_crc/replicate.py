# Part A2: DII Cox model on CN-9 share
import argparse
from pathlib import Path

import pandas as pd

from codex_crc.prep import CD4_CLUSTERS, cn_freq, load, logf, survival_table
from codex_crc.survival import cox1

FEATURE = "CD4+PD-1+|CN-9"


def dii_table(raw):
    df = load(raw / "CRC_clusters_neighborhoods_markers.csv")
    surv = survival_table(raw / "supplementary" / "mmc2.xlsx")
    dii = surv[surv.group == 2]
    return df, dii


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", type=Path)
    ap.add_argument("out", type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    df, dii = dii_table(a.raw)
    freq = cn_freq(df).loc[dii.index]
    cd4 = df["ClusterName"].isin(CD4_CLUSTERS)
    rows = []
    for name, x in [("main: share of all CN-9 cells", freq[FEATURE]),
                    ("CD4+ T-cell clusters only", cn_freq(df, cell_mask=cd4).loc[dii.index, FEATURE])]:
        f = cox1(logf(x), dii.os.values, dii.death.values)
        rows.append({"definition": name, "n": len(dii), "deaths": int(dii.death.sum()), "zero_patients": int((x == 0).sum()), **f})
    # denominator: CD4+ T-cell clusters in CN-9
    d9 = df[(df["neighborhood number final"] == 9) & cd4]
    alt = (d9.groupby("patients")["CD4+PD-1+"].mean()).reindex(dii.index).fillna(0)
    f = cox1(logf(alt), dii.os.values, dii.death.values)
    rows.append({"definition": "share among CD4+ T-cell-cluster cells in CN-9", "n": len(dii), "deaths": int(dii.death.sum()),
                 "zero_patients": int((alt == 0).sum()), **f})
    res = pd.DataFrame(rows)
    res.to_csv(a.out / "replication.csv", index=False)
    freq.join(dii).to_csv(a.out / "dii_features.csv")
    print(res.round(4).to_string(index=False))
    p = res.p_lr.iloc[0]
    print(f"\nreplication target (same direction, LR p within 0.003-0.012): {'met' if 0.003 <= p <= 0.012 and res.beta.iloc[0] < 0 else 'NOT met'}")


if __name__ == "__main__":
    main()
