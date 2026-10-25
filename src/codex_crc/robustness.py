# Part B robustness, 18 DII patients
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from codex_crc.prep import CD4_CLUSTERS, cn_freq, load, logf, spot_counts, survival_table
from codex_crc.replicate import FEATURE
from codex_crc.survival import cox1, icc_oneway, maxstat, westfall_young

N_PERM, N_BOOT = 10000, 2000


def holm(p):
    p = np.asarray(p)
    o = np.argsort(p)
    adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o])
    out = np.empty_like(p)
    out[o] = np.minimum(adj, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", type=Path)
    ap.add_argument("out", type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    df = load(a.raw / "CRC_clusters_neighborhoods_markers.csv")
    surv = survival_table(a.raw / "supplementary" / "mmc2.xlsx")
    dii = surv[surv.group == 2]
    t, e = dii.os.values, dii.death.values
    freq = cn_freq(df).loc[dii.index]

    ## B1 multiplicity
    X = logf(freq.values)
    wy = westfall_young(X, t, e, n_perm=N_PERM)
    fits = [cox1(X[:, j], t, e) for j in range(X.shape[1])]
    fam = pd.DataFrame({"feature": freq.columns, "hr": [f["hr"] for f in fits], "p_lr": [f["p_lr"] for f in fits],
                        "score_z": wy["z"], "p_score": wy["p_raw"], "p_westfall_young": wy["p_adj"]})
    fam["p_holm"] = holm(fam.p_lr.values)
    fam = fam.sort_values("p_lr")
    fam.to_csv(a.out / "robustness_multiplicity.csv", index=False)
    print(fam.head(8).round(4).to_string(index=False))

    ## B2 spot resampling and reliability
    num, den = spot_counts(df, "CD4+PD-1+", 9)
    sp = df.groupby("spot")["patients"].first()
    sp = sp[sp.isin(dii.index)]
    rng = np.random.default_rng(11)
    # numerator + denominator, same spots
    boot = []
    for _ in range(N_BOOT):
        x = {}
        for p, s in sp.groupby(sp).groups.items():
            pick = rng.choice(np.asarray(s), len(s))
            x[p] = num[pick].sum() / den[pick].sum() if den[pick].sum() > 0 else 0.0
        f = cox1(logf(pd.Series(x).loc[dii.index].values), t, e)
        boot.append((f["beta"], f["p_lr"]))
    boot = np.array(boot)
    ok = den[sp.index] > 0
    spot_log = logf(num[sp.index][ok] / den[sp.index][ok])
    b2 = {"share_p_below_0.05": (boot[:, 1] < 0.05).mean(), "beta_2.5": np.percentile(boot[:, 0], 2.5),
          "beta_97.5": np.percentile(boot[:, 0], 97.5), "icc_spot_log_freq": icc_oneway(spot_log, sp[ok].values),
          "spots_with_cn9": int(ok.sum()), "spots": int(len(sp))}

    ## B3 leave one out
    x0 = logf(freq[FEATURE].values)
    loo = pd.DataFrame([{"left_out": p, **cox1(np.delete(x0, i), np.delete(t, i), np.delete(e, i))}
                        for i, p in enumerate(dii.index)])
    loo.to_csv(a.out / "robustness_loo.csv", index=False)

    ## B4 pseudocount and scale
    raw = freq[FEATURE].values
    cd4 = cn_freq(df, cell_mask=df["ClusterName"].isin(CD4_CLUSTERS)).loc[dii.index, FEATURE].values
    sens = pd.DataFrame([{"version": k, **cox1(v, t, e)} for k, v in [
        ("log, pseudocount 1e-3 (paper)", logf(raw)), ("log, pseudocount 1e-4", logf(raw, 1e-4)),
        ("log, pseudocount 1e-2", logf(raw, 1e-2)), ("rank", rankdata(raw)),
        ("CD4+ T-cell clusters only (amendment 1)", logf(cd4))]])
    sens.to_csv(a.out / "robustness_sensitivity.csv", index=False)

    ## B5 best cut point
    ms = maxstat(x0, t, e, n_perm=N_PERM)

    summary = pd.Series({**b2, "loo_max_p": loo.p_lr.max(), "loo_n_above_0.05": int((loo.p_lr > 0.05).sum()),
                         "wy_p_paper_feature": fam.set_index("feature").loc[FEATURE, "p_westfall_young"],
                         "holm_p_paper_feature": fam.set_index("feature").loc[FEATURE, "p_holm"],
                         "rank_paper_feature": int(np.flatnonzero(fam.feature.values == FEATURE)[0] + 1),
                         "cut_naive_p": ms["p_naive"], "cut_adjusted_p": ms["p_adjusted"], "cut_value_log": ms["cut"]})
    summary.to_csv(a.out / "robustness_summary.csv", header=["value"])
    print(sens.round(4).to_string(index=False))
    print(summary.round(4).to_string())
    pred = {"adjusted p > 0.05": summary.wy_p_paper_feature > 0.05,
            "< 80% of spot resamples p < 0.05": summary["share_p_below_0.05"] < 0.8,
            "a leave-one-out fit has p > 0.05": summary["loo_n_above_0.05"] >= 1}
    for k, v in pred.items():
        print(f"prediction '{k}': {'met' if v else 'not met'}")


if __name__ == "__main__":
    main()
