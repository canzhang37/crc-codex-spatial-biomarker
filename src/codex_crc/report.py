# readout figures
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from codex_crc.prep import load, logf, spot_counts, survival_table  # noqa: E402
from codex_crc.replicate import FEATURE  # noqa: E402


def km(t, e):
    order = np.argsort(t)
    t, e = t[order], e[order]
    s, xs, ys = 1.0, [0.0], [1.0]
    for i, tt in enumerate(t):
        if e[i]:
            s *= 1 - 1 / (len(t) - i)
        xs += [tt, tt]
        ys += [ys[-1], s]
    return xs, ys


def multiplicity(res, fig_dir):
    fam = pd.read_csv(res / "robustness_multiplicity.csv")
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    x = np.arange(1, len(fam) + 1)
    ax.scatter(x, -np.log10(fam.p_lr), s=12, label="unadjusted (likelihood ratio)")
    ax.scatter(x, -np.log10(fam.p_westfall_young), s=12, marker="x", label="Westfall-Young adjusted")
    ax.axhline(-np.log10(0.05), color="grey", lw=0.8, ls="--")
    ax.annotate("PD-1+ CD4+ in CN-9", (1, -np.log10(fam.p_lr.iloc[0])), xytext=(12, 2.0), fontsize=8,
                arrowprops={"arrowstyle": "->", "lw": 0.6})
    ax.set_xlabel("72 CN x T-cell-subset features, ranked by unadjusted p")
    ax.set_ylabel("-log10 p")
    ax.set_title("Survival in 18 DII patients (13 deaths)", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "multiplicity.png", dpi=140)
    plt.close(fig)


def spots(raw, res, fig_dir):
    df = load(raw / "CRC_clusters_neighborhoods_markers.csv")
    surv = survival_table(raw / "supplementary" / "mmc2.xlsx")
    dii = surv[surv.group == 2].sort_values("os")
    num, den = spot_counts(df, "CD4+PD-1+", 9)
    sp = df.groupby("spot")["patients"].first()
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for i, p in enumerate(dii.index):
        s = sp[sp == p].index
        v = logf(num[s] / den[s].where(den[s] > 0))
        ax.scatter(np.full(len(v), i), v, s=14, color="firebrick" if dii.death[p] else "steelblue", alpha=0.8)
        pv = logf(num[s].sum() / den[s].sum())
        ax.plot([i - 0.3, i + 0.3], [pv, pv], color="black", lw=1)
    ax.set_xticks(range(len(dii)), [f"{m:.0f}" for m in dii.os], fontsize=7)
    ax.set_xlabel("DII patients by overall survival, months (red: died)")
    ax.set_ylabel("log(0.001 + share), per spot")
    ax.set_title("PD-1+ CD4+ share of CN-9 cells: four spots per patient, bar = pooled", fontsize=9)
    fig.tight_layout()
    fig.savefig(fig_dir / "spots.png", dpi=140)
    plt.close(fig)


def cutpoint(res, fig_dir):
    d = pd.read_csv(res / "dii_features.csv", index_col=0)
    s = pd.read_csv(res / "robustness_summary.csv", index_col=0)["value"]
    x = logf(d[FEATURE].values)
    hi = x > s["cut_value_log"]
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for g, lab, col in [(hi, "above cut", "steelblue"), (~hi, "at or below cut", "firebrick")]:
        xs, ys = km(d.os.values[g], d.death.values[g])
        ax.step(xs, ys, where="post", color=col, label=f"{lab} (n = {g.sum()})")
    ax.set_xlabel("months")
    ax.set_ylabel("overall survival")
    ax.set_title(f"Best split: log-rank p {s['cut_naive_p']:.4f},\ncorrected for the search {s['cut_adjusted_p']:.3f}", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "cutpoint.png", dpi=140)
    plt.close(fig)


def partc(res, fig_dir):
    r = pd.read_csv(res / "composition_vs_spatial.csv")
    fig, ax = plt.subplots(figsize=(6.2, 2.8))
    y = np.arange(len(r))[::-1]
    ax.errorbar(r.auc, y, xerr=[r.auc - r.auc_lo, r.auc_hi - r.auc], fmt="o", capsize=3)
    ax.set_yticks(y, r["set"], fontsize=8)
    ax.axvline(0.5, color="grey", lw=0.8, ls="--")
    ax.set_xlabel("leave-one-patient-out AUC, CLR vs DII (95% bootstrap)")
    fig.tight_layout()
    fig.savefig(fig_dir / "partc_auc.png", dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", type=Path)
    ap.add_argument("results", type=Path)
    ap.add_argument("figures", type=Path)
    a = ap.parse_args()
    a.figures.mkdir(parents=True, exist_ok=True)
    multiplicity(a.results, a.figures)
    spots(a.raw, a.results, a.figures)
    cutpoint(a.results, a.figures)
    if (a.results / "composition_vs_spatial.csv").exists():
        partc(a.results, a.figures)
    print("figures in", a.figures)


if __name__ == "__main__":
    main()
