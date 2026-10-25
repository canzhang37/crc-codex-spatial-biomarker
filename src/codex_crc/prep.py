# patient + spot tables for Parts A-C
import numpy as np
import pandas as pd

T_FLAGS = ["CD4+ICOS+", "CD4+Ki67+", "CD4+PD-1+", "CD8+ICOS+", "CD8+Ki67+", "CD8+PD-1+", "Treg-ICOS+", "Treg-Ki67+", "Treg-PD-1+"]
CN = "neighborhood number final"  # paper's CN-1..9, artifact CN empty
CNS = [1, 2, 3, 4, 6, 7, 8, 9]  # CN-5 (follicle) left out
CD4_CLUSTERS = ["CD4+ T cells", "CD4+ T cells CD45RO+", "CD4+ T cells GATA3+"]


def load(csv):
    use = ["patients", "groups", "spots", "File Name", "ClusterName", "X:X", "Y:Y", "neighborhood10", CN] + T_FLAGS
    df = pd.read_csv(csv, usecols=use)
    df["spot"] = df["File Name"].astype(str) + "|" + df["spots"].astype(str)
    return df


def survival_table(xlsx):
    t = pd.read_excel(xlsx, "A. Patient_data_TMA_annotations")
    t = t[pd.to_numeric(t["Patient"], errors="coerce").notna()].copy()
    t = t.astype({"Patient": int, "Group": int})
    t = t.rename(columns={"Patient": "patient", "Group": "group", "OS": "os", "OS_Censor": "death"})
    return t.set_index("patient")[["group", "os", "death"]].astype({"os": float, "death": int})


def cn_freq(df, flags=T_FLAGS, cns=CNS, by="patients", cell_mask=None):
    """flagged share per CN (0 if no cells)"""
    d = df[df[CN].isin(cns)]
    num = d[flags].copy()
    if cell_mask is not None:
        num = num.mul(cell_mask[d.index].astype(int), axis=0)
    g = num.assign(_u=d[by].values, _cn=d[CN].astype(int).values).groupby(["_u", "_cn"])
    freq = g[flags].sum().div(g.size(), axis=0)
    out = freq.unstack("_cn").reindex(columns=pd.MultiIndex.from_product([flags, cns])).fillna(0.0)
    out.columns = [f"{f}|CN-{c}" for f, c in out.columns]
    out.index.name = by
    return out


def spot_counts(df, flag, cn, cell_mask=None):
    """per spot: flagged and total cells in CN"""
    d = df[df[CN] == cn]
    v = d[flag].astype(int)
    if cell_mask is not None:
        v = v * cell_mask[d.index].astype(int)
    num = v.groupby(d["spot"]).sum().reindex(df["spot"].unique(), fill_value=0)
    den = d.groupby("spot").size().reindex(df["spot"].unique(), fill_value=0)
    return num, den


def logf(x, pseudo=1e-3):
    return np.log(pseudo + np.asarray(x, float))
