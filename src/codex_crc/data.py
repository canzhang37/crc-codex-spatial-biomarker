# CODEX single-cell table + data checks
# column names checked in check 1
import json
from pathlib import Path

import numpy as np
import pandas as pd

COLS = {"x": "X:X", "y": "Y:Y", "cluster": "ClusterName", "patient": "patients", "spot": "spots",
        "group": "groups", "cn": "neighborhood10", "cn_name": "neighborhood name", "file": "File Name"}
GROUPS = {1: "CLR", 2: "DII"}  # checked: 17 CLR, 18 DII


def load_cells(path, usecols=None):
    df = pd.read_csv(path, index_col=0, usecols=usecols)
    missing = [c for c in COLS.values() if c not in df.columns]
    if missing and usecols is None:
        print("columns not found:", missing)
    return df


def spot_key(df):
    # spot ids repeat -> file + spot key
    return df[COLS["file"]].astype(str) + "|" + df[COLS["spot"]].astype(str)


def checks(df, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    c = COLS
    info = {
        "n_cells": int(len(df)),
        "n_patients": int(df[c["patient"]].nunique()),
        "n_spots": int(spot_key(df).nunique()),
        "patients_by_group": df.groupby(c["group"])[c["patient"]].nunique().to_dict(),
        "n_clusters": int(df[c["cluster"]].nunique()),
        "n_cn": int(df[c["cn"]].nunique()),
        "columns": list(map(str, df.columns)),
        "missing_by_column": {k: int(v) for k, v in df.isna().sum().items() if v},
    }
    (out / "checks.json").write_text(json.dumps(info, indent=1, default=str))
    df[c["cluster"]].value_counts().to_csv(out / "cells_by_cluster.csv")
    pd.crosstab(df[c["cn"]], df[c["cn_name"]]).to_csv(out / "cn_numbers_vs_names.csv")
    spot_key(df).value_counts().describe().to_csv(out / "cells_per_spot.csv")
    return info


def patient_table(path):
    """Table S1, one row per patient"""
    t = pd.read_excel(path, sheet_name=0)
    t.columns = [str(x).strip() for x in t.columns]
    return t


def frac(num, den):
    return np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0)


if __name__ == "__main__":
    import sys
    print(json.dumps(checks(load_cells(sys.argv[1]), sys.argv[2]), indent=1, default=str)[:3000])
