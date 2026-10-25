# Part C: CLR vs DII from patient features
# L1 logistic, nested leave-one-patient-out
# AUC w/ patient bootstrap + label permutations
import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

C_GRID = np.logspace(-2, 2, 10)


def loo_probs(X, y, seed=0):
    X, y = np.asarray(X, float), np.asarray(y)
    p = np.empty(len(y))
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        clf = make_pipeline(StandardScaler(), LogisticRegressionCV(
            Cs=C_GRID, penalty="l1", solver="liblinear", scoring="neg_log_loss",
            cv=StratifiedKFold(5, shuffle=True, random_state=seed), max_iter=5000))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            clf.fit(X[tr], y[tr])
        p[i] = clf.predict_proba(X[i:i + 1])[0, 1]
    return p


def boot_auc(y, p, n_boot=2000, seed=1):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n_boot:
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) == 2:
            out.append(roc_auc_score(y[i], p[i]))
    return np.array(out)


def perm_p(X, y, auc_obs, n_perm=1000, seed=2):
    rng = np.random.default_rng(seed)
    null = np.array([roc_auc_score(yp, loo_probs(X, yp)) for yp in (rng.permutation(y) for _ in range(n_perm))])
    return (1 + (null >= auc_obs).sum()) / (n_perm + 1), null
