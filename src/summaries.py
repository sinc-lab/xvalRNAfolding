"""Per-sequence F1 of every method on the test partition of each fold.

Builds, from data/<dataset>/ and predictions/<dataset>/, one table per
cross-validation strategy (index: sequence id; columns: method, fold, f1, rec,
pre; sim_folds adds similarity, fam_folds adds family) and caches it in
results/<dataset>/summary_<strategy>.pkl, which is not part of the repository.
Notebooks call load_summary(); `python src/summaries.py [dataset ...]`
(re)builds every table.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
PREDICTIONS = REPO / "predictions"
RESULTS = REPO / "results"

TRAINED = ["MxFold2", "REDfold", "UFold", "sincFold", "RiNALMo"]
CLASSICAL = ["RNAfold", "RNAstructure", "LinearFoldV", "LinearPartitionV", "ProbKnot", "IPKnot"]
METHODS = TRAINED + CLASSICAL

STRATEGIES = {
    "archiveII": ["random_kfolds", "stratified_kfold", "clustering_folds", "fam_folds", "hl_folds", "sim_folds"],
    "bpRNA": ["random_kfolds", "clustering_folds", "fam_folds", "hl_folds", "sim_folds"],
}


def bp_metrics(ref_bp, pre_bp):
    """F1, recall and precision from base pairs (exact matches only)."""
    if len(ref_bp) == 0 and len(pre_bp) == 0:
        return 1.0, 1.0, 1.0
    ref, pre = {tuple(bp) for bp in ref_bp}, {tuple(bp) for bp in pre_bp}
    tp1 = sum(tuple(bp) in ref for bp in pre_bp)
    tp2 = sum(tuple(bp) in pre for bp in ref_bp)
    fn = len(ref_bp) - tp1
    fp = len(pre_bp) - tp1
    tpr = tp1 / (tp1 + fn) if tp1 + fn > 0 else 0.0
    prec = tp2 / (tp1 + fp) if tp1 + fp > 0 else 0.0
    f1 = 2 * prec * tpr / (prec + tpr) if tpr + prec > 0 else 0.0
    return f1, tpr, prec


def load_dataset(dataset):
    return pd.read_csv(DATA / dataset / f"{dataset}.csv", index_col="id")


def read_split(dataset, strategy):
    """Split table with a `fold` column (fold_number, or fold_name for hl-folds);
    sim-folds add `similarity` (one file per threshold). For bpRNA fam-folds,
    the 20 largest bpRNA-new families (all trained on the bpRNA partition)."""
    d = DATA / dataset
    if strategy == "sim_folds":
        parts = []
        for f in sorted((d / "sim_folds_split").glob("sim*.csv"), key=lambda f: int(f.stem[3:])):
            s = pd.read_csv(f)
            s["similarity"] = f.stem
            parts.append(s)
        split = pd.concat(parts, ignore_index=True)
    elif strategy == "fam_folds" and dataset == "bpRNA":
        split = pd.read_csv(d / "fam_folds_test_split.csv")
    else:
        split = pd.read_csv(d / f"{strategy}_split.csv")
    split["fold"] = split["fold_name"] if strategy == "hl_folds" else split["fold_number"]
    return split


def pred_path(dataset, strategy, method, fold, similarity=None):
    p = PREDICTIONS / dataset
    if strategy == "random_kfolds" and dataset == "archiveII":
        return p / "random_kfolds/standard_kfold" / method / str(fold) / "pred.csv"
    if strategy == "stratified_kfold":
        return p / "random_kfolds/stratified_kfold" / method / str(fold) / "pred.csv"
    if strategy == "hl_folds":
        return p / "hl_folds/F1mean_folds" / method / str(fold) / "pred.csv"
    if strategy == "sim_folds":
        return p / "sim_folds" / method / similarity / str(fold) / "pred.csv"
    return p / strategy / method / str(fold) / "pred.csv"


def _scores(pred, ref):
    m = [bp_metrics(json.loads(ref[i]), json.loads(b)) for i, b in zip(pred.index, pred["base_pairs"])]
    return pd.DataFrame(m, index=pred.index, columns=["f1", "rec", "pre"])


def classical_summary(dataset, rebuild=False):
    """F1 of the classical (untrained) methods on every sequence of the dataset."""
    cache = RESULTS / dataset / "classical_summary.pkl"
    if cache.exists() and not rebuild:
        return pd.read_pickle(cache)
    ref = load_dataset(dataset)["base_pairs"]
    parts = []
    for method in CLASSICAL:
        pred = pd.read_csv(PREDICTIONS / dataset / f"{method}.csv", index_col="id")
        s = _scores(pred, ref)
        s.insert(0, "method", method)
        parts.append(s)
    summary = pd.concat(parts)
    cache.parent.mkdir(parents=True, exist_ok=True)
    summary.to_pickle(cache)
    return summary


def load_summary(dataset, strategy, rebuild=False):
    """Per-sequence scores of all methods on each test fold of a strategy."""
    cache = RESULTS / dataset / f"summary_{strategy}.pkl"
    if cache.exists() and not rebuild:
        return pd.read_pickle(cache)
    ref = load_dataset(dataset)["base_pairs"]
    classical = classical_summary(dataset, rebuild)
    split = read_split(dataset, strategy)
    test = split[split["partition"] == "test"]
    keys = ["similarity", "fold"] if strategy == "sim_folds" else ["fold"]
    groups = [(k if isinstance(k, tuple) else (k,), g) for k, g in test.groupby(keys, sort=False)]
    if strategy == "sim_folds":
        groups.sort(key=lambda kg: (int(kg[0][0][3:]), kg[0][1]))
    else:
        groups.sort(key=lambda kg: int(str(kg[0][0]).replace("hl", "")))

    def tag(s, key):
        if strategy == "sim_folds":
            s.insert(1, "similarity", key[0])
        s.insert(len(keys), "fold", key[-1])
        return s

    # trained methods fold by fold, then classical methods, as in the original summaries
    parts = []
    for key, g in groups:
        for method in TRAINED:
            f = pred_path(dataset, strategy, method, key[-1], key[0] if strategy == "sim_folds" else None)
            if not f.exists():
                continue
            pred = pd.read_csv(f, index_col="id")
            if set(pred.index) != set(g["id"]):
                raise ValueError(f"{f}: ids differ from the test partition of the split")
            s = _scores(pred, ref)
            s.insert(0, "method", method)
            parts.append(tag(s, key))
    for method in CLASSICAL:
        c = classical[classical["method"] == method]
        for key, g in groups:
            parts.append(tag(c.loc[g["id"]].copy(), key))
    summary = pd.concat(parts)
    if strategy == "fam_folds":
        family = split.drop_duplicates("fold_number").set_index("fold_number")["fold_name"]
        summary["family"] = summary["fold"].map(family)
    cache.parent.mkdir(parents=True, exist_ok=True)
    summary.to_pickle(cache)
    return summary


def sim_auc(summary_sim):
    """Area under the sim-fold curve (mean F1 per similarity threshold), plain
    and weighted by difficulty (w = 1 - threshold)."""
    out = {}
    for method, s in summary_sim.groupby("method", sort=False):
        f1 = s.groupby("similarity")["f1"].mean()
        w = 1 - f1.index.str[3:].astype(float).to_numpy() / 100
        out[method] = {"AUC": f1.mean(), "wAUC": (w * f1.to_numpy()).sum() / w.sum()}
    return pd.DataFrame(out).T


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("datasets", nargs="*", default=list(STRATEGIES))
    args = parser.parse_args()
    for dataset in args.datasets:
        classical_summary(dataset, rebuild=True)
        for strategy in STRATEGIES[dataset]:
            s = load_summary(dataset, strategy, rebuild=True)
            print(f"{dataset} {strategy}: {len(s)} rows, {s['method'].nunique()} methods")


if __name__ == "__main__":
    main()
