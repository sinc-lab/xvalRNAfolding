"""Generate a cluster-fold split for archiveII at a given identity threshold.

Usage: python make_clustering_folds.py <th> <tool>
    th: identity threshold (e.g. 0.8 or 80); cdhit supports 0.80-1.00, psi 0.60-0.80
    tool: cdhit (cd-hit-est) or psi (psi-cd-hit, requires BLAST+ in PATH)
Output: data/clustering_folds_split/clustering_<tool>_<th>_folds_split.csv
"""

import argparse
import random
import subprocess
from pathlib import Path

import pandas as pd
from cdhit_reader import read_cdhit

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT_DIR = DATA / "clustering_folds_split"
CDHIT_DIR = OUT_DIR / "cdhit"

N_FOLDS = 5
TEST_FRAC = 0.1
VALID_FRAC = 0.1


def write_fasta(df, path):
    with open(path, "w") as f:
        for seq_id, row in df.iterrows():
            f.write(f">{seq_id}\n{row.sequence}\n")


def run_cdhit(fasta, th, cdhit_bin):
    out = CDHIT_DIR / f"archiveII-cdhit_{th:g}"
    cmd = [
        cdhit_bin, "-i", str(fasta), "-o", str(out),
        "-c", f"{th:g}", "-T", "0", "-d", "100", "-g", "1",
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    with open(f"{out}.log", "w") as f:
        f.write(" ".join(cmd) + "\n\n" + res.stdout)
    return Path(f"{out}.clstr")


def run_psi_cdhit(fasta, th, psi_script):
    out = CDHIT_DIR / f"archiveII-psi_{th:g}"
    clstr = Path(f"{out}.clstr")
    if clstr.exists():
        return clstr
    cmd = [
        "perl", str(psi_script), "-i", str(fasta), "-o", str(out),
        "-c", f"{th:g}", "-G", "1", "-g", "1",
        "-prog", "blastn", "-exec", "local", "-para", "6", "-blp", "3",
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    with open(f"{out}.run.log", "w") as f:
        f.write(" ".join(cmd) + "\n\n" + res.stdout)
    return clstr


def load_clusters(clstr_path, ids):
    id2cluster = {}
    for c in read_cdhit(str(clstr_path)):
        cluster = int(c.name.split()[-1])
        for seq in c.sequences:
            id2cluster[seq.name] = cluster
    missing = set(ids) - set(id2cluster)
    if missing:
        raise ValueError(f"{clstr_path}: {len(missing)} ids missing from clusters")
    return id2cluster


def make_splits(id2cluster, n_seqs):
    members = {}
    for seq_id, cluster in id2cluster.items():
        members.setdefault(cluster, []).append(seq_id)
    clusters = sorted(members)

    rows = []
    for fold in range(N_FOLDS):
        random.shuffle(clusters)
        test, valid, train = [], [], []
        for c in clusters:
            if len(test) <= n_seqs * TEST_FRAC:
                test += members[c]
            elif len(valid) <= n_seqs * VALID_FRAC:
                valid += members[c]
            else:
                train += members[c]
        for partition, ids in [("train", train), ("valid", valid), ("test", test)]:
            rows += [(fold, fold, partition, i) for i in ids]
    return pd.DataFrame(rows, columns=["fold_name", "fold_number", "partition", "id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("th", type=float)
    parser.add_argument("tool", choices=["cdhit", "psi"])
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cdhit-bin", default="cd-hit-est")
    parser.add_argument("--psi-script", default="psi-cd-hit.pl")
    args = parser.parse_args()

    th = args.th / 100 if args.th > 1 else args.th

    df = pd.read_csv(DATA / "archiveII.csv", index_col="id")
    CDHIT_DIR.mkdir(parents=True, exist_ok=True)
    fasta = CDHIT_DIR / "archiveII.fasta"
    write_fasta(df, fasta)

    if args.tool == "cdhit":
        clstr = run_cdhit(fasta, th, args.cdhit_bin)
    else:
        clstr = run_psi_cdhit(fasta, th, args.psi_script)
    id2cluster = load_clusters(clstr, df.index)

    random.seed(args.seed + round(th * 100))
    split = make_splits(id2cluster, len(df))
    out_csv = OUT_DIR / f"clustering_{args.tool}_{round(th * 100)}_folds_split.csv"
    split.to_csv(out_csv, index=False)
    print(out_csv)


if __name__ == "__main__":
    main()
