"""Generate a cluster-fold split for archiveII at a given identity threshold.

Usage: python make_clustering_folds.py <th> <tool> [--reps]
    th: identity threshold (e.g. 0.8 or 80); cdhit supports 0.80-1.00, psi >=0.30
    tool: cdhit (cd-hit-est), psi (psi-cd-hit, needs BLAST+) or mmseqs
          (MMseqs2 --cluster-mode 1, transitive/connected-component clustering)
    --reps: drop all but the representative of each cluster before splitting
Output: data/archiveII/multithreshold_clustering/clustering_<tool>[rep]_<th>_folds_split.csv
    (.clstr in its clstr/ subfolder; mmseqs splits and all working files in
    work/clustering/)
"""

import argparse
import os
import random
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from cdhit_reader import read_cdhit

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "data" / "archiveII"
OUT_DIR = DATA / "multithreshold_clustering"
CLSTR_DIR = OUT_DIR / "clstr"
WORK_DIR = REPO / "work" / "clustering"

N_FOLDS = 5
TEST_FRAC = 0.1
VALID_FRAC = 0.1


def write_fasta(df, path):
    with open(path, "w") as f:
        for seq_id, row in df.iterrows():
            f.write(f">{seq_id}\n{row.sequence}\n")


def run_cdhit(fasta, th, cdhit_bin):
    out = WORK_DIR / f"archiveII-cdhit_{th:g}"
    clstr = CLSTR_DIR / f"{out.name}.clstr"
    if clstr.exists():
        return clstr
    cmd = [
        cdhit_bin, "-i", str(fasta), "-o", str(out),
        "-c", f"{th:g}", "-T", "0", "-d", "100", "-g", "1",
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    with open(f"{out}.log", "w") as f:
        f.write(" ".join(cmd) + "\n\n" + res.stdout)
    shutil.move(f"{out}.clstr", clstr)
    return clstr


def find_psi_script(name):
    """psi-cd-hit.pl is not on PATH in the bioconda cd-hit package: it ships
    under <env>/share/cdhit-src/psi-cd-hit/. Look there before giving up."""
    if os.path.sep in name or Path(name).exists():
        return Path(name)
    found = shutil.which(name)
    if found:
        return Path(found)
    prefix = os.environ.get("CONDA_PREFIX")
    candidates = []
    if prefix:
        candidates.append(Path(prefix) / "share/cdhit-src/psi-cd-hit" / name)
    candidates += list(Path.home().glob(f"*/envs/*/share/cdhit-src/psi-cd-hit/{name}"))
    for c in candidates:
        if c.exists():
            return c
    raise SystemExit(
        f"{name} not found. Activate the cdhit env (conda activate cdhit) so that "
        f"blastn is on PATH, or pass --psi-script /path/to/{name}")


def run_psi_cdhit(fasta, th, psi_script):
    out = WORK_DIR / f"archiveII-psi_{th:g}"
    clstr = CLSTR_DIR / f"{out.name}.clstr"
    if clstr.exists():
        return clstr
    script = find_psi_script(psi_script)
    if not shutil.which("blastn"):
        # psi-cd-hit calls blastn/blastp by name: borrow the bin/ of the env the
        # script was found in, so this also runs from base or sincfold.
        env_bin = script.parents[3] / "bin"      # <env>/share/cdhit-src/psi-cd-hit/
        if not (env_bin / "blastn").exists():
            raise SystemExit("blastn not found on PATH; conda activate cdhit")
        os.environ["PATH"] = f"{env_bin}{os.pathsep}{os.environ['PATH']}"
    # -i must be RELATIVE and the cwd must be its directory: psi-cd-hit builds
    # its worker script path as "./" + the -i value, so an absolute path yields
    # ".//abs/path-bl.pl: No such file", blast never runs, and every sequence
    # silently comes back as its own singleton cluster.
    # NOTE: -para > 1 makes the result depend on the degree of parallelism.
    # psi-cd-hit pre-launches BLAST for a window of candidates against the
    # representative DB as it stood when the batch started, so a sequence in
    # the batch cannot see a representative created earlier in the same batch
    # and both become representatives. At 0.4 this yields 653 clusters with
    # -para 6 vs 643 with -para 1, and splits 5 groups of byte-identical
    # sequences across clusters (see cluster_rationale.md 8.6). The committed
    # splits were generated with -para 6; use -para 1 for a clean clustering.
    cmd = [
        "perl", str(script), "-i", fasta.name, "-o", str(out),
        "-c", f"{th:g}", "-G", "1", "-g", "1",
        "-prog", "blastn", "-exec", "local", "-para", "6", "-blp", "3",
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True,
                         cwd=fasta.parent)
    with open(f"{out}.run.log", "w") as f:
        f.write(" ".join(cmd) + "\n\n" + res.stdout)
    shutil.move(f"{out}.clstr", clstr)
    return clstr


def run_mmseqs(fasta, th, mmseqs_bin):
    out = WORK_DIR / f"archiveII-mmseqs_{th:g}"
    tsv = Path(f"{out}.tsv")
    if tsv.exists():
        return tsv
    db = WORK_DIR / "archiveII-mmseqsdb"
    if not db.exists():
        subprocess.run([mmseqs_bin, "createdb", str(fasta), str(db)],
                       check=True, capture_output=True)
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [
            mmseqs_bin, "cluster", str(db), str(out), tmp,
            "--min-seq-id", f"{th:g}", "--cluster-mode", "1", "--seq-id-mode", "1",
            "-c", "0", "--single-step-clustering", "-s", "7.5",
        ]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        subprocess.run([mmseqs_bin, "createtsv", str(db), str(db), str(out), str(tsv)],
                       check=True, capture_output=True)
    with open(f"{out}.log", "w") as f:
        f.write(" ".join(cmd) + "\n\n" + res.stdout)
    return tsv


def load_clusters(path, ids):
    id2cluster, reps = {}, set()
    if path.suffix == ".tsv":
        rep2num = {}
        for line in open(path):
            rep, member = line.split()
            reps.add(rep)
            id2cluster[member] = rep2num.setdefault(rep, len(rep2num))
    else:
        for c in read_cdhit(str(path)):
            cluster = int(c.name.split()[-1])
            for seq in c.sequences:
                id2cluster[seq.name] = cluster
                if seq.is_ref:
                    reps.add(seq.name)
    missing = set(ids) - set(id2cluster)
    if missing:
        raise ValueError(f"{path}: {len(missing)} ids missing from clusters")
    return id2cluster, reps


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
    parser.add_argument("tool", choices=["cdhit", "psi", "mmseqs"])
    parser.add_argument("--reps", action="store_true",
                        help="keep only the representative of each cluster")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cdhit-bin", default="cd-hit-est")
    parser.add_argument("--psi-script", default="psi-cd-hit.pl")
    parser.add_argument("--mmseqs-bin", default="mmseqs")
    args = parser.parse_args()

    th = args.th / 100 if args.th > 1 else args.th

    df = pd.read_csv(DATA / "archiveII.csv", index_col="id")
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    CLSTR_DIR.mkdir(parents=True, exist_ok=True)
    fasta = WORK_DIR / "archiveII.fasta"
    write_fasta(df, fasta)

    if args.tool == "cdhit":
        clstr = run_cdhit(fasta, th, args.cdhit_bin)
    elif args.tool == "psi":
        clstr = run_psi_cdhit(fasta, th, args.psi_script)
    else:
        clstr = run_mmseqs(fasta, th, args.mmseqs_bin)
    id2cluster, reps = load_clusters(clstr, df.index)

    tool = args.tool
    if args.reps:
        id2cluster = {i: c for i, c in id2cluster.items() if i in reps}
        tool += "rep"

    random.seed(args.seed + round(th * 100))
    split = make_splits(id2cluster, len(id2cluster))
    out_dir = WORK_DIR / "mmseqs_splits" if args.tool == "mmseqs" else OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"clustering_{tool}_{round(th * 100)}_folds_split.csv"
    split.to_csv(out_csv, index=False)
    print(f"{out_csv}  ({len(id2cluster)} sequences, "
          f"{len(set(id2cluster.values()))} clusters)")


if __name__ == "__main__":
    main()
