"""Check test-train similarity leakage of a fold split file.

Usage: python check_split_leakage.py <split_csv> <th> [--top N] [--metric cdhit]
                                    [--no-filter]
    split_csv: any *_folds_split.csv (columns fold_number, partition, id)
    th: identity threshold to test against (e.g. 0.8 or 80)
Reports, per fold, the test sequences having >= th identity with any train
sequence. Metrics: blast (default) = identical bases / shorter sequence length
from all-vs-all blastn (best HSP per pair, cached in
data/archiveII/distances/; --no-filter uses an uncensored table
built without blast's 50% identity / e-value cutoffs); cdhit = cd-hit-est-2d with the original
run flags (cd-hit's own identity; th >= 0.8 only, no pair listing). Note:
cd-hit only guarantees member-to-representative identity, so greedy
clusterings are expected to show leakage (see cluster_rationale.md).
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "archiveII"
WORK_DIR = REPO / "work" / "clustering"
BLAST_TSV = DATA / "distances" / "archiveII-blastn-allvall.tsv"
BLAST_TSV_NOFILTER = DATA / "distances" / "archiveII-blastn-allvall-nofilter.tsv"


def write_fasta(df, path):
    with open(path, "w") as f:
        for seq_id, row in df.iterrows():
            f.write(f">{seq_id}\n{row.sequence}\n")


def blast_tsv(no_filter=False):
    return BLAST_TSV_NOFILTER if no_filter else BLAST_TSV


def run_blast(df, no_filter=False):
    """All-vs-all blastn, cached. The default table keeps blast's own
    -perc_identity 50 / -evalue 1e-6 cutoffs (every pair below ~50% identity is
    missing, i.e. censored); --no-filter drops them so the low-identity tail of
    the distribution is observable, at the cost of a much larger table."""
    tsv = blast_tsv(no_filter)
    if tsv.exists():
        return
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    tsv.parent.mkdir(parents=True, exist_ok=True)
    fasta = WORK_DIR / "archiveII.fasta"
    write_fasta(df, fasta)
    db = WORK_DIR / "archiveII-blastdb"
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
                    "-out", str(db)], check=True, capture_output=True)
    cutoffs = (["-evalue", "10", "-perc_identity", "0"] if no_filter
               else ["-evalue", "1e-6", "-perc_identity", "50"])
    with open(tsv, "w") as out:
        subprocess.run(
            ["blastn", "-task", "blastn", "-query", str(fasta), "-db", str(db),
             "-outfmt", "6 qseqid sseqid nident qlen slen",
             *cutoffs,
             "-max_target_seqs", "5000", "-dust", "no",
             # BLAST+ >= 2.16 rejects -num_threads 0 (the original table was
             # built with an older version that read it as "all cores")
             "-num_threads", str(os.cpu_count() or 1)],
            check=True, stdout=out)


def load_identities(no_filter=False):
    best = {}
    for line in open(blast_tsv(no_filter)):
        q, s, nident, qlen, slen = line.split()
        if q == s:
            continue
        pair = (q, s) if q < s else (s, q)
        ident = int(nident) * 100 / min(int(qlen), int(slen))
        if ident > best.get(pair, 0):
            best[pair] = ident
    return best


def check_cdhit_metric(df, split, th, cdhit_bin):
    print(f"{'fold':<5} {'test':<6} {'leaking':<8} {'%':<6}")
    with tempfile.TemporaryDirectory() as tmp:
        tr_fa, te_fa, out = Path(tmp) / "tr.fa", Path(tmp) / "te.fa", Path(tmp) / "out"
        for fold, part in split.groupby("fold_number"):
            test = part[part.partition == "test"].id
            write_fasta(df.loc[part[part.partition == "train"].id], tr_fa)
            write_fasta(df.loc[test], te_fa)
            subprocess.run(
                [cdhit_bin, "-i", str(tr_fa), "-i2", str(te_fa), "-o", str(out),
                 "-c", f"{th / 100:g}", "-T", "0", "-d", "100", "-g", "1"],
                check=True, capture_output=True)
            kept = sum(1 for line in open(out) if line.startswith(">"))
            leak = len(test) - kept
            print(f"{fold:<5} {len(test):<6} {leak:<8} {leak * 100 / len(test):<6.1f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split_csv", type=Path)
    parser.add_argument("th", type=float)
    parser.add_argument("--top", type=int, default=0,
                        help="print the N highest-identity leaking pairs per fold")
    parser.add_argument("--metric", choices=["blast", "cdhit"], default="blast")
    parser.add_argument("--no-filter", action="store_true",
                        help="use (and build if missing) the uncensored blastn "
                             "table, without blast's 50%% identity / e-value cutoffs")
    parser.add_argument("--cdhit-2d-bin", default="cd-hit-est-2d")
    args = parser.parse_args()

    th = args.th if args.th > 1 else args.th * 100

    df = pd.read_csv(DATA / "archiveII.csv", index_col="id")

    if args.metric == "cdhit":
        split = pd.read_csv(args.split_csv)
        print(f"{args.split_csv.name} @ {th:g}% (cd-hit-est-2d identity)")
        check_cdhit_metric(df, split, th, args.cdhit_2d_bin)
        return

    run_blast(df, args.no_filter)
    best = load_identities(args.no_filter)

    neighbors = {}
    for (a, b), ident in best.items():
        if ident >= th:
            neighbors.setdefault(a, {})[b] = ident
            neighbors.setdefault(b, {})[a] = ident

    split = pd.read_csv(args.split_csv)
    print(f"{args.split_csv.name} @ {th:g}%")
    print(f"{'fold':<5} {'test':<6} {'leaking':<8} {'%':<6} {'pairs':<7} {'max id':<7}")
    for fold, part in split.groupby("fold_number"):
        test = part[part.partition == "test"].id
        train = set(part[part.partition == "train"].id)
        pairs = []
        leaking = set()
        for t in test:
            hits = {n: v for n, v in neighbors.get(t, {}).items() if n in train}
            if hits:
                leaking.add(t)
                pairs += [(v, t, n) for n, v in hits.items()]
        mx = max((p[0] for p in pairs), default=0)
        print(f"{fold:<5} {len(test):<6} {len(leaking):<8} "
              f"{len(leaking) * 100 / len(test):<6.1f} {len(pairs):<7} {mx:<7.1f}")
        for v, t, n in sorted(pairs, reverse=True)[:args.top]:
            print(f"      {v:6.1f}%  {t}  |  {n}")


if __name__ == "__main__":
    main()
