import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import pairwise2
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parents[4] / "data" / "bpRNA"
OUT_PATH = DATA_PATH / "distances" / "chunks"  # merged by compute_sequence_similarity.ipynb

# first row of each chunk (and the number of sequences at the end), obtained in the notebook compute_sequence_similarity
CHUNKS = [0, 553, 1124, 1714, 2326, 2962, 3625, 4319, 5049, 5820, 6641, 7524, 8483, 9544, 10749, 12178, 14040, 18537]

def valid_sequence_ids(df):
    valid_ids = []
    for i in range(len(df)):
        save = True
        for l in df.iloc[i].sequence:
            if l not in ["A", "C", "G", "U"]: 
                save = False
        if save:
            valid_ids.append(df.iloc[i].name)
    print(f"Total IDs: {len(df)} , Valid IDs: {len(valid_ids)}")
    df_valid = df.loc[valid_ids]
    return df_valid

def get_match_dict(symbols, Match, Mismatch): 
    match_dic = {}
    for simbolA in symbols:
        for simbolB in symbols:
            if simbolA==simbolB:
                match_dic[(simbolA, simbolB)] = Match
            else:
                match_dic[(simbolA, simbolB)] = Mismatch
    return match_dic
    
def compute_similarity(col_seq, row_seq, sim_conf):
    alignm = pairwise2.align.localds(col_seq, row_seq, sim_conf[0], sim_conf[1], sim_conf[2], one_alignment_only=True)[0]
    max_len = max(len(alignm.seqB),len(alignm.seqA))
    IDscore = sum(a==b for a,b in zip(alignm.seqA, alignm.seqB)) / max_len
    return IDscore

def expected_num_comparisons(n, start, end):
    return sum(n - i - 1 for i in range(start, end))

def get_sim_conf():
    Match = 1
    Mismatch = -2
    GapOpen=-5
    GapExtend=-2
    simbolos=["A","C","G","U"]
    match_dic = get_match_dict(simbolos, Match, Mismatch)
    sim_conf = match_dic, GapOpen, GapExtend
    return sim_conf

def compute_slice(df, start, end):
    ids = df.index.to_list()
    seqs = df["sequence"].astype(str).to_dict()
    sim_conf = get_sim_conf()

    row_ids = ids[start:end]
    sims = pd.DataFrame(0.0, index=row_ids, columns=ids, dtype=np.float32)

    total = expected_num_comparisons(len(ids), start, end)

    with tqdm(total=total, desc=f"Rows {start}:{end}") as pbar:
        for i, row_id in enumerate(ids[start:end], start=start):
            seq_i = seqs[row_id]

            for col_id in ids[i + 1:]:
                sims.loc[row_id, col_id] = compute_similarity(
                    seq_i,
                    seqs[col_id],
                    sim_conf,
                )
                pbar.update()
    return sims

def main():
    parser = argparse.ArgumentParser() 
    parser.add_argument("--chunk", type=int, required=True)
    args = parser.parse_args()
    OUT_PATH.mkdir(parents=True, exist_ok=True)

    data_file = DATA_PATH / "bpRNA.csv"
    df = pd.read_csv(data_file, index_col="id")
    df = valid_sequence_ids(df)

    start = CHUNKS[args.chunk]
    end = CHUNKS[args.chunk + 1]
    out = compute_slice(df, start, end)
    out_path = OUT_PATH / f"bprna_seqsim_{end}.h5"
    out.to_hdf(
        out_path,
        key="seqsim",
        mode="w",
        format="fixed",
    )

if __name__ == "__main__":
    main()