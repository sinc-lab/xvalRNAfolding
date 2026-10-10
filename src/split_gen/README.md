# Split generation

Code that generated the cross-validation splits in `data/` and the distance matrices used to
build or describe them. Except `sim_folds.ipynb`, paths are set for the machine where each one
was run.

`sim_folds.ipynb` generates the sim-folds of both datasets into `data/<dataset>/sim_folds_split/`,
from the dataset and its sequence similarity matrix (`data/<dataset>/distances/`, see the main
README). ArchiveII sim-folds are 10 random partitions with 20% of the sequences for test each;
bpRNA sim-folds are 10 folds with mutually exclusive test partitions.

## ArchiveII

| Split | Code |
|---|---|
| k-fold, fam-fold | `dir2list_splits_archiveii.ipynb`: converts the partitions of the dl-rna release (Szikszai et al., 2022) to CSV |
| cluster-fold, multi-threshold cluster-fold | `make_clustering_folds.py` (cd-hit-est / psi-cd-hit / MMseqs2) |
| hl-fold | `hl_subsets_20260804.ipynb` (F1mean and F1rnafold versions) |
| stratified k-fold | `stratified_kfold.ipynb` |
| structural distance matrix (`rnadist_f_all.h5`) | `rnadist_all.ipynb` (RNAdistance) |

## bpRNA

| Split | Code |
|---|---|
| k-fold | `1_random_kfold.ipynb` |
| cluster-fold | `2_clustered_kfold.ipynb` (from a cd-hit-est 80% `.clstr`, not included) |
| fam-fold | `3_famfold.ipynb` (bpRNA-new families from an Rfam 14.2 mapping, not included) |
| hl-fold | `hl_subsets_20260804_f1mean.ipynb` |
| sequence similarity matrix | `seqsim/seqsim_bpRNA.py` (by row slices), `seqsim/compute_sequence_similarity.ipynb` |

## Reproducibility

Checked by rerunning each generator's logic on the files in `data/`:

- Reproduced exactly: ArchiveII stratified k-fold, the sim-folds of both datasets (bpRNA files
  also have a `split_type` column), the 30 multi-threshold cluster-folds from their `.clstr` files,
  and bpRNA k-fold.
- hl-folds (both datasets): the test partitions are reproduced; the train/valid partitions are
  not, because the validation set is sampled from a list built from a Python `set`, whose order
  changes between runs.
- Not reproducible: the original ArchiveII cluster-fold (80%) and bpRNA cluster-fold, whose
  clusters were shuffled without a seed; `clustering_psi_45`, whose `.clstr` was not kept.
- Not checked: ArchiveII k-fold and fam-fold (converted from the dl-rna release), bpRNA fam-fold
  (needs an Rfam 14.2 mapping file not included).

## Missing

- ArchiveII: the sequence similarity matrix (`seqsim_f_all.h5`).
- bpRNA: the structural distance matrix.
