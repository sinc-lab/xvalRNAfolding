# Homology-aware cross-validation strategies for generalization assessment in RNA structure prediction

This repository contains the data and source code for the manuscript *Homology-aware cross-validation strategies for generalization assessment in RNA structure prediction*, by L.A. Bugnon, G. Kulemeyer, M. Gerard, L. Di Persia, G. Stegmayer and D.H. Milone, 2026 (under review). Research Institute for Signals, Systems and Computational Intelligence, [sinc(i)](https://sinc.unl.edu.ar/).

In this work we revise existing cross-validation strategies for RNA secondary structure prediction: random k-fold, clustering fold and family fold. 

We explain and analyze in detail the advantages and disadvantages of each one of them, additionally proposing two novel evaluation methodologies for testing: human-learned fold and similarity fold. 

All validation strategies are applied to state-of-the-art methods for RNA secondary structure prediction and comparative results are analyzed.

## Cross-validation methodologies for RNA secondary structure prediction revised

<div align="center">
  
  ![](figures/FIG1.png)

</div>

**a)** Random **k-fold**: (top) the complete dataset is randomly divided into k groups, in a fold a group is used for testing (red) and the rest for training (gray). 

In **cluster-fold** (bottom), the complete dataset is split into clusters of similar sequences, in a fold a subset of them is assigned to the training set, and the rest is used for testing. 

**b) fam-fold**: the illustration has 6 structural families (triangles, lines, etc.). At each fold, one complete family is left out and used only for testing, while all the others are used for training.

**c) hl-fold**: each fold, in training, has all the sequences for which the mean F1 of four classical methods (RNAfold, RNAstructure, LinearFold and LinearPartition) is above a threshold, and the rest of the sequences are used for testing. Several thresholds define the folds. 

**d) sim-fold**: groups of increasing sequence simlarity are built, then inside each group many random train/test folds can be sampled.


## Performance comparison for RNA secondary structure prediction methods according to different cross-validation strategies
<div align="center">
  
  ![](figures/FIG2.png)
  
</div>

**a) k-fold**: (top) distribution of test-to-train structural distances of each fold; (bottom) median F1 for each prediction method (DL methods in blue, classical methods in red). Gray bars show the proportion of test to train sequences. 

**b) cluster-fold**: (top) test-to-train distance distributions; (bottom) structure prediction performance. 

**c) fam-fold**: (top) test-to-train distance distributions; (bottom) structure prediction performance. 

**d)** Analysis of connection distances for tmRNA family used in test set of fam-fold strategy: (top) distribution of distances between bases of canonical pairs in the reference structures of training (gray) and testing (blue) partitions, (bottom) distribution of connection distances distribution for canonical base pairs GU, AU and GC. 

**e) hl-fold**: (top) test to train distance distributions; (bottom) median F1 performance for each prediction method, with the 95% confidence interval shaded behind each trend line. 

**f) sim-fold**: (top) test-to-train distance distributions; (bottom) structure prediction performance. wAUC: area under the performance curve for each method, weighted by the difficulty of the partition.

## Data

`data/` has one folder per dataset, with the sequences and structures (`<dataset>.csv`) and one split file per cross-validation strategy. Each split file lists, for every fold, the partition (`train`, `valid` or `test`) of each sequence `id`; sim-fold has one file per similarity threshold.

The all-vs-all distance matrices can be download into `data/<dataset>/distances/`:

| File | Dataset | Content | Size |
|---|---|---|---|
| [seqsim_f_all.h5](https://drive.google.com/file/d/1qLPI9ePyUksuJ7eu3Q8PV-kIMhC5NXDC/view) | ArchiveII | sequence similarity (IDscore) | 127 MB |
| [rnadist_f_all.h5](https://drive.google.com/file/d/1qVk5IRY8_JgVFIuPVZg4MCsYNVm-kNGe/view) | ArchiveII | structural distance | 127 MB |
| [bprna_seqsim_full_ACGU.h5](https://drive.google.com/file/d/1YjRZgRwjp-MdpLjXdHaIEk62HbGIXS4r/view) | bpRNA | sequence similarity (IDscore) | 1.4 GB |
| [bprna_rnadist_full.h5](https://drive.google.com/file/d/1jPHR_eBI6dB3CCCEjz7FbAVGiNxrujPV/view) | bpRNA | structural distance | 1.4 GB |

The all-vs-all blastn identity tables of ArchiveII (`data/archiveII/distances/archiveII-blastn-allvall[-nofilter].tsv`) are built with `python src/check_split_leakage.py data/archiveII/clustering_folds_split.csv 80 [--no-filter]`.

## Reproducible research
- [This notebook](https://colab.research.google.com/github/sinc-lab/xvalRNAfolding/blob/main/src/Figure_2_dist_dist.ipynb) shows the train-to-test structural distance distributions of each cross-validation strategy. 

- [This notebook](https://colab.research.google.com/github/sinc-lab/xvalRNAfolding/blob/main/src/Figure_2_strips_plot.ipynb) shows performance comparison of prediction methods among different cross-validation strategies. The per-sequence $F_1$ scores it plots are computed from `data/` and `predictions/` by `src/summaries.py` on first use and cached in `results/` (not part of the repository); `python src/summaries.py` builds them all at once.

- [This notebook](https://colab.research.google.com/github/sinc-lab/xvalRNAfolding/blob/main/src/Figure_2_connections_dist.ipynb) plots the distribution of train and test connection distances, and the detail for base pairs GU, AU and GC in tmRNA family.

Supplementary material:

- Figures S1 and S2 (critical difference diagrams, stratified k-fold) and S7 (bpRNA) are made by the Figure 2 performance notebook above.
- [This notebook](https://colab.research.google.com/github/sinc-lab/xvalRNAfolding/blob/main/src/Figure_S3_cluster_analysis.ipynb): Figure S3, multi-threshold cluster-fold.
- [This notebook](https://colab.research.google.com/github/sinc-lab/xvalRNAfolding/blob/main/src/Figure_S5.ipynb): Figure S5, per-sequence correlations between methods.
- [This notebook](https://colab.research.google.com/github/sinc-lab/xvalRNAfolding/blob/main/src/Figure_S6.ipynb): Figure S6, per-sequence $F_1$ against the maximum similarity to the training set at sim90.
- [This notebook](https://colab.research.google.com/github/sinc-lab/xvalRNAfolding/blob/main/src/Tables_S1_S2_S3.ipynb): Tables S1, S2 and S3.

