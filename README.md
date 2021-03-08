# 2021_pms_mt

This repository contain scripts utilized in the comparative analysis of ...

### `introns_identity`
R script utilized to align introns and plot a heatmap of %identity. Intron sequences are in the file `intron_sequences.tar.gz`. Untar this file and `intron_sequences` will be created, containing several subdirectories, one per intron insertion sites. Inside each one, there are individual fasta files of the intron inserted in the respective site. Species acronym is in the fasta file name: `enec`: *E. necator*, `episi`: *E. pisi*, `gcic`: G. cichoracearum, `bgram`: *B. graminis* f. sp. *tritici*.

Untar: `tar -xzvf intron_sequences.tar.gz`, and them run the script `pairwise_alignments.R` by calling snakemake:
```bash=
snakemake -j 1 --use-conda
```
The required R packages: `seqinr`, `Biostrings`, `ggplot2`, `reshape2` and `viridis` will be installed automatically.

The heat map is plotted to `summary_table_identity.pdf`. The file `summary_table_identity.txt` contains the %identity values for each pairwise alignment.

Because introns can have different sizes, the %identity is calculated with in `local-global` type, which means that the subject (sequence 2), must be entirely aligned (global), whereas subject (sequence 1), does not need to be aligned entirely (local). Size of subject is always at most the size of pattern. Alignment strategy can be changed in the `pairwiseAlignment()` within `getIdentity()` function.

###### tags: `readme` `mt`
