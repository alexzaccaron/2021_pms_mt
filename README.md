# 2021_pms_mt

This repository contain scripts utilized in the comparative analysis of ...

### `download_mt_genomes_gc_hist`
Contains a `Snakefile` that downloads with `efetch` >700 fungal mt genomes and makes a GC content histogram for each mt genome. GenBank accession numbers and other info of the mt genomes are in the TSV file `table_mt_genomes.txt`. Go to this directory and run:

```bash=
snakemake -j 1 --use-conda
```

It might take a while to download all of them. You can try to adjust the number of threads (`-j`) to make it faster.

Two folders are produced: `mtgenomes`, which contains all mt genomes in fasta file, and `mtgenomes_gc_hists` with the GC histograms in PDF.


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
