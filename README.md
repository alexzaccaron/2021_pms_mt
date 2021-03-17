# 2021_pms_mt

This repository contain scripts utilized in the comparative analysis of ...


### `general_stats`
In this directory there's a Snakefile that downloads the mt genomes and estimates some basic stats.

Important output files are `{sample}_fasta.stats`, which contains four numbers: length, GC fraction, GC-skew, and AT-skew. The file `{sample}_repeats.txt` contains the estimate total base pairs that are likely repeats. This number is estimated based on self-blastn searches. Other useful files are `{sample}_{sample}_blastn.out`, which contains pairwise blastn searches results.

There's a `scripts` subdirectory that contains R scripts to generate some plots. These plots will be in the sundirectory `plots`. One is a syntheny plot based on the pairwise blastn searches. Another is a boxplot that shows intron size, and the other shows total bases that correspond to exons, introns and intergenic regions. Note the subdirectory `data` that supplies input files to the R scripts to generate the plots.

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


### `small_tree`
In this directory there's a snakefile to construct a small tree with `IQTREE` based on mt genes.

There's a `main_table.csv` with accession numbers of the proteins and info about the species. However, snakemake actually uses `main_table_melt.txt`, which is a melted version of `main_table.csv`. The R script `scripts/melt_main_table.R` can melt `main_table.csv`:
```r=
scripts/melt_main_table.R main_table.csv main_table.csv
```

Snakemake will perform all steps, from downloading the sequences with `efetch`, align sequences with `MAFFT` and construct a tree with `IQTREE`. The concatenated and thee tree files will be in the subdirectory `concatenated_alignment/`.

###### tags: `readme` `mt`
