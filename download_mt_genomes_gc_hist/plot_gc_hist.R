
library(seqinr)
library(knitr)

#GIVE THE SEQUENCE FILE (SINGLE SEQUENCE IN FASTA) AND A NAME TO PLOT
# Rscript plot_gc_hist.R myseq.fasta myseq.pdf

#------
args = commandArgs(trailingOnly=TRUE)
sequence = args[1]
plotfile = args[2]
#------


fasta = read.fasta(sequence)
fastaseq = fasta[[1]]

starts = seq(1, length(fastaseq)-200, by = 200)

n = length(starts)
GCarray = NULL

for(i in 1:n){
  chunk <- fastaseq[starts[i]:(starts[i]+199)]
  GCarray = append(GCarray, GC(chunk))
}

pdf(plotfile, width = 5, height = 4)
hist(GCarray, seq(0,1,by=0.02), freq=FALSE, col="grey80", ylim=c(0,10), xlim=c(0,1),  main=NA, axes=F, xlab="", ylab="")
axis(1)
axis(2, las=2)
box()
lines(density(GCarray), lwd=2.5, col="red")
dev.off()

plot_crop(plotfile)
