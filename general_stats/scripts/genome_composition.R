
library(ggplot2)

args   = commandArgs(trailingOnly=TRUE)
outpdf_introns_length  = args[1]
outpdf_genome_compos   = args[2]
outpdf_genome_shared_uniq = args[3]

introns = read.table("data/introns_length.txt", header = T)

introns$Species = factor(introns$Species, levels = rev(c("Enec", "Episi", "Gcic", "BGT")))

pdf(outpdf_introns_length, width = 7, height = 5)
ggplot(introns, aes(Species, Length)) +
  geom_boxplot(fill="grey80", width=0.9, lwd=0.8) +
  #geom_jitter(size=2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9, binwidth = 150, fill="grey50") +
  coord_flip() +
  theme_classic()
dev.off()

#=======

compos = read.table("data/genome_composition.txt", header = T)
compos$Species = factor(compos$Species, levels = rev( c("Enec", "Episi", "Gcic", "BGT")) )
compos$Category = factor(compos$Category, levels = c("Intergenic", "Intronic", "Exonic"))

pdf(outpdf_genome_compos,width = 9, height = 6)
ggplot(compos, aes(fill=Category, y=Length, x=Species)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  theme_classic()
dev.off()


#=======


compos = read.table("data/genome_shared_uniq.txt", header = T)
compos$Species = factor(compos$Species, levels = c("Enec", "Episi", "Gcic", "BGT"))
compos$mtDNA = factor(compos$mtDNA, levels = c("Unique", "Shared"))

pdf(outpdf_genome_shared_uniq, width = 4, height = 6)
ggplot(compos, aes(fill=mtDNA, y=Bases, x=Species)) + 
  geom_bar(position="fill", stat="identity")
  #coord_flip() +
  #scale_fill_brewer(palette="Paired")
dev.off()
