# libs to read fasta and align sequences
library(seqinr)
library(Biostrings)
# libs to plot heatmap
library(ggplot2)
library(reshape2)
library(viridis)


#===== Reading data ======
# this table must pre-exist. Rows contain pairwise combinations of the species, and columns name of the introns (also name of the directories where intron sequences are)
# only row names and colnames are important. The table will be filled with %identity between pairwise alignments
introns = read.table("summary_table_identity.txt", sep = '\t', header = T)
#========================



#===== Functions =========
# function to return a single string given an array of characters
array2string <- function(arr){
  strg = toupper( c2s(arr) )
  
  return(strg)
}

# function that returns %identity (or NA) given two filenames. It will try to reads the fasta files. If successfull, then align sequences and return %identity. NA otherwise.
getIdentity <- function(sp1_filename, sp2_filename){
  
  # if either file is nonexistent, return NA
  if( identical(sp1_filename, character(0)) || identical(sp2_filename, character(0)) ) {
    ident = NA
  }else{
    #read fasta files
    seq1 = read.fasta(sp1_filename)
    seq2 = read.fasta(sp2_filename)
    
    #get DNA string from character array
    seq1 = array2string(seq1[[1]])
    seq2 = array2string(seq2[[1]])
    
    # make sure s1 is the longest of the two sequences
    if( nchar(seq1) < nchar(seq2) ){
      aux = seq1
      seq1 = seq2
      seq2 = aux
    }
    
    # aligning sequences. Alignment should cover the entire s2 (shortest) sequence. This way,
    # penalty caused by difference in length will be minimized 
    align = pairwiseAlignment(seq1, seq2, scoreOnly = F, type = "local-global")
    ident = round(pid(align), digits = 1)
  }
  
  return(ident)
}
#=========================

# go in the main folder containing all intron sequences
setwd("intron_sequences")

# for each intron. Folder named as the intron name, must exist.
for(intron in colnames(introns) ){
  # go in the folder containing the respective intron sequences
  setwd(intron)
  
  # for each pair, try to find fasta files and perform alignments
  for(pair in rownames(introns)){
    
    # get all files, i.e. all fasta files, individual sequences. There's nothing else inside, other than fasta files
    fil = list.files()
    
    # getting acronym of each species. They match fasta file names. 
    species1 = strsplit(pair, split = '_')[[1]][1]
    species2 = strsplit(pair, split = '_')[[1]][2]
    
    # Find the fasta file of each species. For some there will be one. In this case, `character(0)` is returned instead of the filename.
    sp1_filename = fil[grep(pattern = species1, fil)]
    sp2_filename = fil[grep(pattern = species2, fil)]
    
    # calculate identity between sequences. NA if either file does not exist
    introns[pair, intron] = getIdentity(sp1_filename, sp2_filename)
  }
  
  # go back to previous dir
  setwd("../")
    
}

# go back once more, where this script is located
setwd("../")

#------------------
# at this point, the table with identity values is constructed. Now plot a heatmap
#------------------

# add a column with species pairs (=rownames)
introns$pair=rownames(introns)

# melt table to long format
tab = melt(introns, variable.name =  'pair')

# rename columns
colnames(tab) = c("pair", "intron", "identity")

# set levels so they appear in the desire order in the plot 
tab$pair = factor(tab$pair , levels = rev( c("enec_episi", "enec_gcic", "enec_bgram", "episi_gcic", "episi_bgram", "gcic_bgram")) )

# plot heatmap
pdf("summary_table_identity.pdf", width = 12, height = 2.2)
ggplot(tab, aes(intron, pair, fill= identity)) + 
  geom_tile() +
  #scale_fill_gradientn(colors = colors, breaks = b, labels = format(b))
  scale_fill_viridis(discrete=FALSE, na.value = 'white', option = "E", limits = c(40, 100)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



