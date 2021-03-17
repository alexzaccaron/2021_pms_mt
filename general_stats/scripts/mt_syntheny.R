
# get arguments. Output PDF should be passed
args   = commandArgs(trailingOnly=TRUE)
outpdf = args[1]


#color to draw blast ribbons
colfunc <- colorRampPalette(c("black", "green", "yellow", "red"))
pallet = colfunc(100)



#======= Functions ==============

#function to add features (genes) at specific height (y). The format is gff with an additional column with color
add_features <- function(features, y){
  #draw genes first
  sub_features = subset(features, feature=="gene")
  for(i in 1:nrow(sub_features)){
    color = as.character(sub_features[i,'color'])
    rect(sub_features[i,'start'], y+0.2, sub_features[i,'end'], y+1.8, lwd = 0.4, col = color)
  }
  
  #now draw exons in color
  sub_features = subset(features, feature!="gene")
  for(i in 1:nrow(sub_features)){
    color = as.character(sub_features[i,'color'])
    rect(sub_features[i,'start'], y+0.2, sub_features[i,'end'], y+1.8, lwd=0.1, col = color)
  }
  
}


add_blast_ribbons <- function(blast_out, topy, bottomy, pallet){
  #adding blast ribbons
  #for each HSP...
  for(i in 1:nrow(blast_out)){
    #get HSP start and end in query and subject
    ref_start = blast_out[i,'sstart']
    ref_end   = blast_out[i,'send']
    qry_start = blast_out[i,'qstart']
    qry_end   = blast_out[i,'qend']
    #get color based on the % identity. It needs to be round
    color = pallet[round(blast_out[i,'pident'])]
    
    #draw the ribbon and a polygon
    polygon(x = c(ref_start, qry_start, qry_end, ref_end), y = c(topy, bottomy, bottomy, topy), col = color, border = NA)
  }
}
#=============================



#======== READING DATA ===========
# gff files with an additional column with color
enec_gff = read.table("data/gff_files/enec_genes.gff", comment.char="&", sep = '\t', col.names = 
                        c("scaff", "note", "feature", "start", "end", "score", "strand", "frame", "tags", "color"))
episi_gff = read.table("data/gff_files/episi_genes.gff", comment.char="&", sep = '\t', col.names = 
                        c("scaff", "note", "feature", "start", "end", "score", "strand", "frame", "tags", "color"))
gcic_gff = read.table("data/gff_files/golcic_genes.gff", comment.char="&", sep = '\t', col.names = 
                        c("scaff", "note", "feature", "start", "end", "score", "strand", "frame", "tags", "color"))
bgt_gff = read.table("data/gff_files/bgt_genes.gff", comment.char="&", sep = '\t', col.names = 
                       c("scaff", "note", "feature", "start", "end", "score", "strand", "frame", "tags", "color"))


#blast output table (format 6)
blast1 = read.table("MT880589_MT880588_blastn.out",
                       col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
blast2 = read.table("MT880590_MT880589_blastn.out",
                               col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
blast3 = read.table("MT880591_MT880590_blastn.out",
                               col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

#================================


#Size of the genomes
enec_size   = 188577
episi_size  = 188623
golcic_size = 332165
bgt_size    = 109800


#======= ADDING OFFSET =====
#this is just a shift to the right, aligning both vertically
# offset = enec_size/2 - blgt_size/2
# blgt_gff$start = blgt_gff$start + offset
# blgt_gff$end = blgt_gff$end + offset
# blast_out$qstart = blast_out$qstart + offset
# blast_out$qend = blast_out$qend + offset
#==========================



#======== MAIN =======
pdf(outpdf, width = 14, height = 6)
#making an empty plot
plot(1, type="n", xlim=c(0,350000), ylim=c(0,25), axes=F, xlab="", ylab="")
axis(side=1, pos=3)

#ading segment to represent the genomes
segments(0,20,enec_size,20, lwd=3)
segments(0,15,episi_size,15, lwd=3)
segments(0,10,golcic_size,10, lwd=3)
segments(0,5,bgt_size,5, lwd=3)

#adding features. Adjust height (second argument) if needed
add_features(enec_gff, 19)
add_features(episi_gff, 14)
add_features(gcic_gff, 9)
add_features(bgt_gff, 4)


add_blast_ribbons(blast1, 19, 16, pallet)
add_blast_ribbons(blast2, 14, 11, pallet)
add_blast_ribbons(blast3, 9, 6, pallet)


#= LEGEND =
step = 500
height = 1
x1=285000
y1=20
for(i in 1:100){
  x2=x1+step
  y2=y1+height
  rect(x1, y1, x2, y2, lwd = NA, col=pallet[i])
  x1=x2
}
#=======

dev.off()
#=================

