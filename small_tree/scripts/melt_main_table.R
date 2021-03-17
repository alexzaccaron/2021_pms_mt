#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

main_table_fname        = args[1]
main_table_melted_fname = args[2]

# reshape package to call melt() functiom
if (!require("reshape")) install.packages("reshape", repos = "http://cran.us.r-project.org")

# reading main table with species info and accession numbers
main_table = read.csv(main_table_fname, header = T, fill = T)

# I will select these fields
fields = c("atp6", "nad1","nad2","nad3","nad4","nad4L","nad5","nad6","cox1","cox2","cox3","cob","TaxID","Organism")

# select the fields
main_table = main_table[,fields]

# melt the table
main_table_melt = melt(main_table, id.vars=c("Organism", "TaxID"))

#adding a new column to store filename to output sequences
main_table_melt$filename = paste(main_table_melt[,'TaxID'],main_table_melt[,'variable'], main_table_melt[,'value'], sep = "_")

# giving names to columns
colnames(main_table_melt) = c("organism", "taxID", "gene", "accession", "filename")

# write melted table
write.table(main_table_melt, main_table_melted_fname, sep ='\t', col.names = T, row.names = F, quote = F)
