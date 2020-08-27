
# read gene file from command line
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
genes <- peach1_genes <- scan(filename, what="character", strip.white = TRUE, sep = ",")
genes <- gsub("\"","",genes) # remove the quotations marks in gene names
genes <- gsub("\\\\","\\",genes) # remove the backlashes
genes <- gsub("\\{ | |\\}","",genes) # remove the { } marks in gene names
# read annotation files
annotation <- read.csv("/staton/projects/apricot/peach_genome/Ppersica_298_v2.1.annotation_info.txt", header = T, stringsAsFactors = F, sep="\t",row.names = 1)
annotation <- annotation[!duplicated(annotation[,1]),] # remove the isoforms

# get the gene functions for gene list
genes_annot <- annotation[annotation$locusName %in% genes, c(1,4,9:12)] # only keep the pfam, GO and Arabidopsis hit columns

# output the gene_annot file as csv format
write.csv(genes_annot, paste0(filename,"annot.csv"))