setwd("C:/Users/raymo/Downloads/AS/PSI") # set the working folder
library(plyr)

# read peach, apricot, and cherry transcript raw counts and reciprocal best hit file
<<<<<<< HEAD
peach1_data <- read.csv("data/peach_transcript_counts.csv", header = T, row.names = 1)
apricot1_data <- read.csv("data/apricot_transcript_counts.csv", header = T, row.names = 1)
peach2_data <- read.csv("data/peach2_transcript_counts.csv", header = T, row.names = 1, sep = "\t")
apricot2_data <- read.csv("data/apricot2_transcript_counts.csv", header = T, row.names = 1,sep="\t")

rownames(peach1_data) <- gsub("MSTRG.\\d+;","",rownames(peach1_data))
rownames(apricot1_data) <- gsub("MSTRG.\\d+;","",rownames(apricot1_data))
rownames(peach2_data) <- gsub("MSTRG.\\d+;","",rownames(peach2_data))
rownames(apricot2_data) <- gsub("MSTRG.\\d+;","",rownames(apricot2_data))


# combine transcript counts for all three
common_genes <- Reduce(intersect, list(rownames(apricot1_data), rownames(peach1_data), rownames(apricot2_data), rownames(peach2_data)))
combined_data <- data.frame(row.names = common_genes) # create an empty dataframe with common genes
combined_data <- transform(merge(combined_data, peach1_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, apricot1_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, peach2_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, apricot2_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
=======
peach_data <- read.csv("p_persica_psi_isoform.psi", header = T, row.names = 1, sep="\t")
peach_2_data <- read.csv("p_persica_2_psi_isoform.psi", header = T, row.names = 1, sep="\t")
apricot_data <- read.csv("p_armeniaca_psi_isoform.psi", header = T, row.names = 1, sep="\t")
apricot_2_data <- read.csv("p_armeniaca_2_psi_isoform.psi", header = T, row.names = 1, sep="\t")
RBH <- read.table("reciprocal_best_hits.txt")

#Removing ‘MSTRG.XX;’ in front of the transcript ID
rownames(peach_data) <- gsub("MSTRG.\\d+;","",rownames(peach_data))
rownames(peach_2_data) <- gsub("MSTRG.\\d+;","",rownames(peach_2_data))
rownames(apricot_data) <- gsub("MSTRG.\\d+;","",rownames(apricot_data))
rownames(apricot_2_data) <- gsub("MSTRG.\\d+;","",rownames(apricot_2_data))


# combine transcript counts for all three
common_genes <- Reduce(intersect, list(rownames(peach_data), rownames(peach_2_data), rownames(apricot_2_data), rownames(apricot_data)))
combined_data <- data.frame(row.names = common_genes) # create an empty dataframe with common genes
combined_data <- transform(merge(combined_data, peach_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, peach_2_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, apricot_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, apricot_2_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data[is.na(combined_data)] <- 0
>>>>>>> 5c325ddfcc40be486b2d4d0ac5bd26d31cba00c6

# normalized the gene counts for combined dataset ???

# Calculate JSD, will be using the package philentropy
#install.packages("philentropy") # only need to run once for installation
library(philentropy)
# it is taking a matrix with rows as samples and columns as genes so we transpose the data and make it a matrix
JSD_matrix <- JSD(as.matrix(t(combined_data)), unit = "log2", est.prob = "empirical")
rownames(JSD_matrix) <- colnames(combined_data) # change the column and row names of JSD matrix
colnames(JSD_matrix) <- colnames(combined_data)
normalized_JSD <- 1 - sqrt(JSD_matrix) # follow the method in the science paper
heatmap(normalized_JSD,col = hcl.colors(256), xlab = "1-sqrt(JSD)", na.rm=TRUE)

