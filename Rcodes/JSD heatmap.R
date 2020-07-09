setwd("~/Desktop/Jiali/UTK/2020 Summer/AS/alternative_splicing_project/") # set the working folder
library(plyr)

# read peach, apricot, and cherry transcript raw counts and reciprocal best hit file
peach_data <- read.csv("data/peach_transcript_counts.csv", header = T, row.names = 1)
apricot_data <- read.csv("data/apricot_transcript_counts.csv", header = T, row.names = 1)
cherry_data <- read.csv("data/cherry_transcript_count.csv", header = T, row.names = 1)
RBH <- read.table("data/reciprocal_best_hits.txt")

rownames(cherry_data) <- gsub("MSTRG.\\d+;","",rownames(peach_data))
rownames(peach_data) <- gsub("MSTRG.\\d+;","",rownames(peach_data))
rownames(apricot_data) <- gsub("MSTRG.\\d+;","",rownames(apricot_data))

# Extract the transcripts in RBH from all three species
cherry_filter <- merge(cherry_data, RBH, by.x = "row.names", by.y = "V2") # add the peach gene corresponding to cherry
cherry_filter <- na.omit(cherry_filter) # only keep the cherry genes have peach hits
cherry_filter$V1 <- gsub("\\.p","",cherry_filter$V1) # remove the .p at the end of the peach gene names
rownames(cherry_filter) <- cherry_filter$V1

rownames(peach_data) <- gsub("_v2.0.a1","",rownames(peach_data)) # remove '_v2.0.a1' at the end of peach genes
peach_filter <- peach_data[rownames(peach_data) %in% cherry_filter$V1,] # keep the peach genes in RBH

rownames(apricot_data) <- gsub("_v2.0.a1","",rownames(apricot_data))
apricot_filter <- apricot_data[rownames(apricot_data) %in% cherry_filter$V1,]

# combine transcript counts for all three
common_genes <- Reduce(intersect, list(rownames(apricot_filter), rownames(peach_filter), cherry_filter$V1))
combined_data <- data.frame(row.names = common_genes) # create an empty dataframe with common genes
combined_data <- transform(merge(combined_data, peach_filter, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, apricot_filter, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, cherry_filter[,2:7], by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)

# normalized the gene counts for combined dataset ???


# Calculate JSD, will be using the package philentropy
install.packages("philentropy") # only need to run once for installation
library(philentropy)
# it is taking a matrix with rows as samples and columns as genes so we transpose the data and make it a matrix
JSD_matrix <- JSD(as.matrix(t(combined_data)), unit = "log2", est.prob = "empirical") 
rownames(JSD_matrix) <- colnames(combined_data) # change the column and row names of JSD matrix
colnames(JSD_matrix) <- colnames(combined_data)
normalized_JSD <- 1 - sqrt(JSD_matrix) # follow the method in the science paper
heatmap(normalized_JSD,col = hcl.colors(256), xlab = "1-sqrt(JSD)")

