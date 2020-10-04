setwd("~/Desktop/Jiali/UTK/2020 Summer/AS/alternative_splicing_project/") # set the working folder

# peach genotype1 
peach1_data <- read.csv("data/peach_transcript_counts.csv", header = T, row.names = 1)
meta <- read.csv("data/metadata.csv", header = T, row.names = 1)
meta$labels <- paste0(meta$Species,"-",meta$Tissue)
library(edgeR)
d <- DGEList(counts = peach1_data, group = meta$Tissue[1:6])
d <- calcNormFactors(d, method = "TMM")
d = estimateCommonDisp(d, verbose=TRUE)
normalized_counts <- d$pseudo.counts
write.table(normalized_counts[,1:3], "data/normalized_data/peach1_normalized_flower.txt", sep = "\t", quote = F)
write.table(normalized_counts[,4:6],"data/normalized_data/peach1_normalized_bud.txt", sep = "\t",quote = F)

# apricot genotype1
apricot1_data <- read.csv("data/apricot_transcript_counts.csv", header = T, row.names = 1)
meta <- read.csv("data/metadata.csv", header = T, row.names = 1)
meta$labels <- paste0(meta$Species,"-",meta$Tissue)

d <- DGEList(counts = apricot1_data, group = meta$Tissue[7:12])
d <- calcNormFactors(d, method = "TMM")
d = estimateCommonDisp(d, verbose=TRUE)
normalized_counts <- d$pseudo.counts
write.table(normalized_counts[,c(2,4,6)], "data/normalized_data/apricot1_normalized_flower.txt", sep = "\t", quote = F)
write.table(normalized_counts[,c(1,3,5)],"data/normalized_data/apricot1_normalized_bud.txt", sep = "\t",quote = F)

# peach genotype2 
peach2_data <- read.csv("data/peach2_transcript_counts.csv", header = T, row.names = 1)
meta <- read.csv("data/metadata.csv", header = T, row.names = 1)
meta$labels <- paste0(meta$Species,"-",meta$Tissue)

d <- DGEList(counts = peach2_data, group = meta$Tissue[13:18])
d <- calcNormFactors(d, method = "TMM")
d = estimateCommonDisp(d, verbose=TRUE)
normalized_counts <- d$pseudo.counts
write.table(normalized_counts[,4:6], "data/normalized_data/peach2_normalized_flower.txt", sep = "\t", quote = F)
write.table(normalized_counts[,1:3],"data/normalized_data/peach2_normalized_bud.txt", sep = "\t",quote = F)

# apricot genotype1
apricot2_data <- read.csv("data/apricot2_transcript_counts.csv", header = T, row.names = 1)
meta <- read.csv("data/metadata.csv", header = T, row.names = 1)
meta$labels <- paste0(meta$Species,"-",meta$Tissue)

d <- DGEList(counts = apricot2_data, group = meta$Tissue[19:24])
d <- calcNormFactors(d, method = "TMM")
d = estimateCommonDisp(d, verbose=TRUE)
normalized_counts <- d$pseudo.counts
write.table(normalized_counts[,1:3], "data/normalized_data/apricot2_normalized_flower.txt", sep = "\t", quote = F)
write.table(normalized_counts[,4:6],"data/normalized_data/apricot2_normalized_bud.txt", sep = "\t",quote = F)
