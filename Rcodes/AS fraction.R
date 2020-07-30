setwd("~/Desktop/Jiali/UTK/2020 Summer/AS/alternative_splicing_project/") # set the working folder
library(plyr)

# read peach, apricot, and cherry transcript raw counts and reciprocal best hit file
peach1_data <- read.csv("data/peach_transcript_counts.csv", header = T, row.names = 1)
apricot1_data <- read.csv("data/apricot_transcript_counts.csv", header = T, row.names = 1)
peach2_data <- read.csv("data/peach2_transcript_counts.csv", header = T, row.names = 1, sep = "\t")
apricot2_data <- read.csv("data/apricot2_transcript_counts.csv", header = T, row.names = 1,sep="\t")

peach1_data <- peach1_data[grep("Prupe",rownames(peach1_data)),]
peach2_data <- peach2_data[grep("Prupe",rownames(peach2_data)),]
apricot1_data <- apricot1_data[grep("Prupe",rownames(apricot1_data)),]
apricot2_data <- apricot2_data[grep("Prupe",rownames(apricot2_data)),]

# combine transcript counts for all three
common_genes <- Reduce(intersect, list(rownames(apricot1_data), rownames(peach1_data), rownames(apricot2_data), rownames(peach2_data)))
combined_data <- data.frame(row.names = common_genes) # create an empty dataframe with common genes
combined_data <- transform(merge(combined_data, peach1_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, apricot1_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, peach2_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)
combined_data <- transform(merge(combined_data, apricot2_data, by="row.names", all.x = T), row.names=Row.names, Row.names=NULL)

# load PSI data

# load gene length file
transcript_data <- read.table("data/t_data.ctab", header = T, sep = "\t")
combined_data <- transform(merge(combined_data, transcript_data[,c(6,8)], by.x = "row.names", by.y = "t_name", all.x=T), row.names=Row.names, Row.names=NULL)

# RPKM = (CDS read count * 10^9) / (CDS length * total transcript count)
# compute RPKM
combined_data_tpm <- (combined_data[,-25] * 10^9)/(rep(combined_data$length,24) * as.numeric(colSums(combined_data[,-25])))

# extract each tissue type
peach1_flower <- data.frame("TPM"=rowMeans(combined_data_tpm[,c(22,23,24)]))
# load PSI data
peach1_PSI <- read.table("data/p_armeniaca_2_psi_isoform.psi", header = T, sep = "\t", row.names = 1)
rownames(peach1_PSI) <- gsub("MSTRG.\\d+;","",rownames(peach1_PSI))
peach1_PSI<- peach1_PSI[grep("Prupe",rownames(peach1_PSI)),]

#divide into 20 bins, from low expression to high expression
peach1_flower <- peach1_flower[order(peach1_flower$TPM, decreasing = F), , drop = FALSE]
peach1_PSI <- peach1_PSI[rownames(peach1_flower),] # make psi the same order as expression data
peach1_flower_data <- data.frame("Bin" = character(),
                                 "AS fraction" = numeric()) # create empty dataframe to store results
peach1_flower$logRPKM <- log10(peach1_flower$TPM+1)
bin_size <- round(peach1_flower$logRPKM[length(peach1_flower$logRPKM)] / 20, digits = 2)
PSI_data <- rowMeans(peach1_PSI[,c(4,5,6)])
for (i in c(0:19)) {
  peach1_flower_data[i+1,1] <- paste0(i*bin_size,"-",(i+1)*bin_size)
  genes_range <- which(peach1_flower$logRPKM >= (i)*bin_size & peach1_flower$logRPKM < (i+1)*bin_size)
  PSI_number <- length(which(PSI_data[genes_range] > 0))
  AS_fraction <- PSI_number/length(genes_range)
  peach1_flower_data[i+1,2] <- AS_fraction
}
library(ggplot2)
ggplot(data=peach1_flower_data, aes(x=Bin, y=AS.fraction, group=1)) +
  geom_line()+
  theme_classic(base_size = 12)+ ylim(0,1)+
  labs(title = "apricot2 bud", x= "log10 RPKM", y = "AS Fraction")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


