# calculate conserved AS genes 
setwd("~/Desktop/Jiali/UTK/2020 Summer/AS/alternative_splicing_project/")
# read the RBH
#RBH <- read.table("data/reciprocal_best_hits.txt")

# convert cherry to peach gene IDs
#cherry_genes <- scan("splice_events/avium/events_AL_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")
#cherry_genes <- gsub("\\[|]","",cherry_genes)
#cherry_rbh <- data.frame(row.names = cherry_genes)
#cherry_rbh <- transform(merge(cherry_rbh, RBH, by.x = "row.names",by.y = "V2", all.x = T), row.names=Row.names, Row.names=NULL)
#cherry_rbh_noNA <- na.omit(cherry_rbh)
#cherry_rbh_noNA$V1 <- gsub("\\.\\d+\\.p","",cherry_rbh_noNA$V1)
peach1_genes <- scan("splice_events/persica/events_SE_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")
peach2_genes <- scan("splice_events/persica2/events_SE_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")
apricot1_genes <- scan("splice_events/armeniaca/events_SE_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")
apricot2_genes <- scan("splice_events/armeniaca2/events_SE_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")

# peach1 + peach2
length(intersect(peach1_genes, peach2_genes))
# apricot1 + apricot2
length(intersect(apricot1_genes, apricot2_genes))
# all 
length(Reduce(intersect, list(peach1_genes, peach2_genes,apricot1_genes, apricot2_genes)))

