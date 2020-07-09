# calculate conserved AS genes 
setwd("~/Desktop/Jiali/UTK/2020 Summer/AS/alternative_splicing_project/")
# read the RBH
RBH <- read.table("data/reciprocal_best_hits.txt")

# convert cherry to peach gene IDs
cherry_genes <- scan("splice_events/avium/events_AL_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")
cherry_genes <- gsub("\\[|]","",cherry_genes)
cherry_rbh <- data.frame(row.names = cherry_genes)
cherry_rbh <- transform(merge(cherry_rbh, RBH, by.x = "row.names",by.y = "V2", all.x = T), row.names=Row.names, Row.names=NULL)
cherry_rbh_noNA <- na.omit(cherry_rbh)
cherry_rbh_noNA$V1 <- gsub("\\.\\d+\\.p","",cherry_rbh_noNA$V1)
peach_genes <- scan("splice_events/persica/events_AL_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")
peach_genes <- gsub("\\[|]|_v2.0","", peach_genes)
apricot_genes <- scan("splice_events/armeniaca/events_AL_strict.gene_list.txt", what="character", strip.white = TRUE, sep = ",")
apricot_genes <- gsub("\\[|]|_v2.0","", apricot_genes)
# peach + apricot
length(intersect(peach_genes, apricot_genes))
# peach+ cherry
length(intersect(peach_genes, cherry_rbh_noNA$V1))
# Cherry + apricot
length(intersect(cherry_rbh_noNA$V1, apricot_genes))
# all three
Reduce(intersect, list(apricot_genes, peach_genes, cherry_rbh_noNA$V1))
