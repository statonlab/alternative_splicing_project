setwd("~/Desktop/Jiali/UTK/2020 Summer/AS")

apricot <- read.csv("apricot_transcript_counts.csv", header = T, stringsAsFactors = F)
apricot$bud <- rowMeans(apricot[,c(2,4,6)]) # mean of bud samples
apricot$flower <- rowMeans(apricot[,c(3,5,7)]) # mean of flower samples
apricot$gene_id <- gsub(".\\d+_v2.0.a1|.{0,1}$", "",apricot$transcript_id)

# Summary - The number of transcripts in buds/flowers/both
bud_transcripts_apri <- apricot[apricot$bud > 0 & apricot$flower == 0,]
flower_transcripts_apri <- apricot[apricot$bud == 0 & apricot$flower > 0, ]
both_apri <- apricot[apricot$bud > 0 & apricot$flower> 0,]

length(unique(apricot$gene_id[grep("Prupe",apricot$gene_id)]))
length(bud_transcripts_apri$transcript_id)
length(unique(bud_transcripts_apri$gene_id))
length(unique(bud_transcripts_apri$gene_id[grep("Prupe",bud_transcripts_apri$gene_id)]))
length(flower_transcripts_apri$transcript_id)
length(unique(flower_transcripts_apri$gene_id))
length(unique(flower_transcripts_apri$gene_id[grep("Prupe",flower_transcripts_apri$gene_id)]))
length(both_apri$transcript_id)
length(unique(both_apri$gene_id))
length(unique(both_apri$gene_id[grep("Prupe",both_apri$gene_id)]))

# peach transcripts
peach <- read.csv("peach_transcript_counts.csv", header = T,stringsAsFactors = F)
peach$bud <- rowMeans(peach[,c(2,4,6)]) # mean of bud samples
peach$flower <- rowMeans(peach[,c(3,5,7)]) # mean of flower samples
peach$gene_id <- gsub(".\\d+_v2.0.a1|.{0,1}$", "",peach$transcript_id)

# Summary - The number of transcripts in buds/flowers/both
bud_transcripts_peach <- peach[peach$bud > 0 & peach$flower == 0,]
flower_transcripts_peach <- peach[peach$bud == 0 & peach$flower > 0, ]
both_peach <- peach[peach$bud > 0 & peach$flower> 0,]

length(unique(peach$gene_id))
length(unique(peach$gene_id[grep("Prupe",peach$gene_id)]))
length(bud_transcripts_peach$transcript_id)
length(unique(bud_transcripts_peach$gene_id))
length(unique(bud_transcripts_peach$gene_id[grep("Prupe",bud_transcripts_peach$gene_id)]))
length(flower_transcripts_peach$transcript_id)
length(unique(flower_transcripts_peach$gene_id))
length(unique(flower_transcripts_peach$gene_id[grep("Prupe",flower_transcripts_peach$gene_id)]))
length(both_peach$transcript_id)
length(unique(both_peach$gene_id))
length(unique(both_peach$gene_id[grep("Prupe",both_peach$gene_id)]))

# cherryripts
cherry <- read.csv("cherry_transcript_count.csv", header = T,stringsAsFactors = F)
cherry$bud <- rowMeans(cherry[,c(5,6,7)]) # mean of bud samples
cherry$flower <- rowMeans(cherry[,c(2,3,4)]) # mean of flower samples
cherry$gene_id <- gsub(".\\d+_v2.0.a1|.{0,1}$", "",cherry$transcript_id)

# Summary - The number of transcripts in buds/flowers/both
bud_transcripts_cherry <- cherry[cherry$bud > 0 & cherry$flower == 0,]
flower_transcripts_cherry <- cherry[cherry$bud == 0 & cherry$flower > 0, ]
both_cherry <- cherry[cherry$bud > 0 & peach$flower> 0,]

length(unique(cherry$transcript_id))
length(bud_transcripts_cherry$transcript_id)
length(unique(bud_transcripts_peach$gene_id))
length(flower_transcripts_peach$transcript_id)
length(unique(flower_transcripts_peach$gene_id))
length(both_peach$transcript_id)
length(unique(both_peach$gene_id))


# common in both peach and apricot
length(intersect(peach$transcript_id, apricot$transcript_id))
length(intersect(peach$gene_id, apricot$gene_id))
length(intersect(bud_transcripts_apri$transcript_id, bud_transcripts_peach$transcript_id))
length(intersect(flower_transcripts_apri$transcript_id, flower_transcripts_peach$transcript_id))
length(intersect(flower_transcripts_apri$gene_id, flower_transcripts_peach$gene_id))
length(intersect(bud_transcripts_apri$gene_id, bud_transcripts_peach$gene_id))
length(intersect(both_apri$transcript_id, both_peach$transcript_id))
length(intersect(both_apri$gene_id, both_peach$gene_id))

## read the rbh results
library(plyr)
rbh <- read.table("reciprocal_best_hits.txt", sep = "\t", stringsAsFactors = F)
cherry_peach_hit <- merge(cherry, rbh, by.x = "transcript_id", by.y = "V2", all.x = T)
cherry_rbh <- na.omit(cherry_peach_hit)
cherry_rbh$V1 <- gsub("\\.p", "_v2.0.a1",cherry_rbh$V1)
cherry_rbh_bud <- cherry_rbh[cherry_rbh$bud > 0 & cherry_rbh$flower ==0, ]
cherry_rbh_flower <- cherry_rbh[cherry_rbh$bud == 0 & cherry_rbh$flower >0, ]
cherry_rbh_both <- cherry_rbh[cherry_rbh$bud > 0 & cherry_rbh$flower >0, ]

length(intersect(cherry_rbh$transcript_id, both_cherry$transcript_id))

# compare cherry and peach
length(intersect(cherry_rbh_both$V1, both_peach$transcript_id))


