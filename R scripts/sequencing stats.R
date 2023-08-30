library(tidyverse)#to allow qiime2R to work
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

 
 
rawseqstat <- read.csv("R inputs/Ladpt sequencing stats.csv") #information provided by the sequencing facility
sum(rawseqstat$Number.of.Reads) # 1,924,507
sum(rawseqstat$Number.of.Bases) # 962,253,500
range(rawseqstat$Number.of.Reads) #26,029 129,597
mean(rawseqstat$Number.of.Reads) #96,225.35
range(rawseqstat$Number.of.Bases) #13,014,500 64,798,500
mean(rawseqstat$Number.of.Reads[1:10]) #95,508.7
mean(rawseqstat$Number.of.Reads[11:20]) #96,942

table_prefilter	<- read_qza("R inputs/Lemna_GxGxE_denoisetable_2.qza")$data
feat.tax	<- read_qza("R inputs/Lemna_GxGxE_taxonomy_den.qza")$data ##NOT FILTERED to kept features, or even to rm streptophyta (chloroplast, which are removed from the dataset, along with mitochondrial sequences)
taxlevs <- sapply(1:nrow(feat.tax), function(z) strsplit(as.character(feat.tax$Taxon[z]),split="; "))
tax <- matrix(NA,ncol=7,nrow=nrow(feat.tax))
for(i in 1:nrow(tax)){ tax[i,1:length(taxlevs[[i]])]<- taxlevs[[i]]}
rownames(tax) <- as.character(feat.tax$Feature.ID)


totreads <- colSums(table_prefilter)
mean(totreads[c(1,3,5,7,9,11,13,15,17,19)]) # field, 66850.1
mean(totreads[c(2,4,6,8,10,12,14,16,18,20)]) # inocula 65665

streptophyta <- table_prefilter[names(which(tax[,4]=="o__Streptophyta")),]
mean(colSums(streptophyta[,c(1,3,5,7,9,11,13,15,17,19)])) # field, 29238.3
mean(colSums(streptophyta[,c(2,4,6,8,10,12,14,16,18,20)])) # inocula 14.7

