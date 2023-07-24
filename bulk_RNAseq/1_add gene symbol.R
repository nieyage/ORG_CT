# add gene symbol and QC 
#######change gene ID and get matrix ######
rm(list = ls())
setwd("D:/project/NRCM_TA")
counts<-read.table("./counts/join.count",header=T)
head(counts)
tail(counts)
#rownames(counts) <- counts
library('biomaRt')
library("curl")
library(ggrepel)
library(dplyr)
my_mart <-useMart("ensembl")
list<-listDatasets(my_mart)
#####查看database中有哪些物种的基因####
mart <- useDataset("rnorvegicus_gene_ensembl", 
                   useMart("ensembl"))#####Rat gene:rnorvegicus_gene_ensembl
listAttributes(mart)

my_ensembl_gene_id<-counts$gene
head(my_ensembl_gene_id)
options(timeout = 4000000)
#####use ensembl_transcript_id to trans######
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',
                                 'ensembl_transcript_id',"description"),
                    filters = 'ensembl_gene_id', 
                    values= my_ensembl_gene_id, mart = mart)
mms_symbols$gene<-mms_symbols$ensembl_gene_id
counts<-merge(counts,mms_symbols,by="gene")
str(counts)
counts<-counts[,c(2:7,9)]
counts<-counts[!duplicated(counts$external_gene_name), ]
rownames(counts)<-counts$external_gene_name
head(counts)
counts<-counts[,1:6]
colnames(counts)<-c("AdHmmr-rep1","AdHmmr-rep2","AdHmmr-rep3",
                    "AdVector-rep1","AdVector-rep2","AdVector-rep3")

write.table(counts,"./counts/rawcounts_add_genesymbol.txt")
write.csv(counts,"./counts/rawcounts_add_genesymbol.txt")

# global raw_counts pca
library(ggpubr)
library(ggthemes)
library(gmodels)
library(export)
counts<-read.csv("/Users/fraya/Documents/project/ORG-CT/01_bulkRNA/rawcounts_add_genesymbol.csv",row.names = 1)
head(counts)
pca.info <- fast.prcomp(counts)
head(pca.info)
summary(pca.info) 
head(pca.info$rotation) 
pca.data <- data.frame(sample = rownames(pca.info$rotation), 
                       Type=c(rep("Nor",3),rep("H-R-DMSO",3),rep("H-R-CT",3)),
                       pca.info$rotation)
pdf("ORG_CT_global_pcs.pdf")
p<-ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, ellipse.type="convex")
p
dev.off()
