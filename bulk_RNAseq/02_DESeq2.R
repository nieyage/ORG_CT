library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("/Users/fraya/Documents/project/ORG-CT/01_bulkRNA/")
counts<-read.csv("/Users/fraya/Documents/project/ORG-CT/01_bulkRNA/rawcounts_add_genesymbol.csv",row.names = 1)
head(counts)
str(counts)
condition <- factor(c(rep("Nor",3),rep("H_R_DMSO",3),rep("H_R_CT",3)))
colData <- data.frame(row.names=colnames(counts),condition=condition)
colData
countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData <- na.omit(countData)
dds <- DESeqDataSetFromMatrix(countData,colData,formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts,"./ORG_CT-normalized_counts.csv")
dds <- DESeq(dds)
res <- results(dds)
# condition in HR global PCA
vsd <- vst(dds)
rld <- rlog(dds)
library(ggplot2)
vsd <- varianceStabilizingTransformation(dds)
head(vsd)

###sample corealation heatmap####
rlogMat <- assay(rld)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
library(amap)
library(pheatmap)
pheatmap(pearson_cor,
         cluster_cols = T,cluster_rows = T,
         color = hmcol,
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)


pca_data <- plotPCA(vsd, intgroup=c("condition"),
                    returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca_data, "percentVar"))
head(pca_data)
pdf("./condition in HR global PCA.pdf",width=8,height=6)
ggplot(pca_data, aes(PC1, PC2, color =condition)) +
  geom_point(size=3) +geom_text(label=paste(pca_data$name),colour="black",size=4)+
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+theme(panel.grid.major=element_line(colour=NA))
dev.off()


#3D PCA#
library("FactoMineR")
library("factoextra")
library("scatterplot3d")
library("gmodels")
pca_count<- normalized_counts

pca.info <- fast.prcomp(pca_count)
head(pca.info$rotation)
head(pca.info)
pca.data <- data.frame(sample =rownames(pca.info$rotation),
                       condition=condition,
                       pca.info$rotation)

library(ggsci)
library("scales")
colors=pal_npg("nrc")(10)
show_col(pal_npg("nrc")(10))
colors_pal<-colors[c(3,2,5,1)]
colors <- colors_pal[as.factor(pca.data$condition)]

 pVar <- pca.info$sdev^2/sum(pca.info$sdev^2)
 pVar = round(pVar,digits = 3)

   paste0("PC1 (",as.character(pVar[1] * 100 ),"%)")
   paste0("PC2 (",as.character(pVar[2] * 100 ),"%)")
   paste0("PC3 (",as.character(pVar[3] * 100 ),"%)")

str(pca.data)
s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],pch=15,color = colors,mar=c(5,5,5,5),
                     angle = 60, type="p",cex.symbols = 1,
                     main = "ORG_CT 3D PCA plot",
                     xlab="PC1 ",
                     ylab = "PC2",
                     zlab = "PC3 ") 
legend("topright", legend = c("Nor","H_R_DMSO","H_R_CT"),
       col =colors_pal, pch = 15, bty="n",
       inset = -0.2,xpd = TRUE, horiz = FALSE)

DMSO_Nor <-results(dds,contrast = c("condition","H_R_DMSO","Nor"))
H_R_CT_DMSO <-results(dds,contrast = c("condition","H_R_CT","H_R_DMSO"))
H_R_CT_Nor <-results(dds,contrast = c("condition","H_R_CT","Nor"))
# save the DEG results 
write.csv(DMSO_Nor,"DMSOvsNor-DEG.csv")
write.csv(H_R_CT_DMSO,"H_R_CTvsDMSO-DEG.csv")
write.csv(H_R_CT_Nor,"H_R_CTvsNor-DEG.csv")

DEG<- rownames(subset(DMSO_Nor, padj < 0.05 & abs(log2FoldChange) >1))

rlogMat <- assay(rld)[DEG,]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
pheatmap(pearson_cor,
         cluster_cols = T,cluster_rows = T,
         color = hmcol,
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)

pearson_cor

library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
#DMSO_Nor
#volcano plot#
data<-as.data.frame(DMSO_Nor)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'H_R_DMSO','Nor'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
pdf("DMSOvsNor-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("DMSO vs Control DEG:",DEG_data[1,1],":",DEG_data[1,2]," ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(DMSO_Nor, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),1:6]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("DMSOvsNor-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(DMSO_Nor, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(DMSO_Nor, pvalue < 0.01 & log2FoldChange< -1) )
pdf("DMSOvsCT-up_in_DMSO-GO.pdf")
gene.df <- bitr(up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-CC.csv")
dev.off()
######KEGG##########
pdf("DMSOvsCT-up_in_DMSO-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(ego,"DMSOvsCT-up_in_DMSO-KEGG.csv")
dev.off()

pdf("DMSOvsCT-up_in_Nor-GO.pdf")
gene.df <- bitr(down, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_Nor-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_Nor-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_Nor-CC.csv")
dev.off()
######KEGG##########
pdf("DMSOvsCT-up_in_Nor-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(ego,"DMSOvsCT-up_in_Nor-KEGG.csv")
dev.off()

#H_R_CT_DMSO
#volcano plot#
data<-as.data.frame(H_R_CT_DMSO)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'H_R_CT','H_R_DMSO'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
pdf("H_R_CT_DMSO-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("H_R_CT vs DMSO DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(H_R_CT_DMSO, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),4:9]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("H_R_CT_DMSO-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(H_R_CT_DMSO, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(H_R_CT_DMSO, pvalue < 0.01 & log2FoldChange< -1) )
pdf("H_R_CTvsDMSO-up_in_H_R_CT-GO.pdf")
gene.df <- bitr(up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsDMSO-up_in_H_R_CT-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsDMSO-up_in_H_R_CT-GO-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsDMSO-up_in_H_R_CT-GO-CC.csv")
dev.off()
######KEGG##########
pdf("H_R_CTvsDMSO-up_in_H_R_CT-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(ego,"H_R_CTvsDMSO-up_in_H_R_CT-KEGG.csv")
dev.off()

pdf("H_R_CTvsDMSO-up_in_DMSO-GO.pdf")
gene.df <- bitr(down, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsDMSO-up_in_DMSO-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsDMSO-up_in_DMSO-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsDMSO-up_in_DMSO-CC.csv")
dev.off()
######KEGG##########
pdf("H_R_CTvsDMSO-up_in_DMSO-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(ego,"H_R_CTvsDMSO-up_in_DMSO-KEGG.csv")
dev.off()


#HR_CT_Nor
#volcano plot#
data<-as.data.frame(H_R_CT_Nor)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'R5','DMSO'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
pdf("H_R_CT_Nor-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("H_R_CT vs control DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(H_R_CT_Nor, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),c(1:3,7:9)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("H_R_CT_Nor-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(H_R_CT_Nor, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(H_R_CT_Nor, pvalue < 0.01 & log2FoldChange< -1) )
pdf("H_R_CTvsNor-up_in_H-R-CT-GO.pdf")
gene.df <- bitr(up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsNor-up_in_H-R-CT-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsNor-up_in_H-R-CT-GO-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsNor-up_in_H-R-CT-GO-CC.csv")
dev.off()
######KEGG##########
pdf("H_R_CTvsNor-up_in_H-R-CT-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(ego,"H_R_CTvsNor-up_in_H-R-CT-KEGG.csv")
dev.off()

pdf("H_R_CTvsNor-up_in_Nor-GO.pdf")
gene.df <- bitr(down, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsNor-up_in_Nor-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsNor-up_in_Nor-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R_CTvsNor-up_in_Nor-CC.csv")
dev.off()
######KEGG##########
pdf("H_R_CTvsNor-up_in_Nor-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(ego,"H_R_CTvsNor-up_in_Nor-KEGG.csv")
dev.off()

DMSO_Nor_DEG<- rownames(subset(DMSO_Nor, pvalue < 0.05 & abs(log2FoldChange) >1))
H_R_CT_DMSO_DEG<- rownames(subset(H_R_CT_DMSO, pvalue < 0.05 & abs(log2FoldChange) >1))
H_R_CT_Nor_DEG<- rownames(subset(H_R_CT_Nor, pvalue < 0.05 & abs(log2FoldChange) >1))
union_DEG<- unique(c(DMSO_Nor_DEG,H_R_CT_DMSO_DEG,H_R_CT_Nor_DEG))

# heatmap 
Group = factor(c(rep("Nor",3),rep("DMSO",3),rep("H_R_CT",3)))
anndf <- data.frame(Group)
rownames(anndf) <- colnames(normalized_counts)
#自定义分组颜色条的颜色；
anncol = list(Group=c(Nor="#EDF1D6",DMSO="#609966",H_R_CT="#9DC08B"))

targetcount<-normalized_counts[which(rownames(normalized_counts)%in%DEG),]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#pdf("CT-DMSO-H-R10-heatmap_allDEGin_HR10.pdf",width=6,height=10)
count<-na.omit(count)
p<-pheatmap(count,cluster_cols = F,cluster_rows = T,
            color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
            #cellwidth = 10, cellheight = 10,
            cutree_rows = 4,
            annotation_col=anndf,
            annotation_colors=anncol,
            show_rownames=F,show_colnames=T)
row_cluster <- cutree(p$tree_row,k=4)
table(row_cluster)
annotation_row <- data.frame(
  type = paste("Cluster",row_cluster,sep="")
)
rownames(annotation_row) <- names(row_cluster);
pdf("All-heatmap_DEG_injury.pdf",width=5,height=10)
#bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         cutree_rows = 4,
         #legend_breaks=seq(-2,2,1),
         #breaks=bk,
         annotation_col=anndf,
         annotation_row=annotation_row,
         annotation_colors=anncol,
         show_rownames=F,show_colnames=T)
dev.off()

#GO and KEGG 
for (i in unique(annotation_row$type)){
  print(i);
  gene_in_cluster<-rownames(annotation_row)[which(annotation_row$type==i)];
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Hs.eg.db)
  #GO
  pdf(paste0("heatmap",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    type="BP"
    ego <- enrichGO(gene.df$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Hs.eg.db,
                    ont = type,
                    pAdjustMethod  = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      p_ego_go <- goplot(ego)
      
      print(p);
      print(p_ego_go);
      write.csv(ego,paste0("heatmap_",i,"_GO_",type,".csv"))
    }
  }
  dev.off()
  #KEGG 
  pdf(paste0("heatmap_",i,"_KEGG.pdf",sep=""))
  ego <- enrichKEGG(
    gene = gene.df$ENTREZID,
    keyType = "kegg",
    organism  = "hsa",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
    p<-barplot(ego, showCategory=20);
    print(p);
    edox<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
    write.table(edox,paste0("heatmap_",i,"_KEGG_",type,".csv"))
  }
  dev.off()
}
