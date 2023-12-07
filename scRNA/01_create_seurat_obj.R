## Load required packages
remotes::install_version("Seurat", version = "4.3")
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(DoubletFinder)
##Load datasets and create Seurat objects with the raw (non-normalized data).
CT   <-  CreateSeuratObject(counts = Read10X(data.dir = "/public/home/nieyg/project/ORG_CT_scRNA/00_data/cellranger/CT/outs/filtered_feature_bc_matrix"), 
	project="CT", assay = "RNA")
DMSO <-  CreateSeuratObject(counts = Read10X(data.dir = "/public/home/nieyg/project/ORG_CT_scRNA/00_data/cellranger/DMSO/outs/filtered_feature_bc_matrix"), 
	project="DMSO", assay = "RNA")
NOR <-  CreateSeuratObject(counts = Read10X(data.dir = "/public/home/nieyg/project/ORG_CT_scRNA/00_data/cellranger/NOR/outs/filtered_feature_bc_matrix"), 
	project="NOR1", assay = "RNA")

#####添加sample信息####

##############################
#######pre-processing#########
#######      1_QC      #########
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#注意参考基因组里面线粒体相关基因的名称
CT[["percent.mt"]] <- PercentageFeatureSet(CT, pattern = "^MT-")
DMSO[["percent.mt"]] <- PercentageFeatureSet(DMSO, pattern = "^MT-")
NOR[["percent.mt"]] <- PercentageFeatureSet(NOR, pattern = "^MT-")

# Visualize QC metrics as a violin plot
pdf("./00_QC/samplel_qc_plot_nopoint.pdf")
VlnPlot(DMSO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,)
VlnPlot(NOR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,)
VlnPlot(CT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,)
dev.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("./00_QC/samplel_qc_plot.pdf")
plot1 <- FeatureScatter(CT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot3 <- FeatureScatter(NOR, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(NOR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4

plot3 <- FeatureScatter(DMSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(DMSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
dev.off()
######################################################################################
# 2_DOUBLETFINDER
######################################################################################

# NORMALIZE for DoubletFinder
library(gridExtra)
library(DoubletFinder)
sample <- c("NOR","CT","DMSO")
objList <- list(DMSO,CT,NOR)

for (i in seq_len(length(objList))) {
  # Normalizing
  dataset.name <- sample[i]
  sratDecontx    <- SCTransform(objList[[i]], verbose = F)
  # Run cluster analysis
  sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
  sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindClusters(sratDecontx, verbose = T)
  # PLOT: umap/tsne
  p1 <- DimPlot(sratDecontx, reduction = "umap", label = TRUE)
  p2 <- DimPlot(sratDecontx, reduction = "tsne", label = TRUE)
  g = arrangeGrob(p1,p2, ncol = 2)
  filename <- paste0("./00_QC/1b_filtering_UMAP-tSNE_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(g)
  objList[[i]] <- sratDecontx
    }


objList2<-list()
for (i in seq_len(length(objList))) {
  dataset.name<-sample[i]
  print(dataset.name)
  sratDecontx <- objList[[i]]
  # Compute expected doublet rate
  cellN=nrow(sratDecontx@meta.data)
  expDoubletRate = (cellN*0.0008 + 0.0527)*0.01
  normalizationMethod='SCTransform'
  sweep.res.list_scData <- paramSweep(sratDecontx, 
                                         PCs = 1:50, 
                                         sct = TRUE,num.cores = 4) #num.cores = 4
  sweep.stats_scData <- summarizeSweep(sweep.res.list_scData, GT = FALSE)
  bcmvn_scData <- find.pK(sweep.stats_scData)
  bcmvn_scData$pK <- as.numeric(as.character(bcmvn_scData$pK))
  pK1=bcmvn_scData$pK[bcmvn_scData$BCmetric==max(bcmvn_scData$BCmetric)]
  print(head(pK1))
  # PLOT: pK selection
  p1=ggplot(data=bcmvn_scData, 
            aes(x=pK, y=BCmetric, group=2)) +
    geom_line(color="blue")+
    geom_point()+
    geom_vline(xintercept=pK1, linetype="dashed", color = "red")+
    labs(title="pK Selection",x="pK", y = "BCmvn")+
    theme_classic()
  filename <- paste0("./00_QC/2a_doubletfinder_pkselection_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)
  # More doublet finder
  pK1=as.numeric(as.character( pK1 ))
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- sratDecontx@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)   
  nExp_poi <- round(expDoubletRate*nrow(sratDecontx@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  sratDecontx <- doubletFinder( sratDecontx, PCs = sratDecontx@commands$RunUMAP.SCT.pca$dims,
                                   pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  sratDecontx@meta.data$DoubletFinder =  sratDecontx@meta.data[,grep('DF.classifications', colnames( sratDecontx@meta.data))]
  # PLOT: Doublet Finder graphs
  p2 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'DoubletFinder')
  p3 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'DoubletFinder')
  g = arrangeGrob(p2,p3, ncol = 2)
  filename <- paste0("./00_QC/2b_doubletfinder_ScatterPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(filename, g)
  # PLOT: Violin Plots
  p4 <- VlnPlot(sratDecontx, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = 'DoubletFinder', pt.size = 0)
  filename <- paste0("./00_QC/2c_doubletfinder_DoubletFinder-ViolinPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p4)
  objList2[[i]] <- sratDecontx
}

######################################################################################
# 3_FEATURE FILTERING
######################################################################################

objList<- objList2
# Subset identified singlets
for (i in seq_len(length(objList))) {
    objList[[i]] <- subset(objList[[i]], subset = DoubletFinder == "Singlet")
}


library(ggplot2)

for (i in seq_len(length(objList))) {
    objList[[i]] <- subset(objList[[i]],nFeature_RNA >200 & nFeature_RNA < 8000 & percent.mt < 40 )
}
#objList <- list(DMSO,CT,NOR)
DMSO<- subset(DMSO,cells=colnames(objList[[1]]))
CT<- subset(CT,cells=colnames(objList[[2]]))
NOR<- subset(NOR,cells=colnames(objList[[3]]))

# Simply merge Seurat objects
merged_obj <- merge(x=NOR,y=c(DMSO,CT),add.cell.ids = c("NOR","DMSO","CT"),project = "ORG_CT")
Idents(merged_obj) <- gsub("_.*", "", colnames(merged_obj))
#split the combined object into a list, with each dataset as an element
ORG_CT.list <- SplitObject(merged_obj,split.by = "ident")

pdf("./00_QC/NOR_normalization_plot.pdf")
ORG_CT.list[[1]] <- FindVariableFeatures(ORG_CT.list[[1]], verbose = FALSE)
top10 <- head(VariableFeatures(ORG_CT.list[[1]]),10)
plot1 <- VariableFeaturePlot(ORG_CT.list[[1]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
dev.off();
pdf("./00_QC/DMSO_normalization_plot.pdf")
ORG_CT.list[[2]] <- NormalizeData(ORG_CT.list[[2]],verbose = FALSE)
ORG_CT.list[[2]] <- FindVariableFeatures(ORG_CT.list[[2]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(ORG_CT.list[[2]]),10)
plot1 <- VariableFeaturePlot(ORG_CT.list[[2]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
dev.off();
pdf("./00_QC/CT_normalization_plot.pdf")
ORG_CT.list[[3]] <- NormalizeData(ORG_CT.list[[3]],verbose = FALSE)
ORG_CT.list[[3]] <- FindVariableFeatures(ORG_CT.list[[3]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(ORG_CT.list[[3]]),10)
plot1 <- VariableFeaturePlot(ORG_CT.list[[3]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
dev.off();

ORG_CT.anchors <- FindIntegrationAnchors(object.list = ORG_CT.list,anchor.features = 2000,dims = 1:30)
ORG_CT.integrated <- IntegrateData(anchorset = ORG_CT.anchors, dims = 1:30,features.to.integrate = rownames(ORG_CT.list[[1]]))
head(ORG_CT.integrated@meta.data)
DefaultAssay(ORG_CT.integrated) <- "integrated"

# scale and center features in the dataset
ORG_CT.integrated <- ScaleData(ORG_CT.integrated, features =rownames(ORG_CT.integrated),verbose = FALSE)
# Perform linear dimensional reduction
ORG_CT.integrated <- RunPCA(ORG_CT.integrated, npcs = 50, verbose = FALSE)
# Determine the ‘dimensionality’ of the dataset
# JackStraw 
ORG_CT.integrated <- JackStraw(ORG_CT.integrated, num.replicate = 100, dims =50)
ORG_CT.integrated  <- ScoreJackStraw(ORG_CT.integrated, dims = 1:50)
pdf("./00_QC/seclect_pca.pdf")
JackStrawPlot(ORG_CT.integrated, dims = 1:50)
# ‘Elbow plot’
ElbowPlot(ORG_CT.integrated,ndims=50)
dev.off()

DefaultAssay(ORG_CT.integrated) <- "integrated"
ORG_CT.integrated <- FindNeighbors(object = ORG_CT.integrated, dims = 1:25)
ORG_CT.integrated <- FindClusters(object = ORG_CT.integrated, resolution = 0.6)
ORG_CT.integrated <- RunUMAP(object = ORG_CT.integrated, dims = 1:25)

table(ORG_CT.integrated@active.ident)

pdf("./01_cluster/Unsupervised_cluster_Umap.pdf")
DimPlot(object = ORG_CT.integrated, group.by = "orig.ident",label = T)
DimPlot(object = ORG_CT.integrated, group.by = "seurat_clusters",label = T)
dev.off()
DefaultAssay(ORG_CT.integrated) <- "RNA"
saveRDS(ORG_CT.integrated, file = "./ORG_CT_integrated.rds")
