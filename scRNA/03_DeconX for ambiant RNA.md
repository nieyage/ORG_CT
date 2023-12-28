# DeconX for ambiant RNA 
## filter bam and Ambient RNA removal:
### After cellranger 

3. Step3: DecontX in R 

```
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(scater)
  #library(Signac)
  library(celda)
  library(DropletUtils)
  })
theme_set(theme_cowplot())
```
* NOR

```
1. Use SoupX for removing ambient RNAs by modified bam files 

NOR <-  CreateSeuratObject(counts = NOR_counts, min.cells=3,  project="NOR", assay = "RNA")
raw.matrix  <- Read10X_h5("./raw_feature_bc_matrix.h5",use.names = T)
srat.raw <- CreateSeuratObject(
        counts = raw.matrix$`GeDMSO Expression`,
        assay = "RNA"
      )
    srat    <- SCTransform(srat, verbose = F)
    # Idents(srat) <- srat@meta.data$seurat_clusters
    srat    <- RunPCA(srat, npcs = 50, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
    srat    <- RunTSDMSO(srat, dims = 1:50, verbose = F)
    srat    <- FindDMSOighbors(srat,reduction = 'pca', dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T,resolution =1.2)

    pdf('./decontX_outs/Original_UMAP_tSDMSO_F.pdf',width=15, height=6)
    print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsDMSO", label = TRUE))
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_VlnPlot_F.pdf', width=15, height=8)
    print( VlnPlot(srat, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_FeaturePlot_F.pdf', width=16, height=16)
    print(FeaturePlot(srat, reduction = 'tsDMSO', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
   # remove ambient RNAs
   library(Matrix)
   data_full<-as(as.matrix(NOR_counts),"dgCMatrix")
   #data_full<-Matrix(NOR_counts, sparse = TRUE) 
   data_empty<-raw.matrix$`GeDMSO Expression`
   data_empty<-data_empty[rownames(data_full),]
   data_full_sce <- SingleCellExperiment(assays = list(counts = data_full))
   data_empty_sce <- SingleCellExperiment(assays = list(counts = data_empty))
   data_full_sce <- decontX(x = data_full_sce, background = data_empty_sce)
   sce.decontX<-data_full_sce
    # Cluster labels on UMAP
    umap <- reducedDim(sce.decontX, "decontX_UMAP")

    pdf('./decontX_outs/decontX_contamination_F.pdf',width=15, height=6)
    print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
            plotDecontXContamination(sce.decontX))
    dev.off()
    DropletUtils:::write10xCounts("./decontX_outs/F_decontXcounts", round(decontXcounts(sce.decontX)),
                                  barcodes = rownames(colData(data_full_sce)))  
    #########################################
    ### decontX
      sratDecontx  <- 
      Read10X("./decontX_outs/F_decontXcounts") %>%
      CreateSeuratObject(project = "NOR", min.cells = 3, min.features = 200)   
    sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-") 
    sratDecontx    <- SCTransform(sratDecontx, verbose = F)
    sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
    sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- RunTSDMSO(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindDMSOighbors(sratDecontx, dims = 1:30, reduction = 'pca',verbose = F)
    sratDecontx    <- FindClusters(sratDecontx, verbose = T,resolution =1.2)   
    pdf('./decontX_outs/decontX_UMAP_tSDMSO_F.pdf', width=15, height=6)
    print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsDMSO", label = TRUE) )
    dev.off()

    pdf('./decontX_outs/decontX_coreceptor_VlnPlot_F.pdf',width=15, height=8)
    print( VlnPlot(sratDecontx, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/decontX_coreceptor_FeaturePlot_F.pdf', width=16, height=16)
    print(FeaturePlot(sratDecontx, reduction = 'tsDMSO', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
```

* CT

``` 
CT_data<-read.table("./counts_no_double_umi_001.tsv.gz")
CT_counts = spread(CT_data, V3, V2)
CT_counts[is.na(CT_counts)] <- 0
rownames(CT_counts)<-CT_counts$V1
CT_counts<-CT_counts[,-1]
CT_barcode<-paste(colnames(CT_counts),"-1",sep="")
colnames(CT_counts)<-CT_barcode

srat <-  CreateSeuratObject(counts = CT_counts, min.cells=3,  project="CT", assay = "RNA")

raw.matrix  <- Read10X_h5("./raw_feature_bc_matrix.h5",use.names = T)

srat.raw <- CreateSeuratObject(
        counts = raw.matrix$`GeDMSO Expression`,
        assay = "RNA"
      )

    srat    <- SCTransform(srat, verbose = F)
    # Idents(srat) <- srat@meta.data$seurat_clusters
    srat    <- RunPCA(srat, npcs = 50, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
    srat    <- RunTSDMSO(srat, dims = 1:50, verbose = F)
    srat    <- FindDMSOighbors(srat,reduction = 'pca', dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T,resolution =1.2)

    pdf('./decontX_outs/Original_UMAP_tSDMSO_CT.pdf',width=15, height=6)
    print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsDMSO", label = TRUE))
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_VlnPlot_CT.pdf', width=15, height=8)
    print( VlnPlot(srat, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_FeaturePlot_CT.pdf', width=16, height=16)
    print(FeaturePlot(srat, reduction = 'tsDMSO', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()

   # remove ambient RNAs
   library(Matrix)
   data_full<-as(as.matrix(CT_counts),"dgCMatrix")
   #data_full<-Matrix(CT_counts, sparse = TRUE) 
   data_empty<-raw.matrix$`GeDMSO Expression`
   data_empty<-data_empty[rownames(data_full),]
   data_full_sce <- SingleCellExperiment(assays = list(counts = data_full))
   data_empty_sce <- SingleCellExperiment(assays = list(counts = data_empty))
   data_full_sce <- decontX(x = data_full_sce, background = data_empty_sce)
   sce.decontX<-data_full_sce
   # Cluster labels on UMAP
   umap <- reducedDim(sce.decontX, "decontX_UMAP")
   pdf('./decontX_outs/decontX_contamination_CT.pdf',width=15, height=6)
   print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
           plotDecontXContamination(sce.decontX))
   dev.off()
   DropletUtils:::write10xCounts("./decontX_outs/CT_decontXcounts", round(decontXcounts(sce.decontX)),
                                 barcodes = rownames(colData(data_full_sce)))
       
    #########################################
    ### decontX
    sratDecontx  <- 
    Read10X("./decontX_outs/CT_decontXcounts") %>%
    CreateSeuratObject(project = "CT", min.cells = 3, min.features = 200)
    sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-")
    sratDecontx    <- SCTransform(sratDecontx, verbose = F)
    sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
    sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- RunTSDMSO(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindDMSOighbors(sratDecontx, dims = 1:30, reduction = 'pca',verbose = F)
    sratDecontx    <- FindClusters(sratDecontx, verbose = T,resolution =1.2)
    pdf('./decontX_outs/decontX_UMAP_tSDMSO_CT.pdf', width=15, height=6)
    print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsDMSO", label = TRUE) )
    dev.off()
    pdf('./decontX_outs/decontX_coreceptor_VlnPlot_CT.pdf',width=15, height=8)
    print( VlnPlot(sratDecontx, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    pdf('./decontX_outs/decontX_coreceptor_FeaturePlot_CT.pdf', width=16, height=16)
    print(FeaturePlot(sratDecontx, reduction = 'tsDMSO', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
```

* DMSO

``` 
DMSO_data<-read.table("./counts_no_double_umi_001.tsv.gz")
DMSO_counts = spread(DMSO_data, V3, V2)
DMSO_counts[is.na(DMSO_counts)] <- 0
rownames(DMSO_counts)<-DMSO_counts$V1
DMSO_counts<-DMSO_counts[,-1]
DMSO_barcode<-paste(colnames(DMSO_counts),"-1",sep="")
colnames(DMSO_counts)<-DMSO_barcode

srat <-  CreateSeuratObject(counts = DMSO_counts, min.cells=3,  project="DMSO", assay = "RNA")

raw.matrix  <- Read10X_h5("./raw_feature_bc_matrix.h5",use.names = T)

srat.raw <- CreateSeuratObject(
        counts = raw.matrix$`GeDMSO Expression`,
        assay = "RNA"
      )

    srat    <- SCTransform(srat, verbose = F)
    # Idents(srat) <- srat@meta.data$seurat_clusters
    srat    <- RunPCA(srat, npcs = 50, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
    srat    <- RunTSDMSO(srat, dims = 1:50, verbose = F)
    srat    <- FindDMSOighbors(srat,reduction = 'pca', dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T,resolution =1.2)

    pdf('./decontX_outs/Original_UMAP_tSDMSO_DMSO.pdf',width=15, height=6)
    print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsDMSO", label = TRUE))
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_VlnPlot_DMSO.pdf', width=15, height=8)
    print( VlnPlot(srat, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_FeaturePlot_DMSO.pdf', width=16, height=16)
    print(FeaturePlot(srat, reduction = 'tsDMSO', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()

   # remove ambient RNAs
   library(Matrix)

   data_full<-as(as.matrix(DMSO_counts),"dgCMatrix")
   #data_full<-Matrix(DMSO_counts, sparse = TRUE) 

   data_empty<-raw.matrix$`GeDMSO Expression`
   data_empty<-data_empty[rownames(data_full),]
   data_full_sce <- SingleCellExperiment(assays = list(counts = data_full))
   data_empty_sce <- SingleCellExperiment(assays = list(counts = data_empty))
   data_full_sce <- decontX(x = data_full_sce, background = data_empty_sce)
   sce.decontX<-data_full_sce
    # Cluster labels on UMAP
    umap <- reducedDim(sce.decontX, "decontX_UMAP")

    pdf('./decontX_outs/decontX_contamination_DMSO.pdf',width=15, height=6)
    print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
            plotDecontXContamination(sce.decontX))
    dev.off()

    DropletUtils:::write10xCounts("./decontX_outs/DMSO_decontXcounts", round(decontXcounts(sce.decontX)),
                                  barcodes = rownames(colData(data_full_sce)))
        
    #########################################
    ### decontX
      sratDecontx  <- 
      Read10X("./decontX_outs/DMSO_decontXcounts") %>%
      CreateSeuratObject(project = "DMSO", min.cells = 3, min.features = 200)
    
    sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-")
    

    sratDecontx    <- SCTransform(sratDecontx, verbose = F)
    sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
    sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- RunTSDMSO(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindDMSOighbors(sratDecontx, dims = 1:30, reduction = 'pca',verbose = F)
    sratDecontx    <- FindClusters(sratDecontx, verbose = T,resolution =1.2)
    
    
    pdf('./decontX_outs/decontX_UMAP_tSDMSO_DMSO.pdf', width=15, height=6)
    print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsDMSO", label = TRUE) )
    dev.off()
    
    pdf('./decontX_outs/decontX_coreceptor_VlnPlot_DMSO.pdf',width=15, height=8)
    print( VlnPlot(sratDecontx, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/decontX_coreceptor_FeaturePlot_DMSO.pdf', width=16, height=16)
    print(FeaturePlot(sratDecontx, reduction = 'tsDMSO', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
```



