## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

#cluster annotation
#基因占比图
DefaultAssay(ORG_CT.integrated) <- "RNA"

markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1", # cardiomyocyte
"Ccl3","Clec4a1","Rgs1","Fcgr1","Cd14","Csf1r","Cd163","Cd68","Itgam","Lgals3","Mrc1", # Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Igfbp4","Lat","Itk","Cd3g", # T cell
"Cd22","Cd79a","Cd79b","Mzb1","Ly6d","Ms4a1", # B cell
"Cd74","Cd83","Cd86","Flt3","Cd209a","Ccl17", # DC
"Plac8", # monocyte
"S100a9", "S100a8",#Granulocy
"Nkg7","Gzma", #NK
"Chl1","Kcna6", # Glial
"Myct1","Cdh5","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1","Tie1","Fabp4","Esam","Kdr","Tek", # endothelial
"Dpt","Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra", # fibroblast
"Abcc9","Rgs4","Ano1","Acta2","Myh11", # smoothmuscle
"Msln","Krt19","Plet1","Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8",# epicardial
"Gfpt2","Hk2","Ddr2","Slc1a7","Adamtsl3","Layn","Cd248","Mcam","Fgfr1" # Pericyte
)

label<- c(rep("CM",7),rep("MP",11),rep("T",8),rep("B",6),rep("DC",6),rep("monocyte",1),rep("Granulocy",2),rep("NK",2),rep("Glial",2),rep("EC",13),
	rep("FB",8),rep("SMC",5),rep("Epi",9),rep("Pericyte",9))
DefaultAssay(ORG_CT.integrated)<-"RNA"
pdf("./02_annotation/ORG_cluster-annotation-all_celltype.pdf",width=20,height=8)
p<-DotPlot(ORG_CT.integrated, features = markerGenes,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()
pdf("./02_annotation/ORG_cluster-annotation-all_celltype_cluster.pdf",width=20,height=8)
Clustered_DotPlot(seurat_object = ORG_CT.integrated, flip=T,features = markerGenes,k=11)
dev.off()

# make the tree 

# make the trans dist tree 
object <- ORG_CT.integrated
embeddings <- Embeddings(object = object, reduction = "lsi")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

pdf("./03_all_celltype/AR3-combined-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


Idents(combined)<-combined$seurat_clusters
#####further annotation########
combined <- RenameIdents(
  object = combined,
  '0' = 'MP_DC',
  '1' = 'MP_DC',
  '2' = 'FB',
  '3' = 'MP_DC',
  '4' = 'EC',
  '5' = 'CM',
  '6' = 'CM',
  '7' = 'FB',
  '8' = 'EC',
  '9' = 'CM',
  '10' = 'CM',
  '11' = 'FB',
  '12' = 'FB',
  '13' = 'FB',
  '14' = 'Pericyte',
  '15' = 'FB',
  '16' = 'FB',
  '17' = 'EC',
  '18' = 'Granulocy',
  '19' = 'FB',
  '20' = 'MP_DC',
  '21' = 'Epi',
  '22' = 'EC',
  '23' = 'T',
  '24' = 'FB',
  '25' = 'MP_DC',
  '26' = 'SMC',
  '27' = 'EC',
  '28' = 'MP_DC',
  '29' = 'Epi',
  '30' = 'CM',
  '31' = 'Epi',
  '32' = 'FB',
  '33' = 'EC',
  '34' = 'SMC',
  '35' = 'B',
  '36' = 'T',
  '37' = 'MP_DC',
  '38' = 'Glial',
  '39' = 'EC',
  '40' = 'MP_DC',
  '41' = 'MP_DC'
  )
combined@meta.data$Annotation<-Idents(combined)
table(combined$Annotation,combined$orig.ident)

combined$Annotation<-factor(combined$Annotation,levels=c("CM",'EC','FB','Epi','SMC','Pericyte',"MP_DC","T","B","Gra","Glial"))
Idents(combined)<-combined$Annotation;
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
combined$detail_anno<- paste(combined$seurat_clusters,combined$Annotation,sep=":")

pdf("./03_all_celltype/Annotated_allcelltype_UMAP.pdf",width=6,height=5)
DimPlot(combined, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "Annotation")
dev.off()
pdf("./03_all_celltype/Annotated_allcelltype_UMAP_detail.pdf",width=8.5,height=5)
DimPlot(combined, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "detail_anno")
dev.off()
object <- combined
Idents(object)<- object$detail_anno
embeddings <- Embeddings(object = object, reduction = "lsi")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

pdf("./03_all_celltype/AR3-combined-tree-cosine_detail_anno.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


# Save rds have annotation information 
DefaultAssay(combined) <- "RNA"
saveRDS(combined,"./03_all_celltype/AR3_integrated_all_celltype_annotated.rds")

# Marker gene for different cell types (color on UMAP ) 
combined<- readRDS("./03_all_celltype/AR3_integrated_all_celltype_annotated.rds")
pdf('./03_all_celltype/All_Marker_gene_FeaturePlot_uamp.pdf', width=5.5, height=5)
for (i in 1:length(markerGenes)){
  p<-FeaturePlot(combined,order=F, reduction = 'umap',max.cutoff = 10, features = markerGenes[i], ncol = 1)
  print(p)
}
dev.off()



