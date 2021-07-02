################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")

library("Seurat")
library("ggplot2")
library("cowplot")
library("dplyr")

################################################################################
# SUPPLEMENTARY FIGURE 6
################################################################################
load("FIGURE_DATA/CITESeq/smartseq2_bonemarrow_isoformExpresion.RData")

bm=CreateSeuratObject(counts = salmon_isoform_count, assay = "smartseq2")
bm <- NormalizeData(bm) #normalisation
bm <- FindVariableFeatures(bm) #find variable features
bm <- ScaleData(bm) #scaling
bm <- RunPCA(bm,reduction.name = 'pca')#PCA
bm <- RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'smartseq2', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
bm <- FindNeighbors(bm, reduction = "pca", dims = 1:30, verbose = TRUE,graph.name="smartseq2NB")
bm <- FindClusters(bm, graph.name="smartseq2NB",verbose = TRUE,algorithm = 3, resolution = 2)
bm@meta.data$ClusterSmartSeq2=Idents(bm)
bm@meta.data$Celltype=celltype

cols <- c("#F9844A","#C07E62","#88797A","#5D788D","#959476","#CDB060","#E8C155","#A0A771","#578E8E","#3F7696","#876176",
          "#D04C55","#E8544A","#C47F58","#A0AA66","#A4B65D","#C8A842","#EC9A27","#F68D21","#F58126","#F3742A","#C38145",
          "#869466","#49A787","#46A18B","#49998C","#4D908E")

ps6a <- DimPlot(bm, reduction = 'rna.umap', group.by = 'Celltype', cols = cols, label = TRUE, repel = TRUE, label.size = 5, pt.size=2)+ ggtitle("SMART-SEQ2 (ISOFORM) - CELL TYPES")+ xlab("UMAP1")+ylab("UMAP2") + NoLegend()
ps6b <- DimPlot(bm, reduction = 'rna.umap', group.by = 'ClusterSmartSeq2', cols = cols, label = TRUE, repel = TRUE, label.size = 5,pt.size=2)  + ggtitle("SMART-SEQ2 (ISOFORM) - CLUSTERS")+ xlab("UMAP1")+ylab("UMAP2")+ NoLegend()

################################################################################
# SUPPLEMENTARY FIGURE 6a
################################################################################

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_6A.png", width=2000, height=2000, units = "px", res = 300)  
print(ps6a)
dev.off()

################################################################################
# SUPPLEMENTARY FIGURE 6b
################################################################################
png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_6B.png", width=2000, height=2000, units = "px", res = 300)  
print(ps6b)
dev.off()
################################################################################