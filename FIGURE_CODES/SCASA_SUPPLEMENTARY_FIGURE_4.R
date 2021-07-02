################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")

library("seurat")
library("ggplot2")
library("patchwork")

################################################################################
# SUPPLEMENTARY FIGURE 4
################################################################################
bm <- readRDS("FIGURE_DATA/Bone_Marrow_Seurat_All_Dim_Reductions.RDS")
col1 <- c("#4D9DE0","#698FC5","#8581AA","#A2738F","#BE6574","#DB5759","#E1644D","#E17845","#E18C3C","#E1A034","#E1B42C",
          "#CDBA31","#ADB83F","#8EB74D","#6EB55C","#4EB36A","#3FAC77","#4B9E82","#568F8E","#628199","#6D73A4","#7466AC",
          "#695DA6","#5E54A0","#534C9A","#484394","#3D3B8E")
col2 <- c("#F9C80E","#F8B512","#F8A216","#F88F1A","#F87C1E","#F86923","#F55E29","#F3552F","#F04B36","#ED423C","#EB3843",
          "#DA344F","#C13260","#A83170","#8E3080","#752E91","#63389E","#5C54A8","#556FB2","#4F8ABB","#48A6C5","#43B8CB",
          "#47A9C3","#4B99BB","#4F8AB3","#537AAB","#586BA4")

plotx <- bm
plotx$rna.umap@cell.embeddings[,1]=-plotx$rna.umap@cell.embeddings[,1] #flip

ps41 <- DimPlot(plotx, reduction = 'rna.umap', group.by = 'Celltype', cols = col1, pt.size = 1,label = TRUE, repel = TRUE, label.size = 5)+ ggtitle("SCASA (ISOFORM) CELL-TYPES")+
NoLegend()+ labs(x = NULL, y=NULL)+ theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+panel_border()+
  theme(title = element_text(size = 20, face = "bold"))

ps42 <- DimPlot(plotx, reduction = 'rna.umap', group.by = 'ClusterSCiso', cols = col2, pt.size = 1, label = TRUE, repel = TRUE, label.size = 6)  + ggtitle("SCASA (ISOFORM) CLUSTERS")+
NoLegend()+ labs(x = NULL, y=NULL)+ theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+panel_border()+
  theme(title = element_text(size = 20, face = "bold"))

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_4.png", width=4000, height=2500, units = "px", res = 300)  
print(ps41+ps42)
dev.off()
################################################################################