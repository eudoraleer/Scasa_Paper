################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")

library("Seurat")
library("ggplot2")
library("patchwork")
library("cowplot")
################################################################################
# SUPPLEMENTARY FIGURE 3
################################################################################
################################################################################
# SUPPLEMENTARY FIGURE 3a
################################################################################
source("FIGURE_CODES/FUNCTIONS/LOAD_ALL_DATA.R")
equ = equrow(scasaAle_tgene,scasaAle_gene) 
pick = equ[[1]]>0 | equ[[2]]>0
set.seed(2021)
p = 0.10
sam = sample(1:length(pick), p*length(pick)) ## subsample for fast plotting
pick[-sam] = FALSE;   print(sum(pick))

data <- data.frame(TRUE_COUNT = c(equ[[1]][pick])+1,
                   ESTIMATED_COUNT = c(equ[[2]][pick])+1,
                   DATA = "SCASA (GENE)")

all_gene = c("tcount_gene","cr_mat","klbus_mat",
             "alevin_mat", "kallisto_gene", "salmon_gene")
current_names <- c("SCASA (GENE)","CELLRANGER (GENE)","KALLISTO-BUS (GENE)","ALEVIN (GENE)","KALLISTO (GENE)","SALMON (GENE)")
matlist = mget(all_gene)
equ = equrow_col(matlist)

for (i in 2:6){
  pick = equ[[1]]>0 | equ[[i]]>0
  set.seed(2021)
  p = 0.10
  sam = sample(1:length(pick), p*length(pick))
  pick[-sam] = FALSE 
  print(sum(pick))
  data <- rbind(data, data.frame(TRUE_COUNT = c(equ[[1]][pick])+1,
                                 ESTIMATED_COUNT = c(equ[[i]][pick])+1,
                                 DATA = current_names[i]))
}

data$DATA <- factor(data$DATA, levels = current_names)
colors <- c("#9d4edd","#86BBD8","#F6AE2C","#F26419","#33658A","#69995D")
ps3a <- ggplot(data, aes(TRUE_COUNT, ESTIMATED_COUNT, color = DATA)) +
  geom_point()+scale_x_log10()+scale_y_log10()+ theme_classic()+
  xlab("TRUE COUNT + 1")+ylab("ESTIMATED COUNT + 1")+
  facet_wrap(~DATA, scales = "free")+ scale_color_manual(values = colors)+
  geom_abline(intercept = 0, color = "#968147")+
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 35, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 35, margin=margin(0,10,0,0)),
        legend.title = element_text(size =20, face = "bold"),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.background = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"))

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_3A.png", width=5000, height=2800, units = "px", res = 300)  
print(ps3a)
dev.off()

################################################################################
# SUPPLEMENTARY FIGURE 3b
################################################################################
bm <- readRDS("FIGURE_DATA/Bone_Marrow_Seurat_All_Dim_Reductions.RDS")
cols <- c("#4D9DE0","#698FC5","#8581AA","#A2738F","#BE6574","#DB5759","#E1644D","#E17845","#E18C3C","#E1A034","#E1B42C",
          "#CDBA31","#ADB83F","#8EB74D","#6EB55C","#4EB36A","#3FAC77","#4B9E82","#568F8E","#628199","#6D73A4","#7466AC",
          "#695DA6","#5E54A0","#534C9A","#484394","#3D3B8E")

plotx <- bm
plotx$rna.umap@cell.embeddings[,1]=-plotx$rna.umap@cell.embeddings[,1] #flip
plotx$SCgene.umap@cell.embeddings[,2]=-plotx$SCgene.umap@cell.embeddings[,2] #flip

ps3b1 <- DimPlot(plotx, reduction = 'rna.umap', group.by = 'Celltype', label = TRUE, repel = TRUE, cols = cols, label.size = 4, pt.size = 1)+ ggtitle("SCASA (ISOFORM)")+ xlab("UMAP_1")+ylab("UMAP_2") + NoLegend()+ labs(x = NULL, y=NULL)+ theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+panel_border()
ps3b2 <- DimPlot(plotx, reduction = 'SCgene.umap', group.by = 'Celltype', label = TRUE, repel = TRUE, cols = cols,label.size = 4, pt.size = 1)+ ggtitle("SCASA (GENE)")+ xlab("UMAP_1")+ylab("UMAP_2") + NoLegend()+ labs(x = NULL, y=NULL)+ theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+panel_border()
ps3b3 <- DimPlot(plotx, reduction = 'cr.umap', group.by = 'Celltype', label = TRUE, repel = TRUE, cols = cols,label.size = 4, pt.size = 1)+ ggtitle("CELLRANGER (GENE)")+ xlab("UMAP_1")+ylab("UMAP_2") + NoLegend() + labs(x = NULL, y=NULL) + theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+panel_border()
ps3b4 <- DimPlot(plotx, reduction = 'Alevin.umap', group.by = 'Celltype', label = TRUE, repel = TRUE, cols = cols,label.size = 4, pt.size = 1)+ ggtitle("ALEVIN (GENE)")+ xlab("UMAP_1")+ylab("UMAP_2") + NoLegend() + labs(x = NULL, y=NULL) + theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+panel_border()

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_3B1.png", width=1700, height=1700, units = "px", res = 350)  
print(ps3b1)
dev.off()

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_3B2.png", width=1700, height=1700, units = "px", res = 350)  
print(ps3b2)
dev.off()

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_3B3.png", width=1700, height=1700, units = "px", res = 350)  
print(ps3b3)
dev.off()

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_3B4.png", width=1700, height=1700, units = "px", res = 350)  
print(ps3b4)
dev.off()

################################################################################
# SUPPLEMENTARY FIGURE 3c
################################################################################
ref <- read.csv("FIGURE_DATA/transcripts_to_genes_annotation.txt", header=F, sep="\t")
assays <- c('SCiso','SCgene','CR','Alevin')
umaps <- c('rna.umap','SCgene.umap','cr.umap','Alevin.umap')
titles <- c('SCASA (ISOFORM)','SCASA (GENE)','CELLRANGER (GENE)','ALEVIN (GENE)')
output_names <- c('SCASA_ISOFORM','SCASA_GENE','CELLRANGER_GENE','ALEVIN_GENE')

for(i in 1:length(assays)){
  plotx <- bm
  DefaultAssay(plotx) <- assays[i]
  if(assays[i] == "SCgene"){
    plotx@reductions[[umaps[i]]]@cell.embeddings[,2] <- -plotx@reductions[[umaps[i]]]@cell.embeddings[,2]
  }
  TYROBP <- ref[which(ref$V2 == "TYROBP"),]
  pos <- NULL
  for(k in 1:nrow(TYROBP)){
    pos <- c(pos,grep(gsub("_","-",TYROBP$V1[k]),row.names(plotx)))
  }
  pos <- unique(pos)
  
  if(assays[i] != "SCiso"){
    p <- FeaturePlot(plotx, features = "TYROBP", reduction = umaps[i], pt.size = 1) + 
      ggtitle(titles[i])+
      theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
            plot.title = element_text(size =35, face = "bold", hjust = 0.5))+NoLegend()+ labs(x = NULL, y=NULL)
    
    png(paste("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_3C_",output_names[i],"_TYROBP_UMAP.png",sep = ""), width=3000, height=2800, units = "px", res = 300)  
    print(p)
    dev.off()
  }
}
################################################################################
