################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")

library("Seurat")
library("ggplot2")
library("cowplot")
library("dplyr")
library("reshape2")
library("viridis")

################################################################################
# SUPPLEMENTARY FIGURE 7
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

pick1=which(bm@meta.data$ClusterSmartSeq2 %in% c(19)) #TY32.25
ct1=bm@meta.data$Celltype[pick1]
names(ct1)=rownames(bm@meta.data)[pick1]
sort(table(ct1),decreasing=TRUE)
cellid1=names(ct1)

pick2=which(bm@meta.data$Celltype %in% "Monocytes")
ct2=bm@meta.data$Celltype[pick2]
names(ct2)=rownames(bm@meta.data)[pick2]
ct2=ct2[!(names(ct2) %in% cellid1)]
sort(table(ct2),decreasing=TRUE)
cellid2=names(ct2)

fdrThres=0.05
clMarker_Smartseq2 <- FindMarkers(bm$smartseq2, cells.1 = cellid1, cells.2 = cellid2, min.pct = 0.25,test.use="LR")
reformatIsoforms=gsub("-","_",rownames(clMarker_Smartseq2))
genes=genes.tx.map[reformatIsoforms]
clMarker_Smartseq2$GENEID=genes

sum(clMarker_Smartseq2$p_val_adj < fdrThres)
clMarker_Smartseq2=clMarker_Smartseq2[clMarker_Smartseq2$p_val_adj < fdrThres,]
dim(clMarker_Smartseq2)
head(clMarker_Smartseq2)

tx_TYROBP=c("NM-003332","NM-198125","NM-001173514","NR-033390","NM-001173515")
plotData=bm$smartseq2@data[tx_TYROBP,]
plotData=as.matrix(plotData)
plotCt=bm@meta.data$Celltype
plotIdents=as.character(bm@meta.data$ClusterSmartSeq2)
plotIdents[plotIdents=="19"]="TY32.25 Mono"
plotIdents[plotIdents%in%c("13")]="Monocytes"
clusters=names(table(plotIdents))
meanExp=apply(plotData,1,function(x){
  y=tapply(x,plotIdents,function(z){mean(z)})
})

pctExp=apply(plotData,1,function(x){
  y=tapply(x,plotIdents,function(z){sum(z>0)/length(z)})
})

x=rowSums(meanExp)

#order by total reads
myorder=order(x,decreasing=TRUE)
meanExp=meanExp[myorder,]
pctExp=pctExp[myorder,]

cluster_ct=table(plotCt,plotIdents)
cluster_ct=unclass(cluster_ct)
cluster_ct_id=apply(cluster_ct,2,which.max)
cluster_ct_label=rownames(cluster_ct)[cluster_ct_id]
names(cluster_ct_label)=names(cluster_ct_id)

cluster_ct_label_adj=cluster_ct_label
for (i in 1:length(cluster_ct_label)){
  isDup=which(cluster_ct_label==cluster_ct_label[i])
  if (length(isDup)>1){
    cluster_ct_label_adj[isDup]=paste0(cluster_ct_label[i],"_",c(1:length(isDup)))
  }
}

cluster_ct_label_adj[names(cluster_ct_label_adj)=="TY32.25 Mono"]="TY32.25 Mono"
cluster_ct_label_adj[names(cluster_ct_label_adj)=="Monocytes"]="Monocytes"
clusterLabelMap=data.frame(cluster_ct_id,cluster_ct_label,cluster_ct_label_adj)
clusterLabelMap$plotIdents=rownames(clusterLabelMap)

meanExp_plot=t(meanExp)
pctExp_plot=t(pctExp)
matchID=match(names(cluster_ct_label_adj),colnames(meanExp_plot))
colnames(meanExp_plot)[matchID]=cluster_ct_label_adj
colnames(pctExp_plot)[matchID]=cluster_ct_label_adj

myorder=c(1,5,c(2:4,6:ncol(meanExp_plot)))
meanExp_plot=meanExp_plot[,myorder]
pctExp_plot=pctExp_plot[,myorder]
topk=ncol(meanExp_plot)
meanExp_plot=meanExp_plot[,1:topk]
pctExp_plot=pctExp_plot[,1:topk]

plotx <- melt(meanExp_plot)
colnames(plotx) <- c("ISOFORM","CELL_TYPE","MEAN_EXPR")
current <- melt(pctExp_plot)
colnames(current) <- c("ISOFORM","CELL_TYPE","CELL_PROPORTION")
plotx$CELL_PROPORTION <- current[match(paste(plotx$ISOFORM, plotx$CELL_TYPE, sep = "_"),
                                       paste(current$ISOFORM, current$CELL_TYPE, sep = "_")),"CELL_PROPORTION"]

ps7 <- ggplot(plotx, aes(ISOFORM,CELL_TYPE,color=MEAN_EXPR, size = CELL_PROPORTION)) + 
  geom_point(shape=15, stroke = 7) +
  scale_color_viridis(option = "C", direction = -1, guide_legend(title="MEAN EXPRESSION")) +
  theme_classic()+ ylab("CELL TYPE")+
  guides(size = guide_legend(title="CELL PROPORTION"))+
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1,vjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "right",
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
        legend.title = element_text(size =20, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size =30, face = "bold", hjust = 0.5))

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_7.png", width=4000, height=3000, units = "px", res = 300)  
print(ps7)
dev.off()
################################################################################