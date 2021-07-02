################################################################################
## SCASA: MAIN FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")

library("ggplot2")
library("patchwork")
library("Seurat")
library("SeuratData")
library("cowplot")
library("dplyr")
library("viridis")
library("reshape2")
InstallData("bmcite")
bm_seurat <- LoadData(ds = "bmcite")
memory.limit(size=1e10)
################################################################################
# FIGURE 2
################################################################################
metadata=bm_seurat@meta.data
adtcount=bm_seurat$ADT@counts
seuratGene=bm_seurat$RNA@counts
bm_seurat=NULL

load("FIGURE_DATA/CITESeq/Scasa_isoform_citeseq.RData")
SCASA.isomat=isoform_count
isoform_count=NULL

load("FIGURE_DATA/CITESeq/Scasa_gene_citeseq.RData")
dim(gene_count)
SCASA.genemat=gene_count
gene_count=NULL

load("FIGURE_DATA/CITESeq/Alevin_citeseq.RData")
dim(Alevin_count)

#keep only cells from SCASA with barcodes overlapping with the ones of ADT data
cbnames=sapply(colnames(adtcount), function(x){
  y=unlist(strsplit(x,"_"))[2]
  z=unlist(strsplit(y,"-"))[1]
  return(z)
})
names(cbnames)=NULL

colnames(adtcount)=cbnames
colnames(seuratGene)=cbnames
metadata2=t(metadata)
colnames(metadata2)=cbnames
# shared cell IDs between SCASA and cellranger
mycl=intersect(colnames(adtcount),colnames(SCASA.isomat))
length(mycl)
# keep the overlapping cells
matchID=match(mycl,colnames(adtcount))
adtcount=adtcount[,matchID]
metadata=metadata[matchID,]
seuratGene=seuratGene[,matchID]
pick=rowSums(seuratGene)>0
seuratGene=seuratGene[pick,]
rownames(metadata)=colnames(adtcount)
dim(adtcount)

SCASA.genemat=SCASA.genemat[,match(mycl,colnames(SCASA.genemat))]
pick=rowSums(SCASA.genemat)>0
SCASA.genemat=SCASA.genemat[pick,]
dim(SCASA.genemat)

SCASA.isomat=SCASA.isomat[,match(mycl,colnames(SCASA.isomat))]
pick=rowSums(SCASA.isomat)>0
SCASA.isomat=SCASA.isomat[pick,]
dim(SCASA.isomat)

Alevin_count=Alevin_count[,mycl]
pick=rowSums(Alevin_count)>0
Alevin_count=Alevin_count[pick,]
dim(Alevin_count)

#############################################
# create objects containing gene/isoform expression of different methods
bm=CreateSeuratObject(counts = SCASA.isomat, assay = "SCiso")
SCASA.isomat=NULL
bm$SCgene=CreateAssayObject(counts = SCASA.genemat)
SCASA.genemat=NULL
bm$CR=CreateAssayObject(counts = seuratGene)
seuratGene=NULL
bm$Alevin=CreateAssayObject(counts = Alevin_count)
Alevin_count=NULL

#### normalisation, feature selection and run PCA
#process isoform-level gene expression of scasa 
DefaultAssay(bm) <- 'SCiso'
bm <- NormalizeData(bm) #normalisation
bm <- FindVariableFeatures(bm) #find variable features
bm <- ScaleData(bm) #scaling
bm <- RunPCA(bm,reduction.name = 'pca')#PCA

#process gene-level gene expession of scasa
DefaultAssay(bm) <- 'SCgene'
bm <- NormalizeData(bm) %>% FindVariableFeatures() 
bm <- ScaleData(bm) #scaling
bm <- RunPCA(bm,reduction.name = 'SCgenepca')#PCA

#process gene-level gene expession of cell ranger
DefaultAssay(bm) <- 'CR'
bm <- NormalizeData(bm) %>% FindVariableFeatures() 
bm=ScaleData(bm)
bm=RunPCA(bm,reduction.name = 'crpca')

#process gene-level gene expession of alevin
DefaultAssay(bm) <- 'Alevin'
bm <- NormalizeData(bm) %>% FindVariableFeatures() 
bm=ScaleData(bm)
bm=RunPCA(bm,reduction.name = 'Alevinpca')

DefaultAssay(bm) <- 'SCiso'
DefaultAssay(bm) <- 'SCgene'
DefaultAssay(bm) <- 'CR'
DefaultAssay(bm) <- 'Alevin'

bm <- RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA',
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
bm <- RunUMAP(bm, reduction = 'SCgenepca', dims = 1:30, assay = 'SCgene',
              reduction.name = 'SCgene.umap', reduction.key = 'SCgeneUMAP_')
bm <- RunUMAP(bm, reduction = 'crpca', dims = 1:30, assay = 'CR', 
              reduction.name = 'cr.umap', reduction.key = 'crUMAP_')
bm <- RunUMAP(bm, reduction = 'Alevinpca', dims = 1:30, assay = 'Alevin',
              reduction.name = 'Alevin.umap', reduction.key = 'AlevinUMAP_')

### find clusters
#scasa isoform
bm <- FindNeighbors(bm, reduction = "pca", dims = 1:30, verbose = TRUE,graph.name="SCisoNB")
bm <- FindClusters(bm, graph.name="SCisoNB",verbose = TRUE,algorithm = 3, resolution = 2)
bm@meta.data$ClusterSCiso=Idents(bm)

#scasa gene
bm <- FindNeighbors(bm, reduction = "SCgenepca", dims = 1:30, verbose = TRUE,graph.name="SCgeneNB")
bm <- FindClusters(bm, graph.name="SCgeneNB",verbose = TRUE,algorithm = 3, resolution = 2)
bm@meta.data$ClusterSCgene=Idents(bm)

#cellranger
bm <- FindNeighbors(bm, reduction = "crpca", dims = 1:30, verbose = TRUE,graph.name="crNB")
bm <- FindClusters(bm, graph.name="crNB",verbose = TRUE,algorithm = 3, resolution = 2)
bm@meta.data$ClusterCR=Idents(bm)

#alevin
bm <- FindNeighbors(bm, reduction = "Alevinpca", dims = 1:30, verbose = TRUE,graph.name="AlevinNB")
bm <- FindClusters(bm, graph.name="AlevinNB",verbose = TRUE,algorithm = 3, resolution = 2)
bm@meta.data$ClusterAlevin=Idents(bm)

### add metadata to the object to get cell type information
bm@meta.data=cbind(bm@meta.data, metadata[,c("lane","donor","celltype.l2")])
bm@meta.data$Celltype=bm@meta.data$celltype.l2

# saveRDS(bm, "FIGURE_DATA/Bone_Marrow_Seurat_All_Dim_Reductions.RDS")

plotx <- bm
plotx$rna.umap@cell.embeddings[,1]=-plotx$rna.umap@cell.embeddings[,1] #flip
################################################################################
# FIGURE 2a
################################################################################
cols <- c("#4D9DE0","#698FC5","#8581AA","#A2738F","#BE6574","#DB5759","#E1644D","#E17845","#E18C3C","#E1A034","#E1B42C",
          "#CDBA31","#ADB83F","#8EB74D","#6EB55C","#4EB36A","#3FAC77","#4B9E82","#568F8E","#628199","#6D73A4","#7466AC",
          "#695DA6","#5E54A0","#534C9A","#484394","#3D3B8E")
p2a <- DimPlot(plotx, reduction = 'rna.umap', group.by = 'Celltype', cols = cols,label = TRUE,
               repel = TRUE, label.size = 4, pt.size = 1)+ ggtitle("SCASA (ISOFORM)")+ 
               NoLegend()+ labs(x = NULL, y=NULL)+
               theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+panel_border()

png("FIGURES_AND_TABLES/MAIN_FIGURE_2A.png", width=1500, height=2700, units = "px", res = 350)  
print(p2a)
dev.off()
################################################################################

################################################################################
# FIGURE 2b
################################################################################
### identify biomarkers
pick1=which(bm@meta.data$ClusterSCiso %in% c(15)) #new subgroup of CD14
pick2=which(bm@meta.data$ClusterSCiso %in% c(0,11,17))#main subgroup of CD14
length(pick1)
length(pick2)
#extract cell types of two groups
myct=bm@meta.data$celltype.l2
ct1=bm@meta.data$celltype.l2[pick1]
names(ct1)=rownames(bm@meta.data)[pick1]
ct2=bm@meta.data$celltype.l2[pick2]
names(ct2)=rownames(bm@meta.data)[pick2]
sort(table(ct1),decreasing=TRUE)
sort(table(ct2),decreasing=TRUE)
cellid1=names(ct1)
cellid2=names(ct2)

fdrThres=0.05
#cellranger
clMarker_CR <- FindMarkers(bm$CR, cells.1 = cellid1, cells.2 = cellid2, min.pct = 0.25,test.use="LR")
clMarker_CR=clMarker_CR[clMarker_CR$p_val_adj < fdrThres,]
#alevin
clMarker_Alevin <- FindMarkers(bm$Alevin, cells.1 = cellid1, cells.2 = cellid2, min.pct = 0.25,test.use="LR")
clMarker_Alevin=clMarker_Alevin[clMarker_Alevin$p_val_adj < fdrThres,]
#scasa-gene
clMarker_SCgene <- FindMarkers(bm$SCgene, cells.1 = cellid1, cells.2 = cellid2, min.pct = 0.25,test.use="LR")
clMarker_SCgene=clMarker_SCgene[clMarker_SCgene$p_val_adj < fdrThres,]
#scasa-isoform
DefaultAssay(bm) <- 'SCiso'
clMarker_SCiso <- FindMarkers(bm$SCiso, cells.1 = cellid1, cells.2 = cellid2, min.pct = 0.25,test.use="LR")
reformatIsoforms=gsub("-","_",rownames(clMarker_SCiso))
#get isoform-gene mapping
load("FIGURE_DATA/CITESeq/isoform_groups_UCSC_hg38.RData") 
genes=genes.tx.map.all.final[reformatIsoforms]
clMarker_SCiso$GENEID=genes
pick=which(clMarker_SCiso$p_val_adj < fdrThres)
clMarker_SCiso$GENEID[pick]
clMarker_SCiso=clMarker_SCiso[pick,]

### export biomarkers to files
write.table(clMarker_SCiso,"FIGURES_AND_TABLES/clMarker_SCiso.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(clMarker_SCgene,"FIGURES_AND_TABLES/clMarker_SCgene.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(clMarker_CR,"FIGURES_AND_TABLES/clMarker_CR.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(clMarker_Alevin,"FIGURES_AND_TABLES/clMarker_Alevin.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

### plot heatmap in figure 2E
DEiso_list=rownames(clMarker_SCiso)
pick=rownames(bm$SCiso@data) %in% DEiso_list
plotData=bm$SCiso@data[pick,]
matchID=match(DEiso_list,rownames(plotData))

plotData=plotData[matchID,]
plotData=as.matrix(plotData)
plotCt=bm@meta.data$celltype.l2
plotIdents=as.character(bm@meta.data$ClusterSCiso)
#rename idents
plotIdents[plotIdents=="15"]="TY32.25 Mono"
plotIdents[plotIdents%in%c("0","11","17")]="CD14 Mono"
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
#View(cbind(cluster_ct_id,cluster_ct_label))
names(cluster_ct_label)=names(cluster_ct_id)

cluster_ct_label_adj=cluster_ct_label
for (i in 1:length(cluster_ct_label)){
  isDup=which(cluster_ct_label==cluster_ct_label[i])
  if (length(isDup)>1){
    cluster_ct_label_adj[isDup]=paste0(cluster_ct_label[i],"_",c(1:length(isDup)))
  }
}


cluster_ct_label_adj[names(cluster_ct_label_adj)=="TY32.25 Mono"]="TY32.25 Mono"
cluster_ct_label_adj[names(cluster_ct_label_adj)=="CD14 Mono"]="CD14 Mono"
clusterLabelMap=data.frame(cluster_ct_id,cluster_ct_label,cluster_ct_label_adj)
clusterLabelMap$plotIdents=rownames(clusterLabelMap)

meanExp_plot=t(meanExp)
pctExp_plot=t(pctExp)
matchID=match(names(cluster_ct_label_adj),colnames(meanExp_plot))
colnames(meanExp_plot)[matchID]=cluster_ct_label_adj
colnames(pctExp_plot)[matchID]=cluster_ct_label_adj
#remove too tiny clusters
x=unclass(table(plotIdents))
outliers=names(which(x<50))
ClusterOutliers=clusterLabelMap$cluster_ct_label_adj[clusterLabelMap$plotIdents %in% outliers]

meanExp_plot=meanExp_plot[,!(colnames(meanExp_plot)%in%ClusterOutliers)]
pctExp_plot=pctExp_plot[,!(colnames(pctExp_plot)%in%ClusterOutliers)]

myorder=c(1,3,2,c(4:ncol(meanExp_plot)))
meanExp_plot=meanExp_plot[,myorder]
pctExp_plot=pctExp_plot[,myorder]

# Shorten paralog name
rownames(pctExp_plot)[5]=rownames(meanExp_plot)[5]="NM-0012425(24-25) NM-033554"
plotx <- melt(meanExp_plot)
colnames(plotx) <- c("ISOFORM","CELL_TYPE","MEAN_EXPR")
current <- melt(pctExp_plot)
colnames(current) <- c("ISOFORM","CELL_TYPE","CELL_PROPORTION")
plotx$CELL_PROPORTION <- current[match(paste(plotx$ISOFORM, plotx$CELL_TYPE, sep = "_"),
                                       paste(current$ISOFORM, current$CELL_TYPE, sep = "_")),"CELL_PROPORTION"]

p2b <- ggplot(plotx, aes(ISOFORM,CELL_TYPE,color=MEAN_EXPR, size = CELL_PROPORTION)) + 
  geom_point(shape=15, stroke = 6) +
  scale_color_viridis(option = "D", direction = -1, guide_legend(title="MEAN EXPRESSION")) +
  theme_classic()+ ylab("CELL TYPE")+
  guides(size = guide_legend(title="CELL PROPORTION"))+
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust=1,vjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "right",
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
        legend.title = element_text(size =20, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.8, "cm"),
        strip.text.x = element_text(size = 25, face = "bold", ),
        strip.background = element_blank(),
        plot.title = element_text(size =30, face = "bold", hjust = 0.5))

png("FIGURES_AND_TABLES/MAIN_FIGURE_2B.png", width=3000, height=2200, units = "px", res = 200)  
print(p2b)
dev.off()

################################################################################

################################################################################
# FIGURE 2c
################################################################################
load("FIGURE_DATA/CITESeq/smartseq2_bonemarrow_isoformExpresion.RData")
bmsmartseq2 <- CreateSeuratObject(counts = salmon_isoform_count, assay = "smartseq2")
bmsmartseq2 <- NormalizeData(bmsmartseq2) #normalisation
bmsmartseq2 <- FindVariableFeatures(bmsmartseq2) #find variable features
bmsmartseq2 <- ScaleData(bmsmartseq2) #scaling
bmsmartseq2 <- RunPCA(bmsmartseq2,reduction.name = 'pca')#PCA

### run umap
bmsmartseq2 <- RunUMAP(bmsmartseq2, reduction = 'pca', dims = 1:30, assay = 'smartseq2', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

### find clusters
bmsmartseq2 <- FindNeighbors(bmsmartseq2, reduction = "pca", dims = 1:30, verbose = TRUE,graph.name="smartseq2NB")
bmsmartseq2 <- FindClusters(bmsmartseq2, graph.name="smartseq2NB",verbose = TRUE,algorithm = 3, resolution = 2)
bmsmartseq2@meta.data$ClusterSmartSeq2=Idents(bmsmartseq2)

### add metadata to the object to get cell type information
bmsmartseq2@meta.data$Celltype=celltype

### get biomarkers
pick1=which(bmsmartseq2@meta.data$ClusterSmartSeq2 %in% c(19)) #TY32.25
ct1=bmsmartseq2@meta.data$Celltype[pick1]
names(ct1)=rownames(bmsmartseq2@meta.data)[pick1]
sort(table(ct1),decreasing=TRUE)
cellid1=names(ct1)

#remaining monocytes
pick2=which(bmsmartseq2@meta.data$Celltype %in% "Monocytes")
ct2=bmsmartseq2@meta.data$Celltype[pick2]
names(ct2)=rownames(bmsmartseq2@meta.data)[pick2]
ct2=ct2[!(names(ct2) %in% cellid1)]
sort(table(ct2),decreasing=TRUE)
cellid2=names(ct2)

fdrThres=0.05
#cellranger
clMarker_Smartseq2 <- FindMarkers(bmsmartseq2$smartseq2, cells.1 = cellid1, cells.2 = cellid2, min.pct = 0.25,test.use="LR")
reformatIsoforms=gsub("-","_",rownames(clMarker_Smartseq2))
genes=genes.tx.map[reformatIsoforms]
clMarker_Smartseq2$GENEID=genes

sum(clMarker_Smartseq2$p_val_adj < fdrThres)
clMarker_Smartseq2=clMarker_Smartseq2[clMarker_Smartseq2$p_val_adj < fdrThres,]
dim(clMarker_Smartseq2)
head(clMarker_Smartseq2)

tx_TYROBP=c("NM-003332","NM-198125","NM-001173514","NR-033390","NM-001173515")
plotData=bmsmartseq2$smartseq2@data[tx_TYROBP,]
plotData=as.matrix(plotData)
plotCt=bmsmartseq2@meta.data$Celltype
plotIdents=as.character(bmsmartseq2@meta.data$ClusterSmartSeq2)
#rename idents
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
#View(cbind(cluster_ct_id,cluster_ct_label))
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

x=t(plotData)
pick=plotIdents=="TY32.25 Mono"
x1=x[pick,]
chosen1 <- c("NM-198125","NM-003332")
chosen2 <- c("NM-001173515","NM-001173514","NR-033390")
chosen1title <- paste0("TY32.25 MONO\n(SMART-SEQ2,N=",sum(pick),")")
z1=rowSums(x1[,chosen1,drop=FALSE])
z2=rowSums(x1[,chosen2,drop=FALSE])
fig2c <- data.frame(ISOFORM = paste(chosen1,collapse = "\n"),ISOFORM_EXPRESSION = z1, GROUP = chosen1title)
fig2c <- rbind(fig2c, data.frame(ISOFORM = paste(chosen2,collapse = "\n"),ISOFORM_EXPRESSION = z2, GROUP = chosen1title))

pick=plotIdents=="Monocytes"
pick=setdiff(rownames(x)[plotCt=="Monocytes"],rownames(x1))
x2=x[rownames(x) %in% pick,]
z1=rowSums(x2[,chosen1,drop=FALSE])
z2=rowSums(x2[,chosen2,drop=FALSE])
chosen2title <- paste0("OTHER MONO\n(SMART-SEQ2,N=",nrow(x2),")")
fig2c <- rbind(fig2c, data.frame(ISOFORM = paste(chosen1,collapse = "\n"),ISOFORM_EXPRESSION = z1, GROUP = chosen2title))
fig2c <- rbind(fig2c, data.frame(ISOFORM = paste(chosen2,collapse = "\n"),ISOFORM_EXPRESSION = z2, GROUP = chosen2title))

fig2c$GROUP <- factor(fig2c$GROUP, levels = sort(unique(fig2c$GROUP), decreasing = T))
fig2c$ISOFORM <- factor(fig2c$ISOFORM, levels = sort(unique(fig2c$ISOFORM), decreasing = T))

apecols2 <- c("#6153CC","#85cb33")
p2c <- ggplot(fig2c, aes(x=ISOFORM, y=ISOFORM_EXPRESSION, fill=ISOFORM)) +
  stat_boxplot(color = "black")+ geom_hline(yintercept = 1, color = "#A9431F",linetype="dotted")+
  theme_classic() + 
  scale_fill_manual(values = apecols2) + 
  facet_wrap(~GROUP)+
  ylab("ISOFORM EXPRESSION")+
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1,vjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "none",
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
        legend.title = element_text(size =20, face = "bold"),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 20, face = "bold", ),
        strip.background = element_blank(),
        plot.title = element_text(size =20, face = "bold", hjust = 0.5))

png("FIGURES_AND_TABLES/MAIN_FIGURE_2C.png", width=2000, height=2000, units = "px", res = 300)  
print(p2c)
dev.off()

################################################################################

################################################################################
# FIGURE 2d
################################################################################
ref <- read.csv("FIGURE_DATA/transcripts_to_genes_annotation.txt", header=F, sep="\t")
assays <- c('SCiso','SCgene','CR','Alevin')
umaps <- c('rna.umap','SCgene.umap','cr.umap','Alevin.umap')
output_names <- c('SCASA_ISOFORM','SCASA_GENE','CELLRANGER_GENE','ALEVIN_GENE')

for(i in 1:length(assays)){
  plotx <- bm
  DefaultAssay(plotx) <- assays[i]
  plotx@reductions[[umaps[i]]]@cell.embeddings[,1] <- -plotx@reductions[[umaps[i]]]@cell.embeddings[,1]
  TYROBP <- ref[which(ref$V2 == "TYROBP"),]
  pos <- NULL
  for(k in 1:nrow(TYROBP)){
    pos <- c(pos,grep(gsub("_","-",TYROBP$V1[k]),row.names(plotx)))
  }
  pos <- unique(pos)
  
  if(assays[i] == "SCiso"){
    for(j in 1:length(pos)){
      p <- FeaturePlot(plotx, features = row.names(plotx)[pos[j]], reduction = umaps[i], pt.size = 1) + 
        theme(axis.line = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
              plot.title = element_text(size =30, face = "bold", hjust = 0.5))+NoLegend()+ labs(x = NULL, y=NULL)
      
      png(paste("FIGURES_AND_TABLES/MAIN_FIGURE_2D",j,"_",output_names[i],"_TYROBP_UMAP.png",sep = ""), width=1700, height=2000, units = "px", res = 300)  
      print(p)
      dev.off()
    }
    
  }
}

################################################################################