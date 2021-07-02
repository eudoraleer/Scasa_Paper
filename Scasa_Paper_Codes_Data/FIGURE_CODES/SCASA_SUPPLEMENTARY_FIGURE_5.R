################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")
################################################################################
# SUPPLEMENTARY FIGURE 5
################################################################################
fdrThres=0.05

bm <- readRDS("FIGURE_DATA/Bone_Marrow_Seurat_All_Dim_Reductions.RDS")
clMarker_SCiso <- read.table("FIGURES_AND_TABLES/clMarker_SCiso.txt", header = T, sep = "\t")
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

x=t(plotData)

################################################################################
# SUPPLEMENTARY FIGURE 5a
################################################################################
pick1=which(bm@meta.data$ClusterSCiso %in% c(15)) #new subgroup of CD14
pick2=which(bm@meta.data$ClusterSCiso %in% c(0,11,17))#main subgroup of CD14
length(pick1)
length(pick2)

clMarker_SCiso2=clMarker_SCiso[clMarker_SCiso$p_val_adj < fdrThres,]
png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_5A.png", width=2000, height=2500, units = "px", res = 200)  
par(mfrow=c(2,2))
for (i in 1:4){
  feaName=rownames(clMarker_SCiso2)[i]  
  g1=bm$SCiso@data[feaName,pick1]
  g2=bm$SCiso@data[feaName,pick2]
  
  feaName=gsub("-","_",feaName)
  g=clMarker_SCiso2$GENEID[i]
  boxplot(list("Cluster 15"=g1,"CD14 Mono"=g2),main=paste0(feaName,"\n(",g,")"),frame=F,cex.axis=1.5,cex.main=1.5)
}
dev.off()

################################################################################
# SUPPLEMENTARY FIGURE 5b
################################################################################
pick=plotIdents=="TY32.25 Mono"
table(pick)
x1=x[pick,]
z1=rowSums(x1[,c("NM-198125","NM-003332"),drop=FALSE])
z2=rowSums(x1[,c("NM-001173515","NM-001173514 NR-033390"),drop=FALSE])

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_5B.png", width=1500, height=1500, units = "px", res = 200)  
par(mfrow=c(1,2))
par(mar=c(9,4,3,2)+0.1)
boxplot(list("NM_198125\nNM_003332"=z1,"NM_001173515\nNM_001173514\nNR_033390"=z2), main=paste0("TY32.25 Mono\n(n=",sum(pick),")"),frame=F,ylab="Total expression",col="#91dcea",cex.lab=1.2,cex.main=1.5,las=2,cex.axis=1.2, ylim=c(0,8))

pick=plotIdents=="CD14 Mono"
x2=x[pick,]
z1=rowSums(x2[,c("NM-198125","NM-003332"),drop=FALSE])
z2=rowSums(x2[,c("NM-001173515","NM-001173514 NR-033390"),drop=FALSE])

boxplot(list("NM_198125\nNM_003332"=z1,"NM_001173515\nNM_001173514\nNR_033390"=z2),main=paste0("CD14 Mono\n(n=",sum(pick),")"),ylab="Total expression",frame=F,col="#91dcea",cex.lab=1.2,cex.main=1.5,las=2,cex.axis=1.2,ylim=c(0,8))
dev.off()

################################################################################