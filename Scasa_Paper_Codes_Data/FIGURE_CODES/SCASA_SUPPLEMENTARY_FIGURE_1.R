################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")

library("ggplot2")
library("patchwork")
################################################################################
# SUPPLEMENTARY FIGURE 1
################################################################################
load("FIGURE_DATA/Simulated_Data/EstimatedResults/X_matrix_bulk_UCSC_hg38.RData")
hg38 = sapply(CCRP, colnames)
ncrp38 = sapply(hg38, length)
para38 = unlist(hg38); length(para38)
tx38 = unlist(strsplit(para38,' ')); length(tx38); length(unique(tx38))
npara38 = sapply(strsplit(para38,' '), length)
totrow = sum(sapply(CCRP,nrow))

load("FIGURE_DATA/Simulated_Data/EstimatedResults/Xmatrix_hg38.RData")
x10 = sapply(CCRP, colnames)
ncrp10 = sapply(x10, length)
para10 = unlist(x10)
length(para10)
tx10 = unlist(strsplit(para10,' ')); length(tx10); length(unique(tx10))
excl = which(!(tx10 %in% tx38))
tx10 = tx10[-excl]   ## exclude one outlying tx
npara10= sapply(strsplit(para10,' '), length)

load("FIGURE_DATA/Simulated_Data/EstimatedResults/isoform_groups_UCSC_hg38.RData")

ntx = length(tx38)
cat('Total number of tx =', ntx,'\n')
clean=TRUE
if (clean){
  cat('Bulk-seq\n')
  cat('Num of tx in singleton CRPs', sum(tx38 %in% names(para38)),'\n')
  bulk1 = tx38[tx38 %in% names(para38)]
  cat('Num of tx not in paralogs', sum(tx38 %in% (para38)),'\n')
  cat('Num of tx not in paralogs, belonging to singleton genes', 
      sum(tx38 %in% (singleton_isoform)),'\n')
  cat('Num of tx in paralogs', sum(!(tx38 %in% (para38))),'\n')
  cat('scRNAseq\n')
  cat('Num of tx in singleton CRPs', sum(tx10 %in% names(para10)),'\n')
  sc1 = tx10[tx10 %in% names(para10)]
  cat('Num of tx not in paralogs', sum(tx10 %in% (para10)),'\n')
  cat('Num of tx not in paralogs, belonging to singleton genes', 
      sum(tx10 %in% (singleton_isoform)),'\n')
  cat('Num of tx in paralogs', sum(!(tx10 %in% (para10))),'\n')
}

lim=8
p38 = table(npara38)[1:lim]*c(1:lim)/ntx; 
cat('Proportion of tx in paralog forms: bulk=',(1-p38[1]),'\n')
p10 = table(npara10)[1:lim]*c(1:lim)/ntx; 
cat('Proportion of tx in paralog forms: sc=',(1-p10[1]),'\n')

plotx <- data.frame(rbind(cbind(DATA_TYPE = "Bulk RNA-Seq", PARALOG_SIZE = names(p38), PROPORTION = p38),
                          cbind(DATA_TYPE = "scRNA-Seq", PARALOG_SIZE = names(p10), PROPORTION = p10)))
plotx$PARALOG_SIZE <- as.numeric(as.character(plotx$PARALOG_SIZE))
plotx$PROPORTION <- as.numeric(as.character(plotx$PROPORTION))

colors <- c("#6457A6","#F5CB5C")
ps1 <- ggplot(plotx, aes(PARALOG_SIZE, PROPORTION, fill = DATA_TYPE))+
  geom_bar(stat = "identity", position="dodge", color = "black") + theme_classic()+
  scale_fill_manual(values = colors) +
  xlab("PARALOG SIZE") +
  guides(fill = guide_legend(title = "DATA TYPE")) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
        legend.title = element_text(size =25, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.2, "cm"))

png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_1.png", width=3000, height=2000, units = "px", res = 300)  
print(ps1)
dev.off()

################################################################################