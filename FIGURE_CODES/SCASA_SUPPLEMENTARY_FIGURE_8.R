################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")
################################################################################
# SUPPLEMENTARY FIGURE 8
################################################################################
source("FIGURE_CODES/FUNCTIONS/LOAD_ALL_DATA.R")
load("FIGURE_DATA/Simulated_Data/EstimatedResults/Xmatrix_hg38.RData")

scCRP = CRP; scCCRP = CCRP
txmap <- genes.tx.map.all.final
para = rownames(scasaAle_iso)
totx = rowSums(scasaAle_iso)
paralist = strsplit(para,' ')
ntx = sapply(paralist, length)
loc=c(1880, 4454,4937)
txsam = unlist(paralist[loc[2]])
txmap[txsam]

################################################################################
# SUPPLEMENTARY FIGURE 8a
################################################################################
png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_8A.png", width=4000, height=4000, units = "px", res = 800)  
truetx = rownames(tcount_isoform)[rownames(tcount_isoform) %in% txsam]
salmtx = rownames(salmon_isoform)[rownames(salmon_isoform) %in% txsam]
equ = equrow(tcount_isoform[truetx,],salmon_isoform[salmtx,])
cols = c('#fd6f30','blue','green')
plot(c(equ[[1]])+1, c(equ[[2]])+1, log='xy', type='n', 
     cex=0.6, xlab='TRUE COUNT + 1', ylab='ESTIMATED COUNT + 1', xlim=c(1,100), ylim=c(1,100))
for (i in 1:3) {
  a = equ[[1]][i,]; b = equ[[2]][i,]
  pick = a>0 | b>0
  if (i==1) text(a[pick]+1,b[pick]+1.15, i, cex=.5, col=cols[i])
  if (i==2) text(a[pick]+1.1,b[pick]+1.1, i, cex=.5, col=cols[i])
  if (i==3) text(a[pick]+1,b[pick]+1, 5, cex=.5, col=cols[i])
}
abline(0,1)
title('SALMON (scRNA-SEQ)')
dev.off()

################################################################################
# SUPPLEMENTARY FIGURE 8b
################################################################################
load("FIGURE_DATA/Simulated_Data/GenerateSimdata/RPL13A.RData")
true = RPL13A.true[c(1,2,4),]
salmon = RPL13A.salmon[c(1,2,4),]
cols = c('#fd6f30','blue','green')
lab = c(1,2,5)
png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_8B.png", width=4000, height=4000, units = "px", res = 800)  
plot(c(true)+1,c(salmon)+1, log='xy', type='n', 
     xlim=c(1,500), ylim=c(1,500),
     cex=0.6, xlab='TRUE COUNT + 1', ylab='ESTIMATED COUNT + 1')
for (i in 1:3) text(true[i,]+1,salmon[i,]+1, lab[i], cex=.7, col=cols[i])
abline(0,1)
title('SALMON (BULK RNA-SEQ)')
dev.off()
################################################################################
