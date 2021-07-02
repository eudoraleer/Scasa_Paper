################################################################################
## SCASA: SUPPLEMENTARY FIGURES
################################################################################
setwd("Scasa_Paper_Codes_Data/")
################################################################################
# SUPPLEMENTARY FIGURE 2
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
# SUPPLEMENTARY FIGURE 2b
################################################################################
png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_2B.png", width=3000, height=3000, units = "px", res = 800)  
par(mfrow=c(1,1))
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
title("SALMON (ISOFORM)")
dev.off()

################################################################################
# SUPPLEMENTARY FIGURE 2c
################################################################################
png("FIGURES_AND_TABLES/SUPPLEMENTARY_FIGURE_2C.png", width=3000, height=3000, units = "px", res = 800)  
scasatx = rownames(scasaAle_iso)[loc[2]]
plot(scasaAle_tiso[scasatx,]+1, scasaAle_iso[scasatx,]+1, log='xy', 
     xlab='TRUE COUNT + 1', ylab='ESTIMATED COUNT + 1', xlim=c(1,100), ylim=c(1,100))
abline(0,1)
title('SCASA (PARALOG)')
dev.off()

################################################################################
