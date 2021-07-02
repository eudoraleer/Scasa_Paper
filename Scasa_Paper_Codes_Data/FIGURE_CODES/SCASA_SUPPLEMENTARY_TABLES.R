################################################################################
## SCASA: SUPPLEMENTARY TABLES
################################################################################
setwd("Scasa_Paper_Codes_Data/")

################################################################################
# SUPPLEMENTARY TABLE 1
################################################################################
# Run after "FIGURE_CODES/SCASA_FIGURE_1.R" or retrieve the available FIGURE_1.RData from "FIGURE_DATA/" directory
load("FIGURE_DATA/FIGURE_1.RData")

para.prop= sapply(para.num, function(x) x[2]/sum(x))
para.tot = sapply(para.num, function(x) sum(x))
nonpara.prop= sapply(nonpara.num, function(x) x[2]/sum(x))
nonpara.tot = sapply(nonpara.num, function(x) sum(x))
gene.prop= sapply(gene.num, function(x) x[2]/sum(x))
gene.tot = sapply(gene.num, function(x) sum(x))
sing.prop= sapply(sing.num, function(x) x[2]/sum(x))
sing.tot = sapply(sing.num, function(x) sum(x))

PROP = rbind(All_isoforms= c(0,0,0,para.prop),
             Nonparalog = c(0,0,0,nonpara.prop),
             All_genes = gene.prop, 
             Singletons= sing.prop)
TOT = rbind( c(0,0,0,para.tot),
             c(0,0,0,nonpara.tot),
             gene.tot, 
             sing.tot)
ALL = matrix(0, nrow=8, ncol=ncol(PROP))
colnames(ALL) = labels
rownames(ALL) = rep('N', 8)
ALL[c(1,3,5,7),] = round(PROP,2)
rownames(ALL)[c(1,3,5,7)] = rownames(PROP)
ALL[c(2,4,6,8),] = TOT
print(ALL)

write.table(ALL, "FIGURES_AND_TABLES/SUPPLEMENTARY_TABLE_1.txt",quote = F, row.names=T, sep = "\t")

################################################################################
# SUPPLEMENTARY TABLE 2
################################################################################
# SUPPLEMENTARY TABLE 2 is "clMarker_SCiso.txt" file from "FIGURE_CODES/SCASA_FIGURE_2.R"
# File directory: "FIGURES_AND_TABLES/clMarker_SCiso.txt"
################################################################################

################################################################################
# SUPPLEMENTARY TABLE 4
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

## compare CRP and CCRP
X0= scCRP[[grep(txsam[1],names(CRP))]]
round(svd(X0)$d,3)
X = scCCRP[[grep(txsam[1],names(CRP))]]
round(1000*X0)
round(X,3)
round(X0,3)

write.table(round(1000*X0), "FIGURES_AND_TABLES/SUPPLEMENTARY_TABLE_4A.txt",quote = F, row.names=T, sep = "\t")
write.table(round(X,3), "FIGURES_AND_TABLES/SUPPLEMENTARY_TABLE_4B.txt",quote = F, row.names=T, sep = "\t")

################################################################################
# SUPPLEMENTARY TABLE 5
################################################################################
round(X0,3)
write.table(round(X0,3), "FIGURES_AND_TABLES/SUPPLEMENTARY_TABLE_5A.txt",quote = F, row.names=T, sep = "\t")

load("FIGURE_DATA/Simulated_Data/EstimatedResults/X_matrix_bulk_UCSC_hg38.RData")

Xb=CRP[[grep('NM_001270491', names(CRP))]]
CRP[[grep('NR_003932',names(CRP))]]
Xa = cbind(Xb[,1:2],0,Xb[,3:4])
colnames(Xa)[3]='NR_003932'
Xa = rbind(c(0,0,1,0,0), Xa)
rownames(Xa) = c('00100',"00001","00010", "01000", "01001" ,"10000",
                 "10001", "11001", "11011")
a = Xa; a[c(1:3),] = Xa[c(3,1,2),]
rownames(a)[1:3]= rownames(Xa)[c(3,1,2)]
round(a,3)
round(svd(Xa)$d,3)
write.table(round(a,3), "FIGURES_AND_TABLES/SUPPLEMENTARY_TABLE_5B.txt",quote = F, row.names=T, sep = "\t")

################################################################################
