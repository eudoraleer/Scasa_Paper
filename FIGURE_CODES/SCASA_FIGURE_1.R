################################################################################
## SCASA: MAIN FIGURES
################################################################################
library("ggplot2")
library("patchwork")
setwd("Scasa_Paper_Codes_Data/")
################################################################################
# FIGURE 1
################################################################################

################################################################################
# FIGURE 1b
################################################################################

# Load all data
source("FIGURE_CODES/FUNCTIONS/LOAD_ALL_DATA.R")

# Scasa isoform
scasa_iso = c("scasaAle_iso")
# True counts isoform
scasa_true =  c("scasaAle_tiso")
par(mfrow=c(2,3))
scasa = get(scasa_iso[1])
tcount= get(scasa_true[1]) 
equ = equrow(tcount,scasa) 
pick = equ[[1]]>0 | equ[[2]]>0
set.seed(2021)
p = 0.03
sam = sample(1:length(pick), p*length(pick)) ## subsample for fast plotting
pick[-sam] = FALSE;   print(sum(pick))
data <- data.frame(TRUE_COUNT = c(equ[[1]][pick])+1,
                   ESTIMATED_COUNT = c(equ[[2]][pick])+1,
                   DATA = "SCASA (ISOFORM)")

# Output from other quantification software: Kallisto, Salmon
all_isoform = c("tcount_isoform", "kallisto_isoform", "salmon_isoform")
matlist = mget(all_isoform)
equ = equrow_col(matlist)
current_names <- c("SCASA (ISOFORM)","KALLISTO (ISOFORM)","SALMON (ISOFORM)")

for (i in 2:3){
  pick = equ[[1]]>0 | equ[[i]]>0
  set.seed(2021)
  p = 0.02
  sam = sample(1:length(pick), p*length(pick))
  pick[-sam] = FALSE 
  print(sum(pick))
  data <- rbind(data, data.frame(TRUE_COUNT = c(equ[[1]][pick])+1,
                                 ESTIMATED_COUNT = c(equ[[i]][pick])+1,
                                 DATA = current_names[i]))
}

data$DATA <- factor(data$DATA, levels = c("SCASA (ISOFORM)","KALLISTO (ISOFORM)","SALMON (ISOFORM)"))
# colors <- c("#EA8C55","#242423","#69995D")
colors <- c("#EA8C55","#33658A","#69995D")

p1b <- ggplot(data, aes(TRUE_COUNT, ESTIMATED_COUNT, color = DATA)) +
  geom_point()+scale_x_log10()+scale_y_log10()+ theme_classic()+
  facet_wrap(~DATA, scales = "free")+ scale_color_manual(values = colors)+
  geom_abline(intercept = 0, color = c("#968147","#0C0A06","#334B2E"))+
  xlab("TRUE COUNT + 1")+ylab("ESTIMATED COUNT + 1")+
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 30, margin=margin(0,10,0,0)),
        strip.text.x = element_text(size = 30, face = "bold"),
        strip.background = element_blank(),
        legend.position = "none")

png("FIGURES_AND_TABLES/MAIN_FIGURE_1B.png", width=6000, height=2000, units = "px", res = 300)  
print(p1b)
dev.off()

################################################################################
# FIGURE 1c
################################################################################

# Compare APE at isoform-level
sc.res = c("scasaAle_iso", "kallisto_isoform", "salmon_isoform")
para.APE = vector('list', length(sc.res))

# Calculates APE for Scasa
for (i in 1:1){
  scasa = get(sc.res[i])
  tcount = scasa_truemat(scasa, tcount=tcount_isoform)  
  equ = equrow_col(list(true=tcount, est=scasa))
  APE = abs(equ$est-equ$true)/(equ$est+equ$true+1)
  
  para = rownames(scasa)
  npara = sapply(strsplit(para,' '),length)
  para.rep = rep(para, npara)
  txpara = unlist(strsplit(para,' '))   ##list of tx-level 
  para.APE[[i]] = APE[para.rep,]   ## repeat
  rownames(para.APE[[i]]) = txpara
}

# APE Ratio: Calculates APE for Other Softwares and Compare with Scasa APE
ape_p1 <- NULL
ape_p2 <- NULL
ape_p3 <- NULL

################################################################################
# 1. APE Ratio at Isoform-Level
truelim= 0
txpara = rownames(para.APE[[1]])
# txpara = union(rownames(para.APE[[1]]),rownames(para.APE[[2]]))
for (i in 2:length(sc.res)){
  truepick = rownames(tcount_isoform) %in% txpara
  sc.resi = get(sc.res[i])
  sc.respick = rownames(sc.resi) %in% txpara
  ext = equrow(sc.resi[sc.respick,], tcount_isoform[truepick,])
  est = ext$A; 
  true = ext$B[,colnames(est)]
  para.APE[[i]] = abs(est-true)/(est+true+1)
}

# Equalize rows and cols
para.equAPE = equrow_col(para.APE)
para.relAPE = para.num = vector('list', length(sc.res)-1)
tx = rownames(para.equAPE[[1]])
true = matrix(0, nrow= nrow(para.equAPE[[1]]), ncol= ncol(para.equAPE[[1]]) )
rownames(true) = tx
colnames(true) = colnames(para.equAPE[[1]])
pick = tx %in% rownames(tcount_isoform)
txpick = rownames(true)[pick]
true[txpick,] = tcount_isoform[txpick, colnames(true)]
truelim=0
for (i in 1:(length(para.equAPE)-1)){
  a = para.equAPE[[1]]; b= para.equAPE[[i+1]]
  #pick = ((a>0) | (b>0)) ## compare only at least one has pos APE
  pick = abs(a-b)>0.00
  para.num[[i]] = c(sum(a>b), sum(a<b))
  para.relAPE[[i]] = c(b[pick]+0.1)/c(a[pick]+0.1)  
}
labels= c('ALL ISOFORMS KALLISTO','ALL ISOFORMS SALMON')
names(para.relAPE) = labels
## 5% for plotting
set.seed(2021)
para.relAPE.sam = sapply(para.relAPE, function(x) sample(x, 0.05*length(x)))  

for(i in 1:length(para.relAPE.sam)){
  ape_p1 <- rbind(ape_p1,
                  data.frame(APE_RATIO = para.relAPE.sam[[i]],
                             SOFTWARE = toupper(gsub("ALL ISOFORMS |NON-PARALOGS ","",names(para.relAPE.sam)[i])),
                             DATA_TYPE = "ISOFORM-LEVEL",
                             DATA_SUBTYPE = toupper(gsub(" KALLISTO| SALMON","",names(para.relAPE.sam)[i]))))
}

fn = function(x) {
  tx = rownames(x) %in% nonpara_isoform;
  pick = rowSums(x[tx,]) >0  ## pick only rows with >0 values
  return(x[tx,][pick,])
}

# Non-paralogs: isoforms from multiple-isoform genes that are NOT merged together from CCRP
nonpara.list = lapply(para.APE,fn)
nonpara.equlist = equrow_col(nonpara.list)

# pairwise comparisons
true = true[rownames(nonpara.equlist[[1]]),]
truelim= 0         ## subset of true values
nonpara.relAPE = nonpara.num = vector('list', length(nonpara.equlist)-1)
for (i in 1:(length(nonpara.equlist)-1)){
  a = nonpara.equlist[[1]]; b= nonpara.equlist[[i+1]]
  pick = abs(a-b) >0.00
  nonpara.num[[i]] = c(sum(a>b), sum(a<b))
  nonpara.relAPE[[i]] = c(b[pick]+0.1)/c(a[pick]+0.1)  
}
labels= c('NON-PARALOGS KALLISTO','NON-PARALOGS SALMON')
names(nonpara.relAPE) = labels
for(i in 1:length(nonpara.relAPE)){
  ape_p1 <- rbind(ape_p1,
                  data.frame(APE_RATIO = nonpara.relAPE[[i]],
                             SOFTWARE = toupper(gsub("ALL ISOFORMS |NON-PARALOGS ","",names(nonpara.relAPE)[i])),
                             DATA_TYPE = "ISOFORM-LEVEL",
                             DATA_SUBTYPE = toupper(gsub(" KALLISTO| SALMON","",names(nonpara.relAPE)[i]))))
}

################################################################################
# 2. APE Ratio at Gene-Level
gene.APE = vector('list', 6)
equ = equrow_col(list(true=scasaAle_tgene, est=scasaAle_gene))
APE = abs(equ$est-equ$true)/(equ$est+equ$true+1)
gene.APE[[1]] = APE_extgene(APE)

RES = c("cr_mat","klbus_mat","alevin_mat", "kallisto_gene", "salmon_gene")
i  = 2
for (res in RES){
  equ = equrow_col(list(true=tcount_gene, est=get(res)))
  APE = abs(equ$est-equ$true)/(equ$est+ equ$true+1)
  gene.APE[[i]] = APE
  i= i+1
}
sapply(gene.APE, dim)  ## compare sizes
names(gene.APE) = c('Scasa','CellRanger','KallistoBus',
                    'Alevin', 'Kallisto','Salmon')
# equalize all rows and cols
equgene.APE = equrow_col(gene.APE)  #
sapply(equgene.APE, dim)  ## compare sizes

# Pairwise comparisons: gene level, ALL genes
gene.relAPE = gene.num = vector('list', length(equgene.APE)-1)
for (i in 1:(length(equgene.APE)-1)){
  a = equgene.APE[[1]]; b= equgene.APE[[i+1]]
  pick = (abs(a-b)>0.00) ## compare only the APEs are unequal
  # pick = (a>0)|(b>0)  # compare at least one APE>0
  gene.num[[i]] = c(sum(a>b), sum(a<b))
  gene.relAPE[[i]] = c(b[pick]+0.1)/c(a[pick]+0.1)  
}
labels = c('CellRanger','Kallisto-Bus', 'Alevin', 'Kallisto','Salmon')
names(gene.relAPE) = labels
gene.relAPE.sam= gene.relAPE
gene.relAPE.sam$Kallisto = sample(gene.relAPE$Kallisto, 0.1*length(gene.relAPE$Kallisto))

for(i in 1:length(gene.relAPE.sam)){
  ape_p2 <- rbind(ape_p2,
                  data.frame(APE_RATIO = gene.relAPE.sam[[i]],
                             SOFTWARE = toupper(names(gene.relAPE.sam)[i]),
                             DATA_TYPE = "GENE-LEVEL"))
}

################################################################################
# 3. APE Ratio at Singleton-Level
txmap = genes.tx.map.all.final
genes.tx.list = tapply(names(txmap), txmap, c)
ntx = sapply(genes.tx.list, length)
singles = names(genes.tx.list[ntx==1])  ## check = singleton_isoform
##
sing.relAPE = sing.num = vector('list', length(equgene.APE)-1)
truelim = 0
for (i in 1:(length(equgene.APE)-1)){
  a = equgene.APE[[1]]; b= equgene.APE[[i+1]]
  sgene = rownames(a) %in% txmap[singleton_isoform] 
  a = a[sgene,]
  b = b[sgene,]
  #true = equlist[[1]][sgene,]  ## true values
  pick =  (abs(a-b)>0.00)
  sing.num[[i]] = c(sum(a>b), sum(a<b))
  sing.relAPE[[i]] = c(b[pick]+0.1)/c(a[pick]+0.1)  
}
labels = c('CellRanger','Kallisto-Bus', 'Alevin', 'Kallisto','Salmon')
names(sing.relAPE) = labels

for(i in 1:length(sing.relAPE)){
  ape_p3 <- rbind(ape_p3,
                  data.frame(APE_RATIO = sing.relAPE[[i]],
                             SOFTWARE = toupper(names(sing.relAPE)[i]),
                             DATA_TYPE = "SINGLETON-LEVEL"))
}

apecols1 <- c("#33658A","#69995D")
apecols2 <- c("#f26419","#86bbd8","#33658A","#f6ae2d","#69995D")
p1c1 <- ggplot(ape_p1, aes(x=SOFTWARE, y=APE_RATIO, fill=SOFTWARE)) +
  stat_boxplot(color = "black")+ geom_hline(yintercept = 1, color = "#A9431F",linetype="dotted")+
  theme_classic() + ylab("APE RATIO") + xlab("")+
  ggtitle("ISOFORM-LEVEL") +
  scale_fill_manual(values = apecols1) + facet_wrap(~DATA_SUBTYPE)+
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1,vjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "none",
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
        legend.title = element_text(size =20, face = "bold"),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 20),
        # strip.background = element_blank(),
        plot.title = element_text(size =30, face = "bold", hjust = 0.5))

png("FIGURES_AND_TABLES/MAIN_FIGURE_1C_APE_RATIO_1.png", width=2000, height=2000, units = "px", res = 300)  
print(p1c1)
dev.off()

p1c2 <- ggplot(ape_p2, aes(x=SOFTWARE, y=APE_RATIO, fill=SOFTWARE)) +
  stat_boxplot(color = "black")+ geom_hline(yintercept = 1, color = "#A9431F",linetype="dotted")+
  theme_classic() + ylab("APE RATIO") + xlab("")+
  ggtitle("GENE-LEVEL\n") +
  scale_fill_manual(values = apecols2) + 
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1,vjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "none",
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
        legend.title = element_text(size =20, face = "bold"),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 25, face = "bold", ),
        strip.background = element_blank(),
        plot.title = element_text(size =30, face = "bold", hjust = 0.5))

png("FIGURES_AND_TABLES/MAIN_FIGURE_1C_APE_RATIO_2.png", width=2000, height=2000, units = "px", res = 300)  
print(p1c2)
dev.off()

p1c3 <- ggplot(ape_p3, aes(x=SOFTWARE, y=APE_RATIO, fill=SOFTWARE)) +
  stat_boxplot(color = "black")+ geom_hline(yintercept = 1, color = "#A9431F",linetype="dotted")+
  theme_classic() + 
  ggtitle("SINGLETON-LEVEL\n") +
  scale_fill_manual(values = apecols2) + ylab("APE RATIO") + xlab("")+
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1,vjust = 1),
        axis.text.y = element_text(size = 20), legend.position = "none",
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
        legend.title = element_text(size =20, face = "bold"),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 25, face = "bold", ),
        strip.background = element_blank(),
        plot.title = element_text(size =30, face = "bold", hjust = 0.5))

png("FIGURES_AND_TABLES/MAIN_FIGURE_1C_APE_RATIO_3.png", width=2000, height=2000, units = "px", res = 300)  
print(p1c3)
dev.off()

save.image("FIGURE_DATA/FIGURE_1.RData")
################################################################################