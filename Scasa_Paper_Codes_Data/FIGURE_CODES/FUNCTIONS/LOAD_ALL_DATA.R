################################################################################
## MAIN FIGURE: FUNCTIONS AND DATA
## About: Loading of functions and data for analysis and figure generation
################################################################################

# Load functions
source("FIGURE_CODES/FUNCTIONS/SCASA_FUNCTIONS.R")

################################################################################
# Load all data for analysis and plotting
################################################################################
# 1. Simulated data: true counts
#### Variables: tcount_gene, tcount_isoform, tcount_raw
load("FIGURE_DATA/Simulated_Data/EstimatedResults/TrueCounts_sim.RData")

################################################################################
# 2. Scasa result for simulated data: estimated counts
#### Variables: genes.tx.map.all.final, genes.tx.map.raw, nonpara_isoform, para_isoform, singleton_isoform  
load("FIGURE_DATA/Simulated_Data/EstimatedResults/isoform_groups_UCSC_hg38.RData")

################################################################################
# 3. Scasa final isoform result for simulated data: isoform_count (alignment ran with Alevin)
load("FIGURE_DATA/Simulated_Data/EstimatedResults/Scasa_isoform_sim.RData") ## Variables: isoform_count
# print(isoform_count[1:5,1:5])
rsum = rowSums(isoform_count)
scasaAle_iso = isoform_count[rsum>0,]
scasaAle_tiso = scasa_truemat(scasaAle_iso, tcount=tcount_isoform)
scasaAle_gene = scasa_genemat(scasaAle_iso, txmap=genes.tx.map.all.final)
scasaAle_tgene = scasa_truemat(scasaAle_gene, tcount=tcount_gene)

################################################################################
# Results from other quantification tools  
load("FIGURE_DATA/Simulated_Data/EstimatedResults/Alevin_sim.RData")  ## Variables: alevin_mat
load("FIGURE_DATA/Simulated_Data/EstimatedResults/Cellranger_sim.RData") ## Variables: cr_mat
load("FIGURE_DATA/Simulated_Data/EstimatedResults/KallistoBus_sim.RData") ## Variables: klbus_mat
load("FIGURE_DATA/Simulated_Data/EstimatedResults/kallisto_sim.RData")  ## Variables: kallisto_gene, kallisto_isoform
load("FIGURE_DATA/Simulated_Data/EstimatedResults/Salmon_sim.RData")  ## Variables: salmon_gene, salmon_isoform
################################################################################
