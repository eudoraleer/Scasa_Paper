########################################################################
### generate the simulated dataset for the study
# Step 1 - run in R
# - generate reads from polyester
# - produce trancript-gene map file:txp2gene_hg38.tsv
# NOTE: the object cdnaf of original cdnaf.rda file (located in polyester/data/ of the polyester package) needs to be replaced by the 3 prime bias model described in the manuscript.

library("GenomicFeatures")
library("Biostrings")
library("polyester")

### Necessary files
#revised function for Polyester
source("generate_fragments.R")
#annotation
cdnaFasta <- "refMrna.fa"
gtfFile="hg38.refGene.gtf"
# simulation is based on the output of gene expression data of cellranger for one healthy PBMC human donor using 10X Genomics V3.
per_cell_umi <- read.table(gzfile("melt_cellranger_cell_gene_count.txt.gz"))
per_gene_umi <- read.table(gzfile("melt_cellranger_cell_each_gene_count.txt.gz"), header = T)
#isoform ration of each gene is collected from a smartseq dataset. 
reference_ratio <- read.table(gzfile("SMARTSEQ_Samples_Individual_TPM_Counts_Per_Gene.txt.gz"), header = T)

#create new folder
output_dir=paste(getwd(),"/",sep="")
setwd(output_dir)
fasta_output <- paste(output_dir, "CELL/", sep = "")
dir.create(fasta_output)
TrueCount_output <- paste(output_dir, "TrueCount/", sep = "")
dir.create(TrueCount_output)

#load annotation
anntxdb <- makeTxDbFromGFF(file=gtfFile,
                 format="gtf",
                 dataSource=paste("hg38",sep=""),
                 organism="Homo sapiens")
genes.all = genes(anntxdb, single.strand.genes.only = FALSE )
genes.tx.all = suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("GENEID", "TXNAME"), keytype = "GENEID")))
genes.tx.all.ID=paste(genes.tx.all$GENEID,":",genes.tx.all$TXNAME,sep="")

### start processing and generating simulated reads
# Remove the inconsistences between version 2014 and 2015
reference_ratio.GeneTX=paste(reference_ratio$GENE_SYMBOL,":",reference_ratio$TRANSCRIPT_ID,sep="")
pick=reference_ratio.GeneTX %in% genes.tx.all.ID
reference_ratio=reference_ratio[pick,]

# read transcript fasta
tx.all.fasta <- readDNAStringSet(cdnaFasta)
tx.all.ID=sapply(names(tx.all.fasta),function(x) strsplit(x," ")[[1]][1])
names(tx.all.ID)=NULL

#export txp2gene of the whole transcriptome, do only once
pick=tx.all.ID %in% genes.tx.all$TXNAME
x=tx.all.ID[which(!pick)]#some transcripts in fasta file are not found in gtf file
txp2gene_hg38=genes.tx.all[,c("TXNAME","GENEID")]
y=data.frame(TXNAME=x,GENEID=x)
txp2gene_hg38=rbind(txp2gene_hg38,y)
save(txp2gene_hg38,file="txp2gene_hg38.RData")
write.table(txp2gene_hg38,file="txp2gene_hg38.tsv",quote = FALSE,sep="\t", row.names = FALSE,col.names = FALSE)

#Generate true counts and simulate reads
set.seed(2020)
for(i in 1:nrow(per_cell_umi)){

  cell_name <- per_cell_umi[i,"V1"]
  cell_umi_count <- per_cell_umi[i,"V2"]
  
  current <- per_gene_umi[which(per_gene_umi$CELL_ID == cell_name),c("GENE_SYMBOL","GENE_ID","COUNT")]
  current_reference <- data.frame(TRANSCRIPT_ID = reference_ratio$TRANSCRIPT_ID,
                                  GENE_SYMBOL = reference_ratio$GENE_SYMBOL,
                                  READ_RATIO = reference_ratio[,sample(grep("Cell_", colnames(reference_ratio), ignore.case = T), size = 1, replace = T)])
  current_reference$PER_GENE_UMI_COUNT <- current[match(toupper(current_reference$GENE_SYMBOL), toupper(current$GENE_SYMBOL)),"COUNT"]
  current_reference[which(is.na(current_reference$PER_GENE_UMI_COUNT)),"PER_GENE_UMI_COUNT"] <- 0
  current_reference$FINAL_COUNT <- ceiling(current_reference$READ_RATIO*current_reference$PER_GENE_UMI_COUNT)

# Keep only expressed transcripts for simulation
  current_reference=current_reference[current_reference$FINAL_COUNT > 0,]
  #NOTE: just export cdna fasta of expressed transcripts, if non-expressed transcripts are also included, some noises are popped up without controlled

  #make sure the transcripts having fasta sequence
  pick=current_reference$TRANSCRIPT_ID %in% tx.all.ID
  current_reference=current_reference[pick,]
  write.table(current_reference, paste(TrueCount_output,cell_name,"_REFERENCE_RATIO.txt", sep = ""), quote = F, row.names = F)

  tx.export.fasta=tx.all.fasta
  #export cdna of only trascripts in current_reference
  matchID=match(as.character(current_reference$TRANSCRIPT_ID),tx.all.ID)
  tx.export.fasta=tx.export.fasta[matchID]
  output_dir=paste(getwd(),"/",sep="")

  fasta_chosen_transcript <- paste("Cell_level_tx_temp.fa", sep="") #temporary file
  writeXStringSet(tx.export.fasta, paste(output_dir,fasta_chosen_transcript,sep = ""))

  #set read counts
  readmat_input <- as.matrix(current_reference$FINAL_COUNT)  
  simulate_experiment_countmat(fasta_chosen_transcript,
                               readmat=readmat_input, outdir=output_dir, # allow.nonnarrowing = T,
                               readlen = 91, paired = T, seed = 1000,
                               fraglen = 400,
                               fragsd = 100,
                               allow.nonnarrowing=F,
                               distr = 'normal',
                               strand_specific = TRUE,
                               error_rate = 0 , bias="cdnaf")
  #### NOTE: move to the directory, BUT SWITCH THE NAME; read1 <--> read2
  system(paste("mv ",output_dir,"sample_01_1.fasta ", fasta_output, cell_name,"_2.fasta", sep = ""))
  system(paste("mv ",output_dir,"sample_01_2.fasta ", fasta_output, cell_name, "_1.fasta", sep = ""))
  
  # UMI,12bp
  umi_size <- 12
  print("Generating UMI barcode now..")
  
  current_umi <- do.call(paste0, replicate(umi_size, sample(c("A","T","G","C"), sum(current_reference$FINAL_COUNT), TRUE), FALSE))
  
  print("Finished generating UMI barcode..")
  
  read1_barcodes <- data.frame(READ1 = paste(cell_name, current_umi, sep = ""))
  write.table(read1_barcodes, paste(fasta_output,cell_name, ".barcodes.txt",sep = ""), quote = F, row.names = F, col.names = F)
  
} #end of i


