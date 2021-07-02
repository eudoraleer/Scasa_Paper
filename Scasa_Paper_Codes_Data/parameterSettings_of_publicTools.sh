#This file is to provide the parameter settings of the public tools used in the study

###Parameter settings to run Polyester to generate simulated data
# biasModel="cdnaf" #where the object cdnaf of original cdnaf.rda file (located in polyester/data/ of the polyester package) is replaced by the 3 prime bias model described in the manuscript.
# transcript_fasta_input="cdna.fa" #contain the fasta of transcripts for simulation
# truecount_input is the matrix containing the true counts of the transcripts
# cell_output: folder contains the output reads
# Information of parameter setting for simulate_experiment_countmat() is referred to the manual of the Polyester package.
simulate_experiment_countmat(transcript_fasta_input,
                               readmat=truecount_input, outdir=cell_output,
                               readlen = 91, paired = T, seed = 1000,
                               fraglen = 400,
                               fragsd = 100,
                               allow.nonnarrowing=F,
                               distr = 'normal',
                               strand_specific = TRUE,
                               error_rate = 0 , bias=biasModel)

#### Parameter settings to run public tools for the simulated data
#Input:
read1="Allcells/test_allCells_L001_R1_001.fastq.gz"
read2="Allcells/test_allCells_L001_R2_001.fastq.gz"
onecellRead2="Onecell/test_oneCell_R2.fastq" #example of read2 of one individual cell to run bulk RNA-seq methods
fastaFn="refMrna.fa"
tgMapFile="txp2gene_hg38.tsv"
whiteListFn="3M-february-2018.txt"

#cellranger
#refdata-cellranger-GRCh38-3.0.0 is downloaded and uncompressed from http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
./cellranger-3.1.0/cellranger count --id=cellrangerOut --transcriptome=refdata-cellranger-GRCh38-3.0.0 --fastqs=Allcells --localcores=16


# Alevin
salmon index -t $fastaFn -i salmonIndex

salmon alevin -l ISR -1 $read1 -2 $read2 --chromiumV3 -i salmonIndex -p 16 -o AlevinOut --tgMap $tgMapFile --whitelist $whiteListFn --dumpBfh --dumpFeatures --dumpMtx

#Kallisto Bustools
kallisto index -i kallistoIndex  $fastaFn

kallisto bus -i kallistoIndex -o klbusOut -x 10xv3 -t 16 $read1 $read2
bustools correct -w $whiteListFn -o klbusOut/output.correct.bus klbusOut/output.bus
bustools sort -t 16 -o klbusOut/output.correct.sort.bus klbusOut/output.correct.bus
mkdir klbusGenecount
bustools count -o klbusGenecount/gene -g $tgMapFile -e klbusOut/matrix.ec -t klbusOut/transcripts.txt --genecounts klbusOut/output.correct.sort.bus

# Kallisto in single-end mode
kallisto quant -i kallistoIndex -t 16 --single -l 400 -s 100 -o kallistoOut $onecellRead2

# Salmon in single-end mode
salmon quant -i salmonIndex -l A -r $onecellRead2 -p 16 -o salmonOut


#### Parameter settings to run public tools for the Cite-seq data
read1="citeseq_RNA_1.fastq.gz"
read2="citeseq_RNA_2.fastq.gz"

#Alevin
salmon alevin -l ISR -1 $read1 -2 $read2 --chromium -i salmonIndex -p 16 -o citseqAlevinOuput --tgMap $tgMapFile --dumpBfh --dumpFeatures --dumpMtx --expectCells 35000


#### Parameter settings to run public tools for the Smart-seq2 data
read1="smartseq2_cell_1.fastq.gz"
read2="smartseq2_cell_2.fastq.gz"

salmon quant -i salmonIndex -l ISR -1 $read1 -2 $read2 -p 16 -o smartseq2_Out