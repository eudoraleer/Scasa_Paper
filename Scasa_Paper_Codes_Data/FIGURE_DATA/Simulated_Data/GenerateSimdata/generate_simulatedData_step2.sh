########################################################################
### generate the simulated dataset for the study
# Step 2 - run in shell
# - extract reads for individual cells
# - convert files to fastq format
#NOTE: require bbmap to run reformat.sh

cd CELL

perl -e '@a=`ls *_1.fasta`; chomp @a; foreach $i (@a){$i =~/(.*)_1.fasta/; $name = $1; system("awk '"'"'\{print \> \(NR \% 2 \? \"$name.odd.txt\" \: \"$name.even.txt\"\) \}'"'"' $i"); system("paste -d \\\\n $name.odd.txt $name.barcodes.txt > temp.$i");}'

perl -i -pe 's/mate1:(\d+)-(\d+)/"mate1:$1-".($1 + 28 -1)/ge' temp*_1.fasta
perl -i -pe 's/mate1:(\d+)-(\d+)/"mate1:$1-".($1 + 28 -1)/ge' *_2.fasta
perl -e '@a=`ls temp*_1.fasta`; chomp @a; foreach $i(@a){$i=~/temp\.(.*)/; $name = $1; system("mv $i $name");}'

cat *_1.fasta > simhg38_polyster_sampling_01.fasta
cat *_2.fasta > simhg38_polyster_sampling_02.fasta


reformat.sh in1=simhg38_polyster_sampling_01.fasta in2=simhg38_polyster_sampling_02.fasta out1=simhg38_polyster_sampling_01.fastq out2=simhg38_polyster_sampling_02.fastq overwrite=true

########################################################################
### Done