#!/usr/bin/bash

#####################################################
# Ribo-Seq and RNA-Seq data - Uniquely mapped reads #
#####################################################
       
## Ribo-Seq: to remove rRNA
### Human 5S DNA
#### https://www.ncbi.nlm.nih.gov/nuccore/23898
### Human ribosomal DNA complete repeating unit
#### https://www.ncbi.nlm.nih.gov/nuccore/555853

PATH_FILES=$1 # save the folder of the RiboSeq data
REF_GENOME="ReferenceGenome"
TRIM_PATH="Trimmed_RiboASeq"
NEW_PATH_Cont="RiboSeq_RM_Contaminant"
REFERENCE="FinalIsoformAnnotation" # save the folder for the novel transcriptome
PATHNEWFOLDER="Bowtie_Uniq_RiboSeq" # save the folder of the output
FASTQC_Cleaned="RiboSeq_RM_Contaminant/fastq"

mkdir $NEW_PATH_Cont
mkdir $TRIM_PATH
mkdir $PATHNEWFOLDER
mkdir $FASTQC_Cleaned

echo "Ribo-Seq reads trimming with cutadapt (in trim_galore) and contaminants removal with STAR starting the $(date)." # print an OK message
echo "Ribo-Seq reads trimming with cutadapt (in trim_galore) and contaminants removal with STAR starting the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# Contaminant indexes

STAR --runThreadN 30 --runMode genomeGenerate --genomeDir $REF_GENOME/ContaminantIndexes --genomeFastaFiles $REF_GENOME/Contaminants_rRNA_tRNA.fasta

files=($PATH_FILES/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i++)); # for loop with a counter to go inside each file
do 

	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%.*} # Get only the name of the file without the extension

	echo "Trimming ${file%.*} with cutadapt and removing rRNA and tRNA." # print an OK message

	# cutadapt 
	# -q 20	- Trim low-quality ends from reads in addition to adapter removal - default 20
	# --phred33	- ASCII+33 quality scores as Phred scores - default ON 
	# --length <INT>	- Discard reads that became shorter than length INT because of either quality or adapter trimming - default 20 bp
	# --fastqc	- Run FastQC in the default mode on the FastQ file once trimming is complete.
	# --gzip	- Compress the output file with gzip
	trim_galore --cores 4 -q 20 --output_dir $TRIM_PATH --gzip ${files[$i]}

	# STAR: mapping to the rRNA and tRNA
	STAR --runMode alignReads --runThreadN 16 --genomeDir $REF_GENOME/ContaminantIndexes --readFilesIn $TRIM_PATH/${file%.*}_trimmed.fq.gz --readFilesCommand gunzip -c --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix $NEW_PATH_Cont/${file%.*}_rm_repbase_t_rrna.fastq --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --alignEndsType EndToEnd > $NEW_PATH_Cont/${file%.*}_repbase_rrna_comtam.bam
	
	echo "${file%.*} trimmed with cutadapt and aligned to the reference genome with STAR." # print an OK message
	echo "${file%.*} trimmed with cutadapt and aligned to the reference genome with STAR." >> $REPORT_PATH/PipelineBulkRNASeq.txt # print the end message in the progress .txt
done;

mv $NEW_PATH_Cont/*.fastq > $FASTQC_Cleaned

echo "Ribo-Seq - Uniquely mapped reads using Bowtie $(date)." # print an OK message
echo "Ribo-Seq - Uniquely mapped reads using Bowtie $(date)." >> $PATHREPORT/progress.txt # print a message in the progress .txt

REFERENCE="FinalIsoformAnnotation" # save the folder for the novel transcriptome

# Create indexes
bowtie-build ../$REFERENCE/Monocyte_Isoforms.final_oneline.fasta ../$REFERENCE/Monocyte_Isoforms_RiboSeq # Built the indexes

files=($FASTQC_Cleaned/*) # copy all the file names

echo "Indexes created. Alignment with Bowtie (with modified conditions for increasing sensitivity) and extract uniquely mapped reads starting $(date)." # print an OK message
echo "Indexes created. Alignment with Bowtie (with modified conditions for increasing sensitivity) and extract uniquely mapped reads starting $(date)." >> $PATHREPORT/progress.txt # print a message in the progress .txt

for (( i=0; i<(${#files[@]}); i++)); # for loop with a counter to go inside each file
do 

	file=${files[$i]##*/} # Get only the name of the file (without folder)

	# Alignment using bowtie and the parameters included in the -y: Try as hard as possible to find valid alignments when they exist, including paired-end alignments. 
	bowtie -p 16 -m 1 --best --strata -x ../$REFERENCE/Monocyte_Isoforms_RiboSeq ${files[$i]} -S $PATHNEWFOLDER/${file%_*}.sam
	
	# get uniquely mapped reads
	samtools view -@ 16 -q 255 -h -b -o $PATHNEWFOLDER/${file%_*}_Uniquely.bam $PATHNEWFOLDER/${file%_*}.sam
	samtools sort -@ 16 -o $PATHNEWFOLDER/${file%_*}_Uniquely_Sorted.bam $PATHNEWFOLDER/${file%_*}_Uniquely.bam
	samtools index -@ 16 $PATHNEWFOLDER/${file%_*}_Uniquely_Sorted.bam
	# count reads aligning to each isoform
	samtools idxstats -@ 16 $PATHNEWFOLDER/${file%_*}_Uniquely_Sorted.bam | cut -f 1,3 > $PATHNEWFOLDER/${file%_*}_Alignments.tab
done;

## RNA-Seq

PATHFOLDER=$2 # save the folder of the RNASeq data
REFERENCE="FinalIsoformAnnotation" # save the folder for the novel transcriptome
TRIM_PATH="Trimmed_RNASeq_Thp1"
PATHNEWFOLDER="Bowtie2_Uniq_RNASeq" # save the folder of the output

mkdir "Trimmed_RNASeq_Thp1"
mkdir "Bowtie2_Uniq_RNASeq" # make output folder

echo "RNA-Seq-Uniquely mapped reads using Bowtie2 $(date)." # print an OK message
echo "RNA-Seq-Uniquely mapped reads using Bowtie2 $(date)." >> $PATHREPORT/progress.txt # print a message in the progress .txt

files=($PATHFOLDER/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i+=2)); # for loop with counter for using two files each loop
do 
	
	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%.*} # Get only the name of the file without the extension

	file_2=${files[$i+1]##*/} # Get only the name of the file (without folder)
	file2=${file_2%.*} # Get only the name of the file without the extension

	echo "Trimming ${file%.*} and ${file2%.*} with cutadapt." # print an OK message

		# cutadapt 
	# -q 20	- Trim low-quality ends from reads in addition to adapter removal - default 20
	# --phred33	- ASCII+33 quality scores as Phred scores - default ON 
	# --length <INT>	- Discard reads that became shorter than length INT because of either quality or adapter trimming - default 20 bp
	# --stringency <INT>	- Overlap with adapter sequence required to trim a sequence - Defaults to a very stringent setting of 1
	# --fastqc	- Run FastQC in the default mode on the FastQ file once trimming is complete.
	# --gzip	- Compress the output file with gzip
	trim_galore --cores 4 -q 20 --phred33 --output_dir $TRIM_PATH --fastqc --stringency 1 --gzip --length 20 $PATHFOLDER/${files[$i]} $PATHFOLDER/${files[$i+1]}

	echo "${file%.*} trimmed with cutadapt." # print an OK message
	echo "${file%.*} trimmed with cutadapt." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

done;

files=($TRIM_PATH/*) # copy all the file names

for (( i=0; i<(${#files[@]}); i+=2)); # for loop with counter for using two files each loop
do

	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%.*} # Get only the name of the file without the extension

	file_2=${files[$i+1]##*/} # Get only the name of the file (without folder)
	file2=${file_2%.*} # Get only the name of the file without the extension

	# Alignment using bowtie2 and the parameters included in the Very-sensitive-local mode, changing the N from 0 to 1 and the parameters gbar, mp, rdg and rfg one lower than default
	bowtie2 -p 16 --very-sensitive -x ../$REFERENCE/Monocyte_Isoforms -1 ${files[$i]} -2 ${files[$i+1]} -S $PATHNEWFOLDER/${file%.*}.sam

	# get uniquely mapped reads
	samtools view -@ 16 -h -S -b -q 20 -o $PATHNEWFOLDER/${file%_*}_Uniquely.bam $PATHNEWFOLDER/${file%.*}.sam
	sambamba sort -t 16 -o $PATHNEWFOLDER/${file%_*}_Uniquely_Sorted.bam $PATHNEWFOLDER/${file%_*}_Uniquely.bam
	samtools index -@ 16 $PATHNEWFOLDER/${file%_*}_Uniquely_Sorted.bam
	# count reads aligning to each isoform
	samtools idxstats $PATHNEWFOLDER/${file%_*}_Uniquely_Sorted.bam | cut -f 1,3 > $PATHNEWFOLDER/${file%_*}_Alignments.tab
done;

echo "Uniquely mapped reads for RiboSeq and RNASeq finished on $(date)." # print an OK message
echo "Uniquely mapped reads for RiboSeq and RNASeq finished on $(date)." >> $PATHREPORT/progress.txt # print an start message in the progress .txt
