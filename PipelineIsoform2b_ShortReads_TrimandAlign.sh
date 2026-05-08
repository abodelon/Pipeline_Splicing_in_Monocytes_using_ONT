#!/usr/bin/bash

##########################################
# Single-end Short read QC and alignment #
##########################################

# reference genome: https://www.gencodegenes.org/human/
           
PATH_FILES=$1
REF_GENOME="ReferenceGenome"
NEW_PATH="Alignments_ShortReads_InVivo"
REPORT_PATH="Reports"
TRIM_PATH="Trimmed_BulkRNASeq_InVivo"

mkdir "Trimmed_BulkRNASeq_InVivo"
mkdir "Alignments_ShortReads_InVivo"

trim_galore --version
echo "STAR version:"
STAR --version

echo "Short reads trimming with cutadapt (in trim_galore) and aligning with STAR starting the $(date)." # print an OK message
echo "Short reads trimming with cutadapt (in trim_galore) and aligning with STAR starting the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

STAR --runThreadN 30 --runMode genomeGenerate --sjdbGTFfile $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --genomeDir $REF_GENOME/GRCh38_primary_assembly --genomeFastaFiles $REF_GENOME/GRCh38.primary_assembly.genome.fa --sjdbGTFtagExonParentTranscript Parent

files=($PATH_FILES/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i+=2)); # for loop with counter for using two files each loop
do

	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%_*} # Get only the name of the file without the extension

	file_2=${files[$i+1]##*/} # Get only the name of the file (without folder)
	file2=${file_2%_*} # Get only the name of the file without the extension

	echo "Trimming ${file%.*} and ${file2%.*} with cutadapt and quantification with salmon." # print an OK message

	# cutadapt 
	# -q 20	- Trim low-quality ends from reads in addition to adapter removal - default 20
	# --phred33	- ASCII+33 quality scores as Phred scores - default ON 
	# --length <INT>	- Discard reads that became shorter than length INT because of either quality or adapter trimming - default 20 bp
	# --stringency <INT>	- Overlap with adapter sequence required to trim a sequence - Defaults to a very stringent setting of 1
	# --fastqc	- Run FastQC in the default mode on the FastQ file once trimming is complete.
	# --gzip	- Compress the output file with gzip
	trim_galore --cores 8 -q 20 --phred33 --output_dir $TRIM_PATH --fastqc --stringency 1 --gzip --length 20 --paired ${files[$i]} ${files[$i+1]}

	# STAR
	# --quantMode GeneCounts	- Count number reads per gene while mapping using htseq-count with default parameters
	# --outSAMprimaryFlag AllBestScore - which alignments are considered primary - all others will be marked with 0x100 bit in the FLAG - all alignments with the best score are primary (non-default)
	# --outFilterMultimapScoreRange 0	- score range below the maximum score for multimapping alignments (non-default, default 1) - keep multialignments
	# --alignEndsType EndToEnd	- type of read ends alignment - EndToEnd to force end-to-end read alignment, do not soft-clip (non-default, default is Local with soft-clip)
	# --twopassMode Basic	- basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly
	# STAR --runThreadN 16 --outBAMsortingThreadN 16 --runMode alignReads --twopassMode Basic --quantMode GeneCounts --sjdbGTFfile $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $NEW_PATH/${file%.*} --genomeDir $REF_GENOME/GRCh38_primary_assembly --readFilesCommand gunzip -c --readFilesIn $TRIM_PATH/${file%.*}_val_1.fq.gz $TRIM_PATH/${file2%.*}_val_2.fq.gz
	STAR --runThreadN 16 --outBAMsortingThreadN 16 --runMode alignReads --quantMode GeneCounts --sjdbGTFfile $REF_GENOME/gencode.v47.primary_assembly.annotation.gtf --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $NEW_PATH/${file%.*} --genomeDir $REF_GENOME/GRCh38_primary_assembly --readFilesCommand gunzip -c --readFilesIn $TRIM_PATH/${file%.*}_1_val_1.fq.gz $TRIM_PATH/${file2%.*}_2_val_2.fq.gz
	
	echo "${file%.*} trimmed with cutadapt and aligned to the reference genome with STAR." # print an OK message
	echo "${file%.*} trimmed with cutadapt and aligned to the reference genome with STAR." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
done;

echo "Trimming and quantification of short reads finished the $(date)." # print an OK message
echo "Trimming and quantification of short reads finished the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

