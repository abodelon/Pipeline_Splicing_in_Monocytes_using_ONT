#!/usr/bin/bash

###########################
# Short read re-alignment #
###########################

PATH_FILES=$1
REF_GENOME="ReferenceGenome"
REPORT_PATH="Reports"
TRIM_PATH="Trimmed_BulkRNASeq"
NEW_TRANSC="FinalIsoformAnnotation"
KAL_PATH="/home/alejandra/kallisto-0.51.0/build/src"
Salmon_QUANT="Salmon_Realign_Quant"
Salmon_QUANT_Align="Salmon_Realign_Quant_AlignBased"
REAlignments_ShortReads="REAlignments_ShortReads"

mkdir "Salmon_Realign_Quant"
mkdir "Salmon_Realign_Quant_AlignBased"

echo "Re-alignment and quantification of short reads with salmon starting the $(date)." # print an OK message
echo "PRe-alignment and quantification of short reads with salmon starting the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# Create indexes (STAR, salmon, kallisto)
salmon index -t $REF_GENOME/SalmonRef.fasta -i $REF_GENOME/Monocyte_Isoforms_Salmon --decoys $REF_GENOME/decoys.txt -k 31

files=($PATH_FILES/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i++)); # for loop with a counter to go inside each file
do 

	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%.*} # Get only the name of the file without the extension

	# salmon
	##  quant  quantify abundances of the transcripts
   	##  -l is the library type
	salmon quant -i $REF_GENOME/Monocyte_Isoforms_Salmon -l A -r $TRIM_PATH/${file%.*}_trimmed.fq.gz --validateMappings -o $Salmon_QUANT/${file%.*}_salmon	
	
	echo "${file%.*} re-aligned to the new reference transcriptome and quantified using kallisto." # print an OK message
	echo "${file%.*} re-aligned to the new reference transcriptome and quantified using kallisto." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

done;

echo "Re-alignment and quantification of short reads finished the $(date)." # print an OK message
echo "Re-alignment and quantification of short reads finished the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
