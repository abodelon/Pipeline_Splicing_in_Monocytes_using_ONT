#!/usr/bin/bash

##############################
# Long reads re-count Salmon #
##############################

REPORT_PATH="Reports"
NEW_TRANSC="FinalIsoformAnnotation"
SALMON_QUANT="Salmon_Quantification_ONT"
PATH_FASTQ_ONT=$1
PATH_ALIGN_ONT="ONT_REAlignment"

mkdir "Salmon_Quantification_ONT"
mkdir "ONT_REAlignment"

minimap2 --version
salmon -v

echo "Quantification of long reads with salmon starting the $(date)." # print an OK message
echo "Quantification of long reads with salmon starting the $(date).." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

files=($PATH_FASTQ_ONT/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i++)); # for loop with a counter to go inside each file
do 

    file_1=${files[$i]##*/} # Get only the name of the file (without folder)
    file=${file_1%.*} # Get only the name of the file without the extension

   # align with minimap - equal than flair
   minimap2 -t 30 -ax splice --secondary=no -uf -k14 $NEW_TRANSC/Monocyte_Isoforms.final_oneline.fasta ${files[$i]} > $PATH_ALIGN_ONT/${file%.*}.sam

   # salmon: https://combine-lab.github.io/salmon-tutorials/2021/ont-long-read-quantification/
   ##  quant  quantify abundances of the transcripts
   ##  --ont, is required for long read quantification
   ##  -l is the library type (U for “unstranded”
   ##  -t gives the fasta files of the reference target transcripts
   ##  -a gives the bam files with the alignments
   salmon quant --ont -p 30 -t $NEW_TRANSC/Monocyte_Isoforms.final_oneline.fasta -l U -a $PATH_ALIGN_ONT/${file%.*}.sam -o $SALMON_QUANT/${file%.*}_salmon_quant

   echo "${file%.*} quantification of the new reference transcriptome expression using salmon." # print an OK message
   echo "${file%.*} quantification of the new reference transcriptome expression using salmon." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
done;

echo "Quantification of long reads with salmon finished the $(date)." # print an OK message
echo "Quantification of long reads with salmon starting the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
