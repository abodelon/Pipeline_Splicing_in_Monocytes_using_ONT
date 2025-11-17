#!/usr/bin/bash

#####################################
# Overlap novel annotation with TSS #
#####################################

set -e # exit if there is an error
PATHFOLDER="Alignments" # save the forlder information
PATHREPORT="Reports" # create the path of the reports directory
PATHNEWFOLDER="Alignments/BedFromBam" # create the path of the new forlder directory
REFERENCE="FinalIsoformAnnotation" # save the folder for the novel transcriptome

mkdir $PATHNEWFOLDER # create a new folder for the outputs

echo "Overlap between alitghed reads and annotated isoforms (TSS) starting the $(date)." # print an OK message
echo "Overlap between alitghed reads and annotated isoforms (TSS) starting the $(date)." >> $PATHREPORT/progress.txt # print a message in the progress .txt

files=($( ls $PATHFOLDER/ | grep '.bam' | grep -v '.bam.bai' )) # copy all the file names

for (( i=0; i<(${#files[@]}); i++)); # for loop with a counter to go inside each file
do 
	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%.*} # Get only the name of the file without the extension

	bedtools intersect -wao -a $PATHNEWFOLDER/${file}Conver.bed -b $REFERENCE/Monocyte_Isoforms.final.bed > $PATHNEWFOLDER/Monocyte_IsoformsOverlap_${file}.bed

	echo "Read overlap with annotated TSS created for ${files[$i]##*/}!" # print an OK message
	echo "Read overlap with annotated TSS created for ${files[$i]##*/}." >> $PATHREPORT/progress.txt # print a message in the progress .txt

done;

echo "Overlap between aligned reads and TSS finished the $(date). Overlap stored in $PATHNEWFOLDER forlder!" # print an OK message
echo "Overlap between aligned reads and TSS finished the $(date)." >> $PATHREPORT/progress.txt # print the end message in the progress .txt
