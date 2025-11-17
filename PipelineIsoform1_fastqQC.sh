#!/usr/bin/bash

###################
# Quality Control #
###################

set -e # exit if there is an error
PATHFOLDER=$1 # save the path/folder to the RNA-Seq data
PATHREPORT="Reports" # save the path of the reports
PATHNEWFOLDER=$2 # save the path to the output folder

mkdir $2 # create a new folder for the quality outputs
mkdir "Reports" # create a new folder for the reports

echo "Process starting the $(date)." # print an OK message
echo "FastQC report process starting the $(date)." >> $PATHREPORT/progress.txt # print a message in the progress .txt

fastqc --version # print the version of the program for quality control (fastqc)

files=($PATHFOLDER/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i++));
do 
 	fastqc ${files[$i]} -o $PATHNEWFOLDER # create the fastqc report of each fastq file and save it in the new folder
	echo "FastQC ${files[$i]##*/} report created!" # print an OK message
	echo "FastQC ${files[$i]##*/} report created." >> $PATHREPORT/progress.txt # print a message in the progress .txt
done;

multiqc --version # print the version of the program for mergering all quality control reports (multiqc)

multiqc $PATHNEWFOLDER/ -o $PATHNEWFOLDER/ # merge all quality reports located in the new folder together

echo "Process finished the $(date). Summary of all reports stored in the $PATHNEWFOLDER forlder!" # print an OK message
echo "FastQC reports finished. Date: $(date)." >> $PATHREPORT/progress.txt # print the end message in the progress .txt
