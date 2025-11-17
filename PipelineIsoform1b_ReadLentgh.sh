#!/usr/bin/bash

#######################
# Compute read length #
#######################

set -e # exit if there is an error
PATHFOLDER=$1 # save the path/folder to the nanopore data
PATHNEWFOLDER=$2 # save the path to the output folder

files=($PATHFOLDER/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i++));
do 
 	echo ${files[$i]}
 	zcat ${files[$i]} | awk 'NR % 4 == 2 {print length($0)}' > $PATHNEWFOLDER/${files[$i]##*/}.txt # compute leanth of each read
done;
