#!/usr/bin/bash

##################################################################################################
# FLAIR expression quantification and differentiall splicing analysis for long reads uisng SUPPA #
##################################################################################################

NEW_PATH="FinalIsoformAnnotation"
REPORT_PATH="Reports"
FLAIR_PATH="/home/alejandra/flair-master"
IUsage_PATH="IsoformUsage"
REF_GENOME="ReferenceGenome"

echo "Start flair quantify and differential splicing analysis at $(date)." # print an OK message
echo "Start flair quantify and differential splicing analysis at $(date)."  >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# from genePred to bed
/home/alejandra/genePredToBed $NEW_PATH//Monocyte_Isoforms.final.genePred $NEW_PATH/Monocyte_Isoforms.final.bed

mkdir "FinalIsoformAnnotation/flairquantify"

# one-line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $NEW_PATH/Monocyte_Isoforms.final.fasta > $NEW_PATH/Monocyte_Isoforms.final_oneline.fasta

# prepare fasta and bed (add to the fasta and bed the Isoform.id_gene.id)
Rscript PipelineIsoform5b_DifferentialSplicingEvents.R

# quantify the expression of each isoform using flair quantify (stringent does not seem to do anything because was already done before)
python $FLAIR_PATH/flair.py quantify -t 30 -o $NEW_PATH/Monocyte_Isoforms_quantified --isoform_bed $NEW_PATH/Monocyte_Isoforms.final_DEG.bed --temp_dir "FinalIsoformAnnotation/flairquantify" --isoforms $NEW_PATH/Monocyte_Isoforms.final_DEG.fasta  --reads_manifest $REF_GENOME/metadata/reads_manifestFile.txt --tpm

echo "flair quantify finished at $(date)." # print an OK message
echo "flair quantify finished at $(date)."  >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# prepare the annotations
awk '{if ($3=="gene") print $0}' $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf | cut -f9 | awk '{print $2"\t"$4"\t"$6}' | sed 's/"//g' | sed 's/;//g' > $REF_GENOME/gencode.v46.primary_assembly.gene_annotation.gtf
 
### Productivity prediction

# detect productivity
python $FLAIR_PATH/bin/predictProductivity -i $NEW_PATH/Monocyte_Isoforms.final_DEG.bed -f $REF_GENOME/GRCh38.primary_assembly.genome.fa -g $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --longestORF --append_column > $NEW_PATH/Monocyte_Isoforms_productivity.bed
cut -f4 $NEW_PATH/Monocyte_Isoforms_productivity.bed | awk -F_ '{print $1"\t"$2"\t"$3}' > $NEW_PATH/Monocyte_Isoforms_productivity.tab

echo "Productivity prediction" # print the end message in the progress .txt
echo "Productivity prediction" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
 
### Longest-ORF prediction 

# detect the isoform with the longest ORF
python PipelineIsoform5c_FuntionalityAnalysis_MaxORF.py

echo "Funcitonality prediction" # print the end message in the progress .txt
echo "Funcitonality prediction" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

### Differential Splicing Analysis: SUPPA 

# SUPPA reads the ioi or ioe file generated in the previous step and a transcript expression file with the transcript abundances (TPM units) to calculate the relative abundance (PSI) value per sample for each transcript or each local event.
awk '{sub(/_.*/, "", $1); print $1"\t"$2"\t"$4"\t"$6}' FinalIsoformAnnotation/Monocyte_Isoforms_quantified.tpm.tsv > FinalIsoformAnnotation/Monocyte_Isoforms_quantified_Exvivo.tpm.tsv
awk '{sub(/_.*/, "", $1); print $1"\t"$3"\t"$5"\t"$7}' FinalIsoformAnnotation/Monocyte_Isoforms_quantified.tpm.tsv > FinalIsoformAnnotation/Monocyte_Isoforms_quantified_LPS.tpm.tsv
cut -f1 FinalIsoformAnnotation/Monocyte_Isoforms_quantified_Exvivo.tpm.tsv > FinalIsoformAnnotation/ExpressedIsoforms.txt
grep -f FinalIsoformAnnotation/ExpressedIsoforms.txt FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf > FinalIsoformAnnotation/Monocyte_Isoforms.final.expressed.gtf

suppa.py generateEvents -i FinalIsoformAnnotation/Monocyte_Isoforms.final.expressed.gtf -o FinalIsoformAnnotation/SUPPA/Monocyte_Isoform.events.expressed -e SE SS MX RI FL -f ioe
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' FinalIsoformAnnotation/SUPPA/Monocyte_Isoform.events.expressed*.ioe > FinalIsoformAnnotation/Monocyte_Isoform.events.expressed.ioe

suppa.py psiPerEvent --i FinalIsoformAnnotation/Monocyte_Isoform.events.expressed.ioe -e FinalIsoformAnnotation/Monocyte_Isoforms_quantified_Exvivo.tpm.tsv -o FinalIsoformAnnotation/SUPPA/Monocyte_Isoform_Exvivo.perEvent
suppa.py psiPerEvent --i FinalIsoformAnnotation/Monocyte_Isoform.events.expressed.ioe -e FinalIsoformAnnotation/Monocyte_Isoforms_quantified_LPS.tpm.tsv -o FinalIsoformAnnotation/SUPPA/Monocyte_Isoform_LPS.perEvent

suppa.py diffSplice --method empirical --input FinalIsoformAnnotation/Monocyte_Isoform.events.expressed.ioe --psi FinalIsoformAnnotation/SUPPA/Monocyte_Isoform_Exvivo.perEvent.psi FinalIsoformAnnotation/SUPPA/Monocyte_Isoform_LPS.perEvent.psi --tpm FinalIsoformAnnotation/Monocyte_Isoforms_quantified_Exvivo.tpm.tsv FinalIsoformAnnotation/Monocyte_Isoforms_quantified_LPS.tpm.tsv -pa -gc -s -o Monocyte_Isoform.perEvent_Diff
