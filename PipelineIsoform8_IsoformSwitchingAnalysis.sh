#!/bin/bash

##############
#### PFAM ####
##############
# Pfam : Prediction of protein domains, which can be run either locally (using the pfam_scan.pl script as described in the readme found here or via their webserver.
# PFAM downloaded from https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# gzip -d ReferenceGenome/IsoformSwitchAnalyzeRExtraData/Pfam-A.hmm.gz
# hmmpress ReferenceGenome/IsoformSwitchAnalyzeRExtraData/Pfam-A.hmm # specific format needed for pfam_scan.pl

FILE_AA="ReferenceGenome/IsoformSwitchAnalyzeRExtraData/Sequences/InVitro_RNASeq_Published_IsoformSwitchAnalyzeR_isoform_significant_AA.fasta"
FILE_NUCL="ReferenceGenome/IsoformSwitchAnalyzeRExtraData/Sequences/InVitro_RNASeq_Published_IsoformSwitchAnalyzeR_isoform_significant_nt.fasta"
NEW_Name="InVitro_RNASeq_Published"

#############
### CPAT ####
#############
# CPAT : The Coding-Potential Assessment Tool which is a tool for predicting whether an isoform is coding or not. CPAT can be run either locally or via their webserver. This is an alternative to CPC2.
# CPAT: The Coding-Potential Assessment Tool which is a tool for predicting whether an isoform is coding or not. https://cpat.readthedocs.io/en/latest/#run-cpat-online. 
# hexamer models: https://sourceforge.net/projects/rna-cpat/files/prebuilt_models/
cpat -x ReferenceGenome/IsoformSwitchAnalyzeRExtraData/Human_Hexamer.tsv --antisense -d ReferenceGenome/IsoformSwitchAnalyzeRExtraData/Human_logitModel.RData -g $FILE_NUCL -o ReferenceGenome/IsoformSwitchAnalyzeRExtraData/cpat_all_$NEW_Name.result

##############
### TMHMM ####
##############
## https://services.healthtech.dtu.dk/services/TMHMM-2.0/
perl /home/alejandra/tmhmm-2.0c.Linux/tmhmm-2.0c/bin/tmhmm $FILE_AA > ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.txt
cut -f1,3,4 ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.txt > ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.gff 

# Read through the input file
while IFS= read -r line; do
    # If the line starts with a # (sequence metadata)
    if [[ "$line" =~ ^# ]]; then
        # Check if we're moving to a new sequence (if current_sequence is set, add //)
        if [[ -n "$current_sequence" ]]; then
            echo "//" >> ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.gff3
        fi
        # Output the sequence metadata (e.g., "# seq1 Length: 23")
        echo "$line" >> ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.gff3
        # Extract the sequence name from the metadata line
        current_sequence=$(echo "$line" | grep -oP '#\s*(\S+)' | awk '{print $2}')
    elif [[ "$line" =~ ^[a-zA-Z0-9]+ ]]; then
        # Print the sequence feature (e.g., "seq1 inside 1 152")
        echo "$line" >> ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.gff3
    fi
done < ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.gff

# Add // at the end of the file (after the last sequence block)
echo "//" >> ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.tmhmm.gff3

###########################
### IUPred2A - IUPred3 ####
###########################
# Predicts Intrinsically Disordered Regions (IDR) and Intrinsically Disordered Binding Regions (IDBR) — the parts of a protein protein which does not have a fixed three-dimensional structure (as opposite protein domains). 
# It can be run either locally or via their webserver. IUPred3 is also supported by IsoformSwitchAnalyzeR (as it produces identical result files).
# https://iupred2a.elte.hu/help_new
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $FILE_AA > ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.OneLine.fasta
sed -i '1d' ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.OneLine.fasta

mkdir "ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/"
mkdir "ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Results"

while read -r line1 && read -r line2; do
	echo -e "$line1\n$line2" > ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Sequence_$NEW_Name.fasta
  	python3 /home/alejandra/iupred3/iupred3.py --anchor ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Sequence_$NEW_Name.fasta long > ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Sequence_$NEW_Name.txt
  	grep -v "^#" ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Sequence_$NEW_Name.txt | sed -e "1i# POS\tAMINO ACID\tIUPRED SCORE\tANCHOR SCORE" | sed -e "1i>"${line1//>/}"" | sed -e "1i# Nucleic Acids Research 2018, Submitted" | sed -e "1i# Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi" | sed -e "1i# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding" | sed -e "1i################" > ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Sequence_Final_$NEW_Name.txt
  	cat ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Sequence_Final_$NEW_Name.txt >> ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IUPred3/Results/IUPred3_all.txt
done < ReferenceGenome/IsoformSwitchAnalyzeRExtraData/IsoformSwitchAnalyzeR_isoform_significant_AA_$NEW_Name.OneLine.fasta

#################
#### signalP ####
#################
# SignalP : Prediction of Signal Peptides — a short N-terminal sequence of a peptide indicating where a protein should be membrane bound or secreted.
# Signal peptide prediction model based on a Bert protein language model encoder and a conditional random field (CRF) decoder.
# signalp6 -f summary --format txt --fastafile /path/to/input.fasta --organism eukarya --output_dir path/to/be/saved --format txt --mode fast
# signalP: https://www.nature.com/articles/s41587-021-01156-3
signalp6 --model_dir /home/alejandra/signalp6_slow_sequential/signalp-6-package/models/ --mode slow-sequential --organism eukarya -f none --fastafile $FILE_AA --output_dir ReferenceGenome/IsoformSwitchAnalyzeRExtraData/Monocyte_Isoform_classification.final_aa_$NEW_Name.signalP.result

# prepare extract annotation
## Pfam
pfam_scan.pl -fasta $FILE_AA -dir ReferenceGenome/IsoformSwitchAnalyzeRExtraData/ -outfile ReferenceGenome/IsoformSwitchAnalyzeRExtraData/pfam_all_$NEW_Name.result

echo "Finished"

# Export for PLM modeling the predicted protein sequences
cut -f1,47 FinalIsoformAnnotation/Monocyte_Isoforms_classification.final.txt | awk '$2 != "NA" {print $1","$2}' > FinalIsoformAnnotation/Monocyte_Isoforms_predictedprot.csv