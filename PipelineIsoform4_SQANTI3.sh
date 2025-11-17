#!/usr/bin/bash

###########################
# SQANTI3 Quality Control #
###########################

# More information in SQANTI3 in: https://github.com/ConesaLab/SQANTI3

## Downloads to put on the ReferenceGenome folder:
# tappAS (human) from: https://app.tappas.org/resources/downloads/gffs/
# polyA_peak file, from: https://polyasite.unibas.ch/atlas#2

set -e # exit if there is an error
REF_GENOME="ReferenceGenome"
REPORT_PATH="Reports"
FLAIR_OUT_PATH="FLAIRanalysis"
BulkTrim_PATH="Trimmed_BulkRNASeq"
SQANTI3_PATH="/home/alejandra/SQANTI3-5.2.1"
SQANTI3_OUT_PATH="SQANTI3"
FLAIR_PATH="/home/alejandra/flair-master"
BulkAlign_PATH="Alignments_ShortReads"

echo "Running SQANTI3: quality control, filter and rescue" # print a progress 

# Use the STAR outputs as a short-read support for quality control 
# organize short-reads information
mkdir "$BulkAlign_PATH/SpliceJunction"
mv $BulkAlign_PATH/*SJ.out.tab $BulkAlign_PATH/SpliceJunction/

mkdir "$BulkAlign_PATH/BamFiles"
mv $BulkAlign_PATH/*Aligned.sortedByCoord.out.bam $BulkAlign_PATH/BamFiles/

# --polyA_peak file, downloaded from: https://polyasite.unibas.ch/atlas#2
sed 's/^[[:digit:]]/chr&/' $REF_GENOME/SQANTI3_Extradata/atlas.clusters.2.0.GRCh38.96.bed  | sed 's/^Y/chr&/' | sed 's/^X/chr&/' > $REF_GENOME/SQANTI3_Extradata/atlas.clusters.2.0.GRCh38.96_ok.bed

# sqanti3_qc.py (quality control):
## sqanti3_qc.py [-h] [--min_ref_len MIN_REF_LEN] [--force_id_ignore] [--aligner_choice {minimap2,deSALT,gmap,uLTRA}] [--CAGE_peak CAGE_PEAK] [--polyA_motif_list POLYA_MOTIF_LIST] [--polyA_peak POLYA_PEAK] 
## [--phyloP_bed PHYLOP_BED] [--skipORF] [--is_fusion] [--orf_input ORF_INPUT] [--fasta] [-e EXPRESSION] [-x GMAP_INDEX] [-t CPUS] [-n CHUNKS] [-o OUTPUT] [-d DIR] [-c COVERAGE] [-s SITES] [-w WINDOW] [--genename] 
## [-fl FL_COUNT] [-v] [--saturation] [--report {html,pdf,both,skip}] [--isoAnnotLite] [--gff3 GFF3] [--short_reads SHORT_READS] [--SR_bam SR_BAM] [--isoform_hits] [--ratio_TSS_metric {max,mean,median,3quartile}]  isoforms annotation genome
### --force_id_ignore  Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)
### --aligner_choice   minimap2
### --CAGE_peak CAGE_PEAK   FANTOM5 Cage Peak (BED format, optional) - file used by flair
### --polyA_motif_list POLYA_MOTIF_LIST Ranked list of polyA motifs (text, optional)
### --polyA_peak POLYA_PEAK PolyA Peak (BED format, optional) (downloaded from: https://polyasite.unibas.ch/atlas#2) 
### -fl FL_COUNT, --fl_count FL_COUNT   Full-length PacBio abundance file
### --isoAnnotLite  Run isoAnnot Lite to output a tappAS-compatible gff3 file
### --saturation    Include saturation curves into report
### --short_reads SHORT_READS   File Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.
### --gff3 GFF3 Precomputed tappAS species specific GFF3 file. It will serve as reference to transfer functional attributes (human, downloaded from: https://app.tappas.org/resources/downloads/gffs/)
echo "python -W ignore $SQANTI3_PATH/sqanti3_qc.py $FLAIR_OUT_PATH/flairAnalysis_monocytes_collapsed.isoforms.gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf $REF_GENOME/GRCh38.primary_assembly.genome.fa -t 30 -o SQANTI3_Mon -d $SQANTI3_OUT_PATH/ --aligner_choice minimap2 --polyA_peak $REF_GENOME/SQANTI3_Extradata/atlas.clusters.2.0.GRCh38.96_ok.bed --CAGE_peak $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3_ok.bed --polyA_motif_list $SQANTI3_PATH/data/polyA_motifs/mouse_and_human.polyA_motif.txt --force_id_ignore --isoAnnotLite --report both --fl $FLAIR_OUT_PATH/flairAnalysis_Monocytes_quantified_format.tsv --SR_bam $BulkAlign_PATH/BamFiles/ -c $BulkAlign_PATH/SpliceJunction/ --gff3 $REF_GENOME/SQANTI3_Extradata/Homo_sapiens_GRCh38_Ensembl_86.gff3"
python -W ignore $SQANTI3_PATH/sqanti3_qc.py $FLAIR_OUT_PATH/flairAnalysis_monocytes_collapsed.isoforms.gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf $REF_GENOME/GRCh38.primary_assembly.genome.fa -t 30 -o SQANTI3_Mon -d $SQANTI3_OUT_PATH/ --aligner_choice minimap2 --polyA_peak $REF_GENOME/SQANTI3_Extradata/atlas.clusters.2.0.GRCh38.96_ok.bed --CAGE_peak $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3_ok.bed --polyA_motif_list $SQANTI3_PATH/data/polyA_motifs/mouse_and_human.polyA_motif.txt --force_id_ignore --isoAnnotLite --report both --fl $FLAIR_OUT_PATH/flairAnalysis_Monocytes_quantified_format.tsv --SR_bam $BulkAlign_PATH/BamFiles/ -c $BulkAlign_PATH/SpliceJunction/ --gff3 $REF_GENOME/SQANTI3_Extradata/Homo_sapiens_GRCh38_Ensembl_86.gff3

mkdir "SQANTI3/Filtered"

# sqanti3_filter.py: filter isoforms by quality (default set of rules)
## sqanti3_filter.py rules [-h] [--isoAnnotGFF3 ISOANNOTGFF3] [--isoforms ISOFORMS] [--gtf GTF] [--sam SAM] [--faa FAA][-o OUTPUT] [-d DIR][-v] [--skip_report] [-j JSON_FILTER] sqanti_class
### utilities/filter/filter_default.json $SQANTI3_OUT_PATH/SQANTI3_Mon_classification.txt 
### --isoAnnotGFF3 ISOANNOTGFF3 isoAnnotLite GFF3 file to be filtered
### --faa FAA  ORF prediction faa file to be filtered by SQANTI3
### --isoforms ISOFORMS   fasta/fastq isoform file to be filtered
### --gtf GTF             GTF file to be filtered
### --sam SAM             SAM alignment of the input fasta/fastq
### -j JSON_FILTER, --json_filter JSON_FILTER JSON file where filtering rules are expressed. Default: utilities/filter/filter_default.json
echo "python -W ignore $SQANTI3_PATH/sqanti3_filter.py rules $SQANTI3_OUT_PATH/SQANTI3_Mon_classification.txt --isoAnnotGFF3 $SQANTI3_OUT_PATH/SQANTI3_Mon.gff3 --isoforms $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.fasta --gtf $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.gtf --faa $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.faa -j $SQANTI3_PATH/utilities/filter/filter_default_modified.json -d $SQANTI3_OUT_PATH/Filtered -o "SQANTI3_Mon""
python -W ignore $SQANTI3_PATH/sqanti3_filter.py rules $SQANTI3_OUT_PATH/SQANTI3_Mon_classification.txt --isoAnnotGFF3 $SQANTI3_OUT_PATH/SQANTI3_Mon.gff3 --isoforms $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.fasta --gtf $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.gtf --faa $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.faa -j $SQANTI3_PATH/utilities/filter/filter_default_modified.json -d $SQANTI3_OUT_PATH/Filtered -o "SQANTI3_Mon"

mkdir "SQANTI3/Filtered/Rescue/"

# sqanti3_rescue.py: uses the long read-based evidence provided by discarded isoforms (i.e. artifacts) to recover transcripts in the associated reference transcriptome
## sqanti3_rescue.py rules [-h] [--isoforms ISOFORMS] [--gtf GTF] [-g REFGTF] [-f REFGENOME] [-k REFCLASSIF] [-e {all,fsm,none}] [--mode {automatic,full}] [-o OUTPUT] [-d DIR] [-v] [-j JSON]sqanti_filter_classif
### --isoforms  Long read-defined isoforms
### --gtf Long read-defined transcriptome annotation
### --refGTF  Reference transcriptome annotation
### --refClassif  Reference transcriptome classification file (SQANTI3 QC)
### --mode automatic  Those reference transcripts for which all FSM representatives were removed by the filter are rescued, only high-confidence reference isoforms
### -e fsm mono-exonic transcripts will only be considered for the rescue if they belong to the full-splice match (FSM) category. In this case, mono-exons will be rescued via the automatic rescue strategy
echo "python -W ignore $SQANTI3_PATH/sqanti3_rescue.py rules --isoforms $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.fasta --gtf $SQANTI3_OUT_PATH/Filtered/SQANTI3_Mon.filtered.gtf --refGTF $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --refGenome $REF_GENOME/GRCh38.primary_assembly.genome.fa --refClassif $SQANTI3_OUT_PATH/SQANTI3_Mon_classification.txt --mode automatic -e fsm -j $SQANTI3_PATH/utilities/filter/filter_default_modified.json -d $SQANTI3_OUT_PATH/Filtered/Rescue/ -o "SQANTI3_Mon" $SQANTI3_OUT_PATH/Filtered/SQANTI3_Mon_RulesFilter_result_classification.txt"
python -W ignore $SQANTI3_PATH/sqanti3_rescue.py rules --isoforms $SQANTI3_OUT_PATH/SQANTI3_Mon_corrected.fasta --gtf $SQANTI3_OUT_PATH/Filtered/SQANTI3_Mon.filtered.gtf --refGTF $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --refGenome $REF_GENOME/GRCh38.primary_assembly.genome.fa --refClassif $SQANTI3_OUT_PATH/SQANTI3_Mon_classification.txt --mode automatic -e fsm -j $SQANTI3_PATH/utilities/filter/filter_default_modified.json -d $SQANTI3_OUT_PATH/Filtered/Rescue/ -o "SQANTI3_Mon" $SQANTI3_OUT_PATH/Filtered/SQANTI3_Mon_RulesFilter_result_classification.txt

echo "SQANTI3 rescue finished"
echo "SQANTI3 rescue finished" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

echo "Get of SQANTI3 information by runing SQANTI3_qc in the final annotation"
echo "Get of SQANTI3 information by runing SQANTI3_qc in the final annotation" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

mkdir "$SQANTI3_OUT_PATH/FinalQC" # final QC of the rescued reads

# sqanti3_qc.py (quality control):
## sqanti3_qc.py [-h] [--min_ref_len MIN_REF_LEN] [--force_id_ignore] [--aligner_choice {minimap2,deSALT,gmap,uLTRA}] [--CAGE_peak CAGE_PEAK] [--polyA_motif_list POLYA_MOTIF_LIST] [--polyA_peak POLYA_PEAK] 
## [--phyloP_bed PHYLOP_BED] [--skipORF] [--is_fusion] [--orf_input ORF_INPUT] [--fasta] [-e EXPRESSION] [-x GMAP_INDEX] [-t CPUS] [-n CHUNKS] [-o OUTPUT] [-d DIR] [-c COVERAGE] [-s SITES] [-w WINDOW] [--genename] 
## [-fl FL_COUNT] [-v] [--saturation] [--report {html,pdf,both,skip}] [--isoAnnotLite] [--gff3 GFF3] [--short_reads SHORT_READS] [--SR_bam SR_BAM] [--isoform_hits] [--ratio_TSS_metric {max,mean,median,3quartile}]  isoforms annotation genome
### --force_id_ignore  Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)
### --aligner_choice   minimap2
### --CAGE_peak CAGE_PEAK   FANTOM5 Cage Peak (BED format, optional) - file used by flair
### --polyA_motif_list POLYA_MOTIF_LIST Ranked list of polyA motifs (text, optional)
### --polyA_peak POLYA_PEAK PolyA Peak (BED format, optional) (downloaded from: https://polyasite.unibas.ch/atlas#2) 
### -fl FL_COUNT, --fl_count FL_COUNT   Full-length PacBio abundance file
### --isoAnnotLite  Run isoAnnot Lite to output a tappAS-compatible gff3 file
### --saturation    Include saturation curves into report
### --short_reads SHORT_READS   File Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.
### --gff3 GFF3 Precomputed tappAS species specific GFF3 file. It will serve as reference to transfer functional attributes (human, downloaded from: https://app.tappas.org/resources/downloads/gffs/)
echo "python -W ignore $SQANTI3_PATH/sqanti3_qc.py $SQANTI3_OUT_PATH/Filtered/Rescue/SQANTI3_Mon_rescued.gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf $REF_GENOME/GRCh38.primary_assembly.genome.fa -t 30 -o Monocyte_Isoform -d $SQANTI3_OUT_PATH/FinalQC/ --aligner_choice minimap2 --polyA_peak $REF_GENOME/SQANTI3_Extradata/atlas.clusters.2.0.GRCh38.96_ok.bed --CAGE_peak $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3_ok.bed --polyA_motif_list $SQANTI3_PATH/data/polyA_motifs/mouse_and_human.polyA_motif.txt --force_id_ignore --isoAnnotLite --report both --fl $FLAIR_OUT_PATH/flairAnalysis_Monocytes_quantified_format.tsv --SR_bam $BulkAlign_PATH/BamFiles/ -c $BulkAlign_PATH/SpliceJunction/ --gff3 $REF_GENOME/SQANTI3_Extradata/Homo_sapiens_GRCh38_Ensembl_86.gff3"
python -W ignore $SQANTI3_PATH/sqanti3_qc.py $SQANTI3_OUT_PATH/Filtered/Rescue/SQANTI3_Mon_rescued.gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf $REF_GENOME/GRCh38.primary_assembly.genome.fa -t 30 -o Monocyte_Isoform -d $SQANTI3_OUT_PATH/FinalQC/ --aligner_choice minimap2 --polyA_peak $REF_GENOME/SQANTI3_Extradata/atlas.clusters.2.0.GRCh38.96_ok.bed --CAGE_peak $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3_ok.bed --polyA_motif_list $SQANTI3_PATH/data/polyA_motifs/mouse_and_human.polyA_motif.txt --force_id_ignore --isoAnnotLite --report both --fl $FLAIR_OUT_PATH/flairAnalysis_Monocytes_quantified_format.tsv --SR_bam $BulkAlign_PATH/BamFiles/ -c $BulkAlign_PATH/SpliceJunction/ --gff3 $REF_GENOME/SQANTI3_Extradata/Homo_sapiens_GRCh38_Ensembl_86.gff3

mkdir "FinalIsoformAnnotation" # final annotation folder

mv $SQANTI3_OUT_PATH/FinalQC/Monocyte_Isoform_corrected.gtf FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf # move final files to the final folder
mv $SQANTI3_OUT_PATH/FinalQC/Monocyte_Isoform.gff3 FinalIsoformAnnotation/Monocyte_Isoforms_tappas.final.gff3 # move final files to the final folder
mv $SQANTI3_OUT_PATH/FinalQC/Monocyte_Isoform_classification.txt FinalIsoformAnnotation/Monocyte_Isoforms_classification.final.txt # move final files to the final folder
mv $SQANTI3_OUT_PATH/FinalQC/Monocyte_Isoform_corrected.faa FinalIsoformAnnotation/Monocyte_Isoforms_predictedprot.final.faa # move final files to the final folder
mv $SQANTI3_OUT_PATH/FinalQC/Monocyte_Isoform_corrected.fasta FinalIsoformAnnotation/Monocyte_Isoforms.final.fasta # move final files to the final folder
mv $SQANTI3_OUT_PATH/FinalQC/Monocyte_Isoform_corrected.genePred FinalIsoformAnnotation/Monocyte_Isoforms.final.genePred # move final files to the final folder

echo "Final isoform annotation saved in the FinalIsoformAnnotation forlder"
echo "Final isoform annotation saved in the FinalIsoformAnnotation forlder" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

echo "Quality control, filtering and rescue of isoforms (SQANTI3), and classification by type finished at $(date)."  # print the end message in the progress .txt
echo "Quality control, filtering and rescue of isoforms (SQANTI3), and classification by type finished at $(date)."  >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

mkdir "FinalIsoformAnnotation/SUPPA" # final annotation folder
mkdir "$REF_GENOME/SUPPA" # final annotation folder

# SUPPA events annotation
suppa.py generateEvents -i FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf -o FinalIsoformAnnotation/SUPPA/Monocyte_Isoform.events -e SE SS MX RI FL -f ioe
suppa.py generateEvents -i $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf -o $REF_GENOME/SUPPA/gencode.v46.primary_assembly.annotation.events -e SE SS MX RI FL -f ioe

awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' FinalIsoformAnnotation/SUPPA/Monocyte_Isoform.events*.ioe > FinalIsoformAnnotation/Monocyte_Isoform.events.ioe
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' $REF_GENOME/SUPPA/gencode.v46.primary_assembly.annotation.events*.ioe > $REF_GENOME/gencode.v46.primary_assembly.annotation.events.ioe
