#!/usr/bin/bash

####################################################
# Paired-end Short read alignment and re-alignment #
####################################################
           
PATHFOLDER=$1
REF_GENOME="ReferenceGenome"
NEW_PATH="Alignments_ShortReads_Validation"
REPORT_PATH="Reports"
TRIM_PATH="Trimmed_BulkRNASeq_Validation"

mkdir "Trimmed_BulkRNASeq_Validation"
mkdir "Alignments_ShortReads_Validation"

trim_galore --version
echo "STAR version:"
STAR --version

echo "Short read trimming and Re-alignment and quantification of short reads with salmon starting the $(date)." # print an OK message
echo "Short read trimming and Re-alignment and quantification of short reads with salmon starting the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

STAR --runThreadN 30 --runMode genomeGenerate --sjdbGTFfile $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --genomeDir $REF_GENOME/GRCh38_primary_assembly --genomeFastaFiles $REF_GENOME/GRCh38.primary_assembly.genome.fa --sjdbGTFtagExonParentTranscript Parent

files=($PATHFOLDER/*) # copy all the file names in the PATHFOLDER folder

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

	# salmon
	##  quant  quantify abundances of the transcripts
   	##  -l is the library type
	salmon quant -i $REF_GENOME/Monocyte_Isoforms_Salmon -l A -1 $TRIM_PATH/${file%.*}_1_val_1.fq.gz -2 $TRIM_PATH/${file2%.*}_2_val_2.fq.gz --validateMappings -o $Salmon_QUANT/${file%.*}_salmon

	echo "${file%.*} re-aligned to the new reference transcriptome and quantified using salmon." # print an OK message
	echo "${file%.*} re-aligned to the new reference transcriptome and quantified using salmon." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
done;

echo "Trimming and quantification of short reads finished the $(date)." # print an OK message
echo "Pipeline step 2: trimming and quantification of short reads finished the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# SQUANTI validation with new data set
BulkAlign_PATH="Alignments_ShortReads_Validation"
python -W ignore $SQANTI3_PATH/sqanti3_qc.py FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf $REF_GENOME/GRCh38.primary_assembly.genome.fa -t 30 -o Monocyte_Isoform_NewData -d $SQANTI3_OUT_PATH/FinalQC_NewData/ --aligner_choice minimap2 --polyA_peak $REF_GENOME/SQANTI3_Extradata/atlas.clusters.2.0.GRCh38.96_ok.bed --CAGE_peak $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3_ok.bed --polyA_motif_list $SQANTI3_PATH/data/polyA_motifs/mouse_and_human.polyA_motif.txt --force_id_ignore --isoAnnotLite --report both --fl $FLAIR_OUT_PATH/flairAnalysis_Monocytes_quantified_format.tsv --SR_bam $BulkAlign_PATH/BamFiles/ -c $BulkAlign_PATH/SpliceJunction/ --gff3 $REF_GENOME/SQANTI3_Extradata/Homo_sapiens_GRCh38_Ensembl_86.gff3
