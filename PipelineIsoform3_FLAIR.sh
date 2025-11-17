#!/usr/bin/bash

#############################
# FLAIR database generation #
#############################

# Pipeline following: https://flair.readthedocs.io/en/latest/index.html

## Downloads to put on the ReferenceGenome folder:
# CAGE_peaks: https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/DPI_clustering/ - (permissive set! not robust, used for TSS exploration)
# refTSS: https://reftss.riken.jp/datafiles/current/human/ 

PATH_FASTQ_ONT=$1
REF_GENOME="ReferenceGenome"
REPORT_PATH="Reports"
TRIM_PATH="Trimmed_BulkRNASeq"
FLAIR_PATH="/home/alejandra/flair-master"
NEW_PATH="FLAIRanalysis"
NEW_PATH_CORRECTED="FLAIRanalysis/Corrected"
PATH_ALIGN_ONT="Alignments"
SQANTI3_PATH="/home/alejandra/SQANTI3-5.2.1"

mkdir "Alignments"
mkdir "FLAIRanalysis"
mkdir "FLAIRanalysis/Corrected"

python $FLAIR_PATH/flair.py --version

echo "flair align starting the $(date)." # print an OK message
echo "flair align starting the $(date).." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

files=($PATH_FASTQ_ONT/*) # copy all the file names in the PATHFOLDER folder

for (( i=0; i<(${#files[@]}); i++)); # for loop with a counter to go inside each file
do 

	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%.*} # Get only the name of the file without the extension

	echo "Aligning ${file%.*} to the reference genome with minimap2." # print an OK message

	# flair align (minimap2):
	## flair align -g genome.fa -r <reads.fq>|<reads.fa> [options]
	### -a 	Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
	### -x STR	This option applies multiple options at the same time. It should be applied before other options:
	###	splice 	Long-read spliced alignment (-k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -b0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes). 
	###	In the splice mode, 1) long deletions are taken as introns and represented as the ‘N’ CIGAR operator; 2) long insertions are disabled; 
	###	 3) deletion and insertion gap costs are different during chaining; 4) the computation of the ‘ms’ tag ignores introns to demote hits to pseudogenes. 
	### --secondary=no	Do not output secondary alignments
	### --nvrna	Use native-RNA specific alignment parameters for minimap2 (-u f -k 14)
	#### -k14	-k INT 	Minimizer k-mer length [15]
	#### -uf 	How to find canonical splicing sites GT-AG - f: transcript strand; b: both strands; n: no attempt to match GT-AG [n]
	echo "flair align --genome $REF_GENOME/GRCh38.primary_assembly.genome.fa --reads ${files[$i]} --nvrna --threads 30 --output $PATH_ALIGN_ONT/${file%.*}"
	flair align --genome $REF_GENOME/GRCh38.primary_assembly.genome.fa --reads ${files[$i]} --nvrna --threads 30 --output $PATH_ALIGN_ONT/${file%.*}

	echo "${file%.*} aligned to the reference genome with minimap2" # print an OK message
	echo "${file%.*} aligned to the reference genome with minimap2" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
done;

echo "flair align finished the $(date)." # print an OK message
echo "flair align finished the $(date)." >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

echo "Prepare short reads for correction" # print an OK message
python PipelineIsoform3b_ExtractJuncitonAnnotationFromShortReads.py $PATH_ALIGN_ONT "Alignments_ShortReads" "$REF_GENOME/metadata/metadataJIArnaSeq.txt" "$REF_GENOME/metadata/PerStudySamples.txt" "$REF_GENOME/metadata/PerStudyBatch.txt" "FLAIRanalysis/monocytesOnlyNGSjunctionsFromSam.bed"

echo "Start flair correct and collapse at $(date). Parameters for flair correct include -nvrna (make the strand of a read consistent with the input annotation during correction) and --shortread (use bed format splice junctions from short-read sequencing for correction) and flair collapse with parameters --stringent (all supporting reads need to be full-length (80% coverage and spanning 25 bp of the first and last exons) and --filter nosubset (any isoforms that are a proper set of another isoform are removed)" # print an OK message
echo "Start flair correct and collapse at $(date). Parameters for flair correct include -nvrna (make the strand of a read consistent with the input annotation during correction) and --shortread (use bed format splice junctions from short-read sequencing for correction) and flair collapse with parameters --stringent (all supporting reads need to be full-length (80% coverage and spanning 25 bp of the first and last exons) and --filter nosubset (any isoforms that are a proper set of another isoform are removed)" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

files=($( ls $PATH_ALIGN_ONT/ | grep '.bed' )) # copy all the file names

for (( i=0; i<(${#files[@]}); i++)); # for loop with a counter to go inside each file
do 
	file_1=${files[$i]##*/} # Get only the name of the file (without folder)
	file=${file_1%.*} # Get only the name of the file without the extension

 	echo "Running flair correct for ${file%.*} sample"
	echo "Running flair correct for ${file%.*} sample" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

 	# flair correct:
 	## flair correct -q query.bed12 [-f annotation.gtf]|[-j introns.tab] -g genome.fa [options]
 	### --nvrna	Specify this flag to make the strand of a read consistent with the input annotation during correction.
 	### --shortread	Bed format splice junctions from short-read sequencing. You can generate these from SAM format files using the junctions_from_sam program that comes with Flair.
 	echo "python $FLAIR_PATH/flair.py correct --nvrna -t 30 -q $PATH_ALIGN_ONT/${files[$i]} -o $NEW_PATH_CORRECTED/${file%.*}_FLAIRcorrect --genome $REF_GENOME/GRCh38.primary_assembly.genome.fa --gtf $REF_GENOME/gencode.v45.primary_assembly.annotation.gtf --shortread $NEW_PATH/monocytesOnlyNGSjunctionsFromSam.bed --print_check"
 	python $FLAIR_PATH/flair.py correct --nvrna -t 30 -q $PATH_ALIGN_ONT/${files[$i]} -o $NEW_PATH_CORRECTED/${file%.*}_FLAIRcorrect --genome $REF_GENOME/GRCh38.primary_assembly.genome.fa --gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --shortread $NEW_PATH/monocytesOnlyNGSjunctionsFromSam.bed --print_check 
	
	echo "flair correct for ${file%.*} finished"
	echo "flair correct for ${file%.*} finished" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

done;

cat $NEW_PATH_CORRECTED/*_FLAIRcorrect_all_corrected.bed | sort -k1,1 -k2,2n -k3,3n > $NEW_PATH/flairAnalysis_catALLmonocytes_corrected.bed

cat $PATH_FASTQ_ONT/* > $PATH_FASTQ_ONT/AllFastq.fastq.gz

fastq=($PATH_FASTQ_ONT/*)

cut -f1,2,3,4,5,6 $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3.bed > $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3_ok.bed

cat $REF_GENOME/TSS/refTSS_v4.1_human_coordinate.hg38.bed $REF_GENOME/TSS/hg38.tc.decompose_smoothing_merged.ctssMaxCounts3_ok.bed  > $REF_GENOME/TSS/merged_refTSS-v4.1_CAGE_permissive_hg38.bed

echo "flair correct finished. Running flair collapse at $(date)"
echo "flair correct finished. Running flair collapse at $(date)" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# flair collapse (high-confidence isoforms):
## flair collapse -g genome.fa -q <query.bed> -r <reads.fq>/<reads.fa> [options]
### --stringent	Specify if all supporting reads need to be full-length (80% coverage and spanning 25 bp of the first and last exons).
### --check_splice	Enforce coverage of 4 out of 6 bp around each splice site and no insertions greater than 3 bp at the splice site.
### --filter comprehensive	default set + all subset isoforms
### --promoters	Promoter regions bed file to identify full-length reads.
echo "python $FLAIR_PATH/flair.py collapse -t 30 -o $NEW_PATH/flairAnalysis_monocytes_collapsed --promoters $REF_GENOME/TSS/merged_refTSS-v4.1_CAGE_permissive_hg38.bed --filter comprehensive --check_splice --stringent --genome $REF_GENOME/GRCh38.primary_assembly.genome.fa --gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --reads $PATH_FASTQ_ONT/AllFastq.fastq.gz --query $NEW_PATH/flairAnalysis_catALLmonocytes_corrected.bed --temp_dir "./" --generate_map"
python $FLAIR_PATH/flair.py collapse -t 30 -o $NEW_PATH/flairAnalysis_monocytes_collapsed --promoters $REF_GENOME/TSS/merged_refTSS-v4.1_CAGE_permissive_hg38.bed --filter comprehensive --check_splice --stringent --genome $REF_GENOME/GRCh38.primary_assembly.genome.fa --gtf $REF_GENOME/gencode.v46.primary_assembly.annotation.gtf --reads $PATH_FASTQ_ONT/AllFastq.fastq.gz --query $NEW_PATH/flairAnalysis_catALLmonocytes_corrected.bed --temp_dir "./" --generate_map

echo "flair collapse for all samples finished. Start first flair quantify for collapsed isoforms start at $(date)"
echo "flair collapse for all samples finished. Start first flair quantify for collapsed isoforms start at $(date)" >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# flair quantify (quantify the expression of each isoform using flair quantify):
## flair quantify -r reads_manifest.tsv -i isoforms.fa [options]#
### --isoform_bed   must be specified if --stringent or --check-splice is specified.
### --generate_map    Create read-to-isoform assignment files for each sample.  
echo "$FLAIR_PATH/flair.py quantify -t 30 -o $NEW_PATH/flairAnalysis_Monocytes_quantified --isoform_bed $NEW_PATH/flairAnalysis_monocytes_collapsed.isoforms.bed --temp_dir "./" --isoforms $NEW_PATH/flairAnalysis_monocytes_collapsed.isoforms.fa --reads_manifest $REF_GENOME/metadata/reads_manifestFile.txt"
python $FLAIR_PATH/flair.py quantify -t 30 -o $NEW_PATH/flairAnalysis_Monocytes_quantified --isoform_bed $NEW_PATH/flairAnalysis_monocytes_collapsed.isoforms.bed --temp_dir "./" --isoforms $NEW_PATH/flairAnalysis_monocytes_collapsed.isoforms.fa --reads_manifest $REF_GENOME/metadata/reads_manifestFile.txt

echo "First flair quantify finished." # print an OK message
echo "First flair quantify finished."  >> $REPORT_PATH/progress.txt # print the end message in the progress .txt

# prepare isoform-expression counts for long and short reads
sed 's/ids/id/g' $NEW_PATH/flairAnalysis_Monocytes_quantified.counts.tsv | awk -F"\t" '{gsub(/\_.*/,"",$1); print $0'} | sed 's/ /,/g' | sed 's/\t/,/g' > $NEW_PATH/flairAnalysis_Monocytes_quantified_format.tsv 

echo "flair pipeline (align, correct, collapse and quantify) finished at $(date)."  >> $REPORT_PATH/progress.txt # print the end message in the progress .txt
echo "flair pipeline (align, correct, collapse and quantify) finished at $(date)." 
