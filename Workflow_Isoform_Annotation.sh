#!/usr/bin/bash

NanoporePath=$1
ShortReadsPath=$2
ShortReadsPath_InVivo=$3
ShortReadsPath_Validation=$4
Thp1_RNASeq=$5
Thp1_RiboSeq=$6

# Pipeline 1: quality control of nanopore reads and multiqc quality report
bash PipelineIsoform1_fastqQC.sh $NanoporePath "fastqcReport_nanopore"

# Pipeline 2: quality control of short-reads and multiqc quality report
bash PipelineIsoform1_fastqQC.sh $ShortReadsPath "fastqcReport_shortreads"

# Pipeline 3: compute read length for Nanopore reads
bash PipelineIsoform1b_ReadLentgh.sh $NanoporePath "read_length_nanopore"

# Pipeline 4: single-end short-reads process: alignment using STAR
bash PipelineIsoform2_ShortReads_TrimandAlign.sh $ShortReadsPath

# Pipeline 5: paired-end short-reads process: alignment using STAR (GSE87290)
bash PipelineIsoform2b_ShortReads_TrimandAlign.sh $ShortReadsPath_InVivo

# Pipeline 6: novel isoform detection using flair
bash PipelineIsoform3_FLAIR.sh $NanoporePath

# Pipeline 7: quality control of the new isoforms, filter and rescue with SQANTI3
bash PipelineIsoform4_SQANTI3.sh

# Pipeline 8: comparison with other annotations (TRAILS and FLIBase)
bash PipelineIsoform4b_ComparisonWithOtherAnnotations.sh

# Pipeline 9: FLAIR expression quantification and differentiall splicing analysis uisng SUPPA
# longest-ORF and productivity prediciton
bash PipelineIsoform5_DifferentialSplicingEvents.sh

# Pipeline 10: long-read quantification using salmon
bash PipelineIsoform6_LongReadQuantification_salmon.sh $NanoporePath

# Pipeline 11: single-end short-reads re-alignment using salmon
bash PipelineIsoform7_ShortReads_RE-AlignmentandQuantification.sh $ShortReadsPath

# Pipeline 12: paired-end short-reads (GSE103501) alignment and re-alignment using salmon
bash PipelineIsoform7b_ShortReads_TrimandAlignandRE-Align_Validation.sh $ShortReadsPath_Validation

# Pipeline 13: IsoformSwitchingAnalyseR analysis expansion tools
bash PipelineIsoform8_IsoformSwitchingAnalysis.sh

# Pipeline 14: Ribo-Seq and matching RNA-Seq analysis for Thp1 (GSE208041)
bash PipelineIsoform9_RiboSeqData.sh $Thp1_RiboSeq $Thp1_RNASeq

# Pipeline 15: run MaxQuant for PXD004352
bash PipelineIsoform10_MaxQuant.sh

# Pipeline 16: long-read overlap with annotated TSS
bash PipelineIsoform11_TSSOverlap.sh

# Sashimi plot for TDP2
python ggsashimi.py -b ReferenceGenome/BAMFiles.txt -c chr6:24649979-24654530 -g FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf -S plus -M 0 -s ANTISENSE -o TDP2 --color-factor 1 --palette Palette.txt --alpha 1
