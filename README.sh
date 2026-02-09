Pipeline for detecting and quantifying novel isoform expression in 3h LPS activated and non-activated human monocytes using Oxford Nanopore and short-read RNA sequencing, 
followed by assessment of the functional consequences of splicing changes.

PipelineIsoform1_fastqQC.sh: quality control of sequenced reads and multiqc quality report
PipelineIsoform1b_ReadLentgh.sh: ONT read length description
PipelineIsoform2_ShortReads_TrimandAlign.sh: single-end short-reads trimming and alignment using STAR
PipelineIsoform2b_ShortReads_TrimandAlign.sh: paired-end short-reads trimming and alignment using STAR
PipelineIsoform3_FLAIR.sh: novel isoform detection using flair
PipelineIsoform3b_ExtractJuncitonAnnotationFromShortReads.py: short-read junction extraction for flair annotation correction
PipelineIsoform4_SQANTI3.sh: quality control of the new isoforms, filter and rescue with SQANTI3
PipelineIsoform4b_ComparisonWithOtherAnnotations.sh: comparison of the novel annotation with other annotations (TRAILS and FLIBase)
PipelineIsoform5_DifferentialSplicingEvents.sh: FLAIR expression quantification, differentiall splicing analysis uisng SUPPA and functionality prediction (longest-ORF and productivity prediciton)
PipelineIsoform5b_DifferentialSplicingEvents.R
PipelineIsoform5c_FuntionalityAnalysis_MaxORF.py: longest-ORF and productivity prediciton
PipelineIsoform6_LongReadQuantification_salmon.sh: long-read quantification using salmon
PipelineIsoform7_ShortReads_RE-AlignmentandQuantification.sh: single-end short-reads re-alignment using salmon
PipelineIsoform7b_ShortReads_TrimandAlignandRE-Align_Validation.sh: paired-end short-reads alignment and re-alignment using salmon
PipelineIsoform8_IsoformSwitchingAnalysis.sh: IsoformSwitchingAnalyseR analysis expansion tools
PipelineIsoform9_RiboSeqData.sh: Ribo-Seq and matching RNA-Seq analysis for Thp1 analysis using uniquely aligned reads
PipelineIsoform10_MaxQuant.sh: mass spectrometry analysis using MaxQuant
PipelineIsoform11_TSSOverlap.sh: long-read overlap with annotated TSS
Figure1_Analysis.qmd: gene-level short-read RNA-Seq analysis (Figure 1)
Figure2_Analysis.qmd: long-read RNA-Seq annotation and characterization of novel isoforms (Figure 2).
Figure3_Analysis.qmd: long-read RNA-Seq quantification, assessment of functional changes, and validation using short-read data (Figures 3, 4 and 5)
Figure4_Analysis.qmd: long-read RNA-Seq analysis of functional changes, integration of translation data (Ribo-Seq and mass spectrometry), and illustrative examples (Figures 4 and 5).
Workflow_Splicing_in_Monocytes.sh: master pipeline that sequentially executes all other analysis scripts.
