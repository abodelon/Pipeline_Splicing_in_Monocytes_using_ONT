#!/usr/bin/env bash

###################################
# Comparison with other databases #
###################################

set -e # exit if there is an error
REF_GENOME="ReferenceGenome/"
REPORT_PATH="Reports"

mkdir "$REF_GENOME/TRAILS" # final annotation folder

# compare our annotaiton with the annotation of TRAILS
gffcompare -o $REF_GENOME/TRAILS/Iso_Trails_Validated_NotStrict -r $REF_GENOME/TRAILS_validated_by_PBMC15.gtf FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf
awk '$3=="transcript"' $REF_GENOME/TRAILS/Iso_Trails_Validated_NotStrict.annotated.gtf | cut -f9 > $REF_GENOME/TRAILS/Summary.txt


cut -f2,3,8 $REF_GENOME/FLIBase/transcript_info_v4.txt > $REF_GENOME/transcript_info_v4_col.txt
sed -i '1d' $REF_GENOME/FLIBase/transcript_info_v4_col.txt

awk -F"\t" '
{
    transcript_id = $1;
    gene_id = $2;
    split($3, a, ":");
    chrom = a[1];
    strand = a[length(a)];
    split(a[2], exons, ";");

    # Initialize min_start and max_end with first exon coords
    split(exons[1], coords, "-");
    min_start = coords[1];
    max_end = coords[2];

    # Loop from second exon onward
    for(i=2; i<=length(exons); i++) {
        split(exons[i], coords, "-");
        start = coords[1];
        end = coords[2];
        if (start < min_start) min_start = start;
        if (end > max_end) max_end = end;
    }

    # Print transcript line after min/max properly set
    printf "%s\tcustom\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", \
           chrom, min_start, max_end, strand, gene_id, transcript_id;

    # Print exon lines
    for(i=1; i<=length(exons); i++) {
        split(exons[i], coords, "-");
        start = coords[1];
        end = coords[2];
        printf "%s\tcustom\texon\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\";\n", \
               chrom, start, end, strand, gene_id, transcript_id, i
    }
}
' $REF_GENOME/transcript_info_v4_col.txt > $REF_GENOME/transcript_info_v4.gtf

# compare our annotaiton with the annotation of FLIBase
gffcompare -o $REF_GENOME/FLIBase/transcript_info_v4_NotStrict -r $REF_GENOME/FLIBase/transcript_info_v4.gtf FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf
awk '$3=="transcript"' $REF_GENOME/FLIBase/transcript_info_v4_NotStrict.annotated.gtf | cut -f9 > $REF_GENOME/FLIBase/Summary.txt
