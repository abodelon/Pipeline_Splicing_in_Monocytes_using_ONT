library(data.table)
library(dplyr)

Monocyte_Isoforms_Annot=fread("FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf", header=F)
Monocyte_Isoforms_Annot=unique(data.table(Monocyte_Isoforms_Annot %>% subset(V3=="exon") %>% 
                                     mutate(transcript_id=sub('".*','',sub('.*transcript_id "','',V9))) %>%
                                     mutate(gene_id=sub('".*','',sub('.*gene_id "','',V9))) %>%
                                     select(c(transcript_id,gene_id))))

gencodeprimary_assemblyannotation=fread("ReferenceGenome/gencode.v46.primary_assembly.annotation.gtf", header=F)
gencodeprimary_assemblyannotation=gencodeprimary_assemblyannotation %>% 
  filter(V3=="gene") %>% 
  mutate(gene_id=sub('".*','',sub('.*gene_id "','',V9))) %>%
  mutate(gene_name=sub('".*','',sub('.*gene_name "','',V9))) %>%
  select(c(gene_id,gene_name))

Monocyte_Isoforms_Annot=merge(Monocyte_Isoforms_Annot, gencodeprimary_assemblyannotation, all.x=T)
Monocyte_Isoforms_Annot$gene_id=paste(Monocyte_Isoforms_Annot$gene_id, Monocyte_Isoforms_Annot$gene_name, sep="-")

Monocyte_Isoforms=fread("FinalIsoformAnnotation/Monocyte_Isoforms.final_oneline.fasta", header=F)

Header=grep(">",Monocyte_Isoforms$V1)
Count=1

FinalFasta=data.table()
for(HeaderPos in Header){
  if(Count != HeaderPos){
    print("Not one line fasta")
  }
  FinalFasta=rbind(FinalFasta, 
                   data.table(paste(Monocyte_Isoforms$V1[HeaderPos],
                                    Monocyte_Isoforms_Annot[Monocyte_Isoforms_Annot$transcript_id ==
                                                              sub(">","",Monocyte_Isoforms$V1[HeaderPos]),]$gene_id, sep="_")))
  FinalFasta=rbind(FinalFasta,
                   data.table(Monocyte_Isoforms$V1[HeaderPos+1]))
  Count=Count+2
}

write.table(FinalFasta, "FinalIsoformAnnotation/Monocyte_Isoforms.final_DEG.fasta", quote = F, row.names = F, col.names = F, sep="")

Monocyte_Isoforms=fread("FinalIsoformAnnotation/Monocyte_Isoforms.final.bed", header=F)

for(i in 1:nrow(Monocyte_Isoforms)){
  Monocyte_Isoforms[i,]$V4=paste(Monocyte_Isoforms[i,]$V4, 
                                 Monocyte_Isoforms_Annot[Monocyte_Isoforms_Annot$transcript_id ==
                                                           Monocyte_Isoforms[i,]$V4,]$gene_id, sep="_")
}

write.table(Monocyte_Isoforms, "FinalIsoformAnnotation/Monocyte_Isoforms.final_DEG.bed", quote = F, row.names = F, col.names = F, sep="\t")
