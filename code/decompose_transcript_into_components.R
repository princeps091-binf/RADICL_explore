# Exon composition
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(furrr)
#-------------------------------------------------------------
transcript_annotation_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

transcript_tbl<-vroom(transcript_annotation_file,
                      col_types = list("c",'i','i','c','d','c','i','i','c','d','c','c'),
                      col_names = F,
                      delim="\t")

transctipt_GRanges<-GRanges(seqnames=transcript_tbl$X1,
                            ranges = IRanges(start=transcript_tbl$X2,
                                             end=transcript_tbl$X3
                            ),
                            strand=transcript_tbl$X6
)
mcols(transctipt_GRanges)<-tibble(ID=transcript_tbl$X4)

exon_tbl<-transcript_tbl %>% 
  dplyr::select(X4,X1,X2,X3,X6,X10,X11,X12) %>% 
  mutate(RelexonStarts=map(X12,function(x){
    as.integer(unlist(str_split(x,',')))
  }),
  exonLength=map(X11,function(x){
    as.integer(unlist(str_split(x,',')))
  })) %>% 
  unnest(cols=c(RelexonStarts,exonLength)) %>% 
  mutate(exonStart=X2+RelexonStarts) %>% 
  mutate(exonEnd=exonStart+exonLength) %>% 
  dplyr::select(X4,X1,exonStart,exonEnd,X6,X2,X3)

exon_GRanges<-GRanges(seqnames=exon_tbl$X1,
                      IRanges(
                        start=exon_tbl$exonStart,
                        end=exon_tbl$exonEnd,
                      ),
                      strand=exon_tbl$X6)

mcols(exon_GRanges)<-tibble(name=exon_tbl$X4)

reduced_exon_tbl<-exon_tbl %>% 
  distinct(X1,exonStart,exonEnd,X6)

reduced_exon_GRanges<-GRanges(seqnames=reduced_exon_tbl$X1,
                              IRanges(
                                start=reduced_exon_tbl$exonStart,
                                end=reduced_exon_tbl$exonEnd,
                              ),
                              strand=reduced_exon_tbl$X6)

red_exon_to_transcript_tbl<-as_tibble(findOverlaps(reduced_exon_GRanges,transctipt_GRanges)) %>% 
  mutate(tr.ID=mcols(transctipt_GRanges)$ID[subjectHits]) %>% 
  group_by(queryHits) %>% 
  summarise(tr.ID=list(unique(tr.ID)))
mcols(reduced_exon_GRanges)<-tibble(tr.IDs=red_exon_to_transcript_tbl$tr.ID)

red_exon_to_transcript_tbl %>% 
  mutate(n=map_int(tr.ID,function(x)length(x))) %>% 
  ggplot(.,aes(n))+
  geom_density()+
  scale_x_log10()

#memory limit impose per chromosome processing
plan(multisession,workers=4)
intron_tbl<-future_map_dfr(unique(transctipt_GRanges$ID),function(i){
  tmp_gene_GRange<-transctipt_GRanges[transctipt_GRanges@elementMetadata$ID == i]
  tmp_exon_GRanges<-reduced_exon_GRanges[map_lgl(mcols(reduced_exon_GRanges)$tr.IDs,function(x){
    i %in% x
  })]
  tmp_intron<-setdiff(tmp_gene_GRange,tmp_exon_GRanges)
  tibble(as.data.frame(tmp_intron)) %>% 
    mutate(name=i)
  
})
plan(sequential)
