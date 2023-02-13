library(vroom)
library(GenomicRanges)
library(tidyverse)
library(furrr)
RADICL_intra_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/raw/RADICL_iPSC_intra.bed",
                        col_names = F,delim = "\t")

gene_anno<-vroom("~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed",
                 col_names = F,delim = "\t",col_types = list("c",'i','i','c','d','c','i','i','c','d','c','c'))




DNA_read_GRanges<-GRanges(seqnames=RADICL_intra_tbl$X7,
                          ranges = IRanges(start=RADICL_intra_tbl$X8,
                                           end=RADICL_intra_tbl$X9
                          ),
                          strand=RADICL_intra_tbl$X12
)
mcols(DNA_read_GRanges)<-tibble(ID=RADICL_intra_tbl$X10)

RNA_read_GRanges<-GRanges(seqnames=RADICL_intra_tbl$X1,
                          ranges = IRanges(start=RADICL_intra_tbl$X2,
                                           end=RADICL_intra_tbl$X3
                          ),
                          strand=RADICL_intra_tbl$X6
)
mcols(RNA_read_GRanges)<-tibble(ID=RADICL_intra_tbl$X4)



anno_GRanges<-GRanges(seqnames=gene_anno$X1,
                      ranges = IRanges(start=gene_anno$X2,
                                       end=gene_anno$X3
                      ),
                      strand=gene_anno$X6)
mcols(anno_GRanges)<-tibble(name=gene_anno$X4)
# Get read per transcript
transcript_read_inter_tbl<-as_tibble(findOverlaps(anno_GRanges,RNA_read_GRanges)) %>% 
  dplyr::mutate(transcript=mcols(anno_GRanges)$name[queryHits],
         read.ID=mcols(RNA_read_GRanges)$ID[subjectHits]) %>% 
  dplyr::select(transcript,read.ID)
transcript_DNA_read_GRanges<-DNA_read_GRanges[mcols(DNA_read_GRanges)$ID %in% transcript_read_inter_tbl$read.ID]
#process by chromosome
transcript_read_inter_tbl<-transcript_read_inter_tbl %>% 
  left_join(.,gene_anno %>% 
              dplyr::select(X1,X4),by=c('transcript'='X4'))
chr_set<-unique(transcript_read_inter_tbl$X1)
map(chr_set,function(chromo){
  message(chromo)
  chromo_transcript_read_inter_tbl<-transcript_read_inter_tbl %>% 
    filter(X1==chromo)
  plan(multisession, workers=4)
  dist_l<-future_map(unique(chromo_transcript_read_inter_tbl$transcript),function(i){
    transcript_GRange<-anno_GRanges[mcols(anno_GRanges)$name == i]
    tmp_read_ID<-chromo_transcript_read_inter_tbl %>% 
      filter(transcript==i) %>% 
      distinct(read.ID) %>% 
      unlist
    tmp_DNA_GRange<-transcript_DNA_read_GRanges[mcols(transcript_DNA_read_GRanges)$ID %in% tmp_read_ID]
    strand(tmp_DNA_GRange)<-'*'
    return(IRanges::distance(transcript_GRange,tmp_DNA_GRange))
  })
  plan(sequential)
  tibble(d=unlist(dist_l)) %>% 
    ggplot(.,aes(d+1))+
    geom_density()+
    scale_x_log10()
  return(tibble(chr=chromo,transcript=unique(chromo_transcript_read_inter_tbl$transcript),dna.dist=dist_l))
  
})
# Examine distance distribution of DNA-side reads relative to considered source transcript
