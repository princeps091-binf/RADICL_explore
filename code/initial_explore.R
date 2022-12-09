library(tidyverse)
library(vroom)
library(arrow)
library(dplyr)
library(GenomicRanges)
sign_inter_tbl<-vroom("./data/CHICANE_significant_interaction_3cells.tsv",col_names = T)
inter_anno_tbl<-vroom("./data/interaction_ID.annotation.tsv",col_names = T,col_select = c("interaction_ID","DNA_bin"))
MACS_peak_tbl<-vroom("./data/NA_peaks.narrowPeak",col_names = F)
iPSC_inter_tbl<-sign_inter_tbl %>% 
  filter(sig_CHICANE_iPSC == "yes")
iPSC_inter_tbl<-iPSC_inter_tbl %>% 
  left_join(.,inter_anno_tbl)
iPSC_inter_tbl<-iPSC_inter_tbl %>% 
  mutate(chr=str_split_fixed(DNA_bin,"_",3)[,1],
         start=as.integer(str_split_fixed(DNA_bin,"_",3)[,2]),
         end=as.integer(str_split_fixed(DNA_bin,"_",3)[,3]))
MACS_Grange<-GRanges(seqnames=MACS_peak_tbl$X1,
                    ranges = IRanges(start=MACS_peak_tbl$X2,
                                     end=MACS_peak_tbl$X3)
)

sig_inter_GRange<-GRanges(seqnames=iPSC_inter_tbl$chr,
                          ranges = IRanges(start=iPSC_inter_tbl$start,
                                           end=iPSC_inter_tbl$end)
)

length(unique(subjectHits(findOverlaps(sig_inter_GRange,MACS_Grange))))
