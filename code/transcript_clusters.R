# Transcript clusters
library(tidyverse)
library(GenomicRanges)
library(vroom)
library(igraph)
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

edge_tbl<-as_tibble(findOverlaps(transctipt_GRanges,transctipt_GRanges))
  
transcript_g<-graph_from_data_frame(edge_tbl)
components(transcript_g)
