library(tidyverse)
library(vroom)
library(GenomicRanges)

GRange_from_UCSC_fn<-function(anno_file){
  Refseq_anno<-vroom(anno_file,
                     col_names = F,delim = "\t")
  
  Refseq_GRanges<-GRanges(seqnames=Refseq_anno$X1,
                          ranges = IRanges(start=Refseq_anno$X2,
                                           end=Refseq_anno$X3
                          ))
  return(Refseq_GRanges)
}
refseq_file<-"~/Documents/FANTOM6/data/annotation/hg38_Refseq_all.txt"
mirna_file<-"~/Documents/FANTOM6/data/annotation/sno_mirna_hg38.tsv"
lincrna_file<-"~/Documents/FANTOM6/data/annotation/lincRNA_tucp_hg38.tsv"
fantom_cat_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
Refseq_GRanges<-GRange_from_UCSC_fn(refseq_file)
mirna_GRanges<-GRange_from_UCSC_fn(mirna_file)
lincrna_GRanges<-GRange_from_UCSC_fn(lincrna_file)
fantom_cat_GRanges<-GRange_from_UCSC_fn(fantom_cat_file)

length(unique(queryHits(findOverlaps(Refseq_GRanges,fantom_cat_GRanges))))/length(Refseq_GRanges)
