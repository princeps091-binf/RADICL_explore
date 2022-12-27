library(tidyverse)
library(GenomicRanges)
library(vroom)
library(furrr)
chr_dat<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_iPSC_chr19.bed",
      col_names = F)
peak_dat<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/naive_broad_peaks.broadPeak",
                col_names = F)

DNA_GRange<-GRanges(seqnames=chr_dat$X7,
                                   ranges = IRanges(start=chr_dat$X8,
                                                    end=chr_dat$X9
                                   ))
mcols(DNA_GRange)<-tibble(rna.read=paste(chr_dat$X1,chr_dat$X2,chr_dat$X3,sep="_"))
peak_GRange<-GRanges(seqnames=peak_dat$X1,
                     ranges = IRanges(start=peak_dat$X2,
                                      end=peak_dat$X3
                     ))
mcols(peak_GRange)<-tibble(peak.ID=paste(peak_dat$X1,peak_dat$X2,peak_dat$X3,sep="_"))

inter_tbl<-tibble(rna.read=DNA_GRange$rna.read[queryHits(findOverlaps(DNA_GRange,peak_GRange))],
       peak=peak_GRange$peak.ID[subjectHits(findOverlaps(DNA_GRange,peak_GRange))])
plan(multisession,workers=3)
test<-inter_tbl %>% 
  group_by(peak) %>% 
  summarise(rna.set=list(do.call(c,future_map(rna.read,function(x){
    tmp_coord<-str_split_fixed(grep("_[A-z]",x,value = T,invert = T),"_",3)
    
    return(IRanges::reduce(GRanges(seqnames=tmp_coord[,1],
                            ranges = IRanges(start=as.integer(tmp_coord[,2]),
                                             end=as.integer(tmp_coord[,3])
                            ))))
  }))))
plan(sequential)
test %>% 
  mutate(nrna=map_int(rna.set,function(x){
    length(x)
  }))
inter_tbl %>% 
  group_by(peak) %>% 
  count() %>% 
  arrange(desc(n))
