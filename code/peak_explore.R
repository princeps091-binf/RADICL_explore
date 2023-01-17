# Peak RefSeq gene content
library(tidyverse)
library(GenomicRanges)
library(vroom)
peak_file<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/result/DNA/MACS/peaks/RADICL_DNA_tot_peaks.bed"
read_file<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/processed/DNA/cluster/RADICL_iPSC_DNA_chr1.bed_cluster.txt"
# load the peak file
peak_tbl<-vroom(peak_file,
                col_names = F,
                delim="\t")
# load Read files
DNA_read_tbl<-vroom(read_file,
                col_names = F,
                delim="\t")
#----------------------------------------------------
peak_tbl %>% 
  mutate(width=abs(X3-X2)) %>% 
  ggplot(.,aes(width))+
  geom_density()+
  scale_x_log10()
peak_tbl %>% 
  mutate(width=abs(X3-X2)) %>% 
  arrange(desc(width))
cluster_summary<-DNA_read_tbl %>% 
  group_by(X6) %>% 
  count
read_GRanges<-GRanges(seqnames=DNA_read_tbl$X1,
                      ranges = IRanges(start=DNA_read_tbl$X2,
                                       end=DNA_read_tbl$X3
                      ))
mcols(read_GRanges)<-tibble(ID=DNA_read_tbl$X4,cluster.ID=DNA_read_tbl$X6) %>% 
  left_join(.,cluster_summary,by=c('cluster.ID'='X6'))

peak_GRanges<-GRanges(seqnames=peak_tbl$X1,
                      ranges = IRanges(start=peak_tbl$X2,
                                       end=peak_tbl$X3
                      ))
length(findOverlaps(peak_GRanges,read_GRanges))/length(read_GRanges)
as_tibble(findOverlaps(peak_GRanges,read_GRanges)) %>% 
  mutate(n=mcols(read_GRanges)$n[subjectHits]) %>% 
  mutate(set='peak') %>% 
  dplyr::select(n,set) %>% 
  bind_rows(.,as_tibble(mcols(read_GRanges)) %>% 
              mutate(set='all') %>% 
              dplyr::select(n,set)
  ) %>% 
  ggplot(.,aes(n,color=set))+
  geom_density()+
  scale_x_log10()+
  facet_grid(set~.,scales='free_y')

extreme_peak<-reduce(read_GRanges[subjectHits(findOverlaps(peak_GRanges[which(width(peak_GRanges) > 9e4)],read_GRanges))])
tibble(start=start(extreme_peak),end=end(extreme_peak),rc=countOverlaps(extreme_peak,read_GRanges)) %>% 
  ggplot(.,aes(x=start,xend=end,y=rc,yend=rc,linewidth=log10(rc)))+
  geom_segment()
peak_gaps<-GenomicRanges::setdiff(range(extreme_peak),extreme_peak)
plan(multisession,workers=3)
peak_profile<-future_map_dfr(seq_along(extreme_peak),function(i){
  tibble(pos=seq(start(extreme_peak)[i],end(extreme_peak)[i]),rc=countOverlaps(extreme_peak[i],read_GRanges))

})
plan(sequential)

tibble(start=start(extreme_peak),end=end(extreme_peak),rc=countOverlaps(extreme_peak,read_GRanges))

peak_profile %>% 
  ggplot(.,aes(x=pos,y=rc))+
  geom_line()
