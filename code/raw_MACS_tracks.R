library(tidyverse)
library(vroom)
#library(arrow)
#MACS_track_tbl<- read_delim_arrow("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_DNA_track_chr19.bdg",delim="\t",col_names = F, as_data_frame = F) %>% 
MACS_track_tbl<- vroom("./../data/ATAC/FANTOM6/ff_iPSC_ATAC_sorted_chr19.bdg",delim="\t",col_names = F) %>% 
#MACS_track_tbl<- read_delim_arrow("./../data/ATAC/FANTOM6/ff_iPSC_ATAC_sorted_chr19.bdg",delim="\t",col_names = F, as_data_frame = F) %>% 
  #  filter(f3 > 0) %>% 
#  mutate(log.f3=log10(f3)) %>% 
  collect()
MACS_track_tbl<- vroom("./../data/ATAC/FANTOM6/ff_iPSC_ATAC_sorted_chr19.bdg",delim="\t",col_names = F)
  
MACS_peak_tbl<-vroom("./../data/ATAC/FANTOM6/chr19_ATAC_peaks.broadPeak",col_names=F)
max_coord<-as.integer(max(MACS_track_tbl$X3) +(1e6 - max(MACS_track_tbl$X3) %% 1e6))
peak_band_position<-ceiling(log10(max(MACS_track_tbl$X4)))
MACS_track_tbl %>% 
  mutate(flog=log10(X4 + 1)) %>% 
#  filter(f1 > 1.5e7 & f2 > 1.5e7 & f1 < 1.75e7 & f2 < 1.75e7) %>% 
  ggplot(.,aes(x=X2,xend=X3,y=flog,yend=flog))+
  geom_segment(linewidth=5,color="red")+
  geom_segment(data=MACS_peak_tbl,aes(x=X2,xend=X3,y=peak_band_position,yend=peak_band_position),color="black",linewidth=20)+
  scale_x_continuous(breaks = seq(0, max_coord, by=5e6), 
                     limits=c(0,max_coord), 
                     labels = paste0(seq(0, max_coord, by=5e6)/1e6,"Mb"))+
  xlab("chromosome coordinate")+
  ylab("ATAC")+
  theme_classic()

ggsave("./../ATAC_chr19.png")
MACS_track_tbl %>% 
  collect() %>%   
  summarise(d=f2-f1) %>% 
  ggplot(.,aes(d)) +
  geom_density()+
  scale_x_log10()

