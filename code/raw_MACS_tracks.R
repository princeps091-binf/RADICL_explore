library(tidyverse)
library(vroom)
library(arrow)
MACS_track_tbl<- read_delim_arrow("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_DNA_track_chr19.bdg",delim="\t",col_names = F, as_data_frame = F) %>% 
#MACS_track_tbl<- read_delim_arrow("./../data/ATAC/GM12878_ATAC_ENCFF962FMH_sorted_chr19.bdg",delim="\t",col_names = F, as_data_frame = F) %>% 
  filter(f3 > 0) %>% 
#  mutate(log.f3=log10(f3)) %>% 
  collect()

MACS_track_tbl %>% 
  mutate(flog=log10(f3 + 1)) %>% 
#  filter(f1 > 1e7 & f2 > 1e7 & f1 < 2e7 & f2 < 2e7) %>% 
  ggplot(.,aes(x=f1,xend=f2,y=flog,yend=flog))+
  geom_segment(linewidth=5,color="red")+
  theme_classic()

MACS_track_tbl %>% 
  collect() %>%   
  summarise(d=f2-f1) %>% 
  ggplot(.,aes(d)) +
  geom_density()+
  scale_x_log10()

