library(vroom)
library(dplyr)
library(ggplot2)

MACS_track_tbl<- vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/RADICL_DNA_smooth_chr19.bdg",delim="\t",col_names = F)

MACS_peak_tbl<-vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/RADICL_DNA_chr19_peaks.bed",delim="\t",col_names=F,skip = 1)

max_coord<-as.integer(max(MACS_track_tbl$X3) +(1e6 - max(MACS_track_tbl$X3) %% 1e6))
min_coord<-as.integer(min(MACS_track_tbl$X3) - (min(MACS_track_tbl$X3) %% 1e6))
peak_band_position<-ceiling(log10(max(MACS_track_tbl$X4)))
MACS_track_tbl %>% 
  mutate(flog=log10(X4 + 1)) %>% 
  #  filter(f1 > 1.5e7 & f2 > 1.5e7 & f1 < 1.75e7 & f2 < 1.75e7) %>% 
  ggplot(.,aes(x=X2,xend=X3,y=flog,yend=flog))+
  geom_segment(linewidth=5,color="red")+
  geom_segment(data=MACS_peak_tbl,
               aes(x=X2,xend=X3,y=peak_band_position,yend=peak_band_position),
               color="black",linewidth=5)+
  scale_x_continuous(breaks = seq(min_coord, max_coord, by=5e6), 
                     limits=c(min_coord,max_coord), 
                     labels = paste0(seq(min_coord, max_coord, by=5e6)/1e6,"Mb"))+
  xlab("chromosome coordinate")+
  ylab("ATAC")+
  theme_classic()
#-------------------------------------------------------------------------
#Multi-track
RNA_track_tbl<- vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_RNA_smooth_chr19.bdg",delim="\t",col_names = F)
DNA_track_tbl<- vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_DNA_smooth_chr19.bdg",delim="\t",col_names = F)
#ATAC_track_tbl<-vroom("./../data/ATAC/FANTOM6/ff_iPSC_ATAC_sorted_chr19.bdg",delim="\t",col_names = F)
all_track_val<-c(RNA_track_tbl$X3,DNA_track_tbl$X3)
max_coord<-as.integer(max(all_track_val) +(1e6 - max(all_track_val) %% 1e6))
min_coord<-as.integer(min(all_track_val) - (min(all_track_val) %% 1e6))
RNA_track_tbl %>% 
  mutate(type="RNA") %>% 
  bind_rows(., DNA_track_tbl %>% 
              mutate(type="DNA")) %>% 
#  bind_rows(., ATAC_track_tbl %>% 
#              mutate(type="ATAC")) %>% 
  mutate(flog=log10(X4 + 1)) %>% 
  #  filter(f1 > 1.5e7 & f2 > 1.5e7 & f1 < 1.75e7 & f2 < 1.75e7) %>% 
  ggplot(.,aes(x=X2,xend=X3,y=flog,yend=flog,color=type))+
  geom_segment(linewidth=5)+
  scale_x_continuous(breaks = seq(min_coord, max_coord, by=5e6), 
                     limits=c(min_coord,max_coord), 
                     labels = paste0(seq(min_coord, max_coord, by=5e6)/1e6,"Mb"))+
  xlab("chromosome coordinate")+
  ylab("RADICL")+
  scale_color_brewer(type="qual",palette="Set1")+
  facet_grid(type~.,scales = "free_y")+
  theme_classic()
#-------------------------------------------------------------------------
#Single-track
RADICL_track_tbl<- vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/track_viz/RADICL_DNA_singleton_bg.bdg",delim="\t",col_names = F)
DNA_track<-vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/track_viz/RADICL_DNA_smooth.bdg",delim="\t",col_names = F)
all_track_val<-c(RADICL_track_tbl$X3,DNA_track$X3)
max_coord<-as.integer(max(all_track_val) +(1e6 - max(all_track_val) %% 1e6))
min_coord<-as.integer(min(all_track_val) - (min(all_track_val) %% 1e6))
RADICL_track_tbl %>% 
  mutate(type="ctrl") %>% 
  mutate(flog=log10(X4 + 1)) %>% 
  select(-c(X4)) %>% 
  bind_rows(., DNA_track %>% 
              mutate(type="full")  %>% 
              mutate(flog=log10(X4 + 1)) %>% 
              select(-c(X4))) %>% 
  ggplot(.,aes(x=X2,xend=X3,y=flog,yend=flog,color=type))+
  geom_segment(linewidth=5)+
  scale_x_continuous(breaks = seq(min_coord, max_coord, by=5e6), 
                     limits=c(min_coord,max_coord), 
                     labels = paste0(seq(min_coord, max_coord, by=5e6)/1e6,"Mb"))+
  xlab("chromosome coordinate")+
  ylab("RADICL")+
  scale_color_brewer(type="qual",palette="Set1")+
  theme_classic()
#Single-track
qval_track<-vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/RADICL_DNA_qval_chr19.bdg",delim="\t",col_names = F)
all_track_val<-c(qval_track$X3)
max_coord<-as.integer(max(all_track_val) +(1e6 - max(all_track_val) %% 1e6))
min_coord<-as.integer(min(all_track_val) - (min(all_track_val) %% 1e6))

qval_track %>% 
  ggplot(.,aes(x=X2,xend=X3,y=X4,yend=X4))+
  geom_segment(linewidth=5)+
  scale_x_continuous(breaks = seq(min_coord, max_coord, by=5e6), 
                     limits=c(min_coord,max_coord), 
                     labels = paste0(seq(min_coord, max_coord, by=5e6)/1e6,"Mb"))+
  xlab("chromosome coordinate")+
  ylab("RADICL")+
  scale_color_brewer(type="qual",palette="Set1")+
  theme_classic()
