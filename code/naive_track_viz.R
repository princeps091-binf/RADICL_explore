library(vroom)
library(dplyr)
library(ggplot2)

MACS_t_track_tbl<- vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/track_viz/RADICL_DNA_smooth.bdg",delim="\t",col_names = F)

MACS_peak_tbl<-vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/naive/naive_peaks.narrowPeak",delim="\t",col_names=F,skip = 1)

chr_MACS_tbl<-MACS_peak_tbl %>% 
  filter(X1 == unique(MACS_t_track_tbl$X1))
all_track_val<-c(MACS_t_track_tbl$X3)
max_coord<-as.integer(max(all_track_val) +(1e6 - max(all_track_val) %% 1e6))
min_coord<-as.integer(min(all_track_val) - (min(all_track_val) %% 1e6))

peak_band_position<-ceiling(log10(max(c(MACS_t_track_tbl$X4))))

MACS_t_track_tbl %>% 
  mutate(type="RADICL") %>% 
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
  geom_point(data=chr_MACS_tbl,
               aes(x=X2,xend=X3,y=peak_band_position,yend=peak_band_position),
               color="black",linewidth=5)+
  theme_classic()

chr_MACS_tbl %>% 
  ggplot(.,aes(x=X2,xend=X3,y=1,yend=1))+
  geom_segment(color="red",linewidth=5)+
  theme_classic()

