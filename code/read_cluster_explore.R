library(vroom)
library(tidyverse)
intra_chr_radicl<-vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_iPSC_intra_chr1.bed",col_names=F)

intra_chr_radicl %>% 
  mutate(d=abs(X2-X8)) %>%
  ggplot(.,aes(d))+
  geom_density()
#-----------------------------------------------------
chr_dna_read_cluster<-vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_iPSC_chr19_DNA_cluster.txt",delim = "\t",col_names=F)

read_cluster_summary<-chr_dna_read_cluster %>% 
  group_by(X6) %>% 
  summarise(min=min(X2),
            max=max(X3),
            n=n())

read_cluster_summary %>% 
  mutate(d=max-min) %>% 
  ggplot(.,aes(n,d))+
  geom_point(size=0.1,alpha=0.1) + 
  scale_x_log10()+
  scale_y_log10()

read_cluster_summary %>% 
#  filter(n>1) %>% 
  mutate(d=max-min) %>% 
  ggplot(.,aes(n,d))+
  geom_density_2d() + 
  scale_x_log10()+
  scale_y_log10()

read_cluster_summary %>% 
  filter(n==1) %>% 
  mutate(d=max-min) %>% 
  ggplot(.,aes(d))+
  geom_density()

read_cluster_summary %>% 
  filter(n==1) %>% 
  ggplot(.,aes(x=min,xend=max,y=n,yend=n,color=n))+
  geom_segment(linewidth=5)
#-----------------------------------------------------
chr_dna_read_dist<-vroom("./../data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_iPSC_chr19_DNA_read_neighbour.tsv",delim = "\t",col_names=F)
chr_dna_read_dist %>% 
  ggplot(.,aes(X7))+
  geom_histogram()
#-----------------------------------------------------
chr_atac_read_cluster<-vroom("./../data/ATAC/FANTOM6/ff_iPSC_ATAC_sorted_cluster_chr19.txt",delim = "\t",col_names=F)

read_cluster_summary<-chr_atac_read_cluster %>% 
  group_by(X13) %>% 
  summarise(min=min(X2),
            max=max(X3),
            n=n())
read_cluster_summary %>% 
  mutate(d=max-min) %>% 
  ggplot(.,aes(n,d))+
  geom_point(size=0.1,alpha=0.1) + 
  scale_x_log10()+
  scale_y_log10()

read_cluster_summary %>% 
  filter(n>1) %>% 
  mutate(d=max-min) %>% 
  ggplot(.,aes(n,d))+
  geom_density_2d_filled() + 
  scale_x_log10()+
  scale_y_log10()
read_cluster_summary %>% 
  mutate(one=ifelse(n>1,"many","single")) %>% 
  ggplot(.,aes(x=min,xend=max,y=n,yend=n,color=one))+
  geom_segment(linewidth=30)+
  scale_y_log10()+
  facet_grid(one~.)
  