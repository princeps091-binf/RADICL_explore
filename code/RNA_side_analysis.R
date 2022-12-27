library(tidyverse)
library(vroom)

RNA_cluster_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/cluster/RADICL_iPSC_RNA_chr22.bed_cluster.txt",
      col_names = F,delim = "\t")
RNA_cluster_tbl %>% 
  group_by(X7) %>% 
  summarise(start=min(X2),
            end=max(X3),
            n=n()) %>% 
  mutate(io=ifelse(n<2,"single","multi")) %>% 
  ggplot(.,aes(x=start,xend=end,y=1,yend=1,color=log10(n)))+
  geom_segment(linewidth=20)+
  facet_grid(.~io)+
  theme_classic()

RNA_cluster_tbl %>% 
  group_by(X7) %>% 
  summarise(mid=min(X2)+(max(X3)-min(X2))/2,
            n=n()) %>% 
  arrange(mid) %>% 
  ggplot(.,aes(x=mid,y=n))+
  geom_line(linewidth=0.1)+
  scale_y_log10()+
  theme_classic()

cluster_summary_tbl<-RNA_cluster_tbl %>% 
  group_by(X7) %>% 
  summarise(start=min(X2),
            end=max(X3),
            n=n())
cluster_summary_tbl %>% 
  filter(n > 1) %>% 
  ggplot(.,aes(n))+
  geom_density()+
  scale_x_log10()


cluster_summary_tbl %>% 
  filter(n > 1) %>% 
  mutate(d=abs(end-start)) %>% 
  ggplot(.,aes(d))+
  geom_density()+
  scale_x_log10()
