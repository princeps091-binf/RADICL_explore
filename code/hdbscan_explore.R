library(vroom)
library(tidyverse)
library(dbscan)
RADICL_inter_chr_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/raw/RADICL_iPSC_inter_chr22_chr21.txt",
                        col_names = F,delim = "\t")

RADICL_inter_chr_tbl %>% 
  select(X2,X3,X8,X9) %>% 
  mutate(RNA.mid=X2 + (X3-X2)/2,
         DNA.mid=X8 + (X8-X9)/2) %>% 
  ggplot(.,aes(RNA.mid,DNA.mid))+
  geom_point(size=0.1,alpha=0.1)
# memory intensive!!! -> run on cluster
inter_db_cluster<-RADICL_inter_chr_tbl %>% 
  select(X2,X3,X8,X9) %>% 
  mutate(RNA.mid=X2 + (X3-X2)/2,
         DNA.mid=X8 + (X8-X9)/2) %>% 
  select(RNA.mid,DNA.mid) %>% 
  hdbscan(.,minPts=100)
RADICL_inter_chr_tbl %>% 
  select(X2,X3,X8,X9) %>% 
  mutate(RNA.mid=X2 + (X3-X2)/2,
         DNA.mid=X8 + (X8-X9)/2,
         clu=inter_db_cluster$cluster,
         out.cl=inter_db_cluster$outlier_scores) %>% 
  filter(out.cl<0.99) %>% 
  ggplot(.,aes(RNA.mid,DNA.mid,color=as.factor(clu)))+
  geom_point(size=0.1)+
  guides(color='none')

RADICL_inter_chr_tbl %>% 
  select(X2,X3,X8,X9) %>% 
  mutate(RNA.mid=X2 + (X3-X2)/2,
         DNA.mid=X8 + (X8-X9)/2,
         clu=inter_db_cluster$cluster,
         out.cl=inter_db_cluster$outlier_scores) %>% 
  ggplot(.,aes(RNA.mid,DNA.mid,color=as.factor(clu)))+
  geom_point(size=0.1)+
  guides(color='none')

