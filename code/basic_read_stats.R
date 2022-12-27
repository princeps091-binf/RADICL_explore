library(tidyverse)
library(vroom)
dat_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/"
chromo<-"chr22"
chr_dat<-vroom(paste0(dat_folder,"RADICL_iPSC_chr22.bed"),
               col_names = F)
chr_dat %>% 
  mutate(inter=ifelse(X1==X7,"intra","inter")) %>% 
  group_by(inter) %>% 
  count

chr_set<-str_split_fixed(str_split_fixed(grep("^RADICL",
              list.files("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/"),
              value = T),
              "\\.",
              2)[,1],
              "_",3)[,3] 
inter_vs_intra_tbl<-map_df(chr_set,function(chromo){
  message(chromo)
  return(vroom(paste0(dat_folder,"RADICL_iPSC_",chromo,".bed"),
        col_names = F) %>% 
    mutate(inter=ifelse(X1==X7,"intra","inter")) %>% 
    group_by(inter) %>% 
    count %>% 
      mutate(chr=chromo))
})
inter_vs_intra_tbl %>% 
  ggplot(.,aes(chr,n,fill=inter))+
  geom_bar(stat="identity",position="fill")
#-------------------------------------------------
dat_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/"
radicl_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/"
chromo<-"chr19"
chr_dat<-vroom(paste0(dat_folder,"RADICL_iPSC_",chromo,"_DNA.bed_cluster.txt"),
               col_names = F,delim = "\t")
radicl_dat<-vroom(paste0(radicl_folder,"RADICL_iPSC_",chromo,".bed"),
                  col_names = F,delim = "\t")
single_cluster_tbl<-chr_dat %>% 
  group_by(X6) %>% 
  count %>% 
  filter(n<2)
singleton_tag<-chr_dat %>% 
  filter(X6 %in% single_cluster_tbl$X6) %>% 
  select(X4)
radicl_dat %>% 
  mutate(single= ifelse(X4 %in% singleton_tag$X4,"single","multi")) %>% 
  mutate(inter=ifelse(X1==X7,"intra","inter")) %>% 
  group_by(single,inter) %>% 
  count %>% 
  ggplot(.,aes(single,n,fill=inter))+
  geom_bar(stat='identity')
#----------------------------------
dat_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/"
chr_set<-str_split_fixed(str_split_fixed(grep("^RADICL",
                                              list.files(dat_folder),
                                              value = T),
                                         "\\.",
                                         2)[,1],
                         "_",4)[,4]
single_vs_multi_tbl<-map_df(chr_set,function(chromo){
  message(chromo)
  chr_dat<-vroom(paste0(dat_folder,"RADICL_iPSC_DNA_",chromo,".bed_cluster.txt"),
        col_names = F,delim = "\t")
  single_cluster_tbl<-chr_dat %>% 
    group_by(X6) %>% 
    count %>% 
    filter(n<2)
  return(chr_dat %>% 
    mutate(single= ifelse(X6 %in% single_cluster_tbl$X6,"single","multi")) %>% 
    group_by(single) %>% 
    count %>%
    mutate(chr=chromo))
})
single_vs_multi_tbl %>% 
  ggplot(.,aes(chr,n,fill=single))+
  geom_bar(stat="identity",position="fill")
