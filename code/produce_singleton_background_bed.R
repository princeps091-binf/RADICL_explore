#Produce internal control track
library(vroom)
library(tidyverse)
library(furrr)
cluster_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/"
chr_set<-str_split_fixed(
  str_split_fixed(
    grep("^RAD",list.files(cluster_folder),value = T),
    "\\.",2)[,1],
  "_",4)[,4]

purrr::walk(chr_set,function(chromo){
  message(chromo)
  chr_dna_read_cluster<-vroom(paste0(cluster_folder,"RADICL_iPSC_DNA_",chromo,".bed_cluster.txt"),
                              delim = "\t",
                              col_names=F)
  
  singleton_cluster<-chr_dna_read_cluster %>% 
    group_by(X6) %>% 
    count
  levels<-sort(unique(singleton_cluster$n))+1
  
  chr_dna_read_cluster<-chr_dna_read_cluster %>% 
    left_join(.,singleton_cluster)
  
  cluster_cov_tbl<-map_dfr(levels,function(r){
    tmp_cov<-(chr_dna_read_cluster %>% 
                filter(n<r) %>% 
                nrow)/nrow(chr_dna_read_cluster)
    return(tibble(cluster.size=r-1,cov=tmp_cov))
  })
  
  size_thresh<-cluster_cov_tbl %>% 
    filter(cov > 0.5) %>% 
    slice_min(cluster.size)
  
  bg_radicl_dat<-chr_dna_read_cluster %>% 
    filter(n<=size_thresh$cluster.size) %>% 
    dplyr::select(X1,X2,X3)
  readr::write_tsv(bg_radicl_dat,file = paste0(cluster_folder,"/singleton_data/",chromo,"_DNA_singleton.bed"),col_names = F)
  
})
