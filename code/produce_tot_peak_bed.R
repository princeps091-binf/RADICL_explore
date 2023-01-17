# join all peak tables
library(tidyverse)
library(GenomicRanges)
library(vroom)

peak_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/MACS_peak/"

peak_files<-grep(".bed$",list.files(peak_folder),value = T)
peak_tbl<-map_dfr(peak_files,function(i){
  vroom(paste0(peak_folder,i),skip = 1,col_names = F)
})

write_tsv(peak_tbl,file = paste0(peak_folder,"RADICL_DNA_tot_peaks.bed"),col_names = F)
