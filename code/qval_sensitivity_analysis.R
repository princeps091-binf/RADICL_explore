# Peak RefSeq gene content
library(tidyverse)
library(GenomicRanges)
library(vroom)

peak_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/result/DNA/MACS/peaks/sensitivity_analysis/"
peak_files<-grep(".bed$",list.files(peak_folder),value=T)
read_file<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/processed/DNA/cluster/RADICL_iPSC_DNA_chr19.bed_cluster.txt"
# load the peak file

peak_tbl_l<-map(peak_files,function(peak_file){
  vroom(paste0(peak_folder,peak_file),
                col_names = F,
                delim="\t",
                skip=1)
})
names(peak_tbl_l)<-as.character(10**-as.numeric(str_split_fixed(peak_files,"_",4)[,3]))
# load Read files
DNA_read_tbl<-vroom(read_file,
                    col_names = F,
                    delim="\t")
#----------------------------------------------------
cluster_summary<-DNA_read_tbl %>% 
  group_by(X6) %>% 
  count
read_GRanges<-GRanges(seqnames=DNA_read_tbl$X1,
                      ranges = IRanges(start=DNA_read_tbl$X2,
                                       end=DNA_read_tbl$X3
                      ))
mcols(read_GRanges)<-tibble(ID=DNA_read_tbl$X4,cluster.ID=DNA_read_tbl$X6) %>% 
  left_join(.,cluster_summary,by=c('cluster.ID'='X6'))
walk(peak_tbl_l,function(peak_tbl){
  peak_GRanges<-GRanges(seqnames=peak_tbl$X1,
                        ranges = IRanges(start=peak_tbl$X2,
                                         end=peak_tbl$X3
                        ))
  message(length(findOverlaps(peak_GRanges,read_GRanges))/length(read_GRanges))
  
})

map_dfr(names(peak_tbl_l),function(i){
  peak_tbl<-peak_tbl_l[[i]]
  peak_GRanges<-GRanges(seqnames=peak_tbl$X1,
                        ranges = IRanges(start=peak_tbl$X2,
                                         end=peak_tbl$X3
                        ))
  
  return(as_tibble(findOverlaps(peak_GRanges,read_GRanges)) %>% 
         mutate(n=mcols(read_GRanges)$n[subjectHits]) %>% 
         mutate(set='peak') %>% 
         dplyr::select(n,set) %>% 
         bind_rows(.,as_tibble(mcols(read_GRanges)) %>% 
                     mutate(set='all') %>% 
                     dplyr::select(n,set)
         ) %>% 
           mutate(qval=i))
}) %>% 
  mutate(qval=fct_relevel(qval,as.character(sort(as.numeric(names(peak_tbl_l)))))) %>% 
  ggplot(.,aes(n,color=set))+
  geom_density()+
  scale_x_log10()+
  facet_grid(qval~set,scales='free_y')

peak_read_content_summary<-map_dfr(names(peak_tbl_l),function(i){
  peak_tbl<-peak_tbl_l[[i]]
  peak_GRanges<-GRanges(seqnames=peak_tbl$X1,
                        ranges = IRanges(start=peak_tbl$X2,
                                         end=peak_tbl$X3
                        ))
  
  return(as_tibble(findOverlaps(peak_GRanges,read_GRanges)) %>% 
           mutate(n=mcols(read_GRanges)$n[subjectHits],
                  peak.start=start(peak_GRanges)[queryHits],
                  peak.end=end(peak_GRanges)[queryHits]) %>% 
           mutate(qval=i)) 
})
  
peak_read_content_summary %>% 
  mutate(qval=fct_relevel(qval,as.character(sort(as.numeric(names(peak_tbl_l)))))) %>% 
  group_by(peak.start,peak.end,qval) %>% 
  summarise(noise.prop=sum(n<3)/n(),
            cluster.size=mean(n),
            cluster.min=min(n),
            cluster.max=max(n)) %>% 
  ggplot(.,aes(noise.prop))+
  geom_density()+
  scale_x_sqrt()+
  facet_grid(qval~.,scales='free_y')


