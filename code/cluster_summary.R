library(tidyverse)
library(GenomicRanges)
library(vroom)
read_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/processed/DNA/cluster/"

chr_set<-str_split_fixed(str_split_fixed(list.files(read_folder),'\\.',2)[,1],'_',4)[,4]

cluster_summary_table<-map_dfr(chr_set,function(chromo){
  message(chromo,': loading')
  DNA_read_tbl<-vroom(paste0(read_folder,"RADICL_iPSC_DNA_",chromo,".bed_cluster.txt"),
                      col_names = F,
                      delim="\t")
  message(chromo,': summary')
  return(DNA_read_tbl %>% 
    group_by(X6) %>% 
    summarise(n=n(),
              start=min(c(X2,X3)),
              end=max(c(X2,X3))) %>% 
    mutate(chr=chromo))
  
})
# Fragment size indication?
cluster_summary_table %>% 
  filter(n <2 ) %>% 
  mutate(d=end-start) %>% 
  ggplot(.,aes(d))+
  geom_density()
# Lonely clusters constitute bulk of reads
cluster_summary_table %>% 
  mutate(bg=ifelse(n < 3 ,'noise','other')) %>% 
  group_by(bg,chr) %>%
  summarise(nread=sum(n)) %>% 
  ggplot(.,aes(chr,nread,fill=bg))+
  geom_bar(stat='identity',position='fill')
# Large cluster consitute a distinct sub-population concentrating a noticeable proportion of reads
cluster_summary_table %>% 
  group_by(n) %>% 
  summarise(d=sum(end-start),
            nread=sum(n)/sum(cluster_summary_table$n),
            ncluster=n()) %>%
  mutate(rate=nread/d) %>% 
  ggplot(.,aes(n,nread))+
  geom_line()+
#  scale_x_log10()+
  scale_y_log10()
# Large read clusters do represent regions of high density reads compared to lonely clusters
cluster_summary_table %>% 
  group_by(n) %>% 
  summarise(d=sum(end-start),
            md=mean(end-start),
            nread=sum(n)/sum(cluster_summary_table$n),
            ncluster=n()) %>%
  mutate(rate=nread/d) %>% 
  ggplot(.,aes(n,rate,size=nread))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

