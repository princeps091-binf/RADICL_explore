library(tidyverse)
library(GenomicRanges)
library(vroom)
library(purrr)
library(furrr)
library(igraph)
load( "./data/reduced_transcript_to_transcript_ID_tbl.Rda")
load( "./data/peak_transcript_read_inter_tbl.Rda")

peak_transcript_inter_tbl %>% 
  group_by(read.ID) %>% 
  summarise(out=all(is.na(tr.ID))) %>% 
  ungroup %>% 
  group_by(out) %>% 
  count %>% 
  ggplot(.,aes('reads',n,fill=out))+
  geom_bar(stat='identity',position='fill')

peak_transcript_inter_tbl %>% 
  distinct(tr.ID) %>% 
  mutate(in.peaks='in') %>% 
  full_join(.,tibble(tr.ID=1:length(reduced_transcript_Granges))) %>% 
  mutate(in.peaks=ifelse(is.na(in.peaks),'out',in.peaks)) %>% 
  group_by(in.peaks) %>% 
  count %>% 
  ggplot(.,aes('transcript',n,fill=in.peaks))+
  geom_bar(stat='identity',position='fill')

peak_summary_tbl<-peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(out=sum(is.na(tr.ID))/n(),
            ntr=length(unique(tr.ID)),
            nread=n(),
            w=unique(end-start))

peak_summary_tbl %>% 
  ggplot(.,aes(out))+
  geom_density()

peak_summary_tbl %>% 
  ggplot(.,aes(ntr))+
  geom_density()+
  scale_x_log10()

peak_summary_tbl %>% 
  ggplot(.,aes(w))+
  geom_density()+
  scale_x_log10()


peak_summary_tbl %>% 
  mutate(rr=ntr/w) %>% 
  ggplot(.,aes(out,rr))+
  geom_point(size=0.1, alpha=0.1)+
  geom_density_2d(color='red')+
  scale_y_log10()+
  geom_smooth()

peak_summary_tbl %>% 
  ggplot(.,aes(out,w))+
  geom_point(size=0.1, alpha=0.1)+
  geom_density_2d(color='red',alpha=0.5)+
  scale_y_log10()

peak_transcript_inter_tbl %>% 
  filter(!(is.na(tr.ID))) %>% 
  group_by(tr.ID) %>% 
  count %>% 
  ggplot(.,aes(n))+
  geom_density()+
  scale_x_log10()
