# Peak to reads mapping
library(tidyverse)
library(GenomicRanges)
library(vroom)


peak_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/result/DNA/MACS/peaks/"
DNA_read_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/processed/DNA/cluster/"
transcript_annotation_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
black_list_file<-"~/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
chromo<-'chr19'


# Map reads to peak
peak_tbl<-vroom(paste0(peak_folder,'RADICL_DNA_',chromo,'_peaks.bed'),
                col_names = F,
                skip=1,
                delim="\t")
# load Read files
DNA_read_tbl<-vroom(paste0(DNA_read_folder,'RADICL_iPSC_DNA_',chromo,'.bed_cluster.txt'),
                    col_names = F,
                    delim="\t")

# load Black list
black_list<-vroom(black_list_file,
                  col_names = F,
                  delim="\t")

black_list_GRanges<-GRanges(seqnames=black_list$X1,
                      ranges = IRanges(start=black_list$X2,
                                       end=black_list$X3
                      ))


peak_GRanges<-GRanges(seqnames=peak_tbl$X1,
                      ranges = IRanges(start=peak_tbl$X2,
                                       end=peak_tbl$X3
                      ))

read_GRanges<-GRanges(seqnames=DNA_read_tbl$X1,
                      ranges = IRanges(start=DNA_read_tbl$X2,
                                       end=DNA_read_tbl$X3
                      ))
#---------------------------------------------------
# determine cluster size
cluster_size<-DNA_read_tbl %>% 
  group_by(X6) %>% 
  count %>% 
  dplyr::rename(cluster=X6)
mcols(read_GRanges)<-tibble(ID=DNA_read_tbl$X4,cluster=DNA_read_tbl$X6) %>% 
  left_join(.,cluster_size)

# Filter out black list elements

clean_peak_GRange<-peak_GRanges[-unique(queryHits(findOverlaps(peak_GRanges,black_list_GRanges)))]
clean_read_GRanges<-read_GRanges[-unique(queryHits(findOverlaps(read_GRanges,black_list_GRanges)))]

peak_read_content_tbl<-as_tibble(findOverlaps(clean_peak_GRange,clean_read_GRanges)) %>% 
  mutate(start=start(clean_peak_GRange[queryHits]),
         end=end(clean_peak_GRange[queryHits]),
         chr=chromo,
         read.ID=mcols(clean_read_GRanges)$ID[subjectHits],
         read.cluster.size=mcols(clean_read_GRanges)$n[subjectHits]) %>% 
  dplyr::select(chr,start,end,read.ID,read.cluster.size)

peak_read_content_tbl %>% 
  ggplot(.,aes(read.cluster.size))+
  geom_density()+
  scale_x_log10()
peak_read_content_tbl %>% 
  group_by(chr,start,end) %>% 
  count %>% 
  mutate(width=end-start) %>% 
  ggplot(.,aes(width,n))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()
read_set<-peak_read_content_tbl$read.ID
#------------------------------------------------------
# Map reads to transcript
f <- function(x, pos){
  x %>% 
    filter(X4 %in% read_set)
  
}
RNA_side_tbl<-readr::read_delim_chunked('~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/raw/RNA/RADICL_iPSC_RNA.bed',
                                        delim = "\t",
                                        chunk_size = 1e4,
                                        col_names = F,
                                        DataFrameCallback$new(f))

RNA_read_GRanges<-GRanges(seqnames=RNA_side_tbl$X1,
                      ranges = IRanges(start=RNA_side_tbl$X2,
                                       end=RNA_side_tbl$X3
                      ),
                      strand=RNA_side_tbl$X6
)
mcols(RNA_read_GRanges)<-tibble(read.ID=RNA_side_tbl$X4)

clean_RNA_read_GRange<-RNA_read_GRanges[-unique(queryHits(findOverlaps(RNA_read_GRanges,black_list_GRanges)))]

# Map peak to transcripts
transcript_tbl<-vroom(transcript_annotation_file,
                    col_types = list("c",'i','i','c','d','c','i','i','c','d','c','c'),
                    col_names = F,
                    delim="\t")
transcript_GRanges<-GRanges(seqnames=transcript_tbl$X1,
                          ranges = IRanges(start=transcript_tbl$X2,
                                           end=transcript_tbl$X3
                          ),
                          strand=transcript_tbl$X6
)
mcols(transcript_GRanges)<-tibble(ID=transcript_tbl$X4)
reduced_transcript_GRange<-reduce(transcript_GRanges)
#------------------------------------------------------------
peak_read_content_tbl %>% 
  full_join(.,RNA_side_tbl,by=c('read.ID'='X4')) %>% 
  mutate(inter.chr=ifelse(chr==X1,'intra','inter')) %>% 
  group_by(chr,start,end,inter.chr) %>% 
  count
#------------------------------------------------------------

peak_transcript_inter_tbl<-as_tibble(findOverlaps(clean_RNA_read_GRange,reduced_transcript_GRange)) %>% 
  mutate(read.ID=mcols(clean_RNA_read_GRange)$read.ID[queryHits],
         tr.ID=subjectHits) %>% 
  dplyr::select(read.ID,tr.ID) %>% 
  full_join(.,peak_read_content_tbl)

peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(non.tr=sum(is.na(tr.ID))/n(),
            w=unique(end-start)) %>% 
  ggplot(.,aes(non.tr,w))+
  geom_point(size=0.1)+
  geom_vline(xintercept = sum(is.na(peak_transcript_inter_tbl$tr.ID))/nrow(peak_transcript_inter_tbl))+
  scale_y_log10()
peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(non.tr=sum(is.na(tr.ID))/n(),
            w=unique(end-start),
            m.cl.size=max(read.cluster.size)) %>% 
  ggplot(.,aes(m.cl.size,non.tr))+
  geom_point()+
  scale_x_log10()
peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(ntranscript=length(unique(tr.ID)),
            nread=length(unique(read.ID))) %>% 
  ggplot(.,aes(ntranscript))+
  geom_density()+
  scale_x_log10()

peak_transcript_inter_tbl %>% 
  filter(!(is.na(tr.ID))) %>% 
  group_by(chr,start,end,tr.ID) %>% 
  summarise(tr.count=n(),
            w=unique(end-start),
            m=mean(read.cluster.size)) %>% 
  left_join(.,peak_transcript_inter_tbl %>% 
              group_by(chr,start,end) %>% 
              count) %>% 
  ggplot(.,aes(tr.count/n,m))+
  geom_density_2d_filled()+
#  geom_point()+
#  scale_x_log10()+
  scale_y_log10()

