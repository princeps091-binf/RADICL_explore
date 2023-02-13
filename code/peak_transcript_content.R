library(tidyverse)
library(GenomicRanges)
library(vroom)
library(purrr)
library(furrr)
f <- function(x, pos){
  x %>% 
    filter(X4 %in% read_set)
  
}

peak_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/result/DNA/MACS/peaks/"
DNA_read_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/raw/DNA/chr/"
transcript_annotation_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

chr_set<-str_split_fixed(str_split_fixed(list.files(DNA_read_folder),pattern = '_',5)[,4],'\\.',2)[,1]
# Map reads to peak

peak_read_content_tbl<-map_dfr(chr_set,function(chromo){
  message(chromo)
  peak_tbl<-vroom(paste0(peak_folder,'RADICL_DNA_',chromo,'_peaks.bed'),
                col_names = F,
                skip=1,
                delim="\t")
  # load Read files
  DNA_read_tbl<-vroom(paste0(DNA_read_folder,'RADICL_iPSC_DNA_',chromo,'.bed'),
                    col_names = F,
                    delim="\t")

  peak_GRanges<-GRanges(seqnames=peak_tbl$X1,
                      ranges = IRanges(start=peak_tbl$X2,
                                       end=peak_tbl$X3
                      ))

  read_GRanges<-GRanges(seqnames=DNA_read_tbl$X1,
                      ranges = IRanges(start=DNA_read_tbl$X2,
                                       end=DNA_read_tbl$X3
                      ))
  mcols(read_GRanges)<-tibble(ID=DNA_read_tbl$X4)

  return(as_tibble(findOverlaps(peak_GRanges,read_GRanges)) %>% 
    mutate(start=start(peak_GRanges[queryHits]),
         end=end(peak_GRanges[queryHits]),
         chr=chromo,
         read.ID=mcols(read_GRanges)$ID[subjectHits]) %>% 
    dplyr::select(chr,start,end,read.ID))
})
read_set<-peak_read_content_tbl$read.ID

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
reduced_transcript_Granges<-reduce(transcript_GRanges)

red_tr_to_transcript_tbl<-as_tibble(findOverlaps(reduced_transcript_Granges,transcript_GRanges)) %>% 
  mutate(tr.ID=mcols(transcript_GRanges)$ID[subjectHits]) %>% 
  group_by(queryHits) %>% 
  summarise(tr.IDs=list(unique(tr.ID)))
mcols(reduced_transcript_Granges)<-tibble(tr.IDs=red_tr_to_transcript_tbl$tr.IDs)

peak_transcript_inter_tbl<-as_tibble(findOverlaps(RNA_read_GRanges,reduced_transcript_Granges)) %>% 
  mutate(read.ID=mcols(RNA_read_GRanges)$read.ID[queryHits],
         tr.ID=subjectHits) %>% 
  dplyr::select(read.ID,tr.ID) %>% 
  full_join(.,peak_read_content_tbl)

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

peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(out=sum(is.na(tr.ID))/n(),
            ntr=length(unique(tr.ID)),
            nread=n(),
            w=unique(end-start)) %>% 
  ggplot(.,aes(out))+
  geom_density()

peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(out=sum(is.na(tr.ID))/n(),
            ntr=length(unique(tr.ID)),
            nread=n(),
            w=end-start) %>% 
  ggplot(.,aes(ntr))+
  geom_density()+
  scale_x_log10()

peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(out=sum(is.na(tr.ID))/n(),
            ntr=length(unique(tr.ID)),
            nread=n(),
            w=unique(end-start)) %>% 
  ggplot(.,aes(w))+
  geom_density()+
  scale_x_log10()


peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(out=sum(is.na(tr.ID))/n(),
            ntr=length(unique(tr.ID)),
            nread=n(),
            w=unique(end-start)) %>% 
  mutate(rr=ntr/w) %>% 
  ggplot(.,aes(out,rr))+
  geom_point(size=0.1, alpha=0.1)+
  geom_density_2d(color='red')+
  scale_y_log10()+
  geom_smooth()

peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(out=sum(is.na(tr.ID))/n(),
            ntr=length(unique(tr.ID)),
            nread=n(),
            w=unique(end-start)) %>% 
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

#cumulative transcript curve for peaks sorted by peak content
peak_tr_summary_tbl<-peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(ntr=list(unique(tr.ID)),
            nread=n(),
            w=unique(end-start)) %>% 
  mutate(rr=nread/w) %>% 
  arrange(desc(rr))
plan(multisession,workers=4)
cummulative_vec<-future_map_int(sort(unique(peak_tr_summary_tbl$rr),decreasing = T),function(i){
  peak_tr_summary_tbl %>% 
    ungroup %>% 
    filter(rr>=i) %>% 
    select(ntr) %>% 
    unnest(ntr) %>% 
    filter(!(is.na(ntr))) %>% 
    distinct %>% 
    nrow
})
plan(sequential)
tibble(rr=sort(unique(peak_tr_summary_tbl$rr),decreasing = T),cumcurve=cummulative_vec/length(unique(unlist(peak_tr_summary_tbl$ntr)))) %>% 
  ggplot(.,aes(rev(rr),cumcurve))+
  geom_point()
#-----------------------------------------------
#Build network from this table (transfer to cluster)
peak_transcript_inter_tbl
# produce incidence matrix
tr_set<-peak_transcript_inter_tbl %>% 
  distinct(tr.ID) %>% 
  unlist
distinct_peak_tbl<-peak_transcript_inter_tbl %>% 
  distinct(chr,start,end)
# convert to bi-partite graph
# bi-partite projection
g2 <- graph_from_incidence_matrix(M)
g2$name <- "Event network"
proj2 <- bipartite_projection(g2)