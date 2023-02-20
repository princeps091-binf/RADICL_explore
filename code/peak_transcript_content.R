library(tidyverse)
library(GenomicRanges)
library(vroom)
library(purrr)
library(furrr)
library(igraph)
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
save(red_tr_to_transcript_tbl,file = "./data/reduced_transcript_to_transcript_ID_tbl.Rda")
save(peak_transcript_inter_tbl,file = "./data/peak_transcript_read_inter_tbl.Rda")


#cumulative transcript curve for peaks sorted by peak content
peak_tr_summary_tbl<-peak_transcript_inter_tbl %>% 
  group_by(chr,start,end) %>% 
  summarise(ntr=list(unique(tr.ID)),
            nread=n(),
            w=unique(end-start)) %>% 
  mutate(rr=nread/w) %>% 
  arrange(desc(rr))
plan(multisession,workers=4)
cummulative_vec<-future_map_dfr(sort(unique(peak_tr_summary_tbl$rr),decreasing = T),function(i){

  cum_tr<-peak_tr_summary_tbl %>% 
    ungroup %>% 
    filter(rr>=i) %>% 
    select(ntr) %>% 
    unnest(ntr) %>% 
    filter(!(is.na(ntr))) %>% 
    distinct() %>% 
    nrow
  return(peak_tr_summary_tbl %>% 
           ungroup %>% 
           filter(rr==i) %>% 
           select(chr,start,end) %>% 
           mutate(cummtr=cum_tr,
                  rr=i))
  
})
plan(sequential)
cummulative_vec %>% 
  mutate(rank=min_rank(1/rr),
         cummcurve=cummtr/length(unique(unlist(peak_tr_summary_tbl$ntr)))) %>% 
  ggplot(.,aes(rank,cummcurve))+
  geom_point()

tibble(rr=1:length(unique(peak_tr_summary_tbl$rr)),cumcurve=cummulative_vec/length(unique(unlist(peak_tr_summary_tbl$ntr)))) %>% 
  ggplot(.,aes(rr,cumcurve))+
  geom_point()+
  ylim(c(0,1))

#cumulative transcript curve for peaks sorted by peak content
tr_peak_summary_tbl<-peak_transcript_inter_tbl %>% 
  group_by(tr.ID,chr,start,end) %>% 
  summarise(peak.read=n()) %>% 
  ungroup() %>% 
  filter(!(is.na(tr.ID))) %>% 
  group_by(tr.ID) %>% 
  summarise(nread=sum(peak.read),
            peak.set=list(unique(paste(chr,start,end,sep='_'))),
            npeak=n())

plan(multisession,workers=4)
cummulative_tbl<-future_map_dfr(sort(unique(tr_peak_summary_tbl$npeak),decreasing = T),function(i){
  cum_peak<-tr_peak_summary_tbl %>% 
    ungroup %>% 
    filter(npeak>=i) %>% 
    select(peak.set) %>% 
    unnest(cols=c(peak.set)) %>% 
    distinct() %>% 
    nrow
  return(tr_peak_summary_tbl %>% 
    ungroup %>% 
    filter(npeak==i) %>% 
    select(tr.ID) %>% 
    mutate(cummpeak=cum_peak,
           npeak=i))
})
plan(sequential)
cummulative_tbl %>% 
  mutate(rank=min_rank(1/npeak),
         cummcurve=cummpeak/length(unique(unlist(tr_peak_summary_tbl$peak.set)))) %>% 
  ggplot(.,aes(rank,cummcurve))+
  geom_point()

#-----------------------------------------------
#Build network from this table (transfer to cluster)
# produce incidence matrix
tr_set<-peak_transcript_inter_tbl %>% 
  distinct(tr.ID) %>% 
  filter(!(is.na(tr.ID))) %>% 
  unlist
distinct_peak_tbl<-peak_transcript_inter_tbl %>% 
  distinct(chr,start,end)
nrow(distinct_peak_tbl)
plan(multisession,workers=3)
incidence_matrix<-do.call(cbind,future_map(1:nrow(distinct_peak_tbl),function(i){
  tmp_peak<-distinct_peak_tbl %>% 
    dplyr::slice(i)
  peak_tr<-peak_transcript_inter_tbl %>% 
    filter(chr== tmp_peak$chr & start == tmp_peak$start & end == tmp_peak$end) %>% 
    distinct(tr.ID) %>% 
    unlist
  return(as.integer(tr_set %in% peak_tr))
}))
plan(sequential)
colnames(incidence_matrix)<-distinct_peak_tbl %>% 
  mutate(ID=paste(chr,start,end,sep='_')) %>% 
  dplyr::select(ID) %>% 
  unlist
rownames(incidence_matrix)<-as.character(1:nrow(incidence_matrix))
save(incidence_matrix,file = "./data/peak_transcript_incidence_mat.Rda")

# convert to bi-partite graph
# bi-partite projection
load("./data/peak_transcript_incidence_mat.Rda")
g2 <- graph_from_incidence_matrix(incidence_matrix)
rm(incidence_matrix)
proj2 <- bipartite_projection(g2,which="true",multiplicity = T)
save(proj2,file = './data/peak_incidence_projection.Rda')
load('./data/peak_incidence_projection.Rda')
