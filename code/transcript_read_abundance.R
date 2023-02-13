# Get transcript-read set genome-wide
library(tidyverse)
library(GenomicRanges)
library(vroom)

transcript_annotation_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
black_list_file<-"~/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
RNA_read_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/raw/RNA/chr/"
#------------------------------------------------
black_list<-vroom(black_list_file,
                  col_names = F,
                  delim="\t")

black_list_GRanges<-GRanges(seqnames=black_list$X1,
                            ranges = IRanges(start=black_list$X2,
                                             end=black_list$X3
                            ))
transcript_tbl<-vroom(transcript_annotation_file,
                      col_types = list("c",'i','i','c','d','c','i','i','c','d','c','c'),
                      col_names = F,
                      delim="\t")

transctipt_GRanges<-GRanges(seqnames=transcript_tbl$X1,
                            ranges = IRanges(start=transcript_tbl$X2,
                                             end=transcript_tbl$X3
                            ),
                            strand=transcript_tbl$X6
)
mcols(transctipt_GRanges)<-tibble(ID=transcript_tbl$X4)
clean_transcript_GRanges<-transctipt_GRanges[-unique(queryHits(findOverlaps(transctipt_GRanges,black_list_GRanges)))]

# process chromosome-wide (autosomes)
chr_set<-unique(str_split_fixed(list.files(RNA_read_folder),"_|\\.",5)[,4])
transcript_read_count<-map_dfr(chr_set,function(chromo){
  message(chromo)
  chr_transcript<-clean_transcript_GRanges[seqnames(clean_transcript_GRanges)==chromo]
  RNA_read_tbl<-vroom(paste0(RNA_read_folder,'RADICL_iPSC_RNA_',chromo,'.bed'),
                      col_names = F,
                      delim="\t")
  read_GRanges<-GRanges(seqnames=RNA_read_tbl$X1,
                              ranges = IRanges(start=RNA_read_tbl$X2,
                                               end=RNA_read_tbl$X3
                              ),
                              strand=RNA_read_tbl$X6
  )
  mcols(read_GRanges)<-tibble(ID=RNA_read_tbl$X4)
  
  return(as_tibble(findOverlaps(chr_transcript,read_GRanges)) %>% 
    mutate(transcript.ID=mcols(chr_transcript)$ID[queryHits],
           read.ID=mcols(read_GRanges)$ID[subjectHits]) %>% 
    dplyr::select(transcript.ID,read.ID) %>% 
    mutate(chr=chromo) %>% 
    group_by(transcript.ID) %>% 
    count)
  
})
# table with trancript and read content
transcript_read_count %>% 
  left_join(.,tibble(width=width(clean_transcript_GRanges),transcript.ID=mcols(clean_transcript_GRanges)$ID)) %>% 
  ggplot(.,aes(n,width))+
  geom_point(size=0.1,alpha=0.2)+
#  geom_density_2d_filled()+
  scale_x_log10()+
  scale_y_log10()

