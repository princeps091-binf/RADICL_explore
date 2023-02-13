# DNA-cluster distribution for transcript reads
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(furrr)
#-------------------------------------------------------------
transcript_annotation_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"
black_list_file<-"~/Documents/FANTOM6/data/annotation/hg38-blacklist.v2.bed"
RNA_read_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/raw/RNA/chr/"
DNA_read_cluster_folder<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/processed/DNA/cluster/"
#-------------------------------------------------------------
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
#-------------------------------------------------------------
chr_set<-unique(str_split_fixed(list.files(RNA_read_folder),"_|\\.",5)[,4])
# Map reads to transcript
f <- function(x, pos){
  x %>% 
    filter(X4 %in% read_set)
  
}

transcript_read_content<-map_dfr(chr_set,function(chromo){
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
  
  chr_transcript_read_content<-as_tibble(findOverlaps(chr_transcript,read_GRanges)) %>% 
           mutate(transcript.ID=mcols(chr_transcript)$ID[queryHits],
                  read.ID=mcols(read_GRanges)$ID[subjectHits]) %>% 
           dplyr::select(transcript.ID,read.ID) %>% 
           mutate(chr=chromo)
  read_set<-unique(chr_transcript_read_content$read.ID)
  DNA_side_tbl<-readr::read_delim_chunked(paste0(DNA_read_cluster_folder,"RADICL_iPSC_DNA_",chromo,".bed_cluster.txt"),
                                          delim = "\t",
                                          chunk_size = 1e4,
                                          col_names = F,
                                          DataFrameCallback$new(f))
  
})
