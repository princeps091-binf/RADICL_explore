# Overlap with previously documented significant interactions
library(tidyverse)
library(GenomicRanges)
library(vroom)
library(valr)
library(furrr)
peak_file<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/result/DNA/MACS/peaks/RADICL_DNA_tot_peaks.bed"
significant_interaction_file<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/CHICANE_significant_interaction_3cells.tsv"
interactions_annotation_file<-"~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/interaction_ID.annotation.tsv"
genome_file<-"~/Documents/multires_bhicect/data/hg38.chrom.sizes.txt"
###############################
# Genome-wide features
hg19_coord <- read_delim(genome_file, 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")


sign_inter_tbl<-vroom(significant_interaction_file)
interaction_ID_tbl<-vroom(interactions_annotation_file)
peak_tbl<-vroom(peak_file,col_names = F)

peak_GRanges<-GRanges(seqnames=peak_tbl$X1,
                      ranges = IRanges(start=peak_tbl$X2,
                                       end=peak_tbl$X3
                      ))



bin_spec<-str_split_fixed(interaction_ID_tbl$DNA_bin,'_',3)
DNA_bin_GRanges<-GRanges(seqnames=bin_spec[,1],
                      ranges = IRanges(start=as.integer(bin_spec[,2]),
                                       end=as.integer(bin_spec[,3])-1
                      ))
mcols(DNA_bin_GRanges)<-tibble(ID=interaction_ID_tbl$interaction_ID)


tmp_sig_tbl<-sign_inter_tbl %>% 
  filter(sig_CHICANE_iPSC == 'yes')

sign_DNA_bin_GRanges<-DNA_bin_GRanges[mcols(DNA_bin_GRanges)$ID %in% tmp_sig_tbl$interaction_ID]

length(unique(queryHits(findOverlaps(peak_GRanges,sign_DNA_bin_GRanges))))/length(peak_GRanges)

length(unique(subjectHits(findOverlaps(peak_GRanges,sign_DNA_bin_GRanges))))/length(sign_DNA_bin_GRanges)

sum(width(intersect(peak_GRanges,sign_DNA_bin_GRanges)))/sum(width(peak_GRanges))
sum(width(intersect(peak_GRanges,sign_DNA_bin_GRanges)))/sum(width(reduce(sign_DNA_bin_GRanges)))

#ValR routine
peak_bed_tbl<-peak_tbl %>% 
  dplyr::select(X1,X2,X3) %>% 
  dplyr::rename(chrom=X1,start=X2,end=X3)


plan(multisession,workers=4)
rn_count<-future_map_dbl(1:1000,function(i){
  rn_peak_tbl<-bed_shuffle(peak_bed_tbl,within = T,genome = hg19_coord)
  
  rn_peak_GRanges<-GRanges(seqnames=rn_peak_tbl$chrom,
                           ranges = IRanges(start=rn_peak_tbl$start,
                                            end=rn_peak_tbl$end
                           ))
  
  return(length(unique(subjectHits(findOverlaps(sign_DNA_bin_GRanges,rn_peak_GRanges)))))
  
})
plan(sequential)
tmp_obs<-length(unique(subjectHits(findOverlaps(sign_DNA_bin_GRanges,peak_GRanges))))

tibble(rn=rn_count) %>% 
  ggplot(.,aes(rn))+
  geom_density()+
  geom_vline(xintercept = tmp_obs)

# Interestingly when considering one chromosome at a time, some chromosome don't display an enrichment 
chr_set<-unique(peak_bed_tbl$chrom)

rn_count_tbl<-map_dfr(chr_set,function(chromo){
  message(chromo)
  plan(multisession,workers=4)
  rn_count<-future_map_dbl(1:1000,function(i){
    rn_peak_tbl<-bed_shuffle(peak_bed_tbl %>% filter(chrom==chromo),within = T,genome = hg19_coord)
    
    rn_peak_GRanges<-GRanges(seqnames=rn_peak_tbl$chrom,
                             ranges = IRanges(start=rn_peak_tbl$start,
                                              end=rn_peak_tbl$end
                             ))
    
    return(length(unique(subjectHits(findOverlaps(sign_DNA_bin_GRanges[seqnames(sign_DNA_bin_GRanges)==chromo],rn_peak_GRanges)))))
    
  })
  plan(sequential)
  return(tibble(chrom=chromo,rn=rn_count))
})

obs_count_tbl<-map_dfr(chr_set,function(chromo){
  message(chromo)
  obs_inter<-length(unique(subjectHits(findOverlaps(sign_DNA_bin_GRanges[seqnames(sign_DNA_bin_GRanges)==chromo],peak_GRanges))))
  tibble(chrom=chromo,obs=obs_inter)
})
## Consider the possible depletion in non significant regions of the chromosome
rn_count_tbl %>% 
  ggplot(.,aes(rn,group=chrom))+
  geom_density()+
  geom_vline(data=obs_count_tbl,aes(xintercept = obs,group=chrom))+
  facet_wrap(chrom~.,scales='free')
