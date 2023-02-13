# exon-intron explore
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(furrr)
#-------------------------------------------------------------
transcript_annotation_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

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

exon_tbl<-transcript_tbl %>% 
  dplyr::select(X4,X1,X2,X3,X6,X10,X11,X12) %>% 
  mutate(RelexonStarts=map(X12,function(x){
    as.integer(unlist(str_split(x,',')))
  }),
  exonLength=map(X11,function(x){
    as.integer(unlist(str_split(x,',')))
  })) %>% 
  unnest(cols=c(RelexonStarts,exonLength)) %>% 
  mutate(exonStart=X2+RelexonStarts) %>% 
  mutate(exonEnd=exonStart+exonLength) %>% 
  dplyr::select(X4,X1,exonStart,exonEnd,X6,X2,X3)

exon_GRanges<-GRanges(seqnames=exon_tbl$X1,
                      IRanges(
                        start=exon_tbl$exonStart,
                        end=exon_tbl$exonEnd,
                      ),
                      strand=exon_tbl$X6)

mcols(exon_GRanges)<-tibble(name=exon_tbl$X4)
#------------------------------------------------------
chromo<-'chr19'
RNA_cluster_tbl<-vroom(paste0("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/cluster/RADICL_iPSC_RNA_",chromo,".bed_cluster.txt"),
                       col_names = F,delim = "\t")

chr_exon<-exon_tbl %>% 
  filter(X1 == chromo)
chr_gene<-transcript_tbl %>% 
  filter(X1 == chromo)
chr_transcript_GRanges<-transctipt_GRanges[seqnames(transctipt_GRanges)==chromo]
chr_exon_GRanges<-exon_GRanges[seqnames(exon_GRanges)==chromo]
# intron
## setdiff for each gene 
plan(multisession,workers=4)
intron_tbl<-future_map_dfr(unique(chr_gene$X4),function(i){
  tmp_gene_GRange<-chr_transcript_GRanges[chr_transcript_GRanges@elementMetadata$ID == i]
  tmp_exon_GRanges<-chr_exon_GRanges[chr_exon_GRanges@elementMetadata$name == i]
  tmp_intron<-setdiff(tmp_gene_GRange,tmp_exon_GRanges)
  tibble(as.data.frame(tmp_intron)) %>% 
    mutate(name=i)
  
})
plan(sequential)
#to modify
exon_intron_tbl<-chr_exon %>% 
  dplyr::select(X4,X1,X6,exonStart,exonEnd) %>% 
  dplyr::rename(seqnames=X1,strand=X6,start=exonStart,end=exonEnd,name=X4) %>% 
  dplyr::mutate(type='exon') %>% 
  bind_rows(.,
            intron_tbl %>% 
              dplyr::select(name,seqnames,strand,start,end) %>% 
              dplyr::mutate(type='intron'))

exon_intron_GRanges<-GRanges(seqnames=exon_intron_tbl$seqnames,
                             IRanges(
                               start=exon_intron_tbl$start,
                               end=exon_intron_tbl$end,
                             ),
                             strand=exon_intron_tbl$strand)

mcols(exon_intron_GRanges)<-tibble(name=exon_intron_tbl$name,type=exon_intron_tbl$type)

read_GRanges<-GRanges(seqnames=RNA_cluster_tbl$X1,
                      ranges = IRanges(start=RNA_cluster_tbl$X2,
                                       end=RNA_cluster_tbl$X3
                      ),
                      strand=RNA_cluster_tbl$X6
)

exon_intron_tbl %>% 
  mutate(read.count=countOverlaps(exon_intron_GRanges,read_GRanges)) %>% 
#  filter(read.count>0) %>% 
  mutate(rd=(read.count+0.1)/(end-start)) %>% 
  ggplot(.,aes(rd,color=type))+
  geom_density()+
  scale_x_log10()

map_dfr(c('exon','intron'),function(i){
  read_cov<-length(unique(queryHits(findOverlaps(read_GRanges,exon_intron_GRanges[exon_intron_GRanges@elementMetadata$type==i]))))/length(read_GRanges)
  gen_cov<-sum(width(exon_intron_GRanges[exon_intron_GRanges@elementMetadata$type==i]))/sum(width(exon_intron_GRanges))
  return(tibble(cov=c(read_cov,gen_cov),set=c('read','foot'),type=i))
}) %>% 
  ggplot(.,aes(set,cov,fill=type))+
  geom_bar(stat='identity')

exon_intron_read_count_tbl<-as_tibble(exon_intron_GRanges) %>% 
  mutate(count=countOverlaps(exon_intron_GRanges,read_GRanges),
         rate=count/width)
gene_summary_tbl<-exon_intron_read_count_tbl%>% 
  group_by(name) %>% 
  summarise(g.count=sum(count),
            g.width=sum(width)) %>% 
  mutate(g.rate=g.count/g.width)

exon_intron_read_count_tbl<-exon_intron_read_count_tbl %>% 
  left_join(.,gene_summary_tbl) %>% 
  mutate(p.val=pmap_dbl(list(count,width,g.rate),function(count,width,g.rate){
    poisson.test(count,width,r=g.rate,alternative='greater')$p.value
  }))

exon_intron_read_count_tbl %>% 
  mutate(zero=ifelse(count==0,'zero','interacting')) %>% 
  group_by(type,zero) %>% 
  count() %>% 
  ggplot(.,aes(type,n,fill=zero))+
  geom_bar(stat='identity',position='fill')

exon_intron_read_count_tbl %>%
  filter(count > 0) %>% 
  ggplot(.,aes(p.val,color=type))+
  geom_density()
